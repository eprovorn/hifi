c-----------------------------------------------------------------------
c     file post2fluid.f.
c     post-processes output from sel code for variuos two-fluid 
c     equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. postpn_mod.
c     1. postpn_read.
c     2. postpn_main.
c     3. postpn_dUdt.
c     4. postpn_contour.
c     5. postpn_energy.
c     6. postpn_maxUvsT.
c     7. postpn_UxyT.
c     8. postpn_Uprofile.
c     9. postpn_width.
c     10. postpn_recrate.
c-----------------------------------------------------------------------
c     subprogram 0. postpn_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE postpn_mod
      USE postxmhd_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: source,cylinder=.FALSE.
      CHARACTER(20), PRIVATE :: init_type,eta_case,kappa_case,atom
      REAL(r8), PARAMETER, PRIVATE :: gamma=5._r8/3._r8,qe=1.602e-19,
     $     lnlam=10.,me=9.109e-31,ep0=8.854e-12,
     $     mp=1.673e-27,gsun=2.74e2,gearth=9.81,nol=1.e-6
      REAL(r8), PRIVATE :: b0,L0,n0,eta,mu,mu_sv,nu,beta0,beta_e,
     $     lambda_psi,lambda_phi,lx,ly,delta=0.,deltatau,kappa_par,
     $     kappa_perp,kappa_min=0.,kappa_max=1.e10,Dn,lwidth,llength,
     $     x1,y1,x2,y2,pmin,rhomin,pmax,rhonmin,di_fac,ion_fac,
     $     recomb_fac,cx_fac,rin_fac,kapn_fac,civ_fac,mfp_fac,bplane,
     $     bguide,initpn,initrho,initrhon,r_eta,etavac,cc,equil_flow,
     $     pert_frac,di,r_en_norm,gravity

      INTEGER, PRIVATE :: ncont
      REAL(r8), PRIVATE :: t_old,tc_old,u_old=0.
      REAL(r8), DIMENSION(2), PRIVATE :: u2_old

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. postpn_read.
c     read necessary post-processing parameters from post.in 
c     and input parameters from sel.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_read(indir)

      CHARACTER(*), INTENT(IN) :: indir

      CHARACTER(80) :: infile
      REAL(r8) :: mi=0.,tnorm=0.
      REAL(r8), DIMENSION(4) :: eta_plist,kappa_plist
      REAL(r8), DIMENSION(5) :: norm_plist

      NAMELIST/pn_input/lwidth,llength,x1,y1,x2,y2
      NAMELIST/pn_ext_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,init_type,source,cylinder,eta,
     $     eta_case,r_eta,etavac,cc,mu,mu_sv,nu,kappa_par,kappa_perp,
     $     kappa_case,kappa_min,kappa_max,Dn,beta0,beta_e,rhomin,pmin,
     $     pmax,rhonmin,n0,b0,L0,lx,ly,lambda_psi,lambda_phi,delta,
     $     deltatau,atom,di_fac,ion_fac,recomb_fac,cx_fac,rin_fac,
     $     kapn_fac,civ_fac,mfp_fac,initrhon,initpn,initrho,equil_flow,
     $     pert_frac,bplane,bguide
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="post.in",STATUS="OLD")
      READ(in_unit,NML=pn_input)
      CLOSE(UNIT=in_unit)
      infile=TRIM(indir)//"/sel.in"
      SELECT CASE(job_type)
      CASE("pn_ext")
         CALL system("rm -f Lprofile.bin Wprofile.bin Lprofile.txt 
     $        Wprofile.txt VecSqInt.txt UxyT.txt layer.txt")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=pn_ext_list)
         CLOSE(UNIT=in_unit)
         ncont=8
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF

         IF(atom=="hydrogen")mi = mp
         tnorm=L0*SQRT(mu0*n0*mi)/b0
         di = di_fac*2.28e8*SQRT(mi/(mp*n0))/L0
         r_en_norm = 1.e-19*two*di*n0*L0*SQRT(me/(pi*mi))

         gravity=0.
         SELECT CASE(init_type)
         CASE("RTsun")
            beta_e=half*beta0
            gravity=gsun*tnorm**2/L0
         CASE("RTearth")
            beta_e=half*beta0
            gravity=gearth*tnorm**2/L0
         END SELECT

         WRITE(out_unit,*)"tnorm=",tnorm
         WRITE(out_unit,*)"di=",di
         WRITE(out_unit,*)"n0=",n0
         WRITE(out_unit,*)"L0=",L0
         WRITE(out_unit,*)"me=",me
         WRITE(out_unit,*)"mi=",mi
         WRITE(out_unit,*)"r_en_norm=",r_en_norm

         norm_plist = (/b0,L0,n0,beta_e/beta0,1._r8/)
         eta_plist = (/eta,cc,etavac,r_eta/)
         kappa_plist = (/kappa_par,kappa_perp,kappa_min,kappa_max/)
         CALL postxmhd_init(indir,norm_plist,eta_case,eta_plist,
     $        kappa_case,kappa_plist)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_read
c-----------------------------------------------------------------------
c     subprogram 2. postpn_main.
c     call the desired specialized postpn_mod subroutines for 
c     post-processing plasma-neutral simulation results  
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_main(t,jac,xyw,xyw_kt,uw,uxyw,flag1d,flag2d,
     $     beginner,header,footer,nt_next,ifile,stride,filename,nxs)

      LOGICAL, INTENT(IN) :: flag1d,flag2d,beginner,header,footer
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nt_next,ifile,stride,nxs
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: jac
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: xyw_kt,uxyw
c-----------------------------------------------------------------------
c     call subroutines.
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("pn_ext")
         IF(flag1d)THEN
            CALL postpn_UxyT(t,xyw,xyw_kt,uw)
            CALL postpn_maxUvsT(t,uw,uxyw)
            CALL postpn_dUdt(t,xyw,xyw_kt,uw,beginner)
            CALL postpn_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL postpn_width(nxs,t,xyw,xyw_kt,uw,uxyw)
            CALL postpn_energy(t,xyw,jac,uw,uxyw)
         ENDIF
         IF(flag2d)THEN
            CALL postpn_contour(t,xyw,uw,uxyw,beginner,
     $           header,footer,nt_next,ifile,stride,filename)
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_main
c-----------------------------------------------------------------------
c     subprogram 3. postpn_dUdt.
c     generate du/dt(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_dUdt(t,xyw,xyw_kt,uw,beginner)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: xyw_kt
      LOGICAL, INTENT(IN) :: beginner

      LOGICAL :: err=.FALSE.
      REAL(r8), DIMENSION(1) :: u,value1,value2
      REAL(r8), DIMENSION(2) :: u2
      REAL(r8), DIMENSION(1,SIZE(uw,2),SIZE(uw,3)) :: utemp
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("pn_ext")
         utemp(1,:,:)=uw(2,:,:)
         CALL plotter_UxyT(t,nol,nol,xyw,xyw_kt,utemp,.FALSE.,value1,
     $        err)
         value1(1)=value1(1)+delta
         IF(value1(1) /= 0)THEN
            u2(1)=LOG(ABS(value1(1)))
         ELSE
            u2(1)=0
         ENDIF
         u2(2)=value1(1)
         IF(err)u2=u2_old
         CALL plotter_dUdt(t,t_old,u2,u2_old,beginner)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_dUdt
c-----------------------------------------------------------------------
c     subprogram 4. postpn_contour.
c     generates contours of desired physical quantities
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_contour(t,xyw,uw,uxyw,beginner,header,
     $     footer,nt_next,ifile,stride,filename)

      LOGICAL, INTENT(IN) :: beginner,header,footer
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nt_next,ifile,stride
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: uxyw

      INTEGER :: ix,nyw,cyl_fac
      INTEGER, DIMENSION(3), SAVE :: xcont=0
      REAL(r8), DIMENSION(SIZE(xyw,2),SIZE(xyw,3)) :: r_faci,rho
      REAL(r8), DIMENSION(ncont,SIZE(uw,2),SIZE(uw,3)) :: utemp
c-----------------------------------------------------------------------
c     zero utemp.
c-----------------------------------------------------------------------
      utemp=0.
c-----------------------------------------------------------------------
c     transform variables.
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("pn_ext")
         cyl_fac=0
         IF(cylinder)cyl_fac=1
c-----------------------------------------------------------------------
c     calculate Bx and By
c-----------------------------------------------------------------------
         utemp(1,:,:) = -(cyl_fac*uw(2,:,:)/xyw(2,:,:) + uxyw(:,:,2,2))
         utemp(2,:,:) = uxyw(:,:,1,2)
         WHERE(xyw(2,:,:) < 1.e-12) 
            utemp(1,:,:) = -uxyw(:,:,2,2)
         END WHERE
c-----------------------------------------------------------------------
c     calculate Tp
c-----------------------------------------------------------------------
         utemp(3,:,:)=uw(8,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate Tn
c-----------------------------------------------------------------------
         utemp(4,:,:)=uw(13,:,:)/uw(9,:,:)
c-----------------------------------------------------------------------
c     calculate Vpx
c-----------------------------------------------------------------------
         utemp(5,:,:) = uw(4,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate Vpy
c-----------------------------------------------------------------------
         utemp(6,:,:) = uw(5,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate Vnx
c-----------------------------------------------------------------------
         utemp(7,:,:) = uw(10,:,:)/uw(9,:,:)
c-----------------------------------------------------------------------
c     calculate Vny
c-----------------------------------------------------------------------
         utemp(8,:,:) = uw(11,:,:)/uw(9,:,:)
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
         DO ix=1,ncont
            IF (MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
         ENDDO
         CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $        nt_next,xcont,ifile,stride,filename,"Ucontour")
c-----------------------------------------------------------------------
c     calculate and plot spatially varying resistivity
c-----------------------------------------------------------------------
         CALL postxmhd_ploteta(t,xyw,uw(1,:,:),uw(8,:,:),uw(7,:,:),
     $        header,footer,nt_next,ifile,stride,filename)
c-----------------------------------------------------------------------
c     calculate and plot spatially varying thermal conductivity
c-----------------------------------------------------------------------
         CALL postxmhd_plotkappa(t,xyw,uw(1,:,:),uw(8,:,:),
     $        utemp(1,:,:)**2+utemp(2,:,:)**2+uw(3,:,:)**2,
     $        header,footer,nt_next,ifile,stride,filename)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_contour
c-----------------------------------------------------------------------
c     subprogram 5. postpn_energy.
c     generate a time plot for an integral of total energy
c     over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_energy(t,xyw,jac,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: jac
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: uxyw

      INTEGER :: nx,ny
      REAL(r8) :: cyl_fac
      REAL(r8), DIMENSION(6) :: value
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uw1,uw2,uw3
      REAL(r8), DIMENSION(SIZE(xyw,2),SIZE(xyw,3)) ::
     $     r_fac,r_facs,r_faci
c-----------------------------------------------------------------------
c     prepare auxiliary parameters
c-----------------------------------------------------------------------
      value=zero
      nx=SIZE(xyw,2)
      ny=SIZE(xyw,3)
c-----------------------------------------------------------------------
c     calculate energy integral
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("pn_ext")
c-----------------------------------------------------------------------
c     compute radius factor for cylindrical coordinates
c-----------------------------------------------------------------------
         cyl_fac=0
         r_fac=1.
         r_faci=1.
         r_facs=1.
         IF(cylinder)THEN
            cyl_fac=1.
            r_fac=xyw(2,:,:)
            r_faci=1._r8/r_fac
            r_facs=SQRT(r_fac)
            WHERE(r_fac < 1.e-16)
               r_faci = 0.
            END WHERE
         ENDIF

         ALLOCATE(uw1(3,nx,ny))
         ALLOCATE(uw2(3,nx,ny))
         ALLOCATE(uw3(3,nx,ny))
c-----------------------------------------------------------------------
c     ion kinetic energy components.
c-----------------------------------------------------------------------
         uw1(1,:,:)=r_facs*uw(4,:,:)/SQRT(uw(1,:,:))
         uw2(1,:,:)=r_facs*uw(5,:,:)/SQRT(uw(1,:,:))
         uw3(1,:,:)=r_facs*uw(6,:,:)/SQRT(uw(1,:,:))
         WHERE(uw(1,:,:) <= 0)
            uw1(1,:,:)=0
            uw2(1,:,:)=0
            uw3(1,:,:)=0
         END WHERE
c-----------------------------------------------------------------------
c     magnetic energy components.
c-----------------------------------------------------------------------
         uw1(2,:,:)=-(cyl_fac*r_faci*uw(2,:,:) + uxyw(:,:,2,2))*r_facs
         uw2(2,:,:)=uxyw(:,:,1,2)*r_facs
         uw3(2,:,:)=uw(3,:,:)*r_facs
c-----------------------------------------------------------------------
c     neutral kinetic energy components.
c-----------------------------------------------------------------------
         uw1(3,:,:)=r_facs*uw(10,:,:)/SQRT(uw(9,:,:))
         uw2(3,:,:)=r_facs*uw(11,:,:)/SQRT(uw(9,:,:))
         uw3(3,:,:)=r_facs*uw(12,:,:)/SQRT(uw(9,:,:))
         WHERE(uw(9,:,:) <= 0)
            uw1(3,:,:)=0
            uw2(3,:,:)=0
            uw3(3,:,:)=0
         END WHERE
c-----------------------------------------------------------------------
c     calculate the integral of vector energies over the domain.
c-----------------------------------------------------------------------
         CALL plotter_VecSqInt(t,xyw,uw1,uw2,uw3,.FALSE.,value(1:3))

         uw1=0
         uw2=0
         uw3=0
c-----------------------------------------------------------------------
c        ion + electron & neutral thermal energy.
c-----------------------------------------------------------------------
         uw1(1,:,:)=r_fac*3._r8*uw(8,:,:)
         uw1(2,:,:)=r_fac*3._r8*uw(13,:,:)

         CALL plotter_integral(t,xyw,uw1(1:2,:,:),.FALSE.,value(4:5))
c-----------------------------------------------------------------------
c        total energy.
c-----------------------------------------------------------------------
         value(6)=SUM(value)
         value=0.5*value
         IF(cylinder)value=value*twopi

         DEALLOCATE(uw1,uw2,uw3)
      END SELECT
c-----------------------------------------------------------------------
c     open, write, and close VecSqInt.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      WHERE(value == 0)
         value=1
      END WHERE
      OPEN(UNIT=VecSqInt_unit,FILE="VecSqInt.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(VecSqInt_unit)REAL(t,4),REAL(LOG(value),4)
      CLOSE(UNIT=VecSqInt_unit)
      OPEN(UNIT=VecSqInt_unit,FILE="VecSqInt.txt",STATUS="UNKNOWN",
     $     POSITION="APPEND")
      WRITE(VecSqInt_unit,'(1p,7(e15.5))')t,LOG(value)
      CLOSE(UNIT=VecSqInt_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_energy
c-----------------------------------------------------------------------
c     subprogram 6. postpn_maxUvsT.
c     generate u(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_maxUvsT(t,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
      REAL(r8), DIMENSION(:), ALLOCATABLE :: value
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("pn_ext")
         ALLOCATE(utemp(6,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),value(6))
c-----------------------------------------------------------------------
c     plasma density
c-----------------------------------------------------------------------
         utemp(1,:,:)=uw(1,:,:)
c-----------------------------------------------------------------------
c     y-velocity
c-----------------------------------------------------------------------
         utemp(2,:,:)=uw(5,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     ionization fraction
c-----------------------------------------------------------------------
         utemp(3,:,:)=uw(1,:,:)/(uw(1,:,:)+uw(9,:,:))
c-----------------------------------------------------------------------
c     ion temperature
c-----------------------------------------------------------------------
         utemp(4,:,:)=uw(8,:,:)*(1. - beta_e/beta0)/uw(1,:,:)
c-----------------------------------------------------------------------
c     neutral temperature
c-----------------------------------------------------------------------
         utemp(5,:,:)=uw(13,:,:)/uw(9,:,:)
c-----------------------------------------------------------------------
c     jz
c-----------------------------------------------------------------------
         utemp(6,:,:)=uw(7,:,:)

         CALL plotter_maxUvsT("maxUvsT",t,utemp,value)
         DEALLOCATE(utemp,value)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_maxUvsT
c-----------------------------------------------------------------------
c     subprogram 7. postpn_UxyT.
c     generate U vs. time plot at a point.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_UxyT(t,xyw,xyw_kt,uw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt

      LOGICAL :: err=.FALSE.
      REAL(r8), DIMENSION(5) :: value1,value2
      REAL(r8), DIMENSION(5,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: utemp
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("pn_ext")
c-----------------------------------------------------------------------
c     poloidal magnetic flux
c-----------------------------------------------------------------------
         utemp(1,:,:)=uw(2,:,:)
c-----------------------------------------------------------------------
c     out-of-plane current density
c-----------------------------------------------------------------------
         utemp(2,:,:)=uw(7,:,:)
c-----------------------------------------------------------------------
c     plasma density
c-----------------------------------------------------------------------
         utemp(3,:,:)=uw(1,:,:)
c-----------------------------------------------------------------------
c     plasma pressure
c-----------------------------------------------------------------------
         utemp(4,:,:)=uw(8,:,:)
c-----------------------------------------------------------------------
c     neutral density
c-----------------------------------------------------------------------
         utemp(5,:,:)=uw(9,:,:)
         CALL plotter_UxyT(t,nol,nol,xyw,xyw_kt,utemp,.FALSE.,value1,
     $        err)
         CALL plotter_UxyT(t,x2,y2,xyw,xyw_kt,utemp,.FALSE.,value2,
     $        err)
c-----------------------------------------------------------------------
c     calculate flux difference between (0,0) and (x2,y2)
c-----------------------------------------------------------------------
         value1(1)=value1(1)-value2(1)
         OPEN(UNIT=UxyT_unit,FILE="UxyT.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(UxyT_unit)REAL(t,4),REAL(value1,4)
         CLOSE(UNIT=UxyT_unit)
         OPEN(UNIT=UxyT_unit,FILE="UxyT.txt",STATUS="UNKNOWN",
     $        POSITION="APPEND")
         WRITE(UxyT_unit,'(1p,6(e15.5))')t,value1
         CLOSE(UNIT=UxyT_unit)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_UxyT
c-----------------------------------------------------------------------
c     subprogram 8. postpn_Uprofile.
c     generate slices of U for each time-step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)

      INTEGER, INTENT(IN) :: nxs
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt,uxyw

      CHARACTER(16) :: filename
      INTEGER :: cyl_fac=0
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: uprofile
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      filename="Uprofile"
      SELECT CASE(job_type)
      CASE("pn_ext")
         IF(cylinder)cyl_fac=1
         ALLOCATE(utemp(13,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(14,0:nxs))
c-----------------------------------------------------------------------
c     Bx
c-----------------------------------------------------------------------
         utemp(3,:,:) = -(cyl_fac*uw(2,:,:)/xyw(2,:,:) + uxyw(:,:,2,2))
         WHERE(xyw(2,:,:) < 1.e-12)
            utemp(3,:,:) = -uxyw(:,:,2,2)
         END WHERE
c-----------------------------------------------------------------------
c     By
c-----------------------------------------------------------------------
         utemp(4,:,:) = uxyw(:,:,1,2)
c-----------------------------------------------------------------------
c     Jz*Bx - d(Pp)/dy
c-----------------------------------------------------------------------
         utemp(1,:,:) = uw(7,:,:)*utemp(3,:,:) - uxyw(:,:,2,8)
c-----------------------------------------------------------------------
c     plasma temperature
c-----------------------------------------------------------------------
         utemp(2,:,:) = uw(8,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     Jz
c-----------------------------------------------------------------------
         utemp(5,:,:) = uw(7,:,:)
c-----------------------------------------------------------------------
c     plasma density
c-----------------------------------------------------------------------
         utemp(6,:,:) = uw(1,:,:)
c-----------------------------------------------------------------------
c     neutral density
c-----------------------------------------------------------------------
         utemp(7,:,:) = uw(9,:,:)
c-----------------------------------------------------------------------
c     Vp_x
c-----------------------------------------------------------------------
         utemp(8,:,:) = uw(4,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     Vn_x
c-----------------------------------------------------------------------
         utemp(9,:,:) = uw(10,:,:)/uw(9,:,:)
c-----------------------------------------------------------------------
c     Vp_y
c-----------------------------------------------------------------------
         utemp(10,:,:) = uw(5,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     Vn_y
c-----------------------------------------------------------------------
         utemp(11,:,:) = uw(11,:,:)/uw(9,:,:)
c-----------------------------------------------------------------------
c     Vp_x-Vn_x
c-----------------------------------------------------------------------
         utemp(12,:,:) = utemp(8,:,:) - utemp(9,:,:)
c-----------------------------------------------------------------------
c     Vp_y-Vn_y
c-----------------------------------------------------------------------
         utemp(13,:,:) = utemp(10,:,:) - utemp(11,:,:)

         CALL plotter_Uprofile(nol,nol,nol,y1,nxs,xyw,
     $        xyw_kt,utemp,uprofile,.TRUE.,"Wprofile")
         CALL plotter_Uprofile(nol,nol,x1,nol,nxs,xyw,
     $        xyw_kt,utemp,uprofile,.TRUE.,"Lprofile")
      END SELECT
      DEALLOCATE(utemp,uprofile)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_Uprofile
c-----------------------------------------------------------------------
c     subprogram 9. postpn_width.
c     generate peak width vs. time plot along a profile.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_width(nxs,t,xyw,xyw_kt,uw,uxyw)

      INTEGER, INTENT(IN) :: nxs
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt,uxyw

      LOGICAL :: err=.FALSE.
      CHARACTER(10) :: wtype
      INTEGER, DIMENSION(2) :: mindex
      REAL(r8) :: x0,y0,width1,width2,width3,ratio
      REAL(r8), DIMENSION(1) :: value1,value2,value3,value4,value5,
     $     value6,value7,value8
      REAL(r8), DIMENSION(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: utemp
      REAL(r8), DIMENSION(6,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: ueta
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      CALL postxmhd_geteta(xyw,uw(1,:,:),uw(8,:,:),uw(7,:,:),ueta)

      SELECT CASE(job_type)
      CASE("pn_ext")
c-----------------------------------------------------------------------
c     Jz profile half-width and half-length at half-max
c-----------------------------------------------------------------------
         mindex=MAXLOC(ABS(uw(7,:,:)))
         x0=xyw(1,mindex(1),mindex(2))
         y0=xyw(2,mindex(1),mindex(2))

         wtype="half_max"
         utemp(1,:,:)=uw(7,:,:)
         CALL plotter_width(x0,y0,x0,lwidth,nxs,xyw,xyw_kt,
     $        utemp,wtype,width1)
         CALL plotter_width(nol,nol,llength,nol,nxs,xyw,xyw_kt,
     $        utemp,wtype,width2)
c-----------------------------------------------------------------------
c     current sheet aspect ratio
c-----------------------------------------------------------------------
         ratio=width2/width1
         IF(width1 <= 0)ratio=0
c-----------------------------------------------------------------------
c     half-length based on plasma outflow peak
c-----------------------------------------------------------------------
         mindex=MAXLOC(ABS(uw(4,:,:)/uw(1,:,:)))
         width3=xyw(1,mindex(1),mindex(2))
c-----------------------------------------------------------------------
c     Bx at Jz half-width
c-----------------------------------------------------------------------
         utemp(1,:,:)=ABS(uxyw(:,:,2,2))
         CALL plotter_UxyT
     $        (t,x0,width1,xyw,xyw_kt,utemp,.FALSE.,value1,err)
c-----------------------------------------------------------------------
c     By at Jz half-length
c-----------------------------------------------------------------------
         utemp(1,:,:)=ABS(uxyw(:,:,1,2))
         CALL plotter_UxyT
     $        (t,width2,y0,xyw,xyw_kt,utemp,.FALSE.,value2,err)
c-----------------------------------------------------------------------
c     Vy at Jz half-width
c-----------------------------------------------------------------------
         utemp(1,:,:)=ABS(uw(5,:,:)/uw(1,:,:))
         CALL plotter_UxyT
     $        (t,x0,width1,xyw,xyw_kt,utemp,.FALSE.,value3,err)
c-----------------------------------------------------------------------
c     maximum Vx
cccc     Vx at Jz half-length
c-----------------------------------------------------------------------
         value4(1)=MAXVAL(ABS(uw(4,:,:)/uw(1,:,:)))
ccc         utemp(1,:,:)=ABS(uw(4,:,:)/uw(1,:,:))
ccc         CALL plotter_UxyT
ccc     $        (t,width2,y0,xyw,xyw_kt,utemp,.FALSE.,value4,err)
c-----------------------------------------------------------------------
c     maximum Jz
c-----------------------------------------------------------------------
         value5=MAXVAL(ABS(uw(7,:,:)))
c-----------------------------------------------------------------------
c     plasma density at peak Jz
c-----------------------------------------------------------------------
         utemp(1,:,:)=uw(1,:,:) + uw(9,:,:)
         CALL plotter_UxyT
     $        (t,x0,y0,xyw,xyw_kt,utemp,.FALSE.,value6,err)
c-----------------------------------------------------------------------
c     effective resistivity at peak Jz
c-----------------------------------------------------------------------
         utemp(1,:,:) = ABS(ueta(1,:,:)*uw(7,:,:)
     $        + r_en_norm*SQRT(beta_e/beta0*uw(8,:,:)/uw(1,:,:) 
     $        + uw(13,:,:)/uw(9,:,:))
     $        *(uw(12,:,:) - uw(9,:,:)/uw(1,:,:)
     $        *(uw(6,:,:) - di*uw(7,:,:))))
         CALL plotter_UxyT
     $        (t,x0,y0,xyw,xyw_kt,utemp,.FALSE.,value7,err)
         value7=value7/value5
c-----------------------------------------------------------------------
c     recrate normalized by upstream Bx and c.s. plasma density
c-----------------------------------------------------------------------
         value8 = 0.
         WHERE(value1 > 1.e-6)
            value8=value7*value5/(value1**2/SQRT(value6))
         END WHERE
c-----------------------------------------------------------------------
c     output
c-----------------------------------------------------------------------
         OPEN(UNIT=width_unit,FILE="layer.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(width_unit)REAL(t,4),REAL(width1,4),REAL(width2,4),
     $        REAL(ratio,4),REAL(width3,4),REAL(value1,4),
     $        REAL(value2,4),REAL(value3,4),REAL(value4,4),
     $        REAL(value5,4),REAL(value6,4),REAL(value7,4),
     $        REAL(value8,4)
         CLOSE(UNIT=width_unit)
         OPEN(UNIT=width_unit,FILE="layer.txt",STATUS="UNKNOWN",
     $        POSITION="APPEND")
         WRITE(width_unit,'(1p,13(e15.5))')t,width1,width2,ratio,width3,
     $        value1,value2,value3,value4,value5,value6,value7,value8
         CLOSE(UNIT=width_unit)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_width
c-----------------------------------------------------------------------
c     subprogram 10. postpn_recrate.
c     generates reconnection rate plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postpn_recrate(t,xyw,xyw_kt,uw,first)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt
      LOGICAL, INTENT(IN) :: first

      CHARACTER(16) :: filename
      INTEGER, PARAMETER :: nxs=500
      REAL(r8) :: value
      REAL(r8), DIMENSION(2) :: u2
      REAL(r8), DIMENSION(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: utemp
      REAL(r8), DIMENSION(2,0:nxs) :: uprofile
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      filename=" "
      SELECT CASE(job_type)
      CASE("pn_ext")
         utemp(1,:,:)=uw(2,:,:)
         CALL plotter_Uprofile(nol,nol,lx-nol,nol,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.FALSE.,filename)
         value=MAXVAL(uprofile(2,:))-MINVAL(uprofile(2,:))
         IF(value /= 0)THEN
            u2(1)=LOG(ABS(value))
         ELSE
            u2(1)=0
         ENDIF
         u2(2)=ABS(value)
         CALL plotter_dUdt(t,t_old,u2,u2_old,first)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postpn_recrate
      END MODULE postpn_mod
