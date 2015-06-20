c-----------------------------------------------------------------------
c     file post2fluid.f.
c     post-processes output from sel code for variuos two-fluid 
c     equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. post2fluid_mod.
c     1. post2fluid_read.
c     2. post2fluid_main.
c     3. post2fluid_dUdt.
c     4. post2fluid_contour.
c     5. post2fluid_energy.
c     6. post2fluid_maxUvsT.
c     7. post2fluid_UxyT.
c     8. post2fluid_Uprofile.
c     9. post2fluid_width.
c     10. post2fluid_recrate.
c     11. post2fluid_frc.
c-----------------------------------------------------------------------
c     subprogram 0. post2fluid_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE post2fluid_mod
      USE postxmhd_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: source,eta_spitzer=.TRUE.,if_kinvisc,
     $     ifascii,flux_inflow,cylinder=.FALSE.
      CHARACTER(20), PRIVATE :: dUdt_name,VecSqInt_name,UxyT_name,
     $     init_type,equilfile,coilfile,interp,equil_type,
     $     eta_case,kappa_case,eta_type
      REAL(r8), PARAMETER, PRIVATE :: mass_r=5.44617e-4_r8,nol=1.e-6
      REAL(r8), PRIVATE :: d_inv,b0,eta,mu,nu,di,bound_eps,tau,Bguide,
     $     beta0,skin,sound_rad,L0,lambda_psi,lambda_phi,lx,ly,mach,
     $     alfven,epsilon,kappa_par,kappa_perp,beta_e,Dn,x1,y1,x2,y2,
     $     de_sq,pmin,rhomin,n0,T0,r_eta,pmax,readj,ddiff,lambda,
     $     ion_fac,recomb_fac,mi,ion_efac,recomb_efac,initp,initrho,
     $     initrhon,bx,etavac,zmin,psource,chod_const,rhomax,p0,j_crit,
     $     Te_frac,delx,dely

      INTEGER, PRIVATE :: ncont
      REAL(r8), PRIVATE :: t_old,tc_old,u_old=0.
      REAL(r8), DIMENSION(2), PRIVATE :: u2_old
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: uc_old

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. post2fluid_read.
c     read necessary post-processing parameters from post.in 
c     and input parameters from sel.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_read(indir)

      CHARACTER(*), INTENT(IN) :: indir
      CHARACTER(80) :: infile

      REAL(r8), DIMENSION(4) :: eta_plist,kappa_plist
      REAL(r8), DIMENSION(5) :: norm_plist

      NAMELIST/twofluid_input/x1,y1,x2,y2,dUdt_name,UxyT_name,
     $     VecSqInt_name,ifascii
      NAMELIST/fivefield_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,di,eta,mu,nu,lx,ly,alfven,mach,
     $     lambda_psi,lambda_phi,epsilon,bound_eps,tau,Bguide,source,
     $     init_type
      NAMELIST/HallMHD_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,di,eta,mu,nu,kappa_par,kappa_perp,
     $     eta_spitzer,Dn,lx,ly,mach,alfven,lambda_psi,lambda_phi,
     $     epsilon,Bguide,beta_e,beta0,source,init_type,rhomin,pmin,
     $     equilfile,b0,T0,n0,coilfile,r_eta,interp,pmax,equil_type,bx,
     $     zmin,etavac,eta_case,flux_inflow,L0,cylinder,p0
      NAMELIST/rmhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,gr_curve,eta,mu,nu,de_sq,sound_rad,tau,source,
     $     lx,ly,lambda_psi,lambda_phi,epsilon,bound_eps,mach,alfven,
     $     init_type,eta_type,j_crit,delx,dely
      NAMELIST/emhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,gr_curve,di,eta,nu,lx,ly,alfven,lambda_psi,
     $     source,epsilon,init_type
      NAMELIST/cmhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,eta,mu,nu,lx,ly,lambda_psi,
     $     kappa_par,kappa_perp,epsilon,beta0,source,init_type,
     $     rhomin,pmin,equilfile,b0,T0,n0,L0,coilfile,r_eta,interp,
     $     if_kinvisc,pmax,readj,equil_type,ddiff,mi,bx,etavac,
     $     eta_case,zmin,cylinder,flux_inflow,kappa_case
      NAMELIST/cmhdn_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,eta,mu,lx,ly,lambda_psi,
     $     kappa_par,kappa_perp,epsilon,beta0,source,init_type,
     $     rhomin,pmin,equilfile,b0,T0,n0,L0,coilfile,r_eta,interp,
     $     if_kinvisc,pmax,readj,equil_type,ddiff,mi,bx,etavac,eta_case,
     $     zmin,cylinder,flux_inflow,chod_const,ion_fac,recomb_fac,
     $     ion_efac,recomb_efac,initp,initrho,initrhon,psource,zmin
      NAMELIST/epmhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,gr_curve,eta,mu,de_sq,tau,source,lx,ly,
     $     lambda_psi,lambda_phi,epsilon,bound_eps,mach,alfven,init_type
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="post.in",STATUS="OLD")
      READ(in_unit,NML=twofluid_input)
      CLOSE(UNIT=in_unit)
      infile=TRIM(indir)//"/sel.in"
      SELECT CASE(job_type)
      CASE("fivefield")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=fivefield_list)
         CLOSE(UNIT=in_unit)
         ncont=5
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
      CASE("fourfieldplus")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=fivefield_list)
         CLOSE(UNIT=in_unit)
         ncont=9
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
      CASE("HallMHD")
         CALL system("rm -f frc.bin frc.asc VecSqInt.asc
     $        screw.bin screw.asc")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=HallMHD_list)
         CLOSE(UNIT=in_unit)
         ncont=7
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
         
         norm_plist = (/b0,L0,n0,beta_e/beta0,1._r8/)
         eta_plist = (/eta,0.1_r8,etavac,r_eta/)
         kappa_plist = (/kappa_par,kappa_perp,0._r8,1.e8_r8/)
         IF(kappa_par /= kappa_perp)THEN
            kappa_case = "anisotropic"
         ELSE
            kappa_case="."
         ENDIF
         CALL postxmhd_init(indir,norm_plist,eta_case,eta_plist,
     $        kappa_case,kappa_plist)
      CASE("HallMHDicold")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=HallMHD_list)
         CLOSE(UNIT=in_unit)
      CASE("cmhd")
         CALL system("rm -f frc.bin frc.asc VecSqInt.asc")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=cmhd_list)
         CLOSE(UNIT=in_unit)
         ncont=4
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF

         norm_plist = (/b0,L0,n0,0.5_r8,1._r8/)
         eta_plist = (/eta,0.1_r8,etavac,r_eta/)
         kappa_plist = (/kappa_par,kappa_perp,0._r8,1.e8_r8/)
         CALL postxmhd_init(indir,norm_plist,eta_case,eta_plist,
     $        kappa_case,kappa_plist)
      CASE("cmhdn")
         CALL system("rm -f frc.bin frc.asc VecSqInt.asc")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=cmhdn_list)
         CLOSE(UNIT=in_unit)
         ncont=4
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
      CASE("rmhd")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=rmhd_list)
         CLOSE(UNIT=in_unit)
         ncont=4
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
      CASE("emhd")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=emhd_list)
         CLOSE(UNIT=in_unit)
      CASE("epmhd")
         OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
         READ(in_unit,NML=epmhd_list)
         CLOSE(UNIT=in_unit)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_read
c-----------------------------------------------------------------------
c     subprogram 2. post2fluid_main.
c     call the desired specialized post2fluid_mod subroutines for 
c     post-processing simulation results  
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_main(t,jac,xyw,xyw_kt,uw,uxyw,flag1d,flag2d,
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
      CASE("fivefield")
         IF(flag1d)THEN
            CALL post2fluid_dUdt(t,xyw,xyw_kt,uw,beginner)
            CALL post2fluid_UxyT(t,xyw,xyw_kt,uw)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
         ENDIF
         IF(flag2d)THEN
            CALL post2fluid_contour(t,xyw,uw,uxyw,beginner,
     $           header,footer,nt_next,ifile,stride,filename)
         ENDIF
      CASE("fourfieldplus")
         IF(flag1d)THEN
            CALL post2fluid_recrate(t,xyw,xyw_kt,uw,beginner)
            CALL post2fluid_UxyT(t,xyw,xyw_kt,uw)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL post2fluid_width(nxs,t,xyw,xyw_kt,uw,uxyw)
         ENDIF
         IF(flag2d)THEN
            CALL post2fluid_contour(t,xyw,uw,uxyw,beginner,
     $           header,footer,nt_next,ifile,stride,filename)
         ENDIF
      CASE("HallMHD")
         IF(flag1d)THEN
            CALL post2fluid_dUdt(t,xyw,xyw_kt,uw,beginner)
            CALL post2fluid_energy(t,xyw,jac,uw,uxyw)
            CALL post2fluid_UxyT(t,xyw,xyw_kt,uw)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL post2fluid_frc(t,xyw,uw,uxyw)
            CALL postxmhd_screw(t,xyw,uw,pmin)
         ENDIF
         IF(flag2d)THEN
            CALL post2fluid_contour(t,xyw,uw,uxyw,beginner,
     $           header,footer,nt_next,ifile,stride,filename)
         ENDIF
      CASE("cmhd")
         IF(flag1d)THEN
            CALL post2fluid_dUdt(t,xyw,xyw_kt,uw,beginner)
            CALL post2fluid_energy(t,xyw,jac,uw,uxyw)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL post2fluid_maxUvsT(t,uw,uxyw)
            CALL post2fluid_frc(t,xyw,uw,uxyw)
         ENDIF
         IF(flag2d)THEN
            CALL post2fluid_contour(t,xyw,uw,uxyw,beginner,
     $           header,footer,nt_next,ifile,stride,filename)
         ENDIF
      CASE("cmhdn")
         IF(flag1d)THEN
            CALL post2fluid_energy(t,xyw,jac,uw,uxyw)
            CALL post2fluid_maxUvsT(t,uw,uxyw)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL post2fluid_frc(t,xyw,uw,uxyw)
         ENDIF
         IF(flag2d)THEN
            CALL post2fluid_contour(t,xyw,uw,uxyw,beginner,
     $           header,footer,nt_next,ifile,stride,filename)
         ENDIF
      CASE("rmhd")
         IF(flag1d)THEN
            CALL post2fluid_UxyT(t,xyw,xyw_kt,uw)
            CALL post2fluid_dUdt(t,xyw,xyw_kt,uw,beginner)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL post2fluid_width(nxs,t,xyw,xyw_kt,uw,uxyw)
         ENDIF
         IF(flag2d)THEN
            CALL post2fluid_contour(t,xyw,uw,uxyw,beginner,
     $           header,footer,nt_next,ifile,stride,filename)
         ENDIF
      CASE("emhd")
         IF(flag1d)THEN
            CALL post2fluid_UxyT(t,xyw,xyw_kt,uw)
            CALL post2fluid_dUdt(t,xyw,xyw_kt,uw,beginner)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL post2fluid_width(nxs,t,xyw,xyw_kt,uw,uxyw)
         ENDIF 
      CASE("epmhd")
         IF(flag1d)THEN
            CALL post2fluid_UxyT(t,xyw,xyw_kt,uw)
            CALL post2fluid_dUdt(t,xyw,xyw_kt,uw,beginner)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL post2fluid_width(nxs,t,xyw,xyw_kt,uw,uxyw)
         ENDIF
      CASE("HallMHDicold")
         IF(flag1d)THEN
            CALL post2fluid_dUdt(t,xyw,xyw_kt,uw,beginner)
            CALL post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
            CALL post2fluid_width(nxs,t,xyw,xyw_kt,uw,uxyw)
         ENDIF         
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_main
c-----------------------------------------------------------------------
c     subprogram 3. post2fluid_dUdt.
c     generate du/dt(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_dUdt(t,xyw,xyw_kt,uw,beginner)

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
      IF(dUdt_name=="Flux")THEN
         SELECT CASE(job_type)
         CASE("fivefield","HallMHD","cmhd")
            utemp(1,:,:)=uw(2,:,:)
            CALL plotter_UxyT(t,nol,nol,xyw,xyw_kt,utemp,.FALSE.,value1,
     $           err)
            value1(1)=value1(1)-epsilon
            IF(value1(1) /= 0)THEN
               u2(1)=LOG(ABS(value1(1)))
            ELSE
               u2(1)=0
            ENDIF
            u2(2)=value1(1)
            IF(err)u2=u2_old
            CALL plotter_dUdt(t,t_old,u2,u2_old,beginner)
         CASE("rmhd","HallMHDicold","emhd","fourfieldplus","epmhd")
            utemp(1,:,:)=uw(1,:,:)
            CALL plotter_UxyT(t,nol,nol,xyw,xyw_kt,utemp,.FALSE.,value1,
     $           err)
            value1(1)=value1(1)-epsilon
            IF(value1(1) /= 0)THEN
               u2(1)=LOG(ABS(value1(1)))
            ELSE
               u2(1)=0
            ENDIF
            u2(2)=value1(1)
            IF(err)u2=u2_old
            CALL plotter_dUdt(t,t_old,u2,u2_old,beginner)
         END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_dUdt
c-----------------------------------------------------------------------
c     subprogram 4. post2fluid_contour.
c     generates contours of desired physical quantities
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_contour(t,xyw,uw,uxyw,beginner,header,
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
      CASE("fivefield")
         nyw=SIZE(uw,3)-1
c-----------------------------------------------------------------------
c     calculate Bz
c-----------------------------------------------------------------------
         utemp(1,:,:)=uw(3,:,:)
c-----------------------------------------------------------------------
c     calculate Bz even in y
c-----------------------------------------------------------------------
         utemp(2,:,:)=(utemp(1,:,:)+utemp(1,:,nyw:0:-1))/2
c-----------------------------------------------------------------------
c     calculate Bz odd in y
c-----------------------------------------------------------------------
         utemp(3,:,:)=(utemp(1,:,:)-utemp(1,:,nyw:0:-1))/2
c-----------------------------------------------------------------------
c     calculate Vx
c-----------------------------------------------------------------------
         utemp(4,:,:)=-uxyw(:,:,2,7)
c-----------------------------------------------------------------------
c     calculate Vy
c-----------------------------------------------------------------------
         utemp(5,:,:)=uxyw(:,:,1,7)
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
         DO ix=1,ncont
            IF (MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
         ENDDO
         CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $        nt_next,xcont,ifile,stride,filename,"Ucontour")
      CASE("fourfieldplus")
c-----------------------------------------------------------------------
c     calculate Bx
c-----------------------------------------------------------------------
         utemp(1,:,:)=-uxyw(:,:,2,1)
c-----------------------------------------------------------------------
c     calculate By
c-----------------------------------------------------------------------
         utemp(2,:,:)=uxyw(:,:,1,1)
c-----------------------------------------------------------------------
c     calculate electron Vx
c-----------------------------------------------------------------------
         utemp(3,:,:)=-(di*uxyw(:,:,2,2) + uxyw(:,:,2,6))
c-----------------------------------------------------------------------
c     calculate electron Vy
c-----------------------------------------------------------------------
         utemp(4,:,:)=di*uxyw(:,:,1,2) + uxyw(:,:,1,6)
c-----------------------------------------------------------------------
c     calculate ion Vx
c-----------------------------------------------------------------------
         utemp(5,:,:)=-uxyw(:,:,2,6)
c-----------------------------------------------------------------------
c     calculate ion Vy
c-----------------------------------------------------------------------
         utemp(6,:,:)=uxyw(:,:,1,6)
c-----------------------------------------------------------------------
c     calculate -Vi x B
c-----------------------------------------------------------------------
         utemp(7,:,:)=uxyw(:,:,2,6)*uxyw(:,:,1,1)
     $        - uxyw(:,:,1,6)*uxyw(:,:,2,1)
c-----------------------------------------------------------------------
c     calculate -Ve x B
c-----------------------------------------------------------------------
         utemp(8,:,:)=-utemp(3,:,:)*uxyw(:,:,1,1)
     $        - utemp(4,:,:)*uxyw(:,:,2,1)
c-----------------------------------------------------------------------
c     calculate eta*Jz
c-----------------------------------------------------------------------
         utemp(9,:,:)=eta/di*(uw(4,:,:)-uw(5,:,:))
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
         DO ix=1,ncont
            IF (MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
         ENDDO
         CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $        nt_next,xcont,ifile,stride,filename,"Ucontour")
      CASE("HallMHD")
         SELECT CASE("init_type")
         CASE("frc")
            r_faci = 1._r8/xyw(2,:,:)
            WHERE(xyw(2,:,:) < 1.e-12)
               r_faci = 0.
            END WHERE
         CASE DEFAULT
            r_faci = 0.
         END SELECT
c-----------------------------------------------------------------------
c     Magnetic Pressure 
c-----------------------------------------------------------------------
         utemp(1,:,:)=.5*(uxyw(:,:,1,2)**2
     $        + (r_faci*uw(2,:,:)+uxyw(:,:,2,2))**2 + uw(3,:,:)**2)
c-----------------------------------------------------------------------
c     Temperature
c-----------------------------------------------------------------------
         utemp(2,:,:)=uw(8,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     Total Pressure
c-----------------------------------------------------------------------
         utemp(3,:,:)=utemp(1,:,:)+uw(8,:,:)
c-----------------------------------------------------------------------
c     -V_e x B/ rho
c-----------------------------------------------------------------------
         utemp(4,:,:)=-(uw(4,:,:)*uxyw(:,:,1,2)
     $        + uw(5,:,:)*(r_faci*uw(2,:,:) + uxyw(:,:,2,2))
     $        + di*(uxyw(:,:,1,3)*(r_faci*uw(2,:,:) + uxyw(:,:,2,2))
     $        - (uxyw(:,:,2,3) + r_faci*uw(3,:,:))*uxyw(:,:,1,2)))
     $        /uw(1,:,:)
c-----------------------------------------------------------------------
c     Alfven speed
c-----------------------------------------------------------------------
         utemp(5,:,:)=SQRT(uw(4,:,:)**2+uw(5,:,:)**2)/
     $        SQRT(uw(1,:,:)*2.*utemp(1,:,:))
c-----------------------------------------------------------------------
c     Mach number
c-----------------------------------------------------------------------
         utemp(6,:,:)=SQRT(uw(4,:,:)**2+uw(5,:,:)**2)
     $        /(SQRT(uw(8,:,:))*uw(1,:,:))
c-----------------------------------------------------------------------
c     Ex
c-----------------------------------------------------------------------
         utemp(7,:,:)=-uw(5,:,:)*uw(3,:,:)/uw(1,:,:)
     $        + uw(6,:,:)*uxyw(:,:,1,2)/uw(1,:,:)
     $        + eta*(uxyw(:,:,2,3) + r_faci*uw(3,:,:))
     $        - di/uw(1,:,:)*(uw(3,:,:)*uxyw(:,:,1,3) 
     $        + uxyw(:,:,1,2)*uw(7,:,:) + uxyw(:,:,1,8))
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
         DO ix=1,ncont
            IF (MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
         ENDDO
         CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $        nt_next,xcont,ifile,stride,filename,"Ucontour")
      CASE("rmhd")
c-----------------------------------------------------------------------
c     Ez-field 
c-----------------------------------------------------------------------
         IF(beginner)ALLOCATE(uc_old(0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         IF(header)THEN
            utemp(1,:,:)=0
         ELSE
            utemp(1,:,:)=(uw(1,:,:)-uc_old)/(t-tc_old)
         ENDIF
         uc_old=uw(1,:,:)
         tc_old=t
         IF(footer)DEALLOCATE(uc_old)
c-----------------------------------------------------------------------
c     -VxB
c-----------------------------------------------------------------------
         utemp(2,:,:)=uxyw(:,:,2,4)*uxyw(:,:,1,1) 
     $        - uxyw(:,:,1,4)*uxyw(:,:,2,1)
c-----------------------------------------------------------------------
c     Vx
c-----------------------------------------------------------------------
         utemp(3,:,:)=-uxyw(:,:,2,4)
c-----------------------------------------------------------------------
c     Vy
c-----------------------------------------------------------------------
         utemp(4,:,:)=uxyw(:,:,1,4)
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
         DO ix=1,ncont
            IF (MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
         ENDDO
         CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $        nt_next,xcont,ifile,stride,filename,"Ucontour")
c-----------------------------------------------------------------------
c     CMHD
c-----------------------------------------------------------------------
      CASE("cmhd")
         cyl_fac=0
         IF(cylinder)cyl_fac=1
c-----------------------------------------------------------------------
c     calculate Bz and Br
c-----------------------------------------------------------------------
         utemp(1,:,:) = -(cyl_fac*uw(2,:,:)/xyw(2,:,:) + uxyw(:,:,2,2))
         utemp(2,:,:) = uxyw(:,:,1,2)
         WHERE(xyw(2,:,:) < 1.e-12) 
            utemp(1,:,:) = -uxyw(:,:,2,2)
         END WHERE
c-----------------------------------------------------------------------
c     calculate T
c-----------------------------------------------------------------------
         utemp(3,:,:)=uw(3,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate psi (=-rAphi)
c-----------------------------------------------------------------------
         utemp(4,:,:) = xyw(2,:,:)**cyl_fac*uw(2,:,:)
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
         CALL postxmhd_ploteta(t,xyw,uw(1,:,:),uw(3,:,:),uw(6,:,:),
     $        header,footer,nt_next,ifile,stride,filename)
c-----------------------------------------------------------------------
c     calculate and plot spatially varying thermal conductivity
c-----------------------------------------------------------------------
         CALL postxmhd_plotkappa(t,xyw,uw(1,:,:),uw(3,:,:),
     $        utemp(1,:,:)**2+utemp(2,:,:)**2,header,footer,nt_next,
     $        ifile,stride,filename)
c-----------------------------------------------------------------------
c     CMHDN
c-----------------------------------------------------------------------
      CASE("cmhdn")
         cyl_fac=0
         IF(cylinder)cyl_fac=1
c-----------------------------------------------------------------------
c     calculate Bz and Br
c-----------------------------------------------------------------------
         utemp(1,:,:) = -(cyl_fac*uw(2,:,:)/xyw(2,:,:) + uxyw(:,:,2,2))
         utemp(2,:,:) = uxyw(:,:,1,2)
         WHERE(xyw(2,:,:) < 1.e-12) 
            utemp(1,:,:) = -uxyw(:,:,2,2)
         END WHERE
c-----------------------------------------------------------------------
c     calculate T
c-----------------------------------------------------------------------
         utemp(3,:,:)=uw(3,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate psi (=-rAphi)
c-----------------------------------------------------------------------
         utemp(4,:,:) = xyw(2,:,:)**cyl_fac*uw(2,:,:)
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
         DO ix=1,ncont
            IF (MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
         ENDDO
         CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $        nt_next,xcont,ifile,stride,filename,"Ucontour")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_contour
c-----------------------------------------------------------------------
c     subprogram 5. post2fluid_energy.
c     generate a time plot for an integral of total energy
c     over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_energy(t,xyw,jac,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: jac
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: uxyw

      INTEGER :: nx,ny
      REAL(r8) :: cyl_fac
      REAL(r8), DIMENSION(4) :: value
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
      CASE("cmhd","cmhdn")
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

         ALLOCATE(uw1(2,nx,ny))
         ALLOCATE(uw2(2,nx,ny))
         ALLOCATE(uw3(2,nx,ny))
c-----------------------------------------------------------------------
c     kinetic energy components.
c-----------------------------------------------------------------------
         uw1(1,:,:)=r_facs*uw(4,:,:)/SQRT(uw(1,:,:))
         uw2(1,:,:)=r_facs*uw(5,:,:)/SQRT(uw(1,:,:))
         WHERE(uw(1,:,:) <= 0)
            uw1(1,:,:)=0
            uw2(1,:,:)=0
         ENDWHERE
c-----------------------------------------------------------------------
c     magnetic energy components.
c-----------------------------------------------------------------------
         uw1(2,:,:)=-(cyl_fac*r_faci*uw(2,:,:) + uxyw(:,:,2,2))*r_facs
         uw2(2,:,:)=uxyw(:,:,1,2)*r_facs
c-----------------------------------------------------------------------
c     for FRC, eliminate contributions outside of the separatrix.
c-----------------------------------------------------------------------
         SELECT CASE(init_type)
         CASE("frc","frc_rm")
            WHERE(uw(2,:,:)<0)
               uw1(1,:,:)=0
               uw1(2,:,:)=0
               uw2(1,:,:)=0
               uw2(2,:,:)=0
            ENDWHERE
         END SELECT
c-----------------------------------------------------------------------
c     calculate the integral of vector energies over the domain.
c-----------------------------------------------------------------------
         uw3=0
         CALL plotter_VecSqInt(t,xyw,uw1,uw2,uw3,.FALSE.,value(1:2))
c-----------------------------------------------------------------------
c        thermal energy.
c-----------------------------------------------------------------------
         uw1(1,:,:)=r_fac*3._r8*uw(3,:,:)
c-----------------------------------------------------------------------
c     for FRC, eliminate contributions outside of the separatrix.
c-----------------------------------------------------------------------
         SELECT CASE(init_type)
         CASE("frc","frc_rm")
            WHERE(uw(2,:,:)<0)
               uw1(1,:,:)=0
            END WHERE
         END SELECT

         CALL plotter_integral(t,xyw,uw1(1:1,:,:),.FALSE.,value(3:3))
c-----------------------------------------------------------------------
c        total energy.
c-----------------------------------------------------------------------
         value(4)=SUM(value)
         value=0.5*value
         IF(cylinder)value=value*twopi

         DEALLOCATE(uw1,uw2,uw3)
      CASE("HallMHD")
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

         ALLOCATE(uw1(2,nx,ny))
         ALLOCATE(uw2(2,nx,ny))
         ALLOCATE(uw3(2,nx,ny))
c-----------------------------------------------------------------------
c     kinetic energy components.
c-----------------------------------------------------------------------
         uw1(1,:,:)=r_facs*uw(4,:,:)/SQRT(uw(1,:,:))
         uw2(1,:,:)=r_facs*uw(5,:,:)/SQRT(uw(1,:,:))
         uw3(1,:,:)=r_facs*uw(6,:,:)/SQRT(uw(1,:,:))
         WHERE(uw(1,:,:) <= 0)
            uw1(1,:,:)=0
            uw2(1,:,:)=0
            uw3(1,:,:)=0
         ENDWHERE
c-----------------------------------------------------------------------
c     magnetic energy components.
c-----------------------------------------------------------------------
         uw1(2,:,:)=-(cyl_fac*r_faci*uw(2,:,:) + uxyw(:,:,2,2))*r_facs
         uw2(2,:,:)=uxyw(:,:,1,2)*r_facs
         uw3(2,:,:)=uw(3,:,:)
c-----------------------------------------------------------------------
c     for FRC, eliminate contributions outside of the separatrix.
c-----------------------------------------------------------------------
         IF(init_type=="frc")THEN
            WHERE(uw(2,:,:)<0)
               uw1(1,:,:)=0
               uw2(1,:,:)=0
               uw3(1,:,:)=0
               uw1(2,:,:)=0
               uw2(2,:,:)=0
               uw3(2,:,:)=0
            ENDWHERE
         ENDIF
c-----------------------------------------------------------------------
c     calculate the integral of vector energies over the domain.
c-----------------------------------------------------------------------
         CALL plotter_VecSqInt(t,xyw,uw1,uw2,uw3,.FALSE.,value(1:2))
c-----------------------------------------------------------------------
c     thermal energy.
c-----------------------------------------------------------------------
         uw1(1,:,:)=r_fac*3._r8*uw(8,:,:)
c-----------------------------------------------------------------------
c     for FRC, eliminate contributions outside of the separatrix.
c-----------------------------------------------------------------------
         IF(init_type=="frc")THEN
            WHERE(uw(2,:,:)<0)
               uw1(1,:,:)=0
            ENDWHERE
         ENDIF

         CALL plotter_integral(t,xyw,uw1(1:1,:,:),.FALSE.,value(3:3))
c-----------------------------------------------------------------------
c        total energy.
c-----------------------------------------------------------------------
         value(4)=SUM(value)
         value=0.5*value
         IF(cylinder)value=value*twopi
         IF(init_type(1:3)=="GEM")value=value*4._r8

         DEALLOCATE(uw1,uw2,uw3)
      END SELECT
c-----------------------------------------------------------------------
c     open, write, and close VecSqInt.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      WHERE(value == 0)
         value=TINY(value)
      END WHERE
      OPEN(UNIT=VecSqInt_unit,FILE="VecSqInt.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(VecSqInt_unit)REAL(t,4),REAL(value,4)
      CLOSE(UNIT=VecSqInt_unit)
c-----------------------------------------------------------------------
c     open, write, and close VecSqInt.asc file to store ascii data for 
c     xdraw.
c-----------------------------------------------------------------------
      IF(ifascii)THEN
         OPEN(UNIT=VecSqInt_unit,FILE="VecSqInt.asc",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="FORMATTED")
         WRITE(VecSqInt_unit,*)REAL(t,4),REAL(value,4)
         CLOSE(UNIT=VecSqInt_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_energy
c-----------------------------------------------------------------------
c     subprogram 6. post2fluid_maxUvsT.
c     generate u(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_maxUvsT(t,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
      REAL(r8), DIMENSION(:), ALLOCATABLE :: value
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("cmhd")
         ALLOCATE(utemp(2,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),value(2))
         utemp(1,:,:)=uw(3,:,:)           ! pressure
         utemp(2,:,:)=uw(4,:,:)/uw(1,:,:) ! z velocity
         CALL plotter_maxUvsT("maxUvsT",t,utemp,value)
         DEALLOCATE(utemp,value)
      CASE("cmhdn")
         ALLOCATE(utemp(4,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),value(4))
         utemp(1,:,:)=uw(1,:,:)   ! density
         utemp(2,:,:)=uw(7,:,:)   ! neutral density
         utemp(3,:,:)=uw(1,:,:)/(uw(1,:,:)+uw(7,:,:)) ! ioniz. fract.
         utemp(4,:,:)=uw(3,:,:)/uw(1,:,:)  ! temperature
         CALL plotter_maxUvsT("maxUvsT",t,utemp,value)
         DEALLOCATE(utemp,value)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_maxUvsT
c-----------------------------------------------------------------------
c     subprogram 7. post2fluid_UxyT.
c     generate U vs. time plot at a point.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_UxyT(t,xyw,xyw_kt,uw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt


      LOGICAL :: err=.FALSE.
      REAL(r8), DIMENSION(1) :: value1,value2
      REAL(r8), DIMENSION(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: utemp
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(UxyT_name)
      CASE("Flux")
         SELECT CASE(job_type)
         CASE("rmhd","emhd","epmhd")
            utemp(1,:,:)=uw(1,:,:)
            CALL plotter_UxyT(t,nol,nol,xyw,xyw_kt,utemp,.FALSE.,value1,
     $           err)
            CALL plotter_UxyT(t,x1,y1,xyw,xyw_kt,utemp,.FALSE.,value2,
     $           err)
            value1=value1-value2
            OPEN(UNIT=UxyT_unit,FILE="UxyT.bin",STATUS="UNKNOWN",
     $           POSITION="APPEND",FORM="UNFORMATTED")
            WRITE(UxyT_unit)REAL(t,4),REAL(value1,4)
            CLOSE(UNIT=UxyT_unit)
         CASE("fivefield","HallMHD")
            utemp(1,:,:)=uw(2,:,:)
            CALL plotter_UxyT(t,x1,y1,xyw,xyw_kt,utemp,.FALSE.,value1,
     $           err)
            CALL plotter_UxyT(t,x2,y2,xyw,xyw_kt,utemp,.FALSE.,value2,
     $           err)
            value1=value1-value2
            OPEN(UNIT=UxyT_unit,FILE="UxyT.bin",STATUS="UNKNOWN",
     $           POSITION="APPEND",FORM="UNFORMATTED")
            WRITE(UxyT_unit)REAL(t,4),REAL(value1,4)
            CLOSE(UNIT=UxyT_unit)
         CASE("fourfieldplus")
            utemp(1,:,:)=uw(1,:,:)
            CALL plotter_UxyT(t,x1,y1,xyw,xyw_kt,utemp,.FALSE.,value1,
     $           err)
            CALL plotter_UxyT(t,x2,y2,xyw,xyw_kt,utemp,.FALSE.,value2,
     $           err)
            value1=value1-value2
            OPEN(UNIT=UxyT_unit,FILE="UxyT.bin",STATUS="UNKNOWN",
     $           POSITION="APPEND",FORM="UNFORMATTED")
            WRITE(UxyT_unit)REAL(t,4),REAL(value1,4)
            CLOSE(UNIT=UxyT_unit)
         END SELECT
      CASE("Current")
         SELECT CASE(job_type)
         CASE("rmhd","emhd","epmhd")
            utemp(1,:,:)=uw(3,:,:)
            CALL plotter_UxyT(t,nol,nol,xyw,xyw_kt,utemp,.TRUE.,value1,
     $           err)
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_UxyT
c-----------------------------------------------------------------------
c     subprogram 8. post2fluid_Uprofile.
c     generate slices of U for each time-step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)

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
      CASE("fivefield")
         ALLOCATE(utemp(3,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(4,0:nxs))
         utemp(1,:,:)=uw(5,:,:)
         utemp(2,:,:)=uw(6,:,:)
         utemp(3,:,:)=uw(5,:,:)-uw(6,:,:)
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
      CASE("fourfieldplus")
         ALLOCATE(utemp(3,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(4,0:nxs))
         utemp=0
         IF(y1 == y2)THEN
            utemp(1,:,:)=-uxyw(:,:,2,6)
            utemp(2,:,:)=-di*uxyw(:,:,2,2) - uxyw(:,:,2,6)
         ELSEIF(x1 == x2)THEN
            utemp(1,:,:)=uxyw(:,:,1,6)
            utemp(2,:,:)=di*uxyw(:,:,1,2) + uxyw(:,:,1,6)
         ENDIF
         utemp(3,:,:)=(uw(4,:,:)-uw(5,:,:))/di
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
      CASE("rmhd","epmhd")
         ALLOCATE(utemp(4,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(5,0:nxs))
         utemp(1,:,:)=-uxyw(:,:,2,1)
         utemp(2,:,:)=uxyw(:,:,1,4)
         utemp(3,:,:)=uxyw(:,:,1,1)*uxyw(:,:,2,4)
     $        - uxyw(:,:,2,1)*uxyw(:,:,1,4)
         utemp(4,:,:)=uw(3,:,:)
         filename="Uprofile1"
         CALL plotter_Uprofile(nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
         utemp(1,:,:)=uxyw(:,:,1,1)
         utemp(2,:,:)=-uxyw(:,:,2,4)
         utemp(3,:,:)=uxyw(:,:,1,1)*uxyw(:,:,2,4)
     $        - uxyw(:,:,2,1)*uxyw(:,:,1,4)
         utemp(4,:,:)=uw(3,:,:)
         filename="Uprofile2"
         CALL plotter_Uprofile(nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
      CASE("emhd")
         ALLOCATE(utemp(7,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(8,0:nxs))
         utemp(1,:,:)=-uxyw(:,:,2,1)
         utemp(2,:,:)=uxyw(:,:,1,1)
         utemp(3,:,:)=-di*uxyw(:,:,2,2)
         utemp(4,:,:)=di*uxyw(:,:,1,2)
         utemp(5,:,:)=di*(uxyw(:,:,1,1)*uxyw(:,:,2,2)
     $        - uxyw(:,:,2,1)*uxyw(:,:,1,2))
         utemp(6,:,:)=uw(3,:,:)
         utemp(7,:,:)=uw(2,:,:)
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
      CASE("HallMHD")
         ALLOCATE(utemp(7,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(8,0:nxs))
         utemp(1,:,:)=uw(1,:,:)
         utemp(2,:,:)=SQRT(uw(4,:,:)**2+uw(5,:,:)**2)
     $        /(SQRT(uw(8,:,:))*uw(1,:,:))
         utemp(3,:,:)=uw(4,:,:)/uw(1,:,:)
         utemp(4,:,:)=eta*(uw(6,:,:)-uw(7,:,:))
         utemp(5,:,:)=di*(uxyw(:,:,2,3)*uxyw(:,:,1,2)
     $        - uxyw(:,:,1,3)*uxyw(:,:,2,2))/uw(1,:,:)
         utemp(6,:,:)=-(uw(5,:,:)*uxyw(:,:,2,2)+uw(4,:,:)*uxyw(:,:,1,2))
         utemp(7,:,:)=utemp(4,:,:)+utemp(5,:,:)+utemp(6,:,:)
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
      CASE("HallMHDicold")
         ALLOCATE(utemp(9,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(10,0:nxs))
         utemp(1,:,:)=(uw(5,:,:)-uw(6,:,:))/di
         utemp(2,:,:)=uw(7,:,:)
         IF(eta_spitzer)THEN
            utemp(3,:,:)=(uw(5,:,:)-uw(6,:,:))/di
     $           *eta*(beta0/uw(7,:,:))**1.5
         ELSE
            utemp(3,:,:)=(uw(5,:,:)-uw(6,:,:))/di*eta
         ENDIF
         IF(x1==x2)THEN
            utemp(4,:,:)=uw(4,:,:)
            utemp(5,:,:)=uw(4,:,:)+di*uxyw(:,:,1,2)
            utemp(6,:,:)=-di*uxyw(:,:,1,2)*uxyw(:,:,2,1)
            utemp(7,:,:)=-uw(4,:,:)*uxyw(:,:,2,1)
            utemp(8,:,:)=-uxyw(:,:,2,1)
         ELSE
            utemp(4,:,:)=uw(3,:,:)
            utemp(5,:,:)=uw(3,:,:)-di*uxyw(:,:,2,2)
            utemp(6,:,:)=di*uxyw(:,:,2,2)*uxyw(:,:,1,1)
            utemp(7,:,:)=-uw(3,:,:)*uxyw(:,:,1,1)
            utemp(8,:,:)=uxyw(:,:,1,1)
         ENDIF
         utemp(9,:,:)=uw(1,:,:)
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
      CASE("cmhd","cmhdn")
         IF(cylinder)cyl_fac=1
         ALLOCATE(utemp(9,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(10,0:nxs))
c-----------------------------------------------------------------------
c     utemp(3): Bz
c     utemp(9): Br
c-----------------------------------------------------------------------
         utemp(3,:,:) = -(cyl_fac*uw(2,:,:)/xyw(2,:,:) + uxyw(:,:,2,2))
         WHERE(xyw(2,:,:) < 1.e-12)
            utemp(3,:,:) = -uxyw(:,:,2,2)
         END WHERE
         utemp(9,:,:) = uxyw(:,:,1,2)
         utemp(1,:,:) = uw(6,:,:)*utemp(3,:,:) - uxyw(:,:,2,3)
         utemp(2,:,:) = uw(3,:,:)/uw(1,:,:)
         utemp(4,:,:) = uw(6,:,:)
         utemp(5,:,:) = uw(1,:,:)
         utemp(6,:,:) = uw(5,:,:)
         utemp(7,:,:) = uw(5,:,:)/uw(1,:,:)
         utemp(8,:,:) = uw(3,:,:)
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
      END SELECT
      DEALLOCATE(utemp,uprofile)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_Uprofile
c-----------------------------------------------------------------------
c     subprogram 9. post2fluid_width.
c     generate peak width vs. time plot along a profile.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_width(nxs,t,xyw,xyw_kt,uw,uxyw)

      INTEGER, INTENT(IN) :: nxs
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt,uxyw

      LOGICAL :: err=.FALSE.
      CHARACTER(10) :: wtype
      REAL(r8) :: width1,width2,width3,width4,width5,width6
      REAL(r8), DIMENSION(1) :: value1,value2,value3,value4,value5
      REAL(r8), DIMENSION(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: utemp
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("HallMHDicold")
         wtype="half_max"
         utemp(1,:,:)=(uw(5,:,:)-uw(6,:,:))/di
         CALL plotter_width
     $        (nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width1)
         CALL plotter_width
     $        (nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,wtype,width2)
         OPEN(UNIT=UxyT_unit,FILE="layer.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(UxyT_unit)REAL(t,4),REAL(width1,4),REAL(width2,4)
         CLOSE(UNIT=UxyT_unit)
      CASE("rmhd","epmhd")
         wtype="half_max"
         utemp(1,:,:)=uw(3,:,:)
         CALL plotter_width
     $        (nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width1)
         CALL plotter_width
     $        (nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,wtype,width2)
         wtype="local_max"
         utemp(1,:,:)=uxyw(:,:,1,4)
         CALL plotter_width
     $        (nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width3)
         utemp(1,:,:)=-uxyw(:,:,2,4)
         CALL plotter_width
     $        (nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,wtype,width4)
         utemp(1,:,:)=ABS(uxyw(:,:,2,1))
         CALL plotter_UxyT
     $        (t,x1,width1,xyw,xyw_kt,utemp,.FALSE.,value1,err)
         utemp(1,:,:)=ABS(uxyw(:,:,1,1))
         CALL plotter_UxyT
     $        (t,width2,y2,xyw,xyw_kt,utemp,.FALSE.,value2,err)
         utemp(1,:,:)=ABS(uxyw(:,:,1,4))
         CALL plotter_UxyT
     $        (t,x1,width1,xyw,xyw_kt,utemp,.FALSE.,value3,err)
         utemp(1,:,:)=ABS(uxyw(:,:,2,4))
         CALL plotter_UxyT
     $        (t,width2,y2,xyw,xyw_kt,utemp,.FALSE.,value4,err)
         utemp(1,:,:)=ABS(uw(3,:,:))
         CALL plotter_UxyT
     $        (t,nol,nol,xyw,xyw_kt,utemp,.FALSE.,value5,err)
         OPEN(UNIT=width_unit,FILE="layer.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(width_unit)REAL(t,4),REAL(width1,4),REAL(width2,4),
     $        REAL(width3,4),REAL(width4,4),REAL(value1,4),
     $        REAL(value2,4),REAL(value3,4),REAL(value4,4),
     $        REAL(value5,4)
         CLOSE(UNIT=width_unit)
      CASE("emhd")
         SELECT CASE(init_type)
         CASE("de_instability")
            wtype="local_max"
            utemp(1,:,:)=uw(2,:,:)
            CALL plotter_width
     $           (nol,y1,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width1)
            wtype="half_max"
            CALL plotter_width
     $           (width1,y1,nol,y1,nxs,xyw,xyw_kt,utemp,wtype,width2)
            CALL plotter_width
     $           (width1,y1,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width3)
            value1(1)=width3+width2
            OPEN(UNIT=width_unit,FILE="ef.bin",STATUS="UNKNOWN",
     $           POSITION="APPEND",FORM="UNFORMATTED")
            WRITE(width_unit)REAL(t,4),REAL(value1,4),REAL(width1,4)
            CLOSE(UNIT=width_unit)
         CASE DEFAULT
            wtype="half_max"
            utemp(1,:,:)=uw(3,:,:)
            CALL plotter_width
     $           (nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width1)
            CALL plotter_width
     $           (nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,wtype,width2)
            wtype="local_max"
            utemp(1,:,:)=uxyw(:,:,1,2)
            CALL plotter_width
     $           (nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width3)
            utemp(1,:,:)=-uxyw(:,:,2,2)
            CALL plotter_width
     $           (nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,wtype,width4)
            utemp(1,:,:)=ABS(uxyw(:,:,2,1))
            CALL plotter_UxyT(t,x1,width3,xyw,xyw_kt,utemp,.FALSE.,
     $           value1,err)
            utemp(1,:,:)=ABS(uxyw(:,:,1,1))
            CALL plotter_UxyT(t,width4,y2,xyw,xyw_kt,utemp,.FALSE.,
     $           value2,err)
            utemp(1,:,:)=ABS(uxyw(:,:,1,2))
            CALL plotter_UxyT(t,x1,width3,xyw,xyw_kt,utemp,.FALSE.,
     $           value3,err)
            utemp(1,:,:)=ABS(uxyw(:,:,2,2))
            CALL plotter_UxyT(t,width4,y2,xyw,xyw_kt,utemp,.FALSE.,
     $           value4,err)
            utemp(1,:,:)=ABS(uw(3,:,:))
            CALL plotter_UxyT(t,nol,nol,xyw,xyw_kt,utemp,.FALSE.,
     $           value5,err)
            OPEN(UNIT=width_unit,FILE="layer.bin",STATUS="UNKNOWN",
     $           POSITION="APPEND",FORM="UNFORMATTED")
            WRITE(width_unit)REAL(t,4),REAL(width1,4),REAL(width2,4),
     $           REAL(width3,4),REAL(width4,4),REAL(value1,4),
     $           REAL(value2,4),REAL(value3,4),REAL(value4,4),
     $           REAL(value5,4)
            CLOSE(UNIT=width_unit)
         END SELECT
      CASE("fourfieldplus")
         wtype="half_max"
         utemp(1,:,:)=(uw(4,:,:)-uw(5,:,:))/di
         CALL plotter_width
     $        (nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width1)
         CALL plotter_width
     $        (nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,wtype,width2)
         wtype="local_max"
         utemp(1,:,:)=di*uxyw(:,:,1,2) + uxyw(:,:,1,6)
         CALL plotter_width
     $        (nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width3)
         utemp(1,:,:)=-di*uxyw(:,:,2,2) - uxyw(:,:,2,6)
         CALL plotter_width
     $        (nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,wtype,width4)
         utemp(1,:,:)=uxyw(:,:,1,6)
         CALL plotter_width
     $        (nol,nol,x1,y1,nxs,xyw,xyw_kt,utemp,wtype,width5)
         utemp(1,:,:)=-uxyw(:,:,2,6)
         CALL plotter_width
     $        (nol,nol,x2,y2,nxs,xyw,xyw_kt,utemp,wtype,width6)
         OPEN(UNIT=width_unit,FILE="layer.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(width_unit)REAL(t,4),REAL(width1,4),REAL(width2,4),
     $        REAL(width3,4),REAL(width4,4),REAL(width5,4),
     $        REAL(width6,4)
         CLOSE(UNIT=width_unit)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_width
c-----------------------------------------------------------------------
c     subprogram 10. post2fluid_recrate.
c     generates reconnection rate plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_recrate(t,xyw,xyw_kt,uw,first)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt
      LOGICAL, INTENT(IN) :: first

      CHARACTER(16) :: filename
      INTEGER, PARAMETER :: nxs=500
      REAL(r8), PARAMETER :: nol=1.e-6
      REAL(r8) :: value
      REAL(r8), DIMENSION(2) :: u2
      REAL(r8), DIMENSION(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: utemp
      REAL(r8), DIMENSION(2,0:nxs) :: uprofile
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      filename=" "
      SELECT CASE(job_type)
      CASE("fourfieldplus")
         utemp(1,:,:)=uw(1,:,:)
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
      END SUBROUTINE post2fluid_recrate
c-----------------------------------------------------------------------
c     subprogram 11. post2fluid_frc.
c     plot various frc quantities vs. time.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post2fluid_frc(t,xyw,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw

      INTEGER :: ix,iy,nx,ny
      REAL(r8), DIMENSION(5) :: a
      REAL(r8), DIMENSION(8) :: value
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: p,j
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uw1
      REAL(r8) :: aux,vol
c-----------------------------------------------------------------------
c     Quit if not post-processing "frc" simulation
c-----------------------------------------------------------------------
      IF(init_type /= "frc")RETURN
c-----------------------------------------------------------------------
c     prepare auxiliary parameters
c-----------------------------------------------------------------------      
      value=zero
      aux=zero
      vol=zero
      nx=SIZE(xyw,2)-1
      ny=SIZE(xyw,3)-1
      IF(init_type == "screw") RETURN
      SELECT CASE(job_type)
      CASE("cmhd","cmhdn","HallMHD","omhd")
         ALLOCATE(uw1(4,0:nx,0:ny),p(0:nx,0:ny),j(0:nx,0:ny))
         IF(job_type=="cmhd".OR.job_type=="cmhdn"
     $        .OR.job_type=="omhd")THEN
            p=uw(3,:,:)
            j=uw(6,:,:)
         ELSEIF(job_type=="HallMHD")THEN
            p=uw(8,:,:)
            j=(uw(6,:,:)-uw(7,:,:))/di
         ENDIF
c-----------------------------------------------------------------------
c        find trapped flux (the maximum value of psi) -> value(4).
c-----------------------------------------------------------------------
         value(4)=MAXVAL(two*pi*uw(2,:,:)*xyw(2,:,:))
c-----------------------------------------------------------------------
c        calculate the integral of density inside the frc -> value(1).
c        calculate the axial location of the center of mass -> value(2).
c        calculate the average axial velocity -> value(3).
c        calculate the average beta -> value(7).
c          (note that beta=2p/B0^2 .. but B0=1)
c        calculate the average eta -> value(8).
c        locate the x- and y-coordinates of the O-point -> value(5,6)
c        (include radius factor for cylindrical coordinates.)
c-----------------------------------------------------------------------
         uw1(1,:,:)=uw(1,:,:)*xyw(2,:,:)*twopi
         uw1(2,:,:)=uw(4,:,:)*xyw(2,:,:)/uw(1,:,:)
         uw1(3,:,:)=2._r8*p*xyw(2,:,:)/1.763**2
         WHERE(xyw(2,:,:) < 1.e-12)
            uw1(3,:,:)=0
         ENDWHERE
c-----------------------------------------------------------------------
c        calculate spitzer/chodura resistive heating
c        (only applies for a specific frc problem)
c-----------------------------------------------------------------------
         uw1(4,:,:)=.00283*(uw(1,:,:)/p)**1.5
         uw1(4,:,:)=uw1(4,:,:) + 0.01*0.6111*uw(1,:,:)**(-.5)*
     $        (1.-EXP(-.2667/3.
     $        *ABS(j)/SQRT(1.6667*uw(1,:,:)*p)))
         WHERE(p/uw(1,:,:) <= 0)
            uw1(4,:,:)=1.
         ENDWHERE
         WHERE(uw1(4,:,:) > 1.)
            uw1(4,:,:)=1.
         ENDWHERE
         SELECT CASE(init_type)
         CASE("frc","frc_rm")
c-----------------------------------------------------------------------
c        eliminate contributions outside of the separatrix.
c-----------------------------------------------------------------------
            WHERE(uw(2,:,:)<0)
               uw1(1,:,:)=0
               uw1(2,:,:)=0
               uw1(3,:,:)=0
               uw1(4,:,:)=0
            ENDWHERE
         END SELECT
         DO ix=0,nx-1
            DO iy=0,ny-1
               a(1)=(xyw(2,ix,iy)+xyw(2,ix,iy+1))*
     $              (xyw(1,ix,iy+1)-xyw(1,ix,iy))/2._r8
               a(2)=(xyw(2,ix,iy+1)+xyw(2,ix+1,iy+1))*
     $              (xyw(1,ix+1,iy+1)-xyw(1,ix,iy+1))/2._r8
               a(3)=(xyw(2,ix+1,iy+1)+xyw(2,ix+1,iy))*
     $              (xyw(1,ix+1,iy)-xyw(1,ix+1,iy+1))/2._r8
               a(4)=(xyw(2,ix,iy)+xyw(2,ix+1,iy))*
     $              (xyw(1,ix,iy)-xyw(1,ix+1,iy))/2._r8
               a(5)=ABS(SUM(a(1:4)))
               
               value(1)=value(1)+a(5)*
     $              SUM(uw1(1,ix:ix+1,iy:iy+1))/4

               aux=(xyw(1,ix,iy)+xyw(1,ix+1,iy)
     $              +xyw(1,ix,iy+1)+xyw(1,ix+1,iy+1))/4
               value(2)=value(2)+a(5)*aux*
     $              SUM(uw1(1,ix:ix+1,iy:iy+1))/4

               value(3)=value(3)+a(5)*SUM(uw1(2,ix:ix+1,iy:iy+1))/4
               value(7)=value(7)+a(5)*SUM(uw1(3,ix:ix+1,iy:iy+1))/4
               value(8)=value(8)+a(5)*SUM(uw1(4,ix:ix+1,iy:iy+1))/4
               IF(uw(2,ix,iy)>=0)THEN
                  vol=vol+a(5)*SUM(xyw(2,ix:ix+1,iy:iy+1))/4
               ENDIF

               IF(two*pi*uw(2,ix,iy)==value(4))THEN
                  value(5)=xyw(1,ix,iy)
                  value(6)=xyw(2,ix,iy)
               ENDIF

            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c        axial location of the center of mass is
c        sum(two*pi*n*z*r*dA)/sum(two*pi*n*r*dA)
c-----------------------------------------------------------------------
         value(2)=value(2)/value(1)
c-----------------------------------------------------------------------
c        volume average velocity = sum(vel*r*dA)/sum(r*dA)
c        volume average beta = sum(beta*r*dA)/sum(r*dA)
c        volume average eta = sum(eta*r*dA)/sum(r*dA)
c-----------------------------------------------------------------------
         value(3)=value(3)/MAX(vol,1e-10_r8)
         value(7)=value(7)/MAX(vol,1e-10_r8)

         DEALLOCATE(uw1,p,j)
      END SELECT
c-----------------------------------------------------------------------
c     open, write, and close frc.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      OPEN(UNIT=oned_unit,FILE="frc.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(oned_unit)REAL(t,4),REAL(value,4)
      CLOSE(UNIT=oned_unit)
c-----------------------------------------------------------------------
c     open, write, and close frc.asc file to store ascii data for xdraw.
c-----------------------------------------------------------------------
      IF(ifascii)THEN
         OPEN(UNIT=oned_unit,FILE="frc.asc",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="FORMATTED")
 10      FORMAT(9e15.5)
         WRITE(oned_unit,10)REAL(t,4),REAL(value,4)
         CLOSE(UNIT=oned_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post2fluid_frc
      END MODULE post2fluid_mod
