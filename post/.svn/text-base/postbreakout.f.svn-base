c-----------------------------------------------------------------------
c     file postbreakout.f.
c     post-processes output from SEL code for the breakout CME problem.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. postbreakout_mod.
c     1. postbreakout_read.
c     2. postbreakout_dUdt.
c     4. postbreakout_contour.
c     5. postbreakout_energy.
c     6. postbreakout_maxUvsT.
c     8. postbreakout_Uprofile.
c     10. postbreakout_recrate.
c     11. postbreakout_cstracking.
c     12. postbreakout_rTtracking.
c-----------------------------------------------------------------------
c     subprogram 0. postbreakout_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE postbreakout_mod
      USE postxmhd_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: source,cylinder=.FALSE.
      CHARACTER(20), PRIVATE :: init_type,eta_case,kappa_case,mu_case
      REAL(r8), PARAMETER, PRIVATE :: nol=1.e-6
      REAL(r8), PRIVATE :: b0,eta,mu,mu_min,nu,hyper_eta,di,beta0,lx,ly,
     $     kappa_par,kappa_perp,Dn,x1,y1,x2,y2,n0,T0,r_eta,lambda,
     $     etavac,alpha,Te_frac,Rmax,v_peak,v_period,v_angle,rscale,
     $     phiscale,j_crit
      REAL(r8), DIMENSION(99), PRIVATE :: x0_list,y0_list

      LOGICAL, DIMENSION(99) :: in_check
      INTEGER, PRIVATE :: ncont,listsize
      REAL(r8), PRIVATE :: gravity,time0,tempr0,vel0
      REAL(r8), DIMENSION(2,99), PRIVATE :: u2_old,kt0_old,xy0_old
      REAL(r8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: uc_old

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. postbreakout_read.
c     read necessary post-processing parameters from post.in 
c     and input parameters from sel.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_read(indir)

      CHARACTER(*), INTENT(IN) :: indir
      CHARACTER(80) :: infile

      REAL(r8), DIMENSION(4) :: eta_plist,kappa_plist
      REAL(r8), DIMENSION(5) :: norm_plist

      NAMELIST/breakout_input/x1,y1,x2,y2,listsize,x0_list,y0_list

      NAMELIST/SolarHallMHD_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,init_type,source,cylinder,di,eta,eta_case,
     $     r_eta,etavac,mu_case,mu,mu_min,nu,kappa_case,kappa_par,
     $     kappa_perp,Dn,Te_frac,n0,T0,b0,Rmax,v_peak,v_period,v_angle,
     $     alpha,rscale,phiscale,lambda,beta0,lx,ly,hyper_eta,j_crit
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      x0_list=0.
      y0_list=0.
      OPEN(UNIT=in_unit,FILE="post.in",STATUS="OLD")
      READ(in_unit,NML=breakout_input)
      CLOSE(UNIT=in_unit)
      infile=TRIM(indir)//"/sel.in"
      OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
      READ(in_unit,NML=SolarHallMHD_list)
      CLOSE(UNIT=in_unit)
      ncont=14
      IF(flag1d)CALL system("rm -f Efield.bin breakout.bin track*.bin")
      IF(flag2d .AND. out_type/="hdf5")THEN
         OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $        STATUS="UNKNOWN",FORM="UNFORMATTED")
         WRITE(Ucontour_unit)1,0,ncont
         CLOSE(UNIT=Ucontour_unit)
      ENDIF 

      norm_plist = (/b0,6.96e8_r8,n0,Te_frac,1._r8/)
      eta_plist = (/eta,0.1_r8,etavac,r_eta/)
      kappa_plist = (/kappa_par,kappa_perp,0._r8,1.e8_r8/)
      CALL postxmhd_init(indir,norm_plist,eta_case,eta_plist,
     $     kappa_case,kappa_plist)

      gravity=4.01e-22*n0/b0**2
c-----------------------------------------------------------------------
c     normalizations: time in sec, temperature in K, velocity in km/s
c-----------------------------------------------------------------------
      time0=3.191e-8*SQRT(n0)/b0
      tempr0=5.764e28*b0**2/n0
      vel0=2.181e13*b0/SQRT(n0)
c-----------------------------------------------------------------------
c     initialize variables
c-----------------------------------------------------------------------
      in_check=.TRUE.
      u2_old=0.
      xy0_old(1,:)=x0_list
      xy0_old(2,:)=y0_list
      kt0_old=.5
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_read
c-----------------------------------------------------------------------
c     subprogram 2. postbreakout_dUdt.
c     generate du/dt(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_dUdt(t,xyw,uw,beginner,first,last)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw,uw
      LOGICAL, INTENT(IN) :: beginner,first,last

      REAL(r8), SAVE :: t_old
      REAL(r8), DIMENSION(1) :: value1
      REAL(r8), DIMENSION(1,SIZE(uw,2),SIZE(uw,3)) :: utemp
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("SolarLnHMHD")
         IF(beginner)ALLOCATE(uc_old(SIZE(uw,2),SIZE(uw,3)))
         IF(.NOT. first)THEN
            utemp(1,:,:)=(uw(2,:,:)-uc_old)/(t-t_old)
            CALL plotter_maxUvsT("Efield",0.625*(t+t_old),utemp,value1)
         ENDIF
         uc_old=uw(2,:,:)
         t_old=t
         IF(last)DEALLOCATE(uc_old)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_dUdt
c-----------------------------------------------------------------------
c     subprogram 4. postbreakout_contour.
c     generates contours of desired physical quantities
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_contour(t,xyw,uw,uxyw,header,footer,
     $     nt_next,ifile,stride,filename)

      LOGICAL, INTENT(IN) :: header,footer
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nt_next,ifile,stride
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: uxyw

      INTEGER :: ix,nyw,cyl_fac
      INTEGER, DIMENSION(3), SAVE :: xcont=0
      REAL(r8), PARAMETER :: gamma=5._r8/3._r8
      REAL(r8), DIMENSION(SIZE(xyw,2),SIZE(xyw,3)) :: r_faci,rho
      REAL(r8), DIMENSION(ncont,SIZE(uw,2),SIZE(uw,3)) :: utemp
c-----------------------------------------------------------------------
c     zero utemp.
c-----------------------------------------------------------------------
      utemp=0.
c-----------------------------------------------------------------------
c     transform variables.
c-----------------------------------------------------------------------
      cyl_fac=0
      r_faci=0
      IF(cylinder)THEN
         cyl_fac=1
         r_faci = 1._r8/xyw(2,:,:)
         WHERE(xyw(2,:,:) < 1.e-12)
            r_faci = 0.
         END WHERE
      ENDIF
      IF(job_type=="SolarLnHMHD")THEN
         rho=EXP(uw(1,:,:))
      ELSE
         rho=MAX(uw(1,:,:),1.e-20_r8)
      ENDIF
c-----------------------------------------------------------------------
c     Magnetic Pressure 
c-----------------------------------------------------------------------
      utemp(1,:,:)=.5*(uxyw(:,:,1,2)**2
     $     + (r_faci*uw(2,:,:)+uxyw(:,:,2,2))**2 + uw(3,:,:)**2)
c-----------------------------------------------------------------------
c     Temperature
c-----------------------------------------------------------------------
      utemp(2,:,:)=uw(8,:,:)/rho
c-----------------------------------------------------------------------
c     Bx 
c-----------------------------------------------------------------------
      utemp(3,:,:)=-(r_faci*uw(2,:,:) + uxyw(:,:,2,2))
c-----------------------------------------------------------------------
c     By 
c-----------------------------------------------------------------------
      utemp(4,:,:)=uxyw(:,:,1,2)
c-----------------------------------------------------------------------
c     Mach number
c-----------------------------------------------------------------------
      utemp(5,:,:)=SQRT((uw(4,:,:)**2+uw(5,:,:)**2)
     $     /(gamma*uw(8,:,:)*rho))
c-----------------------------------------------------------------------
c     calculate psi (=-rAphi)
c-----------------------------------------------------------------------
      utemp(6,:,:)=xyw(2,:,:)**cyl_fac*uw(2,:,:)
c-----------------------------------------------------------------------
c     calculate velocity components
c-----------------------------------------------------------------------
      utemp(7,:,:)=uw(4,:,:)/rho
      utemp(8,:,:)=uw(5,:,:)/rho
      utemp(9,:,:)=uw(6,:,:)/rho
c-----------------------------------------------------------------------
c     calculate ideal force balance in radial direction 
c-----------------------------------------------------------------------
      utemp(10,:,:) = utemp(3,:,:)*uw(7,:,:)
      utemp(11,:,:) = -uw(3,:,:)*(r_faci*uw(3,:,:)+uxyw(:,:,2,3))
      utemp(12,:,:) = -uxyw(:,:,2,8)
      utemp(13,:,:) = -gravity*rho*xyw(2,:,:)/(SUM(xyw**2,1))**1.5
      utemp(14,:,:) = SUM(utemp(10:13,:,:),1)
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
      DO ix=1,ncont
         IF(MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
      ENDDO
      CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $     nt_next,xcont,ifile,stride,filename,"Ucontour")
c-----------------------------------------------------------------------
c     calculate and plot spatially varying resistivity
c-----------------------------------------------------------------------
      utemp(1,:,:) = SQRT((uxyw(:,:,2,3) + r_faci*uw(3,:,:))**2
     $     + uxyw(:,:,1,3)**2 + uw(7,:,:)**2)
      CALL postxmhd_ploteta(t,xyw,rho,uw(8,:,:),utemp(1,:,:),
     $        header,footer,nt_next,ifile,stride,filename)
c-----------------------------------------------------------------------
c     calculate and plot spatially varying thermal conductivity
c-----------------------------------------------------------------------
      utemp(1,:,:) = uxyw(:,:,1,2)**2 
     $     + (r_faci*uw(2,:,:)+uxyw(:,:,2,2))**2 + uw(3,:,:)**2
      CALL postxmhd_plotkappa(t,xyw,rho,0.5*uw(8,:,:),
     $     utemp(1,:,:),header,footer,nt_next,ifile,stride,filename)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_contour
c-----------------------------------------------------------------------
c     subprogram 5. postbreakout_energy.
c     generate a time plot for an integral of total energy
c     over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_energy(t,xyw,jac,uw,uxyw)

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
      CASE("SolarHallMHD")
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
c     calculate the integral of vector energies over the domain.
c-----------------------------------------------------------------------
         CALL plotter_VecSqInt(t,xyw,uw1,uw2,uw3,.FALSE.,value(1:2))
c-----------------------------------------------------------------------
c     thermal energy.
c-----------------------------------------------------------------------
         uw1(1,:,:)=r_fac*3._r8*uw(8,:,:)

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
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_energy
c-----------------------------------------------------------------------
c     subprogram 6. postbreakout_maxUvsT.
c     generate u(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_maxUvsT(t,uw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uw

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
      REAL(r8), DIMENSION(:), ALLOCATABLE :: value
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("SolarLnHMHD")
         ALLOCATE(utemp(3,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),value(3))
         utemp(1,:,:)=uw(8,:,:)*EXP(-uw(1,:,:))   ! temperature
         utemp(2,:,:)=uw(5,:,:)*EXP(-uw(1,:,:))   ! r-velocity
         utemp(3,:,:)=uw(7,:,:) ! Jphi.
         CALL plotter_maxUvsT("maxUvsT",t,utemp,value)
         DEALLOCATE(utemp,value)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_maxUvsT
c-----------------------------------------------------------------------
c     subprogram 8. postbreakout_Uprofile.
c     generate slices of U for each time-step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)

      INTEGER, INTENT(IN) :: nxs
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt,uxyw

      CHARACTER(16) :: filename
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: uprofile
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      filename="Uprofile"
      SELECT CASE(job_type)
      CASE("SolarLnHMHD")
         ALLOCATE(utemp(4,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(5,0:nxs))
         utemp(1,:,:)=EXP(uw(1,:,:))
         utemp(2,:,:)=uw(8,:,:)*EXP(-uw(1,:,:))
         utemp(3,:,:)=uw(5,:,:)*EXP(-uw(1,:,:))
         utemp(4,:,:)=uw(7,:,:)
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
         DEALLOCATE(utemp,uprofile)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_Uprofile
c-----------------------------------------------------------------------
c     subprogram 10. postbreakout_recrate.
c     generates reconnection rate plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_recrate(t,xyw,uw,first)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      LOGICAL, INTENT(IN) :: first

      CHARACTER(16) :: filename
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      filename=" "
      SELECT CASE(job_type)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_recrate
c-----------------------------------------------------------------------
c     subprogram 11. postbreakout_cstracking.
c     plot CME current sheet location vs. time.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_cstracking(t,xyw,uw,beginner)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      LOGICAL, INTENT(IN) :: beginner

      INTEGER :: cyl_fac,ix,iy
      REAL(r8) :: flare_flux
      REAL(r8), SAVE :: t_old,temp_old
      INTEGER, DIMENSION(2) :: mindex
      REAL(r8), DIMENSION(2), SAVE :: u_old
      REAL(r8), DIMENSION(2) :: vel
      REAL(r8), DIMENSION(3) :: loc
      REAL(r8), DIMENSION(5) :: value,a
      REAL(r8), DIMENSION(0:SIZE(uw,2)-1,0:SIZE(uw,3)-1,3) :: utemp
c-----------------------------------------------------------------------
c     Quit if not post-processing CME simulation
c-----------------------------------------------------------------------
      IF(init_type /= "half_breakout")RETURN
      IF(beginner)THEN
         temp_old=0
         u_old=0
      ENDIF

      value=0
      vel=0
      cyl_fac=0
      IF(cylinder)cyl_fac=1
      utemp=0
      utemp(:,:,2)=xyw(2,:,:)**cyl_fac*uw(2,:,:)
      utemp(:,:,3)=uw(3,:,:)
c-----------------------------------------------------------------------
c     locate the CME leading CS and calculate RecRate there
c-----------------------------------------------------------------------
      WHERE(xyw(1,:,:) == 0.)
         utemp(:,:,1)=uw(7,:,:)
      END WHERE

      mindex=MAXLOC(utemp(:,:,1))
      loc(1)=xyw(2,mindex(1),mindex(2))

      IF(loc(1) < temp_old)THEN
         WHERE(xyw(2,:,:) < temp_old)
            utemp(:,:,1)=0
         END WHERE
         mindex=MAXLOC(utemp(:,:,1))
         loc(1)=xyw(2,mindex(1),mindex(2))
      ELSE
         temp_old=loc(1)
      ENDIF
      value(1)=eta*uw(7,mindex(1),mindex(2))
      value(4)=loc(1)**cyl_fac*uw(2,mindex(1),mindex(2))
      IF(.NOT. beginner)vel(1)=(loc(1)-u_old(1))/(t-t_old)
      u_old(1) = loc(1)
c-----------------------------------------------------------------------
c     locate the flare reconnection site and calculate RecRate there
c-----------------------------------------------------------------------
      utemp(:,:,1)=0
      WHERE(xyw(1,:,:) == 0.)
         utemp(:,:,1)=uw(7,:,:)
      END WHERE

      mindex=MINLOC(utemp(:,:,1))
      loc(2)=xyw(2,mindex(1),mindex(2))
      value(2)=eta*uw(7,mindex(1),mindex(2))
      flare_flux=loc(2)**cyl_fac*uw(2,mindex(1),mindex(2))
c-----------------------------------------------------------------------
c     locate the CME flux-rope axis and calculate poloidal flux 
c     in the flux-rope
c-----------------------------------------------------------------------
      utemp(:,:,1)=-HUGE(utemp(:,:,1))
      WHERE(xyw(1,:,:) == 0. .AND. xyw(2,:,:) < loc(1) 
     $     .AND. xyw(2,:,:) > loc(2))
         utemp(:,:,1)=xyw(2,:,:)**cyl_fac*uw(2,:,:)
      END WHERE

      mindex=MAXLOC(utemp(:,:,1))
      loc(3)=xyw(2,mindex(1),mindex(2))
      value(3)=loc(3)**cyl_fac*uw(2,mindex(1),mindex(2))-flare_flux
      IF(.NOT. beginner)vel(2)=(loc(3)-u_old(2))/(t-t_old)
      u_old(2) = loc(3)
c-----------------------------------------------------------------------
c     calculate polodal flux between the flux-rope and leading CS
c-----------------------------------------------------------------------
      value(4) = flare_flux-value(4)
c-----------------------------------------------------------------------
c     calculate toroidal B-flux in the CME 
c-----------------------------------------------------------------------
      WHERE(xyw(2,:,:) < loc(2))
         utemp(:,:,3)=0
      END WHERE
      WHERE(utemp(:,:,2) < flare_flux)
         utemp(:,:,3)=0
      END WHERE

      DO ix=0,SIZE(xyw,2)-2
         DO iy=0,SIZE(xyw,3)-2
            a(1)=0.5*(xyw(2,ix,iy)+xyw(2,ix,iy+1))*
     $           (xyw(1,ix,iy+1)-xyw(1,ix,iy))
            a(2)=0.5*(xyw(2,ix,iy+1)+xyw(2,ix+1,iy+1))*
     $           (xyw(1,ix+1,iy+1)-xyw(1,ix,iy+1))
            a(3)=0.5*(xyw(2,ix+1,iy+1)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix+1,iy)-xyw(1,ix+1,iy+1))
            a(4)=0.5*(xyw(2,ix,iy)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix,iy)-xyw(1,ix+1,iy))
            a(5)=ABS(SUM(a(1:4)))
            value(5) = value(5) + 0.25*a(5)*
     $           SUM(utemp(ix:ix+1,iy:iy+1,3))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write binary file for XDRAW
c-----------------------------------------------------------------------
      OPEN(UNIT=dUdt_unit,FILE="breakout.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(dUdt_unit)REAL(time0*t/3.6e3,4),REAL(ABS(value),4),
     $     REAL(loc,4),REAL(time0*(t+t_old)/7.2e3,4),REAL(vel*vel0,4)
      CLOSE(UNIT=dUdt_unit)
      t_old=t
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_cstracking
c-----------------------------------------------------------------------
c     subprogram 12. postbreakout_rTtracking.
c     tracks and plots plasma properties in a given plasma element
c     vs. time.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postbreakout_rTtracking(t,xyw,xyw_kt,uw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt

      CHARACTER(7) :: filename
      INTEGER :: ix0,iy0,iqty
      REAL(r8), SAVE :: time_old=0.
      REAL(r8), DIMENSION(4) :: value
      REAL(r8), DIMENSION(4,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: utemp
c-----------------------------------------------------------------------
c     calculate desired quantities
c-----------------------------------------------------------------------
      IF(listsize <= 0)RETURN
      SELECT CASE(job_type)
      CASE("SolarHallMHD")
         utemp(1,:,:) = 1.e-6*n0*uw(1,:,:)
         utemp(2,:,:) = tempr0*uw(8,:,:)/uw(1,:,:)
         utemp(3,:,:) = uw(4,:,:)/uw(1,:,:)
         utemp(4,:,:) = uw(5,:,:)/uw(1,:,:)
         WHERE(uw(1,:,:) == 0.)
            utemp(2,:,:) = 0.
            utemp(3,:,:) = 0.
            utemp(4,:,:) = 0.
         END WHERE
      CASE("SolarLnHMHD")
         utemp(1,:,:) = 1.e-6*n0*EXP(uw(1,:,:))
         utemp(2,:,:) = tempr0*uw(8,:,:)*EXP(-uw(1,:,:))
         utemp(3,:,:) = uw(4,:,:)*EXP(-uw(1,:,:))
         utemp(4,:,:) = uw(5,:,:)*EXP(-uw(1,:,:))
      END SELECT
      xy0_old = xy0_old + (t-time_old)*u2_old
      time_old=t
c-----------------------------------------------------------------------
c     find values along the fluid element track
c-----------------------------------------------------------------------
      DO iqty=1,listsize
         IF(.NOT. in_check(iqty))CYCLE
         CALL plotter_newton_search(xy0_old(:,iqty),xyw,xyw_kt,
     $        kt0_old(:,iqty),in_check(iqty))         
         IF(.NOT. in_check(iqty))CYCLE
         ix0 = INT(kt0_old(1,iqty)*(SIZE(xyw,2)-1))
         iy0 = INT(kt0_old(2,iqty)*(SIZE(xyw,3)-1))
         value = plotter_evaluate(xy0_old(1,iqty),xy0_old(2,iqty),
     $        xyw(:,ix0:ix0+1,iy0:iy0+1),utemp(:,ix0:ix0+1,iy0:iy0+1))
         u2_old(:,iqty)=value(3:4)
c-----------------------------------------------------------------------
c     write binary file for XDRAW and text file for further analysis
c-----------------------------------------------------------------------
         value(3:4) = vel0*value(3:4)
         WRITE(filename,'(a,i2.2)')"track",iqty
         OPEN(UNIT=UxyT_unit,FILE=filename//".bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(UxyT_unit)REAL(time0*t,4),REAL(xy0_old(:,iqty),4),
     $        REAL(SQRT(SUM(xy0_old(:,iqty)**2))-1,4),REAL(value,4)
         CLOSE(UNIT=UxyT_unit)
         
         OPEN(UNIT=UxyT_unit,FILE=filename//".txt",STATUS="UNKNOWN",
     $        POSITION="APPEND")
         WRITE(UxyT_unit,'(1p,8(e15.5))')time0*t,xy0_old(:,iqty),
     $        SQRT(SUM(xy0_old(:,iqty)**2))-1,value
         CLOSE(UNIT=UxyT_unit)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postbreakout_rTtracking
      END MODULE postbreakout_mod
