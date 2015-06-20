c-----------------------------------------------------------------------
c     file postmast.f.
c     post-processes output from SEL code for the MAST simulations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. postmast_mod.
c     1. postmast_read.
c     2. postmast_Btor.
c     3. postmast_dUdt.
c     4. postmast_contour.
c     5. postmast_energy.
c-----------------------------------------------------------------------
c     subprogram 0. postmast_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE postmast_mod
      USE postxmhd_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: source,cylinder=.FALSE.
      CHARACTER(20), PRIVATE :: init_type,eta_case,kappa_case,nu_case
      REAL(r8), PARAMETER, PRIVATE :: nol=1.e-6
      REAL(r8), PRIVATE :: eta,mu,nu,di,lx,ly,kappa_par,kappa_perp,
     $     kappa_pare,kappa_perpe,Dn,n0,L0,b0,lambda,etavac,r_nu,nuvac,
     $     Bguide,ieheat,beta0,beta_e

      INTEGER, PRIVATE :: ncont
      REAL(r8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: uc_old

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. postmast_read.
c     read necessary post-processing parameters from post.in 
c     and input parameters from sel.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postmast_read(indir)

      CHARACTER(*), INTENT(IN) :: indir
      CHARACTER(80) :: infile

      REAL(r8), DIMENSION(4) :: eta_plist,kappa_plist
      REAL(r8), DIMENSION(5) :: norm_plist

      NAMELIST/MAST_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,init_type,source,gr_curve,di,eta,eta_case,
     $     etavac,mu,nu_case,r_nu,nuvac,nu,kappa_case,kappa_par,
     $     kappa_perp,kappa_pare,kappa_perpe,ieheat,Dn,n0,b0,L0,Bguide,
     $     lambda,lx,ly,beta0,beta_e
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      infile=TRIM(indir)//"/sel.in"
      OPEN(UNIT=in_unit,FILE=infile,STATUS="OLD")
      READ(in_unit,NML=MAST_list)
      CLOSE(UNIT=in_unit)
      ncont=10
      IF(flag1d)CALL system("rm -f mast.bin")
      IF(flag2d .AND. out_type/="hdf5")THEN
         OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $        STATUS="UNKNOWN",FORM="UNFORMATTED")
         WRITE(Ucontour_unit)1,0,ncont
         CLOSE(UNIT=Ucontour_unit)
      ENDIF 

      norm_plist = (/b0,L0,n0,1._r8,2._r8/)
      eta_plist = (/eta,0._r8,etavac,1._r8/)
      kappa_plist = (/kappa_par,kappa_perp,0._r8,1.e8_r8/)
      CALL postxmhd_init(indir,norm_plist,eta_case,eta_plist,
     $     kappa_case,kappa_plist)
c-----------------------------------------------------------------------
c     initialize variables
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("cylFF")
         cylinder=.TRUE.
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postmast_read
c-----------------------------------------------------------------------
c     subprogram 2. postmast_Btor.
c     prescribes MAST toroidal guide field.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION postmast_Btor(Bguide,r) RESULT(Bz)

      REAL(r8), INTENT(IN) :: Bguide
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: r
      REAL(r8), DIMENSION(SIZE(r,1),SIZE(r,2)) :: Bz
c-----------------------------------------------------------------------
c     compute Bz.
c-----------------------------------------------------------------------           
      Bz = -Bguide*0.85_r8/r
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION postmast_Btor
c-----------------------------------------------------------------------
c     subprogram 3. postmast_dUdt.
c     generate du/dt(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postmast_dUdt(t,xyw,uw,beginner,first,last)

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
      CASE(".")
         IF(beginner)ALLOCATE(uc_old(SIZE(uw,2),SIZE(uw,3)))
         IF(.NOT. first)THEN
            utemp(1,:,:)=(uw(2,:,:)-uc_old)/(t-t_old)
            CALL plotter_maxUvsT("Efield",0.5*(t+t_old),utemp,value1)
         ENDIF
         uc_old=uw(2,:,:)
         t_old=t
         IF(last)DEALLOCATE(uc_old)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postmast_dUdt
c-----------------------------------------------------------------------
c     subprogram 4. postmast_contour.
c     generates contours of desired physical quantities
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postmast_contour(t,xyw,uw,uxyw,header,footer,
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
      REAL(r8), DIMENSION(SIZE(xyw,2),SIZE(xyw,3)) :: r_faci,rho,Bz
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
      rho=MAX(uw(1,:,:),1.e-20_r8)

      Bz=uw(3,:,:)
      SELECT CASE(init_type)
      CASE("cylFF")
         Bz = Bz + postmast_Btor(Bguide,xyw(2,:,:))
      END SELECT
c-----------------------------------------------------------------------
c     Magnetic Pressure 
c-----------------------------------------------------------------------
      utemp(1,:,:)=.5*(uxyw(:,:,1,2)**2
     $     + (r_faci*uw(2,:,:)+uxyw(:,:,2,2))**2 + Bz**2)
c-----------------------------------------------------------------------
c     Ion Temperature
c-----------------------------------------------------------------------
      utemp(2,:,:)=uw(8,:,:)/rho
c-----------------------------------------------------------------------
c     Electron Temperature
c-----------------------------------------------------------------------
      IF(SIZE(uw,1) > 9)utemp(3,:,:)=uw(10,:,:)/rho
c-----------------------------------------------------------------------
c     Bx 
c-----------------------------------------------------------------------
      utemp(4,:,:)=-(r_faci*uw(2,:,:) + uxyw(:,:,2,2))
c-----------------------------------------------------------------------
c     By 
c-----------------------------------------------------------------------
      utemp(5,:,:)=uxyw(:,:,1,2)
c-----------------------------------------------------------------------
c     Mach number
c-----------------------------------------------------------------------
      utemp(6,:,:)=SQRT((uw(4,:,:)**2+uw(5,:,:)**2)
     $     /(gamma*uw(8,:,:)*rho))
c-----------------------------------------------------------------------
c     calculate psi (=-rAphi)
c-----------------------------------------------------------------------
      utemp(7,:,:)=xyw(2,:,:)**cyl_fac*uw(2,:,:)
c-----------------------------------------------------------------------
c     calculate velocity components
c-----------------------------------------------------------------------
      utemp(8,:,:)=uw(4,:,:)/rho
      utemp(9,:,:)=uw(5,:,:)/rho
      utemp(10,:,:)=uw(6,:,:)/rho
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

      IF(SIZE(uw,1) > 9)THEN
         CALL postxmhd_ploteta(t,xyw,rho,uw(10,:,:),utemp(1,:,:),
     $        header,footer,nt_next,ifile,stride,filename)
      ELSE
         CALL postxmhd_ploteta(t,xyw,rho,uw(8,:,:),utemp(1,:,:),
     $        header,footer,nt_next,ifile,stride,filename)
      ENDIF
c-----------------------------------------------------------------------
c     calculate and plot spatially varying thermal conductivity
c-----------------------------------------------------------------------
      utemp(1,:,:) = uxyw(:,:,1,2)**2 
     $     + (r_faci*uw(2,:,:)+uxyw(:,:,2,2))**2 + Bz**2
      IF(SIZE(uw,1) > 9)THEN
         CALL postxmhd_plotkappa(t,xyw,rho,uw(10,:,:),
     $        utemp(1,:,:),header,footer,nt_next,ifile,stride,filename)
      ELSE
         CALL postxmhd_plotkappa(t,xyw,rho,uw(8,:,:),
     $        utemp(1,:,:),header,footer,nt_next,ifile,stride,filename)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postmast_contour
c-----------------------------------------------------------------------
c     subprogram 5. postmast_energy.
c     generate a time plot for an integral of total energy
c     over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postmast_energy(t,xyw,jac,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: jac
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: uxyw

      INTEGER :: nx,ny
      REAL(r8) :: cyl_fac
      REAL(r8), DIMENSION(5) :: value
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uw1,uw2,uw3
      REAL(r8), DIMENSION(SIZE(xyw,2),SIZE(xyw,3)) :: Bz,
     $     r_fac,r_facs,r_faci
c-----------------------------------------------------------------------
c     prepare auxiliary parameters
c-----------------------------------------------------------------------
      value=zero
      nx=SIZE(xyw,2)
      ny=SIZE(xyw,3)

      Bz=uw(3,:,:)
      SELECT CASE(init_type)
      CASE("cylFF")
         Bz = Bz + postmast_Btor(Bguide,xyw(2,:,:))
      END SELECT
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
      uw3(2,:,:)=Bz
c-----------------------------------------------------------------------
c     calculate the integral of vector energies over the domain.
c-----------------------------------------------------------------------
      CALL plotter_VecSqInt(t,xyw,uw1,uw2,uw3,.FALSE.,value(1:2))
c-----------------------------------------------------------------------
c     thermal energy and magnetic helicity.
c-----------------------------------------------------------------------
      uw1(1,:,:)=r_fac*3._r8*uw(8,:,:)
      IF(SIZE(uw,1) > 9)uw1(1,:,:)=uw1(1,:,:) + r_fac*3._r8*uw(10,:,:)
      uw1(2,:,:)=r_fac*two*uw(2,:,:)*Bz

      CALL plotter_integral(t,xyw,uw1,.FALSE.,value(3:4))
c-----------------------------------------------------------------------
c     total energy.
c-----------------------------------------------------------------------
      value(1:3)=0.5*value(1:3)
      value(5)=SUM(value(1:3))

      IF(cylinder)value=value*twopi
      DEALLOCATE(uw1,uw2,uw3)
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
      END SUBROUTINE postmast_energy
      END MODULE postmast_mod
