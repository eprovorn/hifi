c-----------------------------------------------------------------------
c     file posteuler.f.
c     post-processes output from SEL code for pseudo-euler problem.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. posteuler_mod.
c     1. posteuler_read.
c     2. posteuler_Uprofile.
c     3. posteuler_contour.
c     4. posteuler_energy.
c-----------------------------------------------------------------------
c     subprogram 0. posteuler_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE posteuler_mod
      USE plotter_mod
      IMPLICIT NONE

      CHARACTER(20), PRIVATE :: init_type
      CHARACTER(50), PRIVATE :: cubit_file
      LOGICAL, PRIVATE :: ifbound_visc,robin_tmp,ifnoslip,ifbound_diff,
     $     ifconstmu,ifconsteddiff
      INTEGER, PRIVATE :: ncont
      INTEGER, DIMENSION(3), PRIVATE :: xcont=0
      REAL(r8), PARAMETER, PRIVATE :: gamma=1.4_r8
      REAL(r8), PRIVATE :: x1,y1,x2,y2,mach,mu,pin,pout,lx,alpha,rhoin,
     $     vin,rhoout,vout,alpha2,ddiff,Lconv,Ldiv,ediff,eps,muy,mux,
     $     Aconv,Adiv,pmax,rhomax,expand,mu_h,rho0,p0,poutmax

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. posteuler_read.
c     read necessary post-processing parameters from post.in 
c     and input parameters from sel.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE posteuler_read(indir,out_type)

      CHARACTER(*), INTENT(IN) :: indir,out_type

      CHARACTER(80) :: filename

      NAMELIST/euler1D_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $  dtmax,tmax,nstep,init_type,mach,mu,pin,pout,lx
      NAMELIST/euler2D_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $  dtmax,tmax,nstep,init_type,cubit_file,mach,mu,pin,pout,lx,alpha
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      filename=TRIM(indir)//"/sel.in"

      SELECT CASE(job_type)
      CASE("euler1D")
         OPEN(UNIT=in_unit,FILE=filename,STATUS="OLD")
         READ(in_unit,NML=euler1D_list)
         CLOSE(UNIT=in_unit)
      CASE("euler2D")
         CALL system("rm -f Uprofile_x0.bin Uprofile_y0.bin") 
         OPEN(UNIT=in_unit,FILE=filename,STATUS="OLD")
         READ(in_unit,NML=euler2D_list)
         CLOSE(UNIT=in_unit)
         ncont=8
         IF(flag2d .AND. out_type /= "hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="Ucontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE posteuler_read
c-----------------------------------------------------------------------
c     subprogram 2. posteuler_Uprofile.
c     generate slices of U for each time-step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE posteuler_Uprofile(nxs,xyw,xyw_kt,uw)

      INTEGER, INTENT(IN) :: nxs
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt

      CHARACTER(16) :: filename
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: uprofile
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
      REAL(r8), DIMENSION(SIZE(xyw,2),SIZE(xyw,3)) :: area,p
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("euler1D")
         filename="Uprofile"
         area = 1.398 + 0.347*TANH(0.8*xyw(1,:,:) - 4.0) 
         ALLOCATE(utemp(5,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(6,0:nxs))
         x1 = 1e-6
         x2 = lx-1e-6
         y1 = 0.5
         y2 = 0.5
         utemp(1,:,:)=uw(1,:,:)/area
         utemp(2,:,:)=uw(2,:,:)/area
         utemp(3,:,:)=uw(3,:,:)/area
         utemp(4,:,:)=uw(2,:,:)/uw(1,:,:)
         utemp(5,:,:)= (gamma-1)*(utemp(3,:,:) 
     $        - 0.5*utemp(2,:,:)**2/utemp(1,:,:))
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
         DEALLOCATE(utemp,uprofile)
      CASE("euler2D")
         ALLOCATE(utemp(ncont,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(ncont+1,0:nxs))
c-----------------------------------------------------------------------
c     density
c-----------------------------------------------------------------------
         utemp(1,:,:) = uw(1,:,:)
c-----------------------------------------------------------------------
c     x-momentum
c-----------------------------------------------------------------------
         utemp(2,:,:) = uw(2,:,:)
c-----------------------------------------------------------------------
c     y-momentum
c-----------------------------------------------------------------------
         utemp(3,:,:) = uw(3,:,:)
c-----------------------------------------------------------------------
c     x-velocity
c-----------------------------------------------------------------------
         utemp(4,:,:) = uw(2,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     energy
c-----------------------------------------------------------------------
         utemp(5,:,:) = uw(4,:,:)
c-----------------------------------------------------------------------
c     pressure
c-----------------------------------------------------------------------
         utemp(6,:,:) = (gamma-1)*(utemp(5,:,:)
     $        - 0.5/utemp(1,:,:)*( utemp(2,:,:)**2 + utemp(3,:,:)**2))
c-----------------------------------------------------------------------
c     mach number
c-----------------------------------------------------------------------
         p = utemp(6,:,:)
         utemp(7,:,:)=SQRT(uw(2,:,:)**2+uw(3,:,:)**2)/
     $        (SQRT(gamma*p/uw(1,:,:))*uw(1,:,:))
c-----------------------------------------------------------------------
c     stagnation pressure
c-----------------------------------------------------------------------
         utemp(8,:,:)=
     $        p*(1.+(gamma-1.)/2.*utemp(7,:,:)**2)**(gamma/(gamma-1.))
         
         filename="Uprofile_y0"
         SELECT CASE(init_type)
         CASE("open_nozzle")
            x1 = -Lconv + 1e-6
            x2 = Ldiv + ((2. - COS(Adiv))
     $           + TAN(Adiv)*(Ldiv - SIN(Adiv)))*TAN(pi/8.) - 1e-6
         CASE("open_nozzle2")
            x1 = -Lconv + 1e-6
            x2 = Ldiv + Ldiv/4.*TAN(pi/8.) - 1e-6
         CASE DEFAULT
            x1 = 1e-6
            x2 = lx-1e-6
         END SELECT
         y1 = 0.0
         y2 = 0.0
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
         
         filename="Uprofile_x0"
         x1 = 0
         x2 = 0
         y1 = -0.5
         y2 = +0.5
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
         
         DEALLOCATE(utemp,uprofile)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE posteuler_Uprofile
c-----------------------------------------------------------------------
c     subprogram 3. posteuler_contour.
c     generates contours of desired physical quantities
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE posteuler_contour(t,xyw,uw,header,
     $     footer,nt_next,ifile,stride,filename)

      LOGICAL, INTENT(IN) :: header,footer
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nt_next,ifile,stride
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw

      INTEGER :: ix
      REAL(r8), DIMENSION(ncont,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) 
     $     :: utemp
      REAL(r8), DIMENSION(0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: p
c-----------------------------------------------------------------------
c     calculate Density
c-----------------------------------------------------------------------
      utemp(1,:,:)=uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate velx
c-----------------------------------------------------------------------
      utemp(2,:,:)=uw(2,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate vely
c-----------------------------------------------------------------------
      utemp(3,:,:)=uw(3,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate Energy
c-----------------------------------------------------------------------
      utemp(4,:,:)=uw(4,:,:) 
c-----------------------------------------------------------------------
c     calculate Mach number.
c-----------------------------------------------------------------------
      p = (gamma-one)*(uw(4,:,:)
     $     - half*(uw(2,:,:)**2+uw(3,:,:)**2)/uw(1,:,:))
      utemp(5,:,:)=SQRT(uw(2,:,:)**2+uw(3,:,:)**2)/
     $     (SQRT(gamma*p/uw(1,:,:))*uw(1,:,:))
c-----------------------------------------------------------------------
c     calculate Temperature
c-----------------------------------------------------------------------
      utemp(6,:,:)=p/uw(1,:,:)
c-----------------------------------------------------------------------
c     use calculated pressure
c-----------------------------------------------------------------------
      utemp(7,:,:)=p
c-----------------------------------------------------------------------
c     calculate stagnation pressure
c-----------------------------------------------------------------------
      utemp(8,:,:)=
     $     p*(1.+(gamma-1.)/2.*utemp(5,:,:)**2)**(gamma/(gamma-1.))
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
      DO ix=1,ncont
         IF(MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
      ENDDO
      CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $     nt_next,xcont,ifile,stride,filename,"Ucontour")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE posteuler_contour
c-----------------------------------------------------------------------
c     subprogram 4. euler_energy.
c     generate a time plot for an integral of total energy
c     over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE posteuler_energy(t,xyw,jac,uw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: jac
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw

      INTEGER :: nx,ny
      REAL(r8), DIMENSION(1) :: value
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uw1,uw2,uw3
c-----------------------------------------------------------------------
c     prepare auxiliary parameters
c-----------------------------------------------------------------------
      value=zero
      nx=SIZE(xyw,2)-1
      ny=SIZE(xyw,3)-1
c-----------------------------------------------------------------------
c     calculate energy integral
c-----------------------------------------------------------------------
      ALLOCATE(uw1(1,0:nx,0:ny))
      ALLOCATE(uw2(1,0:nx,0:ny))
      ALLOCATE(uw3(1,0:nx,0:ny))
c-----------------------------------------------------------------------
c     kinetic energy components.
c-----------------------------------------------------------------------
      uw1(1,:,:)=uw(2,:,:)/SQRT(2*uw(1,:,:))
      uw2(1,:,:)=uw(3,:,:)/SQRT(2*uw(1,:,:))
      uw3(1,:,:)=0.0
c-----------------------------------------------------------------------
c     calculate the integral of vector energies over the domain.
c-----------------------------------------------------------------------
      CALL plotter_VecSqInt(t,xyw,uw1,uw2,uw3,.TRUE.,value)
      DEALLOCATE(uw1,uw2,uw3)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE posteuler_energy
      END MODULE posteuler_mod
