c-----------------------------------------------------------------------
c     file LnNavStok.f.
c     contains specifications for 2D Navier-Stokes equations
c     in cylindrical coordinates using r*ln(rho) and r*ln(T) 
c     as primary dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. LnNavStok_mod.
c-----------------------------------------------------------------------
c     external subprograms.
c-----------------------------------------------------------------------
c     a. physics_input.
c     b. physics_init_parameters.
c     c. physics_init.
c     d. physics_boundary.
c     e. physics_edge_rhs.
c     f. physics_edge_drdu.
c     g. physics_edge_mass.
c     j. physics_rhs.
c     k. physics_drdu.
c     l. physics_mass.
c     m. physics_grid.
c     n. physics_schur.
c     o. physics_dealloc.
c     p. physics_main
c-----------------------------------------------------------------------
c     subprogram 0. LnNavStok_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE LnNavStok_mod
      USE local_mod
      IMPLICIT NONE

      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,eps=1.e-10_r8
      CHARACTER(16) :: init_type="."
      INTEGER :: ntheta=0
      REAL(r8) :: mach=1.,tau=1.,mu=0.,kappa=0.,rho0=1.,rhomax=1.,
     $     p0=1.,pmax=1.,lx=1.0,alpha=0.,gr_curve=0.


      END MODULE LnNavStok_mod
c-----------------------------------------------------------------------
c     subprogram a. physics_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $     nqty_schur,dt,dtmax,tmax,nstep,du_diagnose,physics_type,
     $     exit_flag)
      USE LnNavStok_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: physics_type
      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty,nqty_schur
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep,exit_flag
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax

      INTEGER :: myios
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/NavStok_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,init_type,mach,tau,mu,kappa,rho0,rhomax,lx,
     $     p0,pmax,alpha,ntheta,gr_curve
c-----------------------------------------------------------------------
c     definitions.
c-----------------------------------------------------------------------
c     mach - maximum inflow velocity
c     tau - characteristic flow-drive time
c     mu - flow viscosity
c     kappa - heat conduction
c     rho0 - background density
c     p0 - background pressure
c     rhomax - central density
c     pmax - central pressure
c     lx - domain radius
c     alpha - magnitude of the azimuthal perturbation
c     ntheta - azimuthal mode of the perturbation
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=NavStok_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="LnNavStok"

      nqty=4
      nqty_schur=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_input
c-----------------------------------------------------------------------
c     subprogram b. physics_init_parameters.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_init_parameters(static,ground,adapt_qty)
      USE LnNavStok_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      ! none.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(ntheta,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(tau,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rho0,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(rhomax,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(p0,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(pmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_init_parameters
c-----------------------------------------------------------------------
c     subprogram c. physics_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_init(x,y,u)
      USE LnNavStok_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u
c-----------------------------------------------------------------------
c     initial condition.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("pellet1","pellet2")
         u(1,:,:) = x*LOG(rho0 + half*(rhomax-rho0)*(COS(x*pi/lx)+1))
     $        + x*LOG((one + alpha*COS(ntheta*y)*SIN(x*pi/lx)))
         u(2,:,:) = 0
         u(3,:,:) = 0
         u(4,:,:) = x*LOG(p0 + half*(pmax-p0)*(COS(x*pi/lx)+1))
     $        - u(1,:,:)
      END SELECT

      RETURN
      END SUBROUTINE physics_init
c-----------------------------------------------------------------------
c     subprogram d. physics_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_boundary(left,right,top,bottom,nqty,edge_order)
      USE LnNavStok_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      top%bc_type="robin"
      top%static=.FALSE.
      bottom%bc_type="robin"
      bottom%static=.FALSE.
      left%bc_type="robin"
      left%static=.FALSE.
      right%bc_type="robin"
      right%static=.FALSE.

      SELECT CASE(init_type)
      CASE("pellet1","pellet2")
         left%bc_type(1:4)="robin"
         left%static(1:4)=.TRUE.
         right%bc_type(1:4)="robin"
         right%static(1:4)=.TRUE.
         top%bc_type="periodic"
         bottom%bc_type="periodic"
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_boundary
c-----------------------------------------------------------------------
c     subprogram e. physics_edge_rhs.
c     computes rhs for edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_rhs(lrtb,t,x,y,nhat,u,ux,uy,uxx,uyy,uxy,c)
      USE LnNavStok_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c  = 0
c-----------------------------------------------------------------------
c     boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("pellet1")
         SELECT CASE(lrtb)
         CASE("left")
            c = u
         CASE("right")
            c(1,:,:) = ux(1,:,:) + ux(4,:,:) - (u(1,:,:) + u(4,:,:))/x
            IF(t < tau)THEN
               c(2,:,:) = u(2,:,:) + x*mach*(SIN(t*pi/tau))**2
            ELSE
               c(2,:,:) = u(2,:,:)
            ENDIF
            c(3,:,:) = u(3,:,:)
            c(4,:,:) = u(4,:,:) - x*LOG(p0/rho0)
         END SELECT
      CASE("pellet2")
         SELECT CASE(lrtb)
         CASE("left")
            c = u
         CASE("right")
            c(1,:,:) = u(1,:,:) - x*LOG(rho0)
            IF(t < tau)THEN
               c(2,:,:) = u(2,:,:) 
     $              + x*mach*EXP(u(1,:,:)/x)*(SIN(t*pi/tau))**2
            ELSE
               c(2,:,:) = u(2,:,:)
            ENDIF
            c(3,:,:) = u(3,:,:)
            c(4,:,:) = ux(1,:,:) + ux(4,:,:) - (u(1,:,:) + u(4,:,:))/x
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_edge_rhs
c-----------------------------------------------------------------------
c     subprogram f. physics_edge_drdu.
c     computes drdu for edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_drdu(lrtb,t,x,y,nhat,
     $     u,ux,uy,uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)
      USE LnNavStok_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy
c-----------------------------------------------------------------------
c     zero arrays.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
      c_uxy=0
c-----------------------------------------------------------------------
c     boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("pellet1")
         SELECT CASE(lrtb)
         CASE("left")
            c_u(1,1,:,:) = one
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_u(4,4,:,:) = one
         CASE("right")
            c_u(1,1,:,:) = -one/x
            c_u(1,4,:,:) = -one/x
            c_ux(1,1,:,:) = one
            c_ux(1,4,:,:) = one
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_u(4,4,:,:) = one
         END SELECT
      CASE("pellet2")
         SELECT CASE(lrtb)
         CASE("left")
            c_u(1,1,:,:) = one
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_u(4,4,:,:) = one
         CASE("right")
            c_u(1,1,:,:) = one
            IF(t < tau)
     $           c_u(2,1,:,:) = mach*EXP(u(1,:,:)/x)*(SIN(t*pi/tau))**2
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_u(4,1,:,:) = -one/x
            c_u(4,4,:,:) = -one/x
            c_ux(4,1,:,:) = one
            c_ux(4,4,:,:) = one
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_edge_drdu
c-----------------------------------------------------------------------
c     subprogram g. physics_edge_mass.
c     computes mass matrices for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_mass(lrtb,x,y,nhat,mass,mass_x,mass_y)
      USE LnNavStok_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     set mass,mass_x,mass_y(iqty,jqty,:,:), ... to coupling 
c     mass matrices for du/dt, d(du/dx)/dt, and d(du/dy)/dt
c-----------------------------------------------------------------------
      mass=0
      mass_x=0
      mass_y=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_edge_mass
c-----------------------------------------------------------------------
c     subprogram j. physics_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_rhs(t,x,y,u,ux,uy,fx,fy,s,first)
      USE LnNavStok_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first
c-----------------------------------------------------------------------
c     initialize
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     density equation
c-----------------------------------------------------------------------
      s(1,:,:) = -EXP(-u(1,:,:)/x)*(ux(2,:,:)+uy(3,:,:))
c-----------------------------------------------------------------------
c     radial momentum equation
c-----------------------------------------------------------------------
      fx(2,:,:) = u(2,:,:)**2*EXP(-u(1,:,:)/x)/x 
     $     + x*EXP((u(1,:,:)+u(4,:,:))/x)
     $     - mu*(ux(2,:,:) + u(2,:,:)/x*(u(1,:,:)/x - ux(1,:,:)))
      fy(2,:,:) = u(2,:,:)*u(3,:,:)*EXP(-u(1,:,:)/x)/x
     $     - mu*(uy(2,:,:)/x**2 - two*u(3,:,:)/x 
     $     - uy(1,:,:)*u(2,:,:)/x**3)
      s(2,:,:) = u(3,:,:)**2*EXP(-u(1,:,:)/x) 
     $     + EXP((u(1,:,:) + u(4,:,:))/x)
     $     - mu*(ux(2,:,:)/x - uy(1,:,:)*u(3,:,:)/x**2)
c-----------------------------------------------------------------------
c     azimuthal momentum equation
c-----------------------------------------------------------------------
      fx(3,:,:) = u(2,:,:)*u(3,:,:)*EXP(-u(1,:,:)/x)
     $     - mu*(x*ux(3,:,:) + u(3,:,:)*(u(1,:,:)/x - ux(1,:,:)))
      fy(3,:,:) = u(3,:,:)**2*EXP(-u(1,:,:)/x)
     $     + EXP((u(1,:,:) + u(4,:,:))/x)
     $     - mu*(uy(3,:,:)/x + two*u(2,:,:)/x**2 
     $     - uy(1,:,:)*u(3,:,:)/x**2)
      s(3,:,:) = -(u(2,:,:)*u(3,:,:)*EXP(-u(1,:,:)/x)
     $     + mu*(u(3,:,:) + uy(1,:,:)*u(2,:,:)/x**2))/x
c-----------------------------------------------------------------------
c     temperature equation
c-----------------------------------------------------------------------
      fx(4,:,:) = (gamma-one)*u(2,:,:)*EXP(-u(1,:,:)/x)
     $     - kappa*(gamma-one)*(ux(4,:,:) - u(4,:,:)/x)
      fy(4,:,:) = (gamma-one)*u(3,:,:)*EXP(-u(1,:,:)/x)
     $     - kappa*(gamma-one)*uy(4,:,:)/x**2
      s(4,:,:) = 
     $     -EXP(-u(1,:,:)/x)/x
     $     *(u(2,:,:)*(ux(4,:,:) - u(4,:,:)/x) + u(3,:,:)*uy(4,:,:))
     $     + kappa*(gamma-one)
     $     *((ux(1,:,:) - u(1,:,:)/x)*(ux(4,:,:) - u(4,:,:)/x)/x
     $     + (ux(4,:,:) - u(4,:,:)/x)**2/x
     $     + uy(1,:,:)*uy(4,:,:)/x**3 + uy(4,:,:)**2/x**3)
     $     + mu*(gamma-one)*x*EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *((ux(2,:,:)/x - u(2,:,:)/x**2 - ux(1,:,:)*u(2,:,:)/x**2
     $     + u(1,:,:)*u(2,:,:)/x**3)**2
     $     + (uy(2,:,:)/x**2 - u(3,:,:)/x - uy(1,:,:)*u(2,:,:)/x**3)**2
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/x 
     $     + u(1,:,:)*u(3,:,:)/x**2)**2
     $     + (uy(3,:,:)/x + u(2,:,:)/x**2 - uy(1,:,:)*u(3,:,:)/x**2)**2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_rhs
c-----------------------------------------------------------------------
c     subprogram k. physics_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_drdu(t,x,y,u,ux,uy,fx_u,fx_ux,fx_uy,fy_u,
     $     fy_ux,fy_uy,s_u,s_ux,s_uy)
      USE LnNavStok_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx_u=0
      fx_ux=0
      fx_uy=0
      fy_u=0
      fy_ux=0
      fy_uy=0
      s_u=0
      s_ux=0
      s_uy=0
c-----------------------------------------------------------------------
c     density equation
c-----------------------------------------------------------------------
      s_u(1,1,:,:) = EXP(-u(1,:,:)/x)*(ux(2,:,:)+uy(3,:,:))/x
      s_ux(1,2,:,:) = -EXP(-u(1,:,:)/x)
      s_uy(1,3,:,:) = -EXP(-u(1,:,:)/x)
c-----------------------------------------------------------------------
c     radial momentum equation
c-----------------------------------------------------------------------
      fx_u(2,1,:,:) = EXP((u(1,:,:)+u(4,:,:))/x)
     $     - u(2,:,:)*(u(2,:,:)*EXP(-u(1,:,:)/x) + mu)/x**2
      fx_u(2,2,:,:) = (two*u(2,:,:)*EXP(-u(1,:,:)/x)
     $     - mu*(u(1,:,:)/x - ux(1,:,:)))/x
      fx_u(2,4,:,:) = EXP((u(1,:,:)+u(4,:,:))/x)
      fx_ux(2,1,:,:) = mu*u(2,:,:)/x
      fx_ux(2,2,:,:) = -mu

      fy_u(2,1,:,:) = -u(2,:,:)*u(3,:,:)*EXP(-u(1,:,:)/x)/x**2
      fy_u(2,2,:,:) = u(3,:,:)*EXP(-u(1,:,:)/x)/x + mu*uy(1,:,:)/x**3
      fy_u(2,3,:,:) = (u(2,:,:)*EXP(-u(1,:,:)/x) + mu*two)/x
      fy_uy(2,1,:,:) = mu*u(2,:,:)/x**3
      fy_uy(2,2,:,:) = -mu/x**2

      s_u(2,1,:,:) = -u(3,:,:)**2*EXP(-u(1,:,:)/x)/x 
     $     + EXP((u(1,:,:) + u(4,:,:))/x)/x
      s_u(2,3,:,:) = two*u(3,:,:)*EXP(-u(1,:,:)/x) + mu*uy(1,:,:)/x**2
      s_u(2,4,:,:) =  EXP((u(1,:,:) + u(4,:,:))/x)/x
      s_ux(2,2,:,:) =  -mu/x
      s_uy(2,1,:,:) = mu*u(3,:,:)/x**2
c-----------------------------------------------------------------------
c     azimuthal momentum equation
c-----------------------------------------------------------------------
      fx_u(3,1,:,:) = -u(3,:,:)*(u(2,:,:)*EXP(-u(1,:,:)/x) + mu)/x
      fx_u(3,2,:,:) = u(3,:,:)*EXP(-u(1,:,:)/x)
      fx_u(3,3,:,:) = u(2,:,:)*EXP(-u(1,:,:)/x)
     $     - mu*(u(1,:,:)/x - ux(1,:,:))
      fx_ux(3,1,:,:) = mu*u(3,:,:)
      fx_ux(3,3,:,:) = -mu*x

      fy_u(3,1,:,:) = -u(3,:,:)**2*EXP(-u(1,:,:)/x)/x
     $     + EXP((u(1,:,:) + u(4,:,:))/x)/x
      fy_u(3,2,:,:) = -mu*two/x**2
      fy_u(3,3,:,:) = two*u(3,:,:)*EXP(-u(1,:,:)/x) + mu*uy(1,:,:)/x**2
      fy_u(3,4,:,:) = EXP((u(1,:,:) + u(4,:,:))/x)/x
      fy_uy(3,1,:,:) = mu*u(3,:,:)/x**2
      fy_uy(3,3,:,:) = -mu/x

      s_u(3,1,:,:) = u(2,:,:)*u(3,:,:)*EXP(-u(1,:,:)/x)/x**2
      s_u(3,2,:,:) = -u(3,:,:)*EXP(-u(1,:,:)/x)/x - mu*uy(1,:,:)/x**3
      s_u(3,3,:,:) = -(u(2,:,:)*EXP(-u(1,:,:)/x) + mu)/x
      s_uy(3,1,:,:) = -mu*u(2,:,:)/x**3
c-----------------------------------------------------------------------
c     temperature equation
c-----------------------------------------------------------------------
      fx_u(4,1,:,:) = -(gamma-one)*u(2,:,:)*EXP(-u(1,:,:)/x)/x
      fx_u(4,2,:,:) = (gamma-one)*EXP(-u(1,:,:)/x)
      fx_u(4,4,:,:) = kappa*(gamma-one)/x
      fx_ux(4,4,:,:) = -kappa*(gamma-one)

      fy_u(4,1,:,:) = -(gamma-one)*u(3,:,:)*EXP(-u(1,:,:)/x)/x
      fy_u(4,3,:,:) = (gamma-one)*EXP(-u(1,:,:)/x)
      fy_uy(4,4,:,:) = -kappa*(gamma-one)/x**2

      s_u(4,1,:,:) = 
     $     EXP(-u(1,:,:)/x)/x**2
     $     *(u(2,:,:)*(ux(4,:,:) - u(4,:,:)/x) + u(3,:,:)*uy(4,:,:))
     $     - kappa*(gamma-one)*(ux(4,:,:) - u(4,:,:)/x)/x**2
     $     + mu*two*(gamma-one)*EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *((ux(2,:,:) - u(2,:,:)/x - ux(1,:,:)*u(2,:,:)/x
     $     + u(1,:,:)*u(2,:,:)/x**2)*u(2,:,:)/x**3
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/x 
     $     + u(1,:,:)*u(3,:,:)/x**2)*u(3,:,:)/x
     $     - (ux(2,:,:)/x - u(2,:,:)/x**2 - ux(1,:,:)*u(2,:,:)/x**2
     $     + u(1,:,:)*u(2,:,:)/x**3)**2
     $     - (uy(2,:,:)/x**2 - u(3,:,:)/x - uy(1,:,:)*u(2,:,:)/x**3)**2
     $     - (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/x 
     $     + u(1,:,:)*u(3,:,:)/x**2)**2
     $     - (uy(3,:,:)/x + u(2,:,:)/x**2 - uy(1,:,:)*u(3,:,:)/x**2)**2)
      s_u(4,2,:,:) = 
     $     -EXP(-u(1,:,:)/x)*(ux(4,:,:) - u(4,:,:)/x)/x
     $     + mu*two*(gamma-one)/x**2*EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *(-(ux(2,:,:) - u(2,:,:)/x - ux(1,:,:)*u(2,:,:)/x
     $     + u(1,:,:)*u(2,:,:)/x**2)*(one + ux(1,:,:) - u(1,:,:)/x)
     $     - (uy(2,:,:)/x - u(3,:,:) - uy(1,:,:)*u(2,:,:)/x**2)
     $     *uy(1,:,:)/x
     $     + (uy(3,:,:) + u(2,:,:)/x - uy(1,:,:)*u(3,:,:)/x))
      s_u(4,3,:,:) = 
     $     -EXP(-u(1,:,:)/x)*uy(4,:,:)/x
     $     - mu*two*(gamma-one)*EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *(uy(2,:,:)/x**2 - u(3,:,:)/x
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/x + u(1,:,:)*u(3,:,:)/x**2)
     $     *(ux(1,:,:) -  u(1,:,:)/x)
     $     + (uy(3,:,:) - uy(1,:,:)*u(3,:,:)/x)*uy(1,:,:)/x**2)
      s_u(4,4,:,:) = 
     $     EXP(-u(1,:,:)/x)*u(2,:,:)/x**2
     $     - kappa*(gamma-one)*(ux(1,:,:) - u(1,:,:)/x 
     $     + two*(ux(4,:,:) - u(4,:,:)/x))/x**2
     $     - mu*(gamma-one)*EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *((ux(2,:,:)/x - u(2,:,:)/x**2 - ux(1,:,:)*u(2,:,:)/x**2
     $     + u(1,:,:)*u(2,:,:)/x**3)**2
     $     + (uy(2,:,:)/x**2 - u(3,:,:)/x - uy(1,:,:)*u(2,:,:)/x**3)**2
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/x 
     $     + u(1,:,:)*u(3,:,:)/x**2)**2
     $     + (uy(3,:,:)/x + u(2,:,:)/x**2 - uy(1,:,:)*u(3,:,:)/x**2)**2)
      s_ux(4,1,:,:) = kappa*(gamma-one)*(ux(4,:,:) - u(4,:,:)/x)/x
     $     - mu*two*(gamma-one)*EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *((ux(2,:,:)/x - u(2,:,:)/x**2 - ux(1,:,:)*u(2,:,:)/x**2
     $     + u(1,:,:)*u(2,:,:)/x**3)*u(2,:,:)/x
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/x + u(1,:,:)*u(3,:,:)/x**2)
     $     *u(3,:,:))
      s_ux(4,2,:,:) = mu*two*(gamma-one)/x
     $     *EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *(ux(2,:,:) - u(2,:,:)/x - ux(1,:,:)*u(2,:,:)/x
     $     + u(1,:,:)*u(2,:,:)/x**2)
      s_ux(4,3,:,:) = mu*two*(gamma-one)*x
     $     *EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *(ux(3,:,:) - u(3,:,:)*ux(1,:,:)/x + u(1,:,:)*u(3,:,:)/x**2)
      s_ux(4,4,:,:) = 
     $     -EXP(-u(1,:,:)/x)*u(2,:,:)/x
     $     + kappa*(gamma-one)*(ux(1,:,:) - u(1,:,:)/x
     $     + two*(ux(4,:,:) - u(4,:,:)/x))/x
      s_uy(4,1,:,:) = kappa*(gamma-one)*uy(4,:,:)/x**3
     $     - mu*two*(gamma-one)*EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *((uy(2,:,:) - uy(1,:,:)*u(2,:,:)/x)*u(2,:,:)/x**2
     $     + (uy(3,:,:) - uy(1,:,:)*u(3,:,:)/x)*u(3,:,:))/x**2
      s_uy(4,2,:,:) = mu*two*(gamma-one)
     $     *EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *(uy(2,:,:)/x - u(3,:,:) - uy(1,:,:)*u(2,:,:)/x**2)/x**2
      s_uy(4,3,:,:) = mu*two*(gamma-one)
     $     *EXP(-(u(4,:,:) + two*u(1,:,:))/x)
     $     *(uy(3,:,:) + u(2,:,:)/x - uy(1,:,:)*u(3,:,:)/x)/x
      s_uy(4,4,:,:) = 
     $     -EXP(-u(1,:,:)/x)*u(3,:,:)/x
     $     + kappa*(gamma-one)*(uy(1,:,:) + two*uy(4,:,:))/x**3
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_drdu
c-----------------------------------------------------------------------
c     subprogram l. physics_mass.
c     computes mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_mass(x,y,mass,mass_x,mass_y)
      USE LnNavStok_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify mass matrix for azimuthal momentum.
c-----------------------------------------------------------------------
      mass(3,3,:,:) = x
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_mass
c-----------------------------------------------------------------------
c     subprogram m. physics_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_grid(ksi,eta,x,y)
      USE LnNavStok_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: ksi,eta
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: x,y
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("pellet1","pellet2")
         x = lx*ksi
         y = twopi*eta
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_grid
c-----------------------------------------------------------------------
c     subprogram n. physics_schur.
c     computes Schur complement.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_schur(t,hfac,x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
      USE LnNavStok_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t,hfac
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx_u=0
      fx_ux=0
      fx_uy=0
      fy_u=0
      fy_ux=0
      fy_uy=0
      s_u=0
      s_ux=0
      s_uy=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_schur
c-----------------------------------------------------------------------
c     subprogram o. physics_dealloc.
c     deallocate variables
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_dealloc
      USE LnNavStok_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_dealloc
c-----------------------------------------------------------------------
c     subprogram p. physics_main.
c     trivial main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM physics_main
      USE driver_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     do it.
c-----------------------------------------------------------------------
      CALL driver_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination.")
      END PROGRAM physics_main
