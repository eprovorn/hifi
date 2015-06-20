c-----------------------------------------------------------------------
c     file NavStok.f.
c     contains specifications for 2D Navier-Stokes equations
c     in Cartesian coordinates.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. NavStok_mod.
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
c     subprogram 0. NavStok_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE NavStok_mod
      USE local_mod
      IMPLICIT NONE

      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,eps=1.e-10_r8
      CHARACTER(16) :: init_type="."
      INTEGER :: ntheta=0
      REAL(r8) :: mach=1.,tau=1.,mu=0.,kappa=0.,rho0=1.,rhomax=1.,
     $     p0=1.,pmax=1.,lx=1.0,alpha=0.,gr_curve=0.

      END MODULE NavStok_mod
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
      USE NavStok_mod
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
      physics_type="NavStok"

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
      USE NavStok_mod
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
      USE NavStok_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: rad,theta
c-----------------------------------------------------------------------
c     initial condition.
c-----------------------------------------------------------------------
      rad = SQRT(x**2 + y**2)
      theta=ATAN2(y,x)

      SELECT CASE(init_type)
      CASE("pellet1","pellet2")
         u(1,:,:) = rho0 + half*(rhomax-rho0)*(COS(rad*pi/lx)+1)
         u(2,:,:) = 0
         u(3,:,:) = 0
         u(4,:,:) = p0 + half*(pmax-p0)*(COS(rad*pi/lx)+1)
c-----------------------------------------------------------------------
c     add perturbation
c-----------------------------------------------------------------------
         u(1,:,:) = u(1,:,:)
     $        *(one + alpha*COS(ntheta*theta)*SIN(rad*pi/lx))
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
      USE NavStok_mod
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
         left%bc_type="polar"
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
      USE NavStok_mod
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
         CASE("right")
            c(1,:,:) = u(4,:,:) - p0/rho0*u(1,:,:)
            IF(t < tau)THEN
               c(2,:,:) = u(2,:,:)
     $              + mach*nhat(1,:,:)*(SIN(t*pi/tau))**2
               c(3,:,:) = u(3,:,:)
     $              + mach*nhat(2,:,:)*(SIN(t*pi/tau))**2
            ELSE
               c(2,:,:) = u(2,:,:)
               c(3,:,:) = u(3,:,:)
            ENDIF
            c(4,:,:) = ux(4,:,:)*nhat(1,:,:) + uy(4,:,:)*nhat(2,:,:)
         END SELECT
      CASE("pellet2")
         SELECT CASE(lrtb)
         CASE("right")
            c(1,:,:) = u(1,:,:) - rho0
            IF(t < tau)THEN
               c(2,:,:) = u(2,:,:)
     $              + mach*u(1,:,:)*nhat(1,:,:)*(SIN(t*pi/tau))**2
               c(3,:,:) = u(3,:,:)
     $              + mach*u(1,:,:)*nhat(2,:,:)*(SIN(t*pi/tau))**2
            ELSE
               c(2,:,:) = u(2,:,:)
               c(3,:,:) = u(3,:,:)
            ENDIF
            c(4,:,:) = ux(4,:,:)*nhat(1,:,:) + uy(4,:,:)*nhat(2,:,:)
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
      USE NavStok_mod
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
         CASE("right")
            c_u(1,1,:,:) = -p0/rho0
            c_u(1,4,:,:) = one
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_ux(4,4,:,:) = nhat(1,:,:)
            c_uy(4,4,:,:) = nhat(2,:,:)
         END SELECT
      CASE("pellet2")
         SELECT CASE(lrtb)
         CASE("right")
            c_u(1,1,:,:) = one
            IF(t < tau)THEN
               c_u(2,1,:,:) = mach*nhat(1,:,:)*(SIN(t*pi/tau))**2
               c_u(3,1,:,:) = mach*nhat(2,:,:)*(SIN(t*pi/tau))**2
            ENDIF
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_ux(4,4,:,:) = nhat(1,:,:)
            c_uy(4,4,:,:) = nhat(2,:,:)
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
      USE NavStok_mod
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
      USE NavStok_mod
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
c     flux and source terms
c-----------------------------------------------------------------------
      fx(1,:,:) = u(2,:,:)
      fy(1,:,:) = u(3,:,:)

      fx(2,:,:) = u(2,:,:)**2/u(1,:,:) + u(4,:,:)
     $     - mu*(ux(2,:,:) - u(2,:,:)*ux(1,:,:)/u(1,:,:))
      fy(2,:,:) = u(2,:,:)*u(3,:,:)/u(1,:,:)
     $     - mu*(uy(2,:,:) - u(2,:,:)*uy(1,:,:)/u(1,:,:))

      fx(3,:,:) = u(2,:,:)*u(3,:,:)/u(1,:,:)
     $     - mu*(ux(3,:,:) - u(3,:,:)*ux(1,:,:)/u(1,:,:))
      fy(3,:,:) = u(3,:,:)**2/u(1,:,:) + u(4,:,:)
     $     - mu*(uy(3,:,:) - u(3,:,:)*uy(1,:,:)/u(1,:,:))

      fx(4,:,:) = gamma*u(2,:,:)/u(1,:,:)*u(4,:,:)
     $     - kappa*(gamma-1)*(ux(4,:,:) - u(4,:,:)*ux(1,:,:)/u(1,:,:))
      fy(4,:,:) = gamma*u(3,:,:)/u(1,:,:)*u(4,:,:)
     $     - kappa*(gamma-1)*(uy(4,:,:) - u(4,:,:)*uy(1,:,:)/u(1,:,:))
      s(4,:,:) = (gamma-1)
     $     *((u(2,:,:)*ux(4,:,:) + u(3,:,:)*uy(4,:,:))/u(1,:,:)
     $     + mu*((ux(2,:,:) - u(2,:,:)*ux(1,:,:)/u(1,:,:))**2
     $     + (uy(2,:,:) - u(2,:,:)*uy(1,:,:)/u(1,:,:))**2
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/u(1,:,:))**2
     $     + (uy(3,:,:) - u(3,:,:)*uy(1,:,:)/u(1,:,:))**2)/u(1,:,:))
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
      USE NavStok_mod
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
c     2D Euler flux Jacobian
c-----------------------------------------------------------------------

      fx_u(1,2,:,:) = one
      fy_u(1,3,:,:) = one

      fx_u(2,1,:,:) = -u(2,:,:)**2/u(1,:,:)**2 
     $     - mu*u(2,:,:)*ux(1,:,:)/u(1,:,:)**2
      fx_u(2,2,:,:) = two*u(2,:,:)/u(1,:,:) + mu*ux(1,:,:)/u(1,:,:)
      fx_u(2,4,:,:) = one
      fx_ux(2,1,:,:) = mu*u(2,:,:)/u(1,:,:)
      fx_ux(2,2,:,:) = -mu

      fy_u(2,1,:,:) = -u(2,:,:)*u(3,:,:)/u(1,:,:)**2
     $     - mu*u(2,:,:)*uy(1,:,:)/u(1,:,:)**2
      fy_u(2,2,:,:) = u(3,:,:)/u(1,:,:) + mu*uy(1,:,:)/u(1,:,:)
      fy_u(2,3,:,:) = u(2,:,:)/u(1,:,:)
      fy_uy(2,1,:,:) = mu*u(2,:,:)/u(1,:,:)
      fy_uy(2,2,:,:) = -mu

      fx_u(3,1,:,:) = -u(2,:,:)*u(3,:,:)/u(1,:,:)**2
     $     - mu*u(3,:,:)*ux(1,:,:)/u(1,:,:)**2
      fx_u(3,2,:,:) = u(3,:,:)/u(1,:,:)
      fx_u(3,3,:,:) = u(2,:,:)/u(1,:,:) + mu*ux(1,:,:)/u(1,:,:)
      fx_ux(3,1,:,:) = mu*u(3,:,:)/u(1,:,:)
      fx_ux(3,3,:,:) = -mu

      fy_u(3,1,:,:) = -u(3,:,:)**2/u(1,:,:)**2 
     $     - mu*u(3,:,:)*uy(1,:,:)/u(1,:,:)**2
      fy_u(3,3,:,:) = two*u(3,:,:)/u(1,:,:) + mu*uy(1,:,:)/u(1,:,:)
      fy_u(3,4,:,:) = one
      fy_uy(3,1,:,:) = mu*u(3,:,:)/u(1,:,:)
      fy_uy(3,3,:,:) = -mu

      fx_u(4,1,:,:) = -gamma*u(2,:,:)*u(4,:,:)/u(1,:,:)**2
     $     - kappa*(gamma-1)*u(4,:,:)*ux(1,:,:)/u(1,:,:)**2
      fx_u(4,2,:,:) = gamma/u(1,:,:)*u(4,:,:)
      fx_u(4,4,:,:) = gamma*u(2,:,:)/u(1,:,:)
     $     + kappa*(gamma-1)*ux(1,:,:)/u(1,:,:)
      fx_ux(4,1,:,:) = kappa*(gamma-1)*u(4,:,:)/u(1,:,:)
      fx_ux(4,4,:,:) = -kappa*(gamma-1)

      fy_u(4,1,:,:) = -gamma*u(3,:,:)*u(4,:,:)/u(1,:,:)**2
     $     - kappa*(gamma-1)*u(4,:,:)*uy(1,:,:)/u(1,:,:)**2
      fy_u(4,3,:,:) = gamma/u(1,:,:)*u(4,:,:)
      fy_u(4,4,:,:) = gamma*u(3,:,:)/u(1,:,:)
     $     + kappa*(gamma-1)*uy(1,:,:)/u(1,:,:)
      fy_uy(4,1,:,:) = kappa*(gamma-1)*u(4,:,:)/u(1,:,:)
      fy_uy(4,4,:,:) = -kappa*(gamma-1)

      s_u(4,1,:,:) = (gamma-1)/u(1,:,:)**2
     $     *(-(u(2,:,:)*ux(4,:,:) + u(3,:,:)*uy(4,:,:))
     $     - mu*((ux(2,:,:) - u(2,:,:)*ux(1,:,:)/u(1,:,:))**2
     $     + (uy(2,:,:) - u(2,:,:)*uy(1,:,:)/u(1,:,:))**2
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/u(1,:,:))**2
     $     + (uy(3,:,:) - u(3,:,:)*uy(1,:,:)/u(1,:,:))**2)
     $     + two*mu/u(1,:,:)
     $     *((ux(2,:,:) - u(2,:,:)*ux(1,:,:)/u(1,:,:))
     $     *u(2,:,:)*ux(1,:,:)
     $     + (uy(2,:,:) - u(2,:,:)*uy(1,:,:)/u(1,:,:))
     $     *u(2,:,:)*uy(1,:,:)
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/u(1,:,:))
     $     *u(3,:,:)*ux(1,:,:)
     $     + (uy(3,:,:) - u(3,:,:)*uy(1,:,:)/u(1,:,:))
     $     *u(3,:,:)*uy(1,:,:)))
      s_u(4,2,:,:) = (gamma-1)*(ux(4,:,:)/u(1,:,:)
     $     - two*mu/u(1,:,:)**2
     $     *((ux(2,:,:) - u(2,:,:)*ux(1,:,:)/u(1,:,:))*ux(1,:,:)
     $     + (uy(2,:,:) - u(2,:,:)*uy(1,:,:)/u(1,:,:))*uy(1,:,:)))
      s_u(4,3,:,:) = (gamma-1)*(uy(4,:,:)/u(1,:,:)
     $     - two*mu/u(1,:,:)**2
     $     *((ux(3,:,:) - u(3,:,:)*ux(1,:,:)/u(1,:,:))*ux(1,:,:)
     $     + (uy(3,:,:) - u(3,:,:)*uy(1,:,:)/u(1,:,:))*uy(1,:,:)))

      s_ux(4,1,:,:) = -(gamma-1)*two*mu/u(1,:,:)**2
     $     *((ux(2,:,:) - u(2,:,:)*ux(1,:,:)/u(1,:,:))*u(2,:,:)
     $     + (ux(3,:,:) - u(3,:,:)*ux(1,:,:)/u(1,:,:))*u(3,:,:))
      s_ux(4,2,:,:) = (gamma-1)*two*mu/u(1,:,:)
     $     *(ux(2,:,:) - u(2,:,:)*ux(1,:,:)/u(1,:,:))
      s_ux(4,3,:,:) = (gamma-1)*two*mu/u(1,:,:)
     $     *(ux(3,:,:) - u(3,:,:)*ux(1,:,:)/u(1,:,:))
      s_ux(4,4,:,:) = (gamma-1)*u(2,:,:)/u(1,:,:)

      s_uy(4,1,:,:) = -(gamma-1)*two*mu/u(1,:,:)**2
     $     *((uy(2,:,:) - u(2,:,:)*uy(1,:,:)/u(1,:,:))*u(2,:,:)
     $     + (uy(3,:,:) - u(3,:,:)*uy(1,:,:)/u(1,:,:))*u(3,:,:))
      s_uy(4,2,:,:) = (gamma-1)*two*mu/u(1,:,:)
     $     *(uy(2,:,:) - u(2,:,:)*uy(1,:,:)/u(1,:,:))
      s_uy(4,3,:,:) = (gamma-1)*two*mu/u(1,:,:)
     $     *(uy(3,:,:) - u(3,:,:)*uy(1,:,:)/u(1,:,:))
      s_uy(4,4,:,:) = (gamma-1)*u(3,:,:)/u(1,:,:)
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
      USE NavStok_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     no modifications needed.
c-----------------------------------------------------------------------
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
      USE NavStok_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: ksi,eta
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: x,y
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("pellet1","pellet2")
         x = lx*ksi*COS(twopi*eta)
         y = lx*ksi*SIN(twopi*eta)
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
      USE NavStok_mod
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
      USE NavStok_mod
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
