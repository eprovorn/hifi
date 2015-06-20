c-----------------------------------------------------------------------
c     file euler1D.f.
c     contains specifications for pseudo1D Euler compressible flow.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. euler1D_mod.
c     1. euler1D_input.
c     2. euler1D_init.
c     3. euler1D_init_special.
c     4. euler1D_boundary.
c     5. euler1D_edge_rhs.
c     6. euler1D_edge_drdu.
c     7. euler1D_rhs.
c     8. euler1D_drdu.
c     9. euler1D_grid.
c-----------------------------------------------------------------------
c     subprogram 0. euler1D_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE euler1D_mod
      USE local_mod
      IMPLICIT NONE

      CHARACTER(16), PRIVATE :: init_type=" "
      REAL(r8), PARAMETER, PRIVATE :: gamma=1.4
      REAL(r8), PRIVATE :: mach=1.0,mu=0.0,pin=1.0,pout=1.0,lx=10.0

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. euler1D_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $     dt,dtmax,tmax,nstep)

      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic
      INTEGER, INTENT(OUT) :: nqty
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax

      INTEGER :: myios
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/euler1D_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,init_type,mach,mu,pin,pout,lx

c-----------------------------------------------------------------------
c     definitions.
c-----------------------------------------------------------------------
c
c     mach - input mach number
c     mu - viscosity added to the system
c     pin - pressure at the inflow
c     pout - pressure at the outflow
c

c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=euler1D_list,IOSTAT=myios)

      nqty=3
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler1D_input
c-----------------------------------------------------------------------
c     subprogram 2. euler1D_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_init(xpi,ypi,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui
      
      REAL(r8), DIMENSION(SIZE(xpi,1),SIZE(xpi,2)) :: area

c-----------------------------------------------------------------------
c     initial condition.
c-----------------------------------------------------------------------

      area = 1.398 + 0.347*TANH(0.8*xpi - 4.0)

      SELECT CASE(init_type)
      CASE("supersonic")
         ui(1,:,:) = 1.0*area
         ui(2,:,:) = 1.0*mach*SQRT(gamma*pin/1.0)*area
         ui(3,:,:) = (pin/(gamma-1) 
     $        + 0.5*1.0*(mach*SQRT(gamma*pin/1.0))**2)*area
      CASE("super_subsonic")
         ui(1,:,:) = 1.0*area
         ui(2,:,:) = 1.0*mach*SQRT(gamma*pin/1.0)*area
         ui(3,:,:) = (pin + (pout-pin)*xpi/lx)/(gamma-1)*area 
     $        + 0.5*1.0*ui(2,:,:)**2/area
      CASE("subsonic")
         ui(1,:,:) = 1.0*area
         ui(2,:,:) = 1.0*mach*SQRT(gamma*pin/1.0)*area
         ui(3,:,:) = (pin + (pout-pin)*xpi/lx)/(gamma-1)*area 
     $        + 0.5*1.0*ui(2,:,:)**2/area        
      END SELECT

c-----------------------------------------------------------------------
c     add perturbation
c-----------------------------------------------------------------------
      ! none.

      RETURN
      END SUBROUTINE euler1D_init
c-----------------------------------------------------------------------
c     subprogram 3. euler1D_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_init_special(static)

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static
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
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pout,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler1D_init_special
c-----------------------------------------------------------------------
c     subprogram 4. euler1D_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_boundary(left,right,top,bottom,nqty)

      INTEGER, INTENT(IN) :: nqty
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("supersonic")
         left%bc_type="robin"
         left%static=.TRUE.
         right%bc_type="natural"
         top%bc_type="zeroflux"
         bottom%bc_type="zeroflux"
      CASE("super_subsonic")
         left%bc_type(1)="robin"
         left%bc_type(2)="robin"
         left%bc_type(3)="robin"
         left%static(1:3)=.TRUE.
         right%bc_type(1)="natural"
         right%bc_type(2)="robin"
         right%bc_type(3)="robin"
         right%static(2:3)=.TRUE.
         top%bc_type="zeroflux"
         bottom%bc_type="zeroflux"
      CASE("subsonic")
         left%bc_type(1)="robin"
         left%bc_type(2)="robin"
         left%bc_type(3)="robin"
         left%static(1:3)=.TRUE.
         right%bc_type(1)="robin"
         right%bc_type(2)="robin"
         right%bc_type(3)="robin"
         right%static(1:3)=.TRUE.
         top%bc_type="zeroflux"
         bottom%bc_type="zeroflux"
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler1D_boundary
c-----------------------------------------------------------------------
c     subprogram 5. euler1D_edge_rhs.
c     computes rhs for edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_edge_rhs(edge_type,t,x,y,u,ux,uy,uxx,uyy,c)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy,uxx,uyy
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: c

      REAL(r8), DIMENSION(SIZE(x)) :: area
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
      area = 1.398 + 0.347*TANH(0.8*x - 4.0)

      SELECT CASE(init_type)
      CASE("supersonic")
         SELECT CASE(edge_type)
         CASE("left")
            c(1,:) = u(1,:) - 1.0*area
            c(2,:) = u(2,:) - 1.0*mach*SQRT(gamma*pin/1.0)*area
            c(3,:) = u(3,:) - (pin/(gamma-1) 
     $           + 0.5*1.0*(mach*SQRT(gamma*pin/1.0))**2)*area
         END SELECT
      CASE("super_subsonic")
         SELECT CASE(edge_type)
         CASE("left")
            c(1,:) = u(1,:) - 1.0*area
            c(2,:) = u(2,:) - 1.0*mach*SQRT(gamma*pin/1.0)*area
            c(3,:) = u(3,:) - (pin/(gamma-1) 
     $           + 0.5*1.0*(mach*SQRT(gamma*pin/1.0))**2)*area
         CASE("right")
            c(2,:) = uxx(2,:)/u(1,:)
     $           - 2*ux(2,:)*ux(1,:)/u(1,:)**2
     $           - u(2,:)*uxx(1,:)/u(1,:)**2
     $           + 2*u(2,:)*ux(1,:)**2/u(1,:)**3
            c(3,:) = (gamma-1)*(u(3,:) 
     $           - 0.5*u(2,:)*u(2,:)/u(1,:))/area - pout
         END SELECT
      CASE("subsonic")
         SELECT CASE(edge_type)
         CASE("left")
            c(1,:) = ux(1,:)
            c(2,:) = u(2,:) - 1.0*mach*SQRT(gamma*pin/1.0)*area
            c(3,:) = ux(3,:)
         CASE("right")
            c(1,:) = u(1,:) - 1.0*area
            c(2,:) = ux(2,:)
            c(3,:) = ux(3,:)
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler1D_edge_rhs
c-----------------------------------------------------------------------
c     subprogram 6. euler1D_edge_drdu.
c     computes drdu for edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_edge_drdu(edge_type,x,y,u,ux,uy,uxx,uyy,c_u,
     $                  c_ux,c_uy,c_uxx,c_uyy)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy,uxx,uyy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,c_uxx,
     $     c_uyy

      REAL(r8), DIMENSION(SIZE(x)) :: area
c-----------------------------------------------------------------------
c     zero out output.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
      area = 1.398 + 0.347*TANH(0.8*x - 4.0)

      SELECT CASE(init_type)
      CASE("supersonic")
         SELECT CASE(edge_type)
         CASE("left")
            c_u(1,1,:) = 1.0
            c_u(2,2,:) = 1.0
            c_u(3,3,:) = 1.0
         END SELECT
      CASE("super_subsonic")
         SELECT CASE(edge_type)
         CASE("left")
            c_u(1,1,:) = 1.0
            c_u(2,2,:) = 1.0
            c_u(3,3,:) = 1.0
         CASE("right")
            c_u(2,1,:) = -uxx(2,:)/u(1,:)**2 
     $           + 4*ux(2,:)*ux(1,:)/u(1,:)**3
     $           + 2*u(2,:)*uxx(1,:)/u(1,:)**3
     $           - 6*u(2,:)*ux(1,:)**2/u(1,:)**4
            c_ux(2,1,:) = -2*ux(2,:)/u(1,:)**2
     $           + 4*u(2,:)*ux(1,:)/u(1,:)**3
            c_uxx(2,1,:) = -u(2,:)/u(1,:)**2
            c_u(2,2,:) = -uxx(1,:)/u(1,:)**2
     $           + 2*ux(1,:)**2/u(1,:)**3
            c_ux(2,2,:) = -2*ux(1,:)/u(1,:)**2
            c_uxx(2,2,:) = 1/u(1,:)
            
            c_u(3,1,:) = (gamma-1)/area*0.5*(u(2,:)/u(1,:))**2
            c_u(3,2,:) =-(gamma-1)/area*(u(2,:)/u(1,:)) 
            c_u(3,3,:) = (gamma-1)/area
         END SELECT
      CASE("subsonic")
         SELECT CASE(edge_type)
         CASE("left")
            c_ux(1,1,:) = 1.0
            c_u(2,2,:) = 1.0
            c_ux(3,3,:) = 1.0
         CASE("right")
            c_u(1,1,:) = 1.0
            c_ux(2,2,:) = 1.0
            c_ux(3,3,:) = 1.0
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler1D_edge_drdu
c-----------------------------------------------------------------------
c     subprogram 7. euler1D_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_rhs(x,y,u,ux,uy,fx,fy,s)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: area,pressureA,dadx,
     $     dudx,mul

c-----------------------------------------------------------------------
c     initialize
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
      area=1.398 + 0.347*TANH(0.8*x - 4.0)
      pressureA=(gamma-1)*(u(3,:,:) - 0.5*u(2,:,:)*u(2,:,:)/u(1,:,:))
      dadx=0.2776/(COSH(0.8*x - 4.0))**2
      dudx= -u(2,:,:)*ux(1,:,:)/u(1,:,:)**2 + ux(2,:,:)/u(1,:,:)
      mul = mu
      WHERE (x>=lx-1e-6)
         mul = mu
      END WHERE
c-----------------------------------------------------------------------
c     flux and source terms
c-----------------------------------------------------------------------

      fx(1,:,:) = u(2,:,:)
      
      fx(2,:,:) = u(2,:,:)*u(2,:,:)/u(1,:,:) + pressureA - mul*dudx*area
      s(2,:,:) = (pressureA/area -mul*dudx)*dadx
      
      fx(3,:,:)=u(2,:,:)/u(1,:,:)*(u(3,:,:) + pressureA - mul*dudx*area)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler1D_rhs
c-----------------------------------------------------------------------
c     subprogram 8. euler1D_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_drdu(x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: area,pressureA,
     $     dadx,dudx,mul

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
      area=1.398+0.347*TANH(0.8*x-4.0)
      pressureA=(gamma-1)*(u(3,:,:)-0.5*u(2,:,:)*u(2,:,:)/u(1,:,:))
      dadx=0.2776/(COSH(0.8*x - 4.0))**2
      dudx=-u(2,:,:)*ux(1,:,:)/u(1,:,:)**2 + ux(2,:,:)/u(1,:,:)
      mul = mu
      WHERE (x>=lx-1e-6)
         mul = mu
      END WHERE
c-----------------------------------------------------------------------
c     pseudo1D Euler flux Jacobian
c-----------------------------------------------------------------------
      
      fx_u(1,2,:,:) = 1.0
      
      fx_u(2,1,:,:) = ((gamma-1)*0.5 - 1)*(u(2,:,:)/u(1,:,:))**2
     $     - mul*area*(2*u(2,:,:)*ux(1,:,:)/u(1,:,:)**3 
     $     - ux(2,:,:)/u(1,:,:)**2)
      fx_u(2,2,:,:) = (u(2,:,:)/u(1,:,:))*(3-gamma)
     $     + mul*area*(ux(1,:,:)/u(1,:,:)**2)
      fx_u(2,3,:,:) = gamma-1
      
      fx_ux(2,1,:,:) = mul*area*u(2,:,:)/u(1,:,:)**2
      fx_ux(2,2,:,:) = -mul*area/u(1,:,:)
      
      s_u(2,1,:,:) = ((gamma-1)*0.5*(u(2,:,:)/u(1,:,:))**2/area
     $     + mul*( -2*u(2,:,:)*ux(1,:,:)/u(1,:,:)**3 
     $     + ux(2,:,:)/u(1,:,:)**2 ))*dadx
      s_u(2,2,:,:) = -(gamma-1)*(u(2,:,:)/u(1,:,:))/area*dadx 
     $     + mul*ux(1,:,:)/u(1,:,:)**2*dadx
      s_u(2,3,:,:) = (gamma-1)/area*dadx
      
      s_ux(2,1,:,:) = mul*u(2,:,:)/u(1,:,:)**2*dadx
      s_ux(2,2,:,:) = -mul/u(1,:,:)*dadx
      
      fx_u(3,1,:,:) = (u(2,:,:)/u(1,:,:))*
     $     ((gamma-1)*0.5*(u(2,:,:)/u(1,:,:))**2
     $     - (u(3,:,:) + pressureA)/u(1,:,:)
     $     + mul*area*dudx/u(1,:,:) - mul*area*(-dudx/u(1,:,:) 
     $     + u(2,:,:)*ux(1,:,:)/u(1,:,:)**3))
      fx_u(3,2,:,:) = ((u(3,:,:) + pressureA)/u(1,:,:) 
     $     - (u(2,:,:)/u(1,:,:))**2*(gamma-1)) 
     $     + mul*area*u(2,:,:)/u(1,:,:)*ux(1,:,:)/u(1,:,:)**2
     $     - mul*area*dudx/u(1,:,:)
      fx_u(3,3,:,:) = (u(2,:,:)/u(1,:,:))*gamma
      
      fx_ux(3,1,:,:) = mul*area*u(2,:,:)**2/u(1,:,:)**3
      fx_ux(3,2,:,:) = -mul*area*u(2,:,:)/u(1,:,:)**2

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler1D_drdu
c-----------------------------------------------------------------------
c     subprogram 9. euler1D_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler1D_grid(x,y,ksi,eta)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("supersonic","super_subsonic","subsonic")
         ksi=x*lx
         eta=y
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler1D_grid
      END MODULE euler1D_mod
