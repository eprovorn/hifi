c-----------------------------------------------------------------------
c     file fourfieldplus.f.
c     contains specs for four-field model plus electron 
c     viscosity/inertia.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. fourfieldplus_mod.
c     1. fourfieldplus_equil.
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
c     p. physics_main.
c-----------------------------------------------------------------------
c     subprogram 0. fourfieldplus_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE fourfieldplus_mod
      USE local_mod
      IMPLICIT NONE
      
      LOGICAL :: source=.FALSE.
      CHARACTER(16) :: init_type="."
      REAL(r8), PARAMETER :: mass_r=5.44617e-4
      REAL(r8) :: di=0,eta=0,mu=0,nu=0,lx=0,ly=0,alfven=1,
     $     mach=0,lambda_psi=0,lambda_phi=0,epsilon=0,bound_eps=0,tau=0,
     $     kx=0,ky=0,de_sq=0,gr_curve=0.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fourfieldplus_equil.
c     computes initial equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfieldplus_equil(x,y,u,ux,uy,deriv)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: deriv

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: coshy_psi,coshy_phi,
     $     sinhy_psi,sinhy_phi,coshx_psi,coshx_phi,sinhx_psi,sinhx_phi
c-----------------------------------------------------------------------
c     compute cosh.
c-----------------------------------------------------------------------
      coshy_psi=COSH(y/lambda_psi)
      coshy_phi=COSH(y/lambda_phi)
      coshx_psi=COSH(x/lambda_psi)
      coshx_phi=COSH(x/lambda_phi)
c-----------------------------------------------------------------------
c        compute sinh.
c-----------------------------------------------------------------------
      sinhx_psi=SINH(x/lambda_psi)
      sinhx_phi=SINH(x/lambda_phi)
      sinhy_psi=SINH(y/lambda_psi)
      sinhy_phi=SINH(y/lambda_phi)
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0

      SELECT CASE(init_type)
      CASE("Alan")
         u(1,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(3,:,:)=mach/(lambda_phi*coshy_phi**2)
         u(5,:,:)=-di*alfven/(lambda_psi*coshy_psi**2)
         u(6,:,:)=mach*lambda_phi*log(coshy_phi)
         u(7,:,:)=u(3,:,:)
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=alfven*TANH(y/lambda_psi)
         uy(3,:,:)=-u(3,:,:)*(two*TANH(y/lambda_phi)/lambda_phi)
         uy(5,:,:)=-u(5,:,:)*(two*TANH(y/lambda_psi)/lambda_psi)
         uy(6,:,:)=mach*TANH(y/lambda_phi)
         uy(7,:,:)=uy(3,:,:)
      CASE("GEM")
         u(1,:,:)=-lambda_psi*LOG(coshy_psi)
         u(5,:,:)=di/(lambda_psi*coshy_psi**2)
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=-TANH(y/lambda_psi)
         uy(5,:,:)=-u(5,:,:)*(two*TANH(y/lambda_psi)/lambda_psi)
      CASE("Fitzpatrick","Bhimsen")
         u(1,:,:)=x**2*half
         u(5,:,:)=-di
      CASE("Harris_local")
         u(1,:,:)=lambda_psi*log(coshy_psi)
         u(5,:,:)=-di/(lambda_psi*coshy_psi**2)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=sinhy_psi/coshy_psi
         uy(5,:,:)=2*di/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfieldplus_equil
      END MODULE fourfieldplus_mod
c-----------------------------------------------------------------------
c     subprogram a. physics_input.
c     sets up input constants.
c     this subroutine is called only by the mpi_rank=0 processor
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $     nqty_schur,dt,dtmax,tmax,nstep,du_diagnose,physics_type,
     $     exit_flag)
      USE fourfieldplus_mod
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
      NAMELIST/fivefield_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,du_diagnose,gr_curve,di,eta,mu,nu,lx,ly,
     $     alfven,mach,lambda_psi,lambda_phi,epsilon,bound_eps,tau,
     $     source,init_type
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=fivefield_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="fourfieldplus"
c-----------------------------------------------------------------------
c     set number of dependent variables.
c-----------------------------------------------------------------------
      nqty=7
c-----------------------------------------------------------------------
c     set number of variables in the Schur complement
c     (set to 0, if not known)
c-----------------------------------------------------------------------
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
      USE fourfieldplus_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set ground(nqty) and static(nqty) flags for each dependent 
c     variable; 
c     ground(iqty)=.TRUE. if u(iqty) is a potential that
c     needs to grounded when the system is doubly periodic
c     static(iqty)=.TRUE. if iqty-equation is a static equation
c-----------------------------------------------------------------------
      ground(6)=.TRUE.
      static(5:7)=.TRUE.
c-----------------------------------------------------------------------
c     broadcast module variables defined in physics_input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_phi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bound_eps,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(tau,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alfven,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      de_sq=mass_r*di**2

      SELECT CASE(init_type)
      CASE("Fitzpatrick","Bhimsen")
         ky=twopi/lx
      CASE("Alan")
         lambda_psi=.2_r8
         lambda_phi=.2_r8
         ly=one
         lx=4._r8
         kx=twopi/lx
         ky=pi/ly
      CASE("GEM")
         lambda_psi=.5_r8
         lx=25.6_r8
         ly=12.8_r8
         di=one
         kx=twopi/lx
         ky=pi/ly
      END SELECT
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
      USE fourfieldplus_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8) :: ksq
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: cosx,cosy
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     output variable u(nqty,:,:) has the first dimension corresponding
c     to the number of dependent variables and the last two index
c     over the 2D spatial grid
c-----------------------------------------------------------------------
c     calculate an initial equilibrium
c-----------------------------------------------------------------------
      CALL fourfieldplus_equil(x,y,u,ux,uy,.FALSE.)
c-----------------------------------------------------------------------
c     add perturbation to an equilibrium
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Alan")
         cosx=COS(kx*x)
         cosy=COS(ky*y)
         ksq = kx**2 + ky**2
         u(1,:,:) = u(1,:,:)-epsilon*cosx*cosy*alfven/ksq
         u(6,:,:) = u(6,:,:)-epsilon*cosx*cosy*mach/ksq
      CASE("GEM")
         u(1,:,:) = u(1,:,:)+epsilon*COS(kx*x)*COS(ky*y)
      CASE("Harris_local")
         u(1,:,:) = u(1,:,:) + epsilon
     $        *EXP(-x**2/(2*lambda_psi)**2)
     $        *EXP(-y**2/(lambda_psi/2)**2)
      CASE("Fitzpatrick","Bhimsen")
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
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
      USE fourfieldplus_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c
c     set edge_type fields to appropriate values:
c     %bc_type(nqty) -- specify the type of b.c.
c     %static(nqty) -- set .TRUE. if left hand side of "robin" b.c. 
c                      equation is zero.
c
c     set the order for imposing boundary conditions on edges using:
c     edge_order
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Fitzpatrick")
         left%bc_type = "robin"
         right%bc_type = "robin"
         top%bc_type = "robin"
         bottom%bc_type = "robin"
         left%static=.FALSE.
         right%static=.FALSE.
         top%static=.FALSE.
         bottom%static=.FALSE.
      CASE("Bhimsen")
         left%bc_type = "robin"
         right%bc_type = "robin"
         top%bc_type = "robin"
         bottom%bc_type = "robin"

         left%static(1:3)=.FALSE.
         left%static(4:5)=.TRUE.
         left%static(6:7)=.FALSE.
         right%static(1:3)=.FALSE.
         right%static(4:5)=.TRUE.
         right%static(6:7)=.FALSE.
         top%static=.FALSE.
         bottom%static=.FALSE.
      CASE("Alan")
         top%bc_type = "robin"
         bottom%bc_type = "robin"

         top%static(1)=.FALSE.
         top%static(2:5)=.TRUE.
         top%static(6)=.FALSE.
         top%static(7)=.TRUE.
         bottom%static(1)=.FALSE.
         bottom%static(2:5)=.TRUE.
         bottom%static(6)=.FALSE.
         bottom%static(7)=.TRUE.
      CASE("GEM")
         top%bc_type = "robin"
         bottom%bc_type = "robin"

         top%static(1)=.FALSE.
         top%static(2:7)=.TRUE.
         bottom%static(1)=.FALSE.
         bottom%static(2:7)=.TRUE.
      CASE("Harris_local")
         left%bc_type = "robin"
         right%bc_type = "robin"
         top%bc_type = "robin"
         bottom%bc_type = "robin"

         left%static=.TRUE.
         right%static=.TRUE.
         bottom%static=.TRUE.

         top%static(1)=.FALSE.
         top%static(2:7)=.TRUE.
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_boundary
c-----------------------------------------------------------------------
c     subprogram e. physics_edge_rhs.
c     computes rhs for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_rhs(lrtb,t,x,y,nhat,u,ux,uy,uxx,uyy,uxy,c)
      USE fourfieldplus_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      REAL(r8) :: dksi_dt
c-----------------------------------------------------------------------
c     set c(nqty,:,:) to the RHS for "robin" b.c. 
c     or to normal outward flux for "normflux" b.c.
c-----------------------------------------------------------------------
      c=0

      SELECT CASE(init_type)
      CASE("Fitzpatrick")
         dksi_dt=bound_eps/tau**2*EXP(-t/tau)
         SELECT CASE(lrtb)
         CASE("left")
            c(1,:,:)=-dksi_dt*t*SIN(ky*y)
            c(6,:,:)=-1/ky*dksi_dt*(1-t/tau)*COS(ky*y)
         CASE("right")
            c(1,:,:)=-dksi_dt*t*SIN(ky*y)
            c(6,:,:)=1/ky*dksi_dt*(1-t/tau)*COS(ky*y)
         END SELECT
      CASE("Bhimsen")
         dksi_dt=bound_eps/tau**2*EXP(-t/tau)
         SELECT CASE(lrtb)
         CASE("left")
            c(1,:,:)=-dksi_dt*t*SIN(ky*y)
            c(4,:,:)=ux(4,:,:)
            c(5,:,:)=ux(5,:,:)
            c(6,:,:)=-1/ky*dksi_dt*(1-t/tau)*COS(ky*y)
         CASE("right")
            c(1,:,:)=-dksi_dt*t*SIN(ky*y)
            c(4,:,:)=ux(4,:,:)
            c(5,:,:)=ux(5,:,:)
            c(6,:,:)=1/ky*dksi_dt*(1-t/tau)*COS(ky*y)
         END SELECT
      CASE("Alan")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=uy(5,:,:)
            c(7,:,:)=u(7,:,:)
         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(2,:,:)=uy(2,:,:)
            c(3:7,:,:)=u(3:7,:,:)
         END SELECT
      CASE("Harris_local")
         SELECT CASE(lrtb)
         CASE("top")
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4:6,:,:)=uy(4:6,:,:)
            c(7,:,:)=u(7,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=uy(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=u(7,:,:)
         CASE("left","right")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=ux(4,:,:)
            c(5,:,:)=ux(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=u(7,:,:)
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_edge_rhs
c-----------------------------------------------------------------------
c     subprogram f. physics_edge_drdu.
c     computes drdu for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_drdu(lrtb,t,x,y,nhat,
     $     u,ux,uy,uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)
      USE fourfieldplus_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy
c-----------------------------------------------------------------------
c     set c_xxx(iqty,jqty), to the derivatives of c(iqty) 
c     given in physics_edge_rhs with respect to 
c     u(jqty),ux(jqty),uy(jqty),uxx(jqty),uyy(jqty),uxy(jqty).
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
      c_uxy=0

      SELECT CASE(init_type)
      CASE("Bhimsen")
         SELECT CASE(lrtb)
         CASE("left","right")
            c_ux(4,4,:,:)=one
            c_ux(5,5,:,:)=one
         END SELECT
      CASE("Alan")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_uy(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_uy(5,5,:,:)=one
            c_u(7,7,:,:)=one
         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_uy(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
            c_u(7,7,:,:)=one
         END SELECT
      CASE("Harris_local")
         SELECT CASE(lrtb)
         CASE("top")
            c_uy(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_uy(5,5,:,:)=one
            c_uy(6,6,:,:)=one
            c_u(7,7,:,:)=one
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_uy(5,5,:,:)=one
            c_u(6,6,:,:)=one
            c_u(7,7,:,:)=one
         CASE("left","right")
            c_ux(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_ux(4,4,:,:)=one
            c_ux(5,5,:,:)=one
            c_u(6,6,:,:)=one
            c_u(7,7,:,:)=one
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
      USE fourfieldplus_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y

      INTEGER :: iqty
c-----------------------------------------------------------------------
c     set mass,mass_x,mass_y(iqty,jqty,:,:), ... to coupling 
c     mass matrices for du/dt, d(du/dx)/dt, and d(du/dy)/dt terms
c     of the general "robin" b.c. equations
c-----------------------------------------------------------------------
      mass=0
      mass_x=0
      mass_y=0

      DO iqty=1,SIZE(mass,1)
         mass(iqty,iqty,:,:)=one
      ENDDO

      SELECT CASE(lrtb)
      CASE("top")
         SELECT CASE(init_type)
         CASE("Harris_local")
            mass(1,1,:,:) = 0
            mass_x(1,1,:,:) = nhat(1,:,:)
            mass_y(1,1,:,:) = nhat(2,:,:)
         END SELECT
      END SELECT
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
      RECURSIVE SUBROUTINE physics_rhs(t,x,y,u,ux,uy,fx,fy,s,first)
      USE fourfieldplus_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: psix,psiy,bzx,bzy,
     $     wxi,wyi,vzx,vzy
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     fluxes fx(nqty,:,:), fy(nqty,:,:) and sources s(nqty,:,:) for 
c     interior PDEs as an arbitrary function of time t, 
c     physical coordinates x(:,:),y(:,:), and dependent variables and
c     their gradients u(nqty,:,:), ux(nqty,:,:), uy(nqty,:,:)
c-----------------------------------------------------------------------
      psix = ux(1,:,:) + mass_r*di*ux(5,:,:)
      psiy = uy(1,:,:) + mass_r*di*uy(5,:,:)
      bzx = ux(2,:,:) - mass_r*di*ux(7,:,:)
      bzy = uy(2,:,:) - mass_r*di*uy(7,:,:)
      wxi=ux(3,:,:) + mass_r*ux(7,:,:)
      wyi=uy(3,:,:) + mass_r*uy(7,:,:)
      vzx=ux(4,:,:) + mass_r*ux(5,:,:)
      vzy=uy(4,:,:) + mass_r*uy(5,:,:)

      fx(1,:,:) = -eta*ux(1,:,:) - di*nu*ux(5,:,:)
      fy(1,:,:) = -eta*uy(1,:,:) - di*nu*uy(5,:,:)
      s(1,:,:) = psix*(uy(6,:,:) + di*uy(2,:,:)) 
     $     - psiy*(ux(6,:,:) + di*ux(2,:,:))

      fx(2,:,:) = -eta*ux(2,:,:) + di*nu*ux(7,:,:)
      fy(2,:,:) = -eta*uy(2,:,:) + di*nu*uy(7,:,:)
      s(2,:,:) = bzx*uy(6,:,:) - bzy*ux(6,:,:) 
     $     + uy(5,:,:)*ux(1,:,:) - ux(5,:,:)*uy(1,:,:)
     $     + de_sq*(uy(7,:,:)*ux(2,:,:) - ux(7,:,:)*uy(2,:,:))

      fx(3,:,:) = -(mu*ux(3,:,:) + nu*ux(7,:,:))
      fy(3,:,:) = -(mu*uy(3,:,:) + nu*uy(7,:,:))
      s(3,:,:) = wxi*uy(6,:,:) - wyi*ux(6,:,:) 
     $     + mass_r*di*(ux(7,:,:)*uy(2,:,:) - uy(7,:,:)*ux(2,:,:))
     $     + ((uy(4,:,:)-uy(5,:,:))*ux(1,:,:)
     $     - (ux(4,:,:)-ux(5,:,:))*uy(1,:,:))/MAX(di,min_eps)

      fx(4,:,:) = -(mu*ux(4,:,:) + nu*ux(5,:,:))
      fy(4,:,:) = -(mu*uy(4,:,:) + nu*uy(5,:,:))
      s(4,:,:) = vzx*uy(6,:,:) - vzy*ux(6,:,:) 
     $     + psix*uy(2,:,:) - psiy*ux(2,:,:)

      fx(5,:,:) = di*ux(1,:,:)
      fy(5,:,:) = di*uy(1,:,:)
      s(5,:,:) = u(4,:,:)-u(5,:,:)

      fx(6,:,:) = ux(6,:,:)
      fy(6,:,:) = uy(6,:,:)
      s(6,:,:) = u(3,:,:)

      fx(7,:,:) = di*ux(2,:,:)
      fy(7,:,:) = di*uy(2,:,:)
      s(7,:,:) = u(7,:,:)-u(3,:,:)
c-----------------------------------------------------------------------
c     source acting to support the initial equilibrium is assumed to be
c     specified in physics_templ_equil 
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("GEM","Alan","Harris_local")
            CALL fourfieldplus_equil(x,y,u0,u0x,u0y,.TRUE.)
            first=.FALSE.
            CALL physics_rhs(t,x,y,u0,u0x,u0y,fx0,fy0,s0,first)
            fx=fx-fx0
            fy=fy-fy0
            s=s-s0
         CASE DEFAULT
            CALL program_stop("No sources for init_type = "//init_type)
         END SELECT
      ENDIF
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
      SUBROUTINE physics_drdu(t,x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
      USE fourfieldplus_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: psix,psiy,bzx,bzy,
     $     wxi,wyi,vzx,vzy
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
c     set fx_xx(iqty,jqty,:,:), fy_xx(iqty,jqty,:,:) and 
c     s_xx(iqty,jqty,:,:) to the derivatives of fx(iqty,:,:), 
c     fy(iqty,:,:), and s(iqty,:,:) given in physics_rhs 
c     with respect to u(jqty,:,:),ux(jqty,:,:),uy(jqty,:,:).
c-----------------------------------------------------------------------
      psix = ux(1,:,:) + mass_r*di*ux(5,:,:)
      psiy = uy(1,:,:) + mass_r*di*uy(5,:,:)
      bzx = ux(2,:,:) - mass_r*di*ux(7,:,:)
      bzy = uy(2,:,:) - mass_r*di*uy(7,:,:)
      wxi = ux(3,:,:) + mass_r*ux(7,:,:)
      wyi = uy(3,:,:) + mass_r*uy(7,:,:)
      vzx = ux(4,:,:) + mass_r*ux(5,:,:)
      vzy = uy(4,:,:) + mass_r*uy(5,:,:)

      fx_ux(1,1,:,:)=-eta
      fx_ux(1,5,:,:)=-di*nu
      fy_uy(1,1,:,:)=-eta
      fy_uy(1,5,:,:)=-di*nu
      
      s_ux(1,1,:,:) = uy(6,:,:) + di*uy(2,:,:)
      s_ux(1,2,:,:) = -di*psiy
      s_ux(1,5,:,:) = mass_r*di*(uy(6,:,:) + di*uy(2,:,:))
      s_ux(1,6,:,:) = -psiy
      s_uy(1,1,:,:) = -(ux(6,:,:) + di*ux(2,:,:))
      s_uy(1,2,:,:) = di*psix
      s_uy(1,5,:,:) = -mass_r*di*(ux(6,:,:) + di*ux(2,:,:))
      s_uy(1,6,:,:) = psix

      fx_ux(2,2,:,:)=-eta
      fx_ux(2,7,:,:)=di*nu
      fy_uy(2,2,:,:)=-eta
      fy_uy(2,7,:,:)=di*nu

      s_ux(2,1,:,:) = uy(5,:,:)
      s_ux(2,2,:,:) = uy(6,:,:) + de_sq*uy(7,:,:)
      s_ux(2,5,:,:) = -uy(1,:,:)
      s_ux(2,6,:,:) = -bzy
      s_ux(2,7,:,:) = -mass_r*di*uy(6,:,:) - de_sq*uy(2,:,:)
      s_uy(2,1,:,:) = -ux(5,:,:)
      s_uy(2,2,:,:) = -ux(6,:,:) - de_sq*ux(7,:,:)
      s_uy(2,5,:,:) = ux(1,:,:)
      s_uy(2,6,:,:) = bzx
      s_uy(2,7,:,:) = mass_r*di*ux(6,:,:) + de_sq*ux(2,:,:)

      fx_ux(3,3,:,:)=-mu
      fx_ux(3,7,:,:)=-nu
      fy_uy(3,3,:,:)=-mu
      fy_uy(3,7,:,:)=-nu
      
      s_ux(3,1,:,:) = (uy(4,:,:)-uy(5,:,:))/MAX(di,min_eps)
      s_ux(3,2,:,:) = -mass_r*di*uy(7,:,:)
      s_ux(3,3,:,:) = uy(6,:,:)
      s_ux(3,4,:,:) = -uy(1,:,:)/MAX(di,min_eps)
      s_ux(3,5,:,:) = uy(1,:,:)/MAX(di,min_eps)
      s_ux(3,6,:,:) = -wyi
      s_ux(3,7,:,:) = mass_r*(uy(6,:,:) + di*uy(2,:,:)) 
      s_uy(3,1,:,:) = -(ux(4,:,:)-ux(5,:,:))/MAX(di,min_eps)
      s_uy(3,2,:,:) = mass_r*di*ux(7,:,:)
      s_uy(3,3,:,:) = -ux(6,:,:)
      s_uy(3,4,:,:) = ux(1,:,:)/MAX(di,min_eps)
      s_uy(3,5,:,:) = -ux(1,:,:)/MAX(di,min_eps)
      s_uy(3,6,:,:) = wxi
      s_uy(3,7,:,:) = -mass_r*(ux(6,:,:) + di*ux(2,:,:)) 

      fx_ux(4,4,:,:)=-mu
      fx_ux(4,5,:,:)=-nu
      fy_uy(4,4,:,:)=-mu
      fy_uy(4,5,:,:)=-nu

      s_ux(4,1,:,:) = uy(2,:,:)
      s_ux(4,2,:,:) = -psiy
      s_ux(4,4,:,:) = uy(6,:,:)
      s_ux(4,5,:,:) = mass_r*(uy(6,:,:) + di*uy(2,:,:))
      s_ux(4,6,:,:) = -vzy
      s_uy(4,1,:,:) = -ux(2,:,:)
      s_uy(4,2,:,:) = psix
      s_uy(4,4,:,:) = -ux(6,:,:)
      s_uy(4,5,:,:) = -mass_r*(ux(6,:,:) + di*ux(2,:,:))
      s_uy(4,6,:,:) = vzx
      
      fx_ux(5,1,:,:)=di
      fy_uy(5,1,:,:)=di
      s_u(5,4,:,:)=one
      s_u(5,5,:,:)=-one
 
      fx_ux(6,6,:,:)=one
      fy_uy(6,6,:,:)=one
      s_u(6,3,:,:)=one
      
      fx_ux(7,2,:,:)=di
      fy_uy(7,2,:,:)=di
      s_u(7,3,:,:)=-one
      s_u(7,7,:,:)=one
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
      USE fourfieldplus_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify interior PDE's mass matrix
c-----------------------------------------------------------------------
      mass(1,5,:,:) = mass_r*di
      mass(2,7,:,:) = -mass_r*di
      mass(3,7,:,:) = mass_r
      mass(4,5,:,:) = mass_r
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
      SUBROUTINE physics_grid(ksi,phi,x,y)
      USE fourfieldplus_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: ksi,phi
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: x,y
c-----------------------------------------------------------------------
c     analytically initialize physical cordinates x,y in terms of 
c     logical coordinates ksi,eta
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Fitzpatrick","Bhimsen")
         x = two*(ksi-half)
         y = lx*phi
      CASE("Alan","GEM")
         x = lx*(ksi-half)
         y = ly*(phi-half)
      CASE("Harris_local")
         x=lx*ksi*(gr_curve*ksi**3 + 1/gr_curve)/(gr_curve + 1/gr_curve)
         y=ly*(TANH(gr_curve*(phi-1))+TANH(gr_curve))/TANH(gr_curve)
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
      USE fourfieldplus_mod
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
c     deallocates private module objects allocated above.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_dealloc
      USE fourfieldplus_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     terminate
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
