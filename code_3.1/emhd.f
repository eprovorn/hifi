c-----------------------------------------------------------------------
c     file emhd.f.
c     contains specifications for electron MHD.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. emhd_mod.
c     1. emhd_equil.
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
c     subprogram 0. emhd_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE emhd_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL :: source=.FALSE.
      CHARACTER(16) :: init_type="."
      REAL(r8) :: gr_curve=0,di=0,eta=0,nu=0,lx=0,ly=0,alfven=1,
     $     lambda_psi=0,epsilon=0,kx=0,ky=0,de_sq=0,mass_r=5.44617e-4

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. emhd_equil.
c     computes equilibria for EMHD.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE emhd_equil(x,y,u,ux,uy,derivs)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: coshy_psi,
     $     sinhy_psi,coshx_psi,sinhx_psi
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      SELECT CASE(init_type)
      CASE("Harris")
         coshy_psi=COSH(y/lambda_psi)
         u(1,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(3,:,:)=alfven/(lambda_psi*coshy_psi**2)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhy_psi=SINH(y/lambda_psi)
         uy(1,:,:)=alfven*sinhy_psi/coshy_psi
         uy(3,:,:)=-2*alfven/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE("Harris_open")
         coshx_psi=COSH(x/(lambda_psi*lx/ly))
         coshy_psi=COSH(y/lambda_psi)
         u(1,:,:)=lambda_psi*(log(coshy_psi)-alfven*log(coshx_psi))
         u(3,:,:)=one/lambda_psi*
     $        (one/coshy_psi**2 - alfven/(lx/ly*coshx_psi)**2)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhx_psi=SINH(x/(lambda_psi*lx/ly))
         sinhy_psi=SINH(y/lambda_psi)
         ux(1,:,:)=-alfven*sinhx_psi/(lx/ly*coshx_psi)
         ux(3,:,:)=2*alfven/lambda_psi**2*sinhx_psi/(lx/ly*coshx_psi)**3
         uy(1,:,:)=sinhy_psi/coshy_psi
         uy(3,:,:)=-2/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE("Harris_local")
         coshy_psi=COSH(y/lambda_psi)
         u(1,:,:)=lambda_psi*log(coshy_psi)
         u(3,:,:)=one/(lambda_psi*coshy_psi**2)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhy_psi=SINH(y/lambda_psi)
         uy(1,:,:)=sinhy_psi/coshy_psi
         uy(3,:,:)=-2/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE("de_instability")
         coshx_psi=COSH((x-lx)/lambda_psi)
         u(1,:,:)=lambda_psi*LOG(coshx_psi) + x
         u(3,:,:)=one/(lambda_psi*coshx_psi**2)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhx_psi=SINH((x-lx)/lambda_psi)
         ux(1,:,:)=sinhx_psi/coshx_psi + one
         ux(3,:,:)=-two/lambda_psi**2*sinhx_psi/coshx_psi**3
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE emhd_equil
      END MODULE emhd_mod
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
      USE emhd_mod
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
      NAMELIST/emhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,gr_curve,di,eta,nu,lx,ly,alfven,lambda_psi,
     $     source,epsilon,init_type
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=emhd_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="emhd"

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
      USE emhd_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(3:4)=.TRUE.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alfven,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      de_sq=mass_r*di**2
      nu=nu*di**2
      SELECT CASE(init_type)
      CASE("Harris")
         kx=twopi/lx
         ky=pi/ly
      CASE("de_instability")
         ky=twopi/ly
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
      USE emhd_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     initial equilibrium.
c-----------------------------------------------------------------------
      CALL emhd_equil(x,y,u,ux,uy,.FALSE.)
c-----------------------------------------------------------------------
c     add perturbation
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Harris")
         u(1,:,:) = u(1,:,:)-epsilon*COS(kx*x)*COS(ky*y)
      CASE("Harris_local")
         u(1,:,:) = u(1,:,:) + epsilon
     $        *EXP(-x**2/(two*lambda_psi)**2)
     $        *EXP(-y**2/(half*lambda_psi)**2)
      CASE("de_instability")
         u(1,:,:) = u(1,:,:)
     $        + epsilon*lambda_psi*EXP(-(x-lx)**2/lambda_psi**2)
     $        * COS(ky*y)
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
      USE emhd_mod
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
      CASE("Harris")
         top%static(2:4)=.TRUE.
         bottom%static(2:4)=.TRUE.
      CASE("Harris_open")
         left%static=.TRUE.
         right%static=.TRUE.
         top%static=.TRUE.
         bottom%static=.TRUE.
      CASE("Harris_local")
         left%static=.TRUE.
         right%static=.TRUE.
         top%static(2:4)=.TRUE.
         bottom%static=.TRUE.
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
      USE emhd_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
      SELECT CASE(init_type)
      CASE("Harris")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(2,:,:)
         END SELECT
      CASE("Harris_open")
         SELECT CASE(lrtb)
         CASE("top")
            c(1,:,:)=uyy(1,:,:)
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=uyy(3,:,:)
            c(4,:,:)=uy(4,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=u(4,:,:)
         CASE("right")
            c(1,:,:)=uxx(1,:,:)
            c(2,:,:)=ux(2,:,:)
            c(3,:,:)=uxx(3,:,:)
            c(4,:,:)=ux(4,:,:)
         CASE("left")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=ux(3,:,:)
            c(4,:,:)=u(4,:,:)
         END SELECT
      CASE("Harris_local")
         SELECT CASE(lrtb)
         CASE("top")
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=uyy(3,:,:)
            c(4,:,:)=u(4,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=u(4,:,:)
         CASE("left","right")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=ux(3,:,:)
            c(4,:,:)=u(4,:,:)
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
      SUBROUTINE physics_edge_drdu(lrtb,t,x,y,nhat,u,ux,uy,uxx,uyy,uxy,
     $     c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)
      USE emhd_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy
c-----------------------------------------------------------------------
c     zero out output.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
      c_uxy=0

      SELECT CASE(init_type)
      CASE("Harris")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_u(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uy(4,2,:,:)=one
         END SELECT
      CASE("Harris_open")
         SELECT CASE(lrtb)
         CASE("top")
            c_uyy(1,1,:,:)=one
            c_uy(2,2,:,:)=one
            c_uyy(3,3,:,:)=one
            c_uy(4,4,:,:)=one
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_u(4,4,:,:)=one
         CASE("right")
            c_uxx(1,1,:,:)=one
            c_ux(2,2,:,:)=one
            c_uxx(3,3,:,:)=one
            c_ux(4,4,:,:)=one
         CASE("left")
            c_ux(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_ux(3,3,:,:)=one
            c_u(4,4,:,:)=one
         END SELECT
      CASE("Harris_local")
         SELECT CASE(lrtb)
         CASE("top")
            c_uy(2,2,:,:)=one
            c_uyy(3,3,:,:)=one
            c_u(4,4,:,:)=one
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_u(4,4,:,:)=one
         CASE("left","right")
            c_ux(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_ux(3,3,:,:)=one
            c_u(4,4,:,:)=one
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
      USE emhd_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y

      INTEGER :: iqty
c-----------------------------------------------------------------------
c     set mass,mass_x,mass_y(iqty,jqty,:,:), ... to coupling 
c     mass matrices for du/dt, d(du/dx)/dt, and d(du/dy)/dt
c-----------------------------------------------------------------------
      mass=0
      mass_x=0
      mass_y=0

      SELECT CASE(init_type)
      CASE("Harris")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            mass(1,1,:,:)=one
         END SELECT
      CASE("Harris_local")
         SELECT CASE(lrtb)
         CASE("top")
            mass_x(1,1,:,:)=nhat(1,:,:)
            mass_y(1,1,:,:)=nhat(2,:,:)
         END SELECT
      CASE("de_instability")
         SELECT CASE(lrtb)
         CASE("left","right")
            DO iqty=1,SIZE(mass,1)
               mass(iqty,iqty,:,:)=one
            ENDDO
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
      USE emhd_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0y,u0x,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
      u0=0
      u0x=0
      u0y=0
c-----------------------------------------------------------------------
c     EMHD
c-----------------------------------------------------------------------
      fx(1,:,:)=-di*uy(2,:,:)*(u(1,:,:) - de_sq*u(3,:,:)) + nu*ux(3,:,:)
      
      fy(1,:,:)=di*ux(2,:,:)*(u(1,:,:) - de_sq*u(3,:,:)) + nu*uy(3,:,:)
      
      s(1,:,:)=eta*u(3,:,:)

      fx(2,:,:)=-di*(uy(1,:,:)*u(3,:,:) - de_sq*uy(2,:,:)*u(4,:,:))
     $     + nu*ux(4,:,:)
      
      fy(2,:,:)=di*(ux(1,:,:)*u(3,:,:) - de_sq*ux(2,:,:)*u(4,:,:))
     $     + nu*uy(4,:,:)

      s(2,:,:)=eta*u(4,:,:)
      
      fx(3,:,:)=ux(1,:,:)
      fy(3,:,:)=uy(1,:,:)
      s(3,:,:)=u(3,:,:)
      
      fx(4,:,:)=ux(2,:,:)
      fy(4,:,:)=uy(2,:,:)
      s(4,:,:)=u(4,:,:)

      IF(source .AND. first)THEN
         CALL emhd_equil(x,y,u0,u0x,u0y,.TRUE.)
         first=.FALSE.
         CALL physics_rhs(t,x,y,u0,u0x,u0y,fx0,fy0,s0,first)
         fx=fx-fx0
         fy=fy-fy0
         s=s-s0
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
      USE emhd_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
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
c     EMHD.
c-----------------------------------------------------------------------
      fx_u(1,1,:,:)=-di*uy(2,:,:)
      fx_u(1,3,:,:)=di*de_sq*uy(2,:,:)
      fx_ux(1,3,:,:)=nu
      fx_uy(1,2,:,:)=-di*(u(1,:,:) - de_sq*u(3,:,:))

      fy_u(1,1,:,:)=di*ux(2,:,:)
      fy_u(1,3,:,:)=-di*de_sq*ux(2,:,:)
      fy_ux(1,2,:,:)=di*(u(1,:,:) - de_sq*u(3,:,:))
      fy_uy(1,3,:,:)=nu

      s_u(1,3,:,:)=eta

      fx_u(2,3,:,:)=-di*uy(1,:,:)
      fx_u(2,4,:,:)=di*de_sq*uy(2,:,:)
      fx_ux(2,4,:,:)=nu
      fx_uy(2,1,:,:)=-di*u(3,:,:)
      fx_uy(2,2,:,:)=di*de_sq*u(4,:,:)

      fy_u(2,3,:,:)=di*ux(1,:,:)
      fy_u(2,4,:,:)=-di*de_sq*ux(2,:,:)
      fy_ux(2,1,:,:)=di*u(3,:,:)
      fy_ux(2,2,:,:)=-di*de_sq*u(4,:,:)
      fy_uy(2,4,:,:)=nu

      s_u(2,4,:,:)=eta
      
      fx_ux(3,1,:,:)=one
      fy_uy(3,1,:,:)=one
      s_u(3,3,:,:)=one
      
      fx_ux(4,2,:,:)=one
      fy_uy(4,2,:,:)=one
      s_u(4,4,:,:)=one
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_drdu
c-----------------------------------------------------------------------
c     subprogram l. physics_mass.
c     computes mass matrix couplings.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_mass(x,y,mass,mass_x,mass_y)
      USE emhd_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      mass(1,3,:,:)=-de_sq
      mass(2,4,:,:)=-de_sq
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
      SUBROUTINE physics_grid(x,y,ksi,phi)
      USE emhd_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,phi
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Harris")
         ksi=lx*(x-.5)
         phi=ly*(y-.5)*((y-.5)**2+.125)/.375
      CASE("Harris_open")
         ksi=lx*(TANH(3.*x-1.5)+TANH(1.5))/(2.*TANH(1.5))
         phi=ly*(TANH(4.*y-2.5)+TANH(2.5))/(TANH(1.5)+TANH(2.5))
      CASE("Harris_local")
         ksi=lx*x*(gr_curve*x**3 + 1./gr_curve)/(gr_curve + 1./gr_curve)
ccc         ksi=lx*(TANH(gr_curve*(x-1.))+TANH(gr_curve))/TANH(gr_curve)
         phi=ly*(TANH(gr_curve*(y-1.))+TANH(gr_curve))/TANH(gr_curve)
      CASE("de_instability")
         ksi=two*lx*x
         phi=ly*y
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_grid
c-----------------------------------------------------------------------
c     subprogram n. physics_schur.
c     computes schur complement.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_schur(t,hfac,x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
      USE emhd_mod
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
      USE emhd_mod
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
