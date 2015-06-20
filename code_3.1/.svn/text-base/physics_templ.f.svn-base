c-----------------------------------------------------------------------
c     file physics_templ.f.
c     contains a template for the desired physics problem.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. physics_templ_mod.
c     1. physics_templ_equil.
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
c     subprogram 0. physics_templ_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE physics_templ_mod
      USE local_mod
      IMPLICIT NONE
      
      LOGICAL :: source=.FALSE.
      CHARACTER(16) :: init_type="."

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. physics_templ_equil.
c     computes initial equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_templ_equil(x,y,u,ux,uy,derivs)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0

      SELECT CASE(init_type)
      CASE(" ")

         u(1,:,:)=
         IF(.NOT. derivs)RETURN
c-----------------------------------------------------------------------
c     have to compute gradients of equilibrium quantities in order 
c     to use equilibrium source term option.
c-----------------------------------------------------------------------
         ux(1,:,:)=
         uy(1,:,:)=

      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_templ_equil
      END MODULE physics_templ_mod
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
      USE physics_templ_mod
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
      NAMELIST/template_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,du_diagnose,source,init_type
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=template_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="physics_templ"
c-----------------------------------------------------------------------
c     set number of dependent variables.
c-----------------------------------------------------------------------
      nqty=
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
c     sets various initialization parameters.
c     subroutine is called once per processor before dependent variables 
c     are initialized in physics_init 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_init_parameters(static,ground,adapt_qty)
      USE physics_templ_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set ground(nqty) and static(nqty) flags for each dependent 
c     variable; 
c     ground(iqty)=.TRUE. if u(iqty) is a potential that
c     needs to grounded when the system is doubly periodic
c     static(iqty)=.TRUE. if iqty-equation is a static equation
c     adapt_qty(iqty)=.FALSE. if iqty-variable should not be included
c     in griderr calculation for grid adaptation
c-----------------------------------------------------------------------
      ground()=.TRUE.
      static()=.TRUE.
      adapt_qty()=.FALSE.
c-----------------------------------------------------------------------
c     broadcast module variables defined in physics_input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE(" ")

      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_init_parameters
c-----------------------------------------------------------------------
c     subprogram c. physics_init.
c     initializes dependent variables on the computational grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_init(x,y,u)
      USE physics_templ_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     output variable u(nqty,:,:) has the first dimension corresponding
c     to the number of dependent variables and the last two index
c     over the 2D spatial grid
c-----------------------------------------------------------------------
c     calculate an initial equilibrium
c-----------------------------------------------------------------------
      CALL physics_templ_equil(x,y,u,ux,uy,.FALSE.)
c-----------------------------------------------------------------------
c     add perturbation to an equilibrium
c-----------------------------------------------------------------------

      SELECT CASE(init_type)
      CASE(" ")

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
      USE physics_templ_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c     set edge_type fields to appropriate values:
c     %bc_type(nqty) -- specify the type of b.c.
c     %static(nqty) -- set .TRUE. if left hand side of "robin" b.c. 
c                      equation is zero.
c
c     set the order for imposing boundary conditions on edges using:
c     edge_order
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE(" ")

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
      USE physics_templ_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     set c(nqty,:,:) to the RHS for "robin" b.c. 
c     or to normal outward flux for "normflux" b.c.
c-----------------------------------------------------------------------
      c=0

      SELECT CASE(lrtb)
      CASE("left")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("right")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("top")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("bottom")
         SELECT CASE(init_type)
         CASE(" ")

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
      USE physics_templ_mod
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

      SELECT CASE(lrtb)
      CASE("left")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("right")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("top")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("bottom")
         SELECT CASE(init_type)
         CASE(" ")

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
      USE physics_templ_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     set mass,mass_x,mass_y(iqty,jqty,:,:), ... to coupling 
c     mass matrices for du/dt, d(du/dx)/dt, and d(du/dy)/dt terms
c     of the general "robin" b.c. equations
c-----------------------------------------------------------------------
      mass=0
      mass_x=0
      mass_y=0

      SELECT CASE(lrtb)
      CASE("left")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("right")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("top")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("bottom")
         SELECT CASE(init_type)
         CASE(" ")

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
      USE physics_templ_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

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
      fx(1,:,:)=
      fy(1,:,:)=
      s(1,:,:)=

      fx(2,:,:)=
      fy(2,:,:)=
      s(2,:,:)=

c-----------------------------------------------------------------------
c     source acting to support the initial equilibrium is assumed to be
c     specified in physics_templ_equil 
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         IF(init_type = " ")THEN
            CALL physics_templ_equil(x,y,u0,u0x,u0y,.TRUE.)
            first=.FALSE.
            CALL physics_rhs(t,x,y,u0,u0x,u0y,fx0,fy0,s0,first)
            fx=fx-fx0
            fy=fy-fy0
            s=s-s0
         ELSE
            CALL program_stop("No source term speified in physics_rhs!")
         ENDIF
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
      USE physics_templ_mod
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
c     set fx_xx(iqty,jqty,:,:), fy_xx(iqty,jqty,:,:) and 
c     s_xx(iqty,jqty,:,:) to the derivatives of fx(iqty,:,:), 
c     fy(iqty,:,:), and s(iqty,:,:) given in physics_rhs 
c     with respect to u(jqty,:,:),ux(jqty,:,:),uy(jqty,:,:).
c-----------------------------------------------------------------------
      fx_u(1,1,:,:)=
      fx_u(1,2,:,:)=

      fx_ux(1,1,:,:)=
      fx_ux(1,2,:,:)=

      fx_uy(1,1,:,:)=
      fx_uy(1,2,:,:)=

      fy_u(1,1,:,:)=
      fy_u(1,2,:,:)=

      fy_ux(1,1,:,:)=
      fy_ux(1,2,:,:)=

      fy_uy(1,1,:,:)=
      fy_uy(1,2,:,:)=

      s_u(1,1,:,:)=
      s_u(1,2,:,:)=

      s_ux(1,1,:,:)=
      s_ux(1,2,:,:)=

      s_uy(1,1,:,:)=
      s_uy(1,2,:,:)=

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
      USE physics_templ_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify interior PDE's mass matrix
c-----------------------------------------------------------------------
      mass(1,1,:,:)=
      mass(2,2,:,:)=
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
      USE physics_templ_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: ksi,eta
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: x,y
c-----------------------------------------------------------------------
c     analytically initialize physical cordinates x,y in terms of 
c     logical coordinates ksi,eta
c-----------------------------------------------------------------------
      x = 
      y =
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
      USE physics_templ_mod
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
      USE physics_templ_mod
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
