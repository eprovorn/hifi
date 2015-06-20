c-----------------------------------------------------------------------
c     file rmhd.f.
c     contains specifications for reduced MHD + electron inertia
c     + sound radius + hyper-resistivity model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. rmhd_mod.
c     1. rmhd_equil.
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
c     subprogram 0. rmhd_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE rmhd_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL :: source=.FALSE.
      CHARACTER(16) :: init_type=".",eta_type="."
      REAL(r8) :: eta=0,mu=0,nu=0,de_sq=0,sound_rad=0,
     $     lambda_psi=0,lambda_phi=0,epsilon=0,mach=0,alfven=1,lx=0,
     $     ly=0,kx=0,ky=0,ksq=0,gr_curve=1.,j_crit=1.e20

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. rmhd_equil.
c     computes equilibrium for reduced MHD.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE rmhd_equil(x,y,u,ux,uy,derivs)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: coshy_psi,coshy_phi,
     $     sinhy_psi,sinhy_phi,coshx_psi,coshx_phi,sinhx_psi,
     $     cosx_psi,sinx_psi
c-----------------------------------------------------------------------
c     compute cosh.
c-----------------------------------------------------------------------
      coshy_psi=COSH(y/lambda_psi)
      coshy_phi=COSH(y/lambda_phi)
      coshx_psi=COSH(x/lambda_psi)
      coshx_phi=COSH(x/lambda_phi)
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      SELECT CASE(init_type)
      CASE("Harris")
         u(1,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(2,:,:)=mach/(lambda_phi*coshy_phi**2)
         u(3,:,:)=alfven/(lambda_psi*coshy_psi**2)
         u(4,:,:)=mach*lambda_phi*log(coshy_phi)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhy_psi=SINH(y/lambda_psi)
         sinhy_phi=SINH(y/lambda_phi)
         uy(1,:,:)=alfven*sinhy_psi/coshy_psi
         uy(2,:,:)=-2*mach/lambda_phi**2*sinhy_phi/coshy_phi**3
         uy(3,:,:)=-2*alfven/lambda_psi**2*sinhy_psi/coshy_psi**3
         uy(4,:,:)=mach*sinhy_phi/coshy_phi
      CASE("Harris_open")
         coshx_psi=COSH(x/(lambda_psi*lx/ly))
         u(1,:,:)=lambda_psi*(log(coshy_psi)-alfven*log(coshx_psi))
         u(3,:,:)=one/lambda_psi*
     $        (one/coshy_psi**2 - alfven/(lx/ly*coshx_psi)**2)
c-----------------------------------------------------------------------
c           compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhx_psi=SINH(x/(lambda_psi*lx/ly))
         sinhy_psi=SINH(y/lambda_psi)
         ux(1,:,:)=-alfven*sinhx_psi/(lx/ly*coshx_psi)
         ux(3,:,:)=2*alfven/lambda_psi**2*sinhx_psi/(lx/ly*coshx_psi)**3
         uy(1,:,:)=sinhy_psi/coshy_psi
         uy(3,:,:)=-2/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE("Harris_local","Harris_res","anom_res","anom_hres")
         u(1,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(3,:,:)=alfven/(lambda_psi*coshy_psi**2)
c-----------------------------------------------------------------------
c           compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhy_psi=SINH(y/lambda_psi)
         uy(1,:,:)=alfven*sinhy_psi/coshy_psi
         uy(3,:,:)=-2*alfven/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE("plasmoids")
         CALL RANDOM_NUMBER(u)
         u = 1.e-16*u

         u(1,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(3,:,:)=alfven/(lambda_psi*coshy_psi**2)
c-----------------------------------------------------------------------
c           compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhy_psi=SINH(y/lambda_psi)
         uy(1,:,:)=sinhy_psi/coshy_psi
         uy(3,:,:)=-two/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE("islands")
         cosx_psi=COS(x/lambda_psi)
         u(1,:,:)=-lambda_psi*LOG(coshy_psi+alfven*cosx_psi)
         u(3,:,:)=(alfven**2-one)/lambda_psi
     $        /(coshy_psi+alfven*cosx_psi)**2
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinx_psi=SIN(x/lambda_psi)
         sinhy_psi=SINH(y/lambda_psi)
         ux(1,:,:)=alfven*sinx_psi/(coshy_psi+alfven*cosx_psi)
         uy(1,:,:)=-sinhy_psi/(coshy_psi+alfven*cosx_psi)
         ux(3,:,:)=u(3,:,:)*two*alfven/lambda_psi*sinx_psi
     $        /(coshy_psi+alfven*cosx_psi)
         uy(3,:,:)=-u(3,:,:)*two/lambda_psi*sinhy_psi
     $        /(coshy_psi+alfven*cosx_psi)
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rmhd_equil
      END MODULE rmhd_mod
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
      USE rmhd_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: physics_type
      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty,nqty_schur
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep,exit_flag
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax

      INTEGER :: myios
c-----------------------------------------------------------------------
c     physics namelist.
c-----------------------------------------------------------------------
      NAMELIST/rmhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,gr_curve,eta,mu,nu,de_sq,sound_rad,source,lx,ly,
     $     lambda_psi,lambda_phi,epsilon,mach,alfven,init_type,eta_type,
     $     j_crit
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=rmhd_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="rmhd"

      nqty=4
      nqty_schur=0

      SELECT CASE(init_type)
      CASE("Harris")
         xperiodic=.TRUE.
      END SELECT
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
      USE rmhd_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      ground(4)=.TRUE.
      static(3:4)=.TRUE.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(eta_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(de_sq,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(sound_rad,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_phi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alfven,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(j_crit,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("islands")
         lambda_psi=one/twopi
         kx=pi/lx
         ky=half*pi
      CASE("Harris")
         kx=twopi/lx
         ky=pi/ly
      CASE("Harris_local","Harris_res","anom_res","anom_hres",
     $        "plasmoids")
         kx=pi/lx
         ky=half*pi/ly
      END SELECT
      ksq=kx**2+ky**2
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
      SUBROUTINE physics_init(xpi,ypi,ui)
      USE rmhd_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: ix,nmodes=1
      REAL(r8), DIMENSION(SIZE(ui,1),0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) ::
     $     ux,uy
      REAL(r8), DIMENSION(0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) :: cosx,sinx,
     $     cosy,siny
c-----------------------------------------------------------------------
c     initial equilibrium.
c-----------------------------------------------------------------------
      CALL rmhd_equil(xpi,ypi,ui,ux,uy,.FALSE.)
      cosx=COS(kx*xpi)
      sinx=SIN(kx*xpi)
      cosy=COS(ky*ypi)
      siny=SIN(ky*ypi)
c-----------------------------------------------------------------------
c     add perturbation
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("islands")
         ui(1,:,:)=ui(1,:,:)+epsilon*cosx*cosy
      CASE("Harris")
         ui(1,:,:)=ui(1,:,:)+epsilon*cosx*cosy*alfven/ksq
         ui(2,:,:)=ui(2,:,:)-epsilon*cosx*cosy*mach
         ui(3,:,:)=ui(3,:,:)-epsilon*cosx*cosy*alfven
         ui(4,:,:)=ui(4,:,:)+epsilon*cosx*cosy*mach/ksq
      CASE("Harris_local","anom_res","anom_hres","plasmoids")
         ui(1,:,:)=ui(1,:,:) + epsilon
     $        *EXP(-xpi**2/(2*lambda_psi)**2)
     $        *EXP(-ypi**2/(lambda_psi/2)**2)
      CASE("Harris_res")
         ui(1,:,:)=ui(1,:,:) + epsilon
     $        *EXP(-xpi**2/(2*lambda_psi)**2)*COS(ky*ypi)
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
      USE rmhd_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("islands")
         left%bc_type="robin"
         left%static=.FALSE.
         right%bc_type="robin"
         right%static=.FALSE.
         top%bc_type="robin"
         top%static=.FALSE.
         bottom%bc_type="robin"
         bottom%static=.FALSE.
      CASE("Harris")
         top%bc_type="robin"
         top%static(1)=.FALSE.
         top%static(2:4)=.TRUE.
         bottom%bc_type="robin"
         bottom%static(1)=.FALSE.
         bottom%static(2:4)=.TRUE.
      CASE("Harris_open")
         left%bc_type="robin"
         left%static=.TRUE.
         right%bc_type="robin"
         right%static(1)=.FALSE.
         right%static(2:4)=.TRUE.
         top%bc_type="robin"
         top%static(1)=.FALSE.
         top%static(2:4)=.TRUE.
         bottom%bc_type="robin"
         bottom%static=.TRUE.
      CASE("Harris_local","anom_res","anom_hres","plasmoids")
         edge_order=(/2,4,1,3/)
         left%bc_type="robin"
         left%static=.TRUE.
         left%bc_type(1)="zeroflux"
         left%bc_type(3)="zeroflux"

         right%bc_type="robin"
         right%static=.TRUE.
         right%bc_type(1)="zeroflux"
         right%bc_type(3)="zeroflux"

         bottom%bc_type="robin"
         bottom%static=.TRUE.
         bottom%bc_type(1)="zeroflux"
         bottom%bc_type(3)="zeroflux"

         top%bc_type="robin"
         top%static=.TRUE.
      CASE("Harris_res")
         left%bc_type="robin"
         left%static=.TRUE.
         right%bc_type="robin"
         right%static=.TRUE.
         top%bc_type="robin"
         top%static(1)=.FALSE.
         top%static(2:4)=.TRUE.
         bottom%bc_type="robin"
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
      USE rmhd_mod
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
            c(4,:,:)=u(4,:,:)
         END SELECT
      CASE("Harris_open")
         SELECT CASE(lrtb)
         CASE("top")
            c(2,:,:)=u(2,:,:)-uxx(4,:,:)
            c(3,:,:)=u(3,:,:)-uxx(1,:,:)-uyy(1,:,:)
            c(4,:,:)=uy(4,:,:)
         CASE("bottom","left")
            c(1,:,:)=ux(1,:,:)*nhat(1,:,:) + uy(1,:,:)*nhat(2,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=ux(3,:,:)*nhat(1,:,:) + uy(3,:,:)*nhat(2,:,:)
            c(4,:,:)=u(4,:,:)
         CASE("right")
            c(2,:,:)=u(2,:,:)-uyy(4,:,:)
            c(3,:,:)=u(3,:,:)-uxx(1,:,:)-uyy(1,:,:)
            c(4,:,:)=ux(4,:,:)
         END SELECT
      CASE("Harris_local","anom_res","anom_hres","plasmoids")
         SELECT CASE(lrtb)
         CASE("top")
            c(1,:,:)=uy(1,:,:)-alfven*TANH(y/lambda_psi)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
         CASE("left","right","bottom")
            c(2,:,:)=u(2,:,:)
            c(4,:,:)=u(4,:,:)
         END SELECT
      CASE("Harris_res")
         SELECT CASE(lrtb)
         CASE("top")
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=u(3,:,:) - (uxx(1,:,:) + uyy(1,:,:))
            c(4,:,:)=u(4,:,:)
         CASE("left","right","bottom")
            c(1,:,:)=ux(1,:,:)*nhat(1,:,:) + uy(1,:,:)*nhat(2,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=ux(3,:,:)*nhat(1,:,:) + uy(3,:,:)*nhat(2,:,:)
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
      USE rmhd_mod
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
            c_u(4,4,:,:)=one
         END SELECT
      CASE("Harris_open")
         SELECT CASE(lrtb)
         CASE("top")
            c_u(2,2,:,:)=one
            c_uxx(2,4,:,:)=-one
            c_u(3,3,:,:)=one
            c_uxx(3,1,:,:)=-one
            c_uyy(3,1,:,:)=-one
            c_uy(4,4,:,:)=one
         CASE("bottom","left")
            c_ux(1,1,:,:)=nhat(1,:,:)
            c_uy(1,1,:,:)=nhat(2,:,:)
            c_u(2,2,:,:)=one
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)
            c_u(4,4,:,:)=one
         CASE("right")
            c_u(2,2,:,:)=one
            c_uyy(2,4,:,:)=-one
            c_u(3,3,:,:)=one
            c_uxx(3,1,:,:)=-one
            c_uyy(3,1,:,:)=-one
            c_ux(4,4,:,:)=one
         END SELECT
      CASE("Harris_local","anom_res","anom_hres","plasmoids")
         SELECT CASE(lrtb)
         CASE("top")
            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uy(4,4,:,:)=one
         CASE("left","right","bottom")
            c_u(2,2,:,:)=one
            c_u(4,4,:,:)=one
         END SELECT
      CASE("Harris_res")
         SELECT CASE(lrtb)
         CASE("top")
            c_uy(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uxx(3,1,:,:)=-one
            c_uyy(3,1,:,:)=-one
            c_u(4,4,:,:)=one
         CASE("left","right","bottom")
            c_ux(1,1,:,:)=nhat(1,:,:)
            c_uy(1,1,:,:)=nhat(2,:,:)
            c_u(2,2,:,:)=one
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)
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
      USE rmhd_mod
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

      SELECT CASE(init_type)
      CASE("islands")
         SELECT CASE(lrtb)
         CASE("top")
            mass(1:4,1:4,:,:)=one
         CASE("left","right","bottom")
            mass_x(1,1,:,:)=nhat(1,:,:)
            mass_y(1,1,:,:)=nhat(2,:,:)
            mass(2,2,:,:)=one
            mass_x(3,3,:,:)=nhat(1,:,:)
            mass_y(3,3,:,:)=nhat(2,:,:)
            mass(4,4,:,:)=one
         END SELECT
      CASE("Harris")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            mass(1,1,:,:)=one
         END SELECT
      CASE("Harris_open")
         SELECT CASE(lrtb)
         CASE("right","top")
            mass_x(1,1,:,:)=nhat(1,:,:)
            mass_y(1,1,:,:)=nhat(2,:,:)
         END SELECT
      CASE("Harris_res")
         SELECT CASE(lrtb)
         CASE("top")
            mass_y(1,1,:,:)=one
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
      USE rmhd_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8) :: d_eta,d2
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: eta_local,nu_local
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
      eta_local=eta
      nu_local=nu
      IF(init_type=="Harris_open")THEN
         d_eta=10*eta
         d2=10
         eta_local=eta + d_eta
     $        *(EXP(-(d2*(x-lx)/lx)**4) + EXP(-(d2*(x+lx)/lx)**4))
      ELSEIF(init_type=="Harris_local")THEN
         WHERE(u(3,:,:) > j_crit)
            eta_local = eta*(one + (u(3,:,:)-j_crit)**2/j_crit**2)
         END WHERE
      ELSEIF(init_type=="anom_res" .OR. eta_type=="a1")THEN
         WHERE(u(3,:,:)**2 >= j_crit**2)
            eta_local=one/min_eps
         ELSEWHERE
            eta_local = eta*half*(one 
     $           + (one - (u(3,:,:)/j_crit)**2)**(-0.5))
         END WHERE
      ELSEIF(init_type=="anom_hres")THEN
         nu_local = nu*SQRT(half + half/(one-(u(3,:,:)/j_crit)**2)**2)
      ENDIF
c-----------------------------------------------------------------------
c     reduced MHD + electron inertia + sound radius model.
c-----------------------------------------------------------------------
      fx(1,:,:)=-uy(4,:,:)*(u(1,:,:)-de_sq*u(3,:,:))
     $     + sound_rad*u(1,:,:)*uy(2,:,:) + nu_local*ux(3,:,:)
      
      fy(1,:,:)=ux(4,:,:)*(u(1,:,:)-de_sq*u(3,:,:))
     $     - sound_rad*u(1,:,:)*ux(2,:,:) + nu_local*uy(3,:,:)
      
      s(1,:,:)=eta_local*u(3,:,:)

      fx(2,:,:)=-uy(4,:,:)*u(2,:,:)+u(3,:,:)*uy(1,:,:) - mu*ux(2,:,:)
      
      fy(2,:,:)= ux(4,:,:)*u(2,:,:)-u(3,:,:)*ux(1,:,:) - mu*uy(2,:,:)
      
      fx(3,:,:)=ux(1,:,:)
      fy(3,:,:)=uy(1,:,:)
      s(3,:,:)=u(3,:,:)
      
      fx(4,:,:)=ux(4,:,:)
      fy(4,:,:)=uy(4,:,:)
      s(4,:,:)=u(2,:,:)

      IF(source .AND. first)THEN
         CALL rmhd_equil(x,y,u0,u0x,u0y,.TRUE.)
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
      USE rmhd_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8) :: d_eta,d2
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: eta_local,eta_u3,
     $     nu_local,nu_u3
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
      eta_local=eta
      eta_u3=0
      nu_local=nu
      nu_u3=0
      IF(init_type=="Harris_open")THEN
         d_eta=10*eta
         d2=10
         eta_local=eta + d_eta
     $        *(EXP(-(d2*(x-lx)/lx)**4) + EXP(-(d2*(x+lx)/lx)**4))
      ELSEIF(init_type=="Harris_local")THEN
         WHERE(u(3,:,:) > j_crit)
            eta_local = eta*(one + (u(3,:,:)-j_crit)**2/j_crit**2)
            eta_u3 = two*eta*(u(3,:,:)-j_crit)/j_crit**2
         END WHERE
      ELSEIF(init_type=="anom_res" .OR. eta_type=="a1")THEN
         WHERE(u(3,:,:)**2 >= j_crit**2)
            eta_local=one/min_eps
            eta_u3=0
         ELSEWHERE
            eta_local = eta*half*(one 
     $           + (one - (u(3,:,:)/j_crit)**2)**(-0.5))
            eta_u3 = eta*u(3,:,:)
     $           /(two*j_crit**2*(one-(u(3,:,:)/j_crit)**2)**1.5)
         ENDWHERE
      ELSEIF(init_type=="anom_hres")THEN         
         nu_local = nu*SQRT(half + half/(one-(u(3,:,:)/j_crit)**2)**2)
         nu_u3 = nu*u(3,:,:)
     $        /(SQRT(half + half/(one-(u(3,:,:)/j_crit)**2)**2)
     $        *j_crit**2*(one-(u(3,:,:)/j_crit)**2)**3)
      ENDIF
c-----------------------------------------------------------------------
c     reduced MHD + electron inertia + sound radius model.
c-----------------------------------------------------------------------
      fx_u(1,1,:,:)=-uy(4,:,:)+sound_rad*uy(2,:,:)
      fx_u(1,3,:,:)=de_sq*uy(4,:,:)
      fx_ux(1,3,:,:)=nu_local
      fx_uy(1,2,:,:)=sound_rad*u(1,:,:)
      fx_uy(1,4,:,:)=-(u(1,:,:)-de_sq*u(3,:,:))

      fy_u(1,1,:,:)=ux(4,:,:)-sound_rad*ux(2,:,:)
      fy_u(1,3,:,:)=-de_sq*ux(4,:,:)
      fy_uy(1,3,:,:)=nu_local
      fy_ux(1,2,:,:)=-sound_rad*u(1,:,:)
      fy_ux(1,4,:,:)=u(1,:,:)-de_sq*u(3,:,:)

      s_u(1,3,:,:)=eta_local + eta_u3*u(3,:,:)
      
      fx_u(2,2,:,:)=-uy(4,:,:)
      fx_u(2,3,:,:)=uy(1,:,:)
      fx_ux(2,2,:,:)=-mu
      fx_uy(2,1,:,:)=u(3,:,:)
      fx_uy(2,4,:,:)=-u(2,:,:)

      fy_u(2,2,:,:)=ux(4,:,:)
      fy_u(2,3,:,:)=-ux(1,:,:)
      fy_ux(2,1,:,:)=-u(3,:,:)
      fy_ux(2,4,:,:)=u(2,:,:)
      fy_uy(2,2,:,:)=-mu

      fx_ux(3,1,:,:)=one
      fy_uy(3,1,:,:)=one
      s_u(3,3,:,:)=one
      
      fx_ux(4,4,:,:)=one
      fy_uy(4,4,:,:)=one
      s_u(4,2,:,:)=one
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
      USE rmhd_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      mass(1,3,:,:)=-de_sq
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
      USE rmhd_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,phi
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("islands")
         ksi=(x**2+.4*x)/1.4
         phi=(y**2+y)/2
      CASE("Harris")
         ksi=lx*(x-.5)
         phi=ly*(y-.5)
      CASE("Harris_open")
         ksi=lx*x
         phi=ly*(TANH(4*y-2.)+TANH(2.))/(2*TANH(2.))
      CASE("Harris_local","Harris_res","anom_hres","plasmoids")
         ksi=lx*x
ccc         ksi=lx*(TANH(gr_curve*(x-1))+TANH(gr_curve))/TANH(gr_curve)
         phi=ly*(TANH(gr_curve*(y-1))+TANH(gr_curve))/TANH(gr_curve)
      CASE("anom_res")
         ksi=lx*(TANH(gr_curve*(x-1))+TANH(gr_curve))/TANH(gr_curve)
         phi=ly*(TANH(gr_curve*(y-1))+TANH(gr_curve))/TANH(gr_curve)
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
      USE rmhd_mod
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
c     deallocates private module objects.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_dealloc
      USE rmhd_mod
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
