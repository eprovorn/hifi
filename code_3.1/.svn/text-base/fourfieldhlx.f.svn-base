c-----------------------------------------------------------------------
c     file fourfieldhlx.f.
c     contains specifications for RMHD equations with helical symmetry.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. fourfieldhlx_mod.
c     1. fourfieldhlx_equil.
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
c     subprogram 0. fourfieldhlx_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE fourfieldhlx_mod
      USE local_mod
      IMPLICIT NONE
      
      LOGICAL :: source=.FALSE.
      CHARACTER(16) :: init_type="."
      REAL(r8) :: gr_curve=0,di=0,eta=0,mu=0,nu=0,eta_0=0,mu_0=0,nu_0=0,
     $     t0_init=0,t0_decay=0,epsilon=0,Tconst=0,beta_e=0,beta0=0,
     $     r_s=0,Rinv=0

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fourfieldhlx_equil.
c     computes equilibrium for Helical RMHD model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fourfieldhlx_equil(x,y,u,ux,uy,derivs)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs

      INTEGER :: ix1,ix2
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      gfac=1/(1+(x*Rinv)**2)
      SELECT CASE(init_type)
      CASE("m1_park")
         u(1,:,:)=x**6*Rinv**2/12 + x**4*(1-r_s**2*Rinv**2)/8 
     $        - x**2*r_s**2/4
         u(2,:,:)=1/Rinv
         u(5,:,:)=(2*x**2-r_s**2)/gfac - 2*gfac
         IF(.NOT. derivs)RETURN
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         ux(1,:,:)=x**5*Rinv**2/2 + x**3*(1 - r_s**2*Rinv**2)/2 
     $        - x*r_s**2/2
         ux(5,:,:)=2*x*((2*x**2-r_s**2)*Rinv**2 + 2/gfac 
     $        + 2*(gfac*Rinv)**2)
      CASE("m1_axial","m1_axial_new")
         IF(r_s <= 0)CALL program_stop("r_s has to be positive!")
         u(1,:,:)=(x**4*Rinv**2/3. + x**2*(1. - (r_s*Rinv)**2)/2. 
     $        - r_s**2)*x**2/4.
         u(2,:,:)=1/(gfac*Rinv) - Rinv*x**2/(2*gfac)*(x**2-r_s**2)
         u(5,:,:)=(2*x**2-r_s**2)/gfac + (Rinv*x)**2*(x**2-r_s**2) - 2
         DO ix1 = 1,SIZE(x,1)
            DO ix2 = 1,SIZE(x,2)
               IF(x(ix1,ix2) <= r_s)THEN
                  u(1,ix1,ix2)=u(1,ix1,ix2) - half/(8*r_s**4)
     $                 *(x(ix1,ix2)**2 - r_s**2)**4
                  u(5,ix1,ix2)=u(5,ix1,ix2) - one/r_s**4
     $                 *(x(ix1,ix2)**2 - r_s**2)**2
     $                 *(gfac(ix1,ix2)*(x(ix1,ix2)**2 - r_s**2) 
     $                 + 3*x(ix1,ix2)**2)
               ENDIF
            ENDDO
         ENDDO
      CASE("m1_zero","m1_zero_new")
         IF(r_s <= 0)CALL program_stop("r_s has to be positive!")
         u(1,:,:) = (half*x**2 - r_s**2)*x**2/4._r8
         u(2,:,:) = one/(gfac*Rinv) - half*Rinv*x**2*(x**2-r_s**2)
         u(5,:,:) = two*x**2 - r_s**2 - two
         DO ix1 = 1,SIZE(x,1)
            DO ix2 = 1,SIZE(x,2)
               IF(x(ix1,ix2) <= r_s)THEN
                  u(1,ix1,ix2)=u(1,ix1,ix2) - one/(16._r8*r_s**4)
     $                 *(x(ix1,ix2)**2 - r_s**2)**4
                  u(5,ix1,ix2)=u(5,ix1,ix2) - one/r_s**4
     $                 *(x(ix1,ix2)**2 - r_s**2)**2
     $                 *(gfac(ix1,ix2)*(x(ix1,ix2)**2 - r_s**2) 
     $                 + 3._r8*x(ix1,ix2)**2)
               ENDIF
            ENDDO
         ENDDO
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourfieldhlx_equil
      END MODULE fourfieldhlx_mod
c-----------------------------------------------------------------------
c     subprogram a. physics_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,
     $     nqty,nqty_schur,dt,dtmax,tmax,nstep,du_diagnose,physics_type,
     $     exit_flag)
      USE fourfieldhlx_mod
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
      NAMELIST/kink_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,di,eta,mu,nu,eta_0,mu_0,nu_0,
     $     t0_init,t0_decay,beta0,beta_e,Tconst,r_s,Rinv,epsilon,
     $     source,init_type
c-----------------------------------------------------------------------
c     Sample namelist.
c-----------------------------------------------------------------------
ccc&kink_list
ccc
ccc	init_type="m1_axial"
ccc	source=f
ccc     gr_curve=1.25
ccc
ccc	beta0=2.
ccc	beta_e=0.5
ccc	r_s=0.5
ccc	Rinv=0.1
ccc	Tconst=0.1
ccc
ccc	di=2.e-2
ccc	eta=1.e-5
ccc	mu=5.e-5
ccc	nu=2.e-8
ccc
ccc	t0_decay=0.
ccc	t0_init=0.
ccc	eta_0=0.
ccc	mu_0=0.
ccc	nu_0=0.
ccc
ccc	epsilon=1.e-3
ccc/
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=kink_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="fourfieldhlx"

      nqty=6
      nqty_schur=0

      xperiodic=.FALSE.
      yperiodic=.TRUE.
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
      USE fourfieldhlx_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      ground(6)=.TRUE.
      static(5:6)=.TRUE.
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
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta_e,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Rinv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(r_s,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(t0_decay,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(t0_init,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta_0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu_0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu_0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Tconst,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("m1_axial","m1_axial_new","m1_zero","m1_zero_new")
         IF(eta == 0)THEN
            CALL program_stop("Have to have non-zero resistivity !")
         ELSEIF(source)THEN
            CALL program_stop("Do not run m1_*** with source term !")
         ENDIF
      END SELECT
      IF((init_type=="m1_zero" .OR. init_type == "m1_zero_new")
     $     .AND. Tconst > (r_s**2+two)/(two*r_s**2))
     $     CALL program_stop("Resistivity is negative on axis !")
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
      USE fourfieldhlx_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: ui

      REAL(r8), DIMENSION(SIZE(ui,1),SIZE(ui,2),SIZE(ui,3)) :: 
     $     ux,uy
c-----------------------------------------------------------------------
c     add perturbation to an equilibrium
c-----------------------------------------------------------------------
      CALL fourfieldhlx_equil(xpi,ypi,ui,ux,uy,.FALSE.)

      SELECT CASE(init_type)
      CASE("m1_park","m1_axial","m1_axial_new","m1_zero","m1_zero_new")
         ui(2,:,:)=ui(2,:,:) + epsilon*COS(ypi)
     $        *EXP(-((xpi-r_s)/gr_curve)**2)
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
      USE fourfieldhlx_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      left%bc_type="polar"
      SELECT CASE(init_type)
      CASE("m1_park")
         right%bc_type="robin"
         right%static(1)=.FALSE.
         right%static(2:4)=.TRUE.
         right%static(5)=.FALSE.
         right%static(6)=.TRUE.
      CASE("m1_axial","m1_axial_new","m1_zero","m1_zero_new")
         right%bc_type="robin"
         right%static(1:2)=.FALSE.
         right%static(3:4)=.TRUE.
         right%static(5)=.FALSE.
         right%static(6)=.TRUE.
      END SELECT
      top%bc_type = "periodic"
      bottom%bc_type = "periodic"
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
      USE fourfieldhlx_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      REAL(r8) :: gfac,gfac_s,Eb,eta_local,nu_t
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
      gfac=1/(1+Rinv**2)
      gfac_s=1/(1+(r_s*Rinv)**2)

      SELECT CASE(init_type)
      CASE("m1_park")
         Eb=eta*(r_s**2/gfac_s - 2*gfac_s)
         eta_local=ABS(Eb/((2 - r_s**2)/gfac - 2*gfac))
         IF(lrtb == "right")THEN
            IF(.NOT. source)c(1,:,:)=Eb
            IF(eta /= 0 .OR. di /= 0)THEN
               IF(source)THEN
                  c(2,:,:)=eta*ux(2,:,:) + di*u(2,:,:)*uy(2,:,:)
     $                 + 2*nu*di*Rinv*gfac*(di*ux(5,:,:) - ux(4,:,:)
     $                 - 2*Rinv**2*gfac*(di*u(5,:,:) - u(4,:,:)))
               ELSE
                  c(2,:,:)=eta_local*ux(2,:,:) + di*u(2,:,:)*uy(2,:,:)
     $                 + 2*nu*di*Rinv*gfac*(di*ux(5,:,:) - ux(4,:,:)
     $                 - 2*Rinv**2*gfac*(di*u(5,:,:) - u(4,:,:)))
               ENDIF
            ELSE
               c(2,:,:)=ux(2,:,:)
            ENDIF
            c(3,:,:)=u(3,:,:) + 2*Rinv*gfac*u(4,:,:) - gfac*ux(6,:,:)
            c(4,:,:)=ux(4,:,:) - gfac*Rinv**2*u(4,:,:)
            c(6,:,:)=u(6,:,:)
         ENDIF
      CASE("m1_axial","m1_axial_new","m1_zero","m1_zero_new")
         IF(init_type == "m1_axial")THEN
            Eb=eta*(r_s**2/gfac_s - 2)
            eta_local=ABS(Eb/((1+2*Rinv**2)*(1-r_s**2) + 1/gfac - 2))
            nu_t = nu
         ELSEIF(init_type == "m1_axial_new")THEN
            Eb=((eta-eta_0)*EXP(-((t-t0_init)/t0_decay)**2) + eta_0)
     $           *(r_s**2/gfac_s - 2)
            eta_local=ABS(Eb/((1+2*Rinv**2)*(1-r_s**2) + 1/gfac - 2))
            nu_t=(nu-nu_0)*EXP(-((t-t0_init)/t0_decay)**2) + nu_0
         ELSEIF(init_type == "m1_zero" .OR. init_type == "m1_zero_new")
     $           THEN
            Eb=eta*(r_s**2 - 2)
            eta_local=eta*(2 - r_s**2)/r_s**2
            nu_t = nu
         ENDIF
         IF(lrtb == "right")THEN
            c(1,:,:)=Eb/gfac 
     $           + Rinv*(di*u(2,:,:)*uy(2,:,:) + eta_local*ux(2,:,:)
     $           + 2*nu_t*di*Rinv*gfac*(di*ux(5,:,:) - ux(4,:,:)
     $           - 2*Rinv**2*gfac*(di*u(5,:,:) - u(4,:,:))))
            c(2,:,:)=-eta_local*ux(2,:,:) - di*u(2,:,:)*uy(2,:,:)
     $           - 2*nu_t*di*Rinv*gfac*(di*ux(5,:,:) - ux(4,:,:)
     $           - 2*Rinv**2*gfac*(di*u(5,:,:) - u(4,:,:)))
            c(3,:,:)=u(3,:,:) + 2*Rinv*gfac*u(4,:,:) - gfac*ux(6,:,:)
            c(4,:,:)=ux(4,:,:) - gfac*Rinv**2*u(4,:,:)
            c(6,:,:)=u(6,:,:)
         ENDIF
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
      USE fourfieldhlx_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy

      REAL(r8) :: gfac,gfac_s,Eb,eta_local,nu_t
c-----------------------------------------------------------------------
c     zero out output.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
      c_uxy=0

      gfac=1/(1+Rinv**2)
      gfac_s=1/(1+(r_s*Rinv)**2)

      SELECT CASE(init_type)
      CASE("m1_park")
         Eb=eta*(r_s**2/gfac_s - 2*gfac_s)
         eta_local=ABS(Eb/((2 - r_s**2)/gfac - 2*gfac))
         IF(lrtb == "right")THEN
            IF(eta /= 0 .OR. di /= 0)THEN
               IF(source)THEN
                  c_ux(2,2,:,:)=eta
               ELSE
                  c_ux(2,2,:,:)=eta_local
               ENDIF
               c_u(2,2,:,:)=di*uy(2,:,:)
               c_u(2,4,:,:)=4*nu*di*Rinv**3*gfac**2
               c_u(2,5,:,:)=-4*nu*di**2*Rinv**3*gfac**2
               c_ux(2,4,:,:)=-2*nu*di*Rinv*gfac
               c_ux(2,5,:,:)=2*nu*di**2*Rinv*gfac
               c_uy(2,2,:,:)=di*u(2,:,:)
            ELSE
               c_ux(2,2,:,:)=one
            ENDIF
            c_u(3,3,:,:)=one
            c_u(3,4,:,:)=2*Rinv*gfac
            c_ux(3,6,:,:)=-gfac
            c_u(4,4,:,:)=-gfac*Rinv**2
            c_ux(4,4,:,:)=one
            c_u(6,6,:,:)=one
         ENDIF
      CASE("m1_axial","m1_axial_new","m1_zero","m1_zero_new")
         IF(init_type == "m1_axial")THEN
            Eb=eta*(r_s**2/gfac_s - 2)
            eta_local=ABS(Eb/((1+2*Rinv**2)*(1-r_s**2) + 1/gfac - 2))
            nu_t = nu
         ELSEIF(init_type == "m1_axial_new")THEN
            Eb=((eta-eta_0)*EXP(-((t-t0_init)/t0_decay)**2) + eta_0)
     $           *(r_s**2/gfac_s - 2)
            eta_local=ABS(Eb/((1+2*Rinv**2)*(1-r_s**2) + 1/gfac - 2))
            nu_t=(nu-nu_0)*EXP(-((t-t0_init)/t0_decay)**2) + nu_0
         ELSEIF(init_type == "m1_zero" .OR. init_type == "m1_zero_new")
     $           THEN
            eta_local=eta*(two - r_s**2)/r_s**2
            nu_t = nu
         ENDIF
         IF(lrtb == "right")THEN
            c_u(1,2,:,:)=Rinv*di*uy(2,:,:)
            c_u(1,4,:,:)=4*nu_t*di*Rinv**4*gfac**2
            c_u(1,5,:,:)=-4*nu_t*di**2*Rinv**4*gfac**2
            c_ux(1,2,:,:)=Rinv*eta_local
            c_ux(1,4,:,:)=-2*nu_t*di*Rinv**2*gfac
            c_ux(1,5,:,:)=2*nu_t*di**2*Rinv**2*gfac
            c_uy(1,2,:,:)=Rinv*di*u(2,:,:)
            c_u(2,2,:,:)=-di*uy(2,:,:)
            c_u(2,4,:,:)=-4*nu_t*di*Rinv**3*gfac**2
            c_u(2,5,:,:)=4*nu_t*di**2*Rinv**3*gfac**2
            c_ux(2,2,:,:)=-eta_local
            c_ux(2,4,:,:)=2*nu_t*di*Rinv*gfac
            c_ux(2,5,:,:)=-2*nu_t*di**2*Rinv*gfac
            c_uy(2,2,:,:)=-di*u(2,:,:)
            c_u(3,3,:,:)=one
            c_u(3,4,:,:)=2*Rinv*gfac
            c_ux(3,6,:,:)=-gfac
            c_u(4,4,:,:)=-gfac*Rinv**2
            c_ux(4,4,:,:)=one
            c_u(6,6,:,:)=one
         ENDIF
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
      USE fourfieldhlx_mod
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
      CASE("m1_park")
         SELECT CASE(lrtb)
         CASE("right")
            mass(1,1,:,:)=one
            mass(5,5,:,:)=one
         END SELECT
      CASE("m1_axial","m1_axial_new","m1_zero","m1_zero_new")
         SELECT CASE(lrtb)
         CASE("right")
            mass(1,1,:,:)=one
            mass(2,1,:,:)=Rinv
            mass(5,5,:,:)=one
            mass_x(5,2,:,:)=-Rinv*nhat(1,:,:)
            mass_y(5,2,:,:)=-Rinv*nhat(2,:,:)
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
      USE fourfieldhlx_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8) :: gfac_s,Eb,mu_t,nu_t
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac,Ve,Vex,Vey,
     $     eta_local
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     calculate auxiliary variables.
c-----------------------------------------------------------------------
      gfac=1/(1+(x*Rinv)**2)
      Ve=u(4,:,:)-di*u(5,:,:)
      Vex=ux(4,:,:)-di*ux(5,:,:)
      Vey=uy(4,:,:)-di*uy(5,:,:)
c-----------------------------------------------------------------------
c     modify resistivity for boundary applied E-field(no current drive).
c-----------------------------------------------------------------------
      eta_local=eta
      mu_t=mu
      nu_t=nu
      gfac_s=1/(1+(r_s*Rinv)**2)
      SELECT CASE(init_type)
      CASE("m1_park")
         IF(.NOT. source)
     $        eta_local=eta*ABS((r_s**2/gfac_s - 2*gfac_s)
     $        /((2*x**2-r_s**2)/gfac - 2*gfac))
      CASE("m1_axial")
         Eb=eta*(r_s**2/gfac_s - 2)
         eta_local=ABS(Eb/((1+2*(x*Rinv)**2)*(x**2-r_s**2) 
     $        + x**2/gfac - 2))
      CASE("m1_axial_new")
         Eb=((eta-eta_0)*EXP(-((t-t0_init)/t0_decay)**2) + eta_0)
     $        *(r_s**2/gfac_s - 2)
         eta_local=ABS(Eb/((1+2*(x*Rinv)**2)*(x**2-r_s**2) 
     $        + x**2/gfac - 2))
         mu_t=(mu - mu_0)*EXP(-((t-t0_init)/t0_decay)**2) + mu_0
         nu_t=(nu - nu_0)*EXP(-((t-t0_init)/t0_decay)**2) + nu_0
      CASE("m1_zero")
         eta_local=eta*(two - r_s**2)/(two + r_s**2 - two*x**2)
         WHERE(x < r_s)
            eta_local=eta*(one + Tconst*two*(x**2 - r_s**2)
     $           /(r_s**2 + two*(one - x**2)))
         END WHERE
      CASE("m1_zero_new")
         eta_local=eta*(two - r_s**2)/(two + r_s**2 - two*x**2)
         WHERE(x < r_s)
            eta_local = eta*(one + Tconst
     $           *((one-eta_0)*EXP(-(t/t0_decay)**2) + eta_0)
     $           *two*(x**2 - r_s**2)/(r_s**2 + two*(one - x**2)))
         END WHERE
      END SELECT
c-----------------------------------------------------------------------
c     Helical RMHD model.
c-----------------------------------------------------------------------
      fx(1,:,:)=-di*nu_t*x*gfac*Vex
      fy(1,:,:)=-di*nu_t/x*Vey

      s(1,:,:)=gfac*(ux(1,:,:)*uy(6,:,:)-uy(1,:,:)*ux(6,:,:)
     $     + eta_local*x*u(5,:,:) 
     $     + di*(ux(1,:,:)*uy(2,:,:) - uy(1,:,:)*ux(2,:,:))
     $     - 4*di*nu_t*Rinv**2*x*gfac**2*Ve)

      fx(2,:,:)=gfac*(-u(2,:,:)*(uy(6,:,:)+di*uy(2,:,:)) 
     $     + uy(1,:,:)*Ve - eta_local*x*ux(2,:,:)
     $     - 2*di*nu_t*Rinv*x*gfac*Vex
     $     + 4*di*nu_t*Rinv**3*x**2*gfac**2*Ve)
      
      fy(2,:,:)=gfac*(u(2,:,:)*(ux(6,:,:) + di*ux(2,:,:)) 
     $     - ux(1,:,:)*Ve - 2*di*nu_t*Rinv/x*Vey)
     $     - eta_local/x*uy(2,:,:)

      fx(3,:,:)=gfac*(uy(1,:,:)*(u(5,:,:) + 2*Rinv*gfac*u(2,:,:))
     $     - uy(6,:,:)*(u(3,:,:) + 2*Rinv*gfac*u(4,:,:)) 
     $     - mu_t*x*ux(3,:,:) + 2*nu_t*Rinv*x*gfac*Vex
     $     - 4*nu_t*Rinv**3*x**2*gfac**2*Ve)
      
      fy(3,:,:)=gfac*(x*gfac*Rinv**2*(u(4,:,:)**2 - u(2,:,:)**2) 
     $     + ux(6,:,:)*(u(3,:,:) + 2*Rinv*gfac*u(4,:,:))
     $     - ux(1,:,:)*(u(5,:,:) + 2*Rinv*gfac*u(2,:,:))
     $     + 2*nu_t*Rinv/x*Vey) - mu_t/x*uy(3,:,:)
      
      s(3,:,:)=-8*Rinv**3*x*gfac**3
     $     *(u(2,:,:)*uy(1,:,:) - u(4,:,:)*uy(6,:,:))
      
      fx(4,:,:)=-x*gfac*(mu_t*ux(4,:,:) + nu_t*Vex)
      fy(4,:,:)=-1/x*(mu_t*uy(4,:,:) + nu_t*Vey)

      s(4,:,:)=gfac*(ux(4,:,:)*uy(6,:,:) - uy(4,:,:)*ux(6,:,:)
     $     + ux(1,:,:)*uy(2,:,:) - uy(1,:,:)*ux(2,:,:)
     $     + 2*mu_t*Rinv*x*gfac*u(3,:,:)
     $     - 4*nu_t*Rinv**2*x*gfac**2*Ve)
      
      fx(5,:,:)=x*gfac*ux(1,:,:)
      fy(5,:,:)=uy(1,:,:)/x
      s(5,:,:)=x*gfac*(u(5,:,:) + 2*Rinv*gfac*u(2,:,:))
      
      fx(6,:,:)=x*gfac*ux(6,:,:)
      fy(6,:,:)=uy(6,:,:)/x
      s(6,:,:)=x*gfac*(u(3,:,:) + 2*Rinv*gfac*u(4,:,:))

      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("m1_park")
            CALL fourfieldhlx_equil(x,y,u0,u0x,u0y,.TRUE.)
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
      USE fourfieldhlx_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8) :: gfac_s,Eb,mu_t,nu_t
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac,hfac,kfac,Ve,
     $     eta_local
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
c     calculate auxiliary variables.
c-----------------------------------------------------------------------
      gfac=1/(1+(x*Rinv)**2)
      hfac=2*Rinv*gfac**2
      kfac=2*Rinv*x*gfac*hfac
      Ve=u(4,:,:)-di*u(5,:,:)
c-----------------------------------------------------------------------
c     modify resistivity for boundary applied E-field(no current drive).
c-----------------------------------------------------------------------
      eta_local=eta
      mu_t=mu
      nu_t=nu
      gfac_s=1/(1+(r_s*Rinv)**2)
      SELECT CASE(init_type)
      CASE("m1_park")
         IF(.NOT. source)
     $        eta_local=eta*ABS((r_s**2/gfac_s - 2*gfac_s)
     $        /((2*x**2-r_s**2)/gfac - 2*gfac))
      CASE("m1_axial")
         Eb=eta*(r_s**2/gfac_s - 2)
         eta_local=ABS(Eb/((1+2*(x*Rinv)**2)*(x**2-r_s**2) 
     $        + x**2/gfac - 2))
      CASE("m1_axial_new")
         Eb=((eta-eta_0)*EXP(-((t-t0_init)/t0_decay)**2) + eta_0)
     $        *(r_s**2/gfac_s - 2)
         eta_local=ABS(Eb/((1+2*(x*Rinv)**2)*(x**2-r_s**2) 
     $        + x**2/gfac - 2))
         mu_t=(mu - mu_0)*EXP(-((t-t0_init)/t0_decay)**2) + mu_0
         nu_t=(nu - nu_0)*EXP(-((t-t0_init)/t0_decay)**2) + nu_0
      CASE("m1_zero")
         eta_local=eta*(two - r_s**2)/(two + r_s**2 - two*x**2)
         WHERE(x < r_s)
            eta_local=eta*(one + Tconst*two*(x**2 - r_s**2)
     $           /(r_s**2 + two*(one - x**2)))
         END WHERE
      CASE("m1_zero_new")
         eta_local=eta*(two - r_s**2)/(two + r_s**2 - two*x**2)
         WHERE(x < r_s)
            eta_local = eta*(one + Tconst
     $           *((one-eta_0)*EXP(-(t/t0_decay)**2) + eta_0)
     $           *two*(x**2 - r_s**2)/(r_s**2 + two*(one - x**2)))
         END WHERE
      END SELECT
c-----------------------------------------------------------------------
c     Helical RMHD model.
c-----------------------------------------------------------------------
      fx_ux(1,4,:,:)=-di*nu_t*x*gfac
      fx_ux(1,5,:,:)=di**2*nu_t*x*gfac
      fy_uy(1,4,:,:)=-di*nu_t/x
      fy_uy(1,5,:,:)=di**2*nu_t/x

      s_u(1,4,:,:)=-di*nu_t*kfac
      s_u(1,5,:,:)=eta_local*x*gfac + di**2*nu_t*kfac
      s_ux(1,1,:,:)=gfac*(uy(6,:,:) + di*uy(2,:,:))
      s_ux(1,2,:,:)=-gfac*di*uy(1,:,:)
      s_ux(1,6,:,:)=-gfac*uy(1,:,:)
      s_uy(1,1,:,:)=-gfac*(ux(6,:,:) + di*ux(2,:,:))
      s_uy(1,2,:,:)=gfac*di*ux(1,:,:)
      s_uy(1,6,:,:)=gfac*ux(1,:,:)
      
      fx_u(2,2,:,:)=-gfac*(uy(6,:,:) + di*uy(2,:,:))
      fx_u(2,4,:,:)=gfac*uy(1,:,:) + di*nu_t*Rinv*x*kfac
      fx_u(2,5,:,:)=-di*(gfac*uy(1,:,:) + di*nu_t*Rinv*x*kfac)
      fx_ux(2,2,:,:)=-eta_local*x*gfac
      fx_ux(2,4,:,:)=-di*nu_t*x*hfac
      fx_ux(2,5,:,:)=di**2*nu_t*x*hfac
      fx_uy(2,1,:,:)=gfac*Ve
      fx_uy(2,2,:,:)=-di*gfac*u(2,:,:)
      fx_uy(2,6,:,:)=-gfac*u(2,:,:)
            
      fy_u(2,2,:,:)=gfac*(ux(6,:,:) + di*ux(2,:,:))
      fy_u(2,4,:,:)=-gfac*ux(1,:,:)
      fy_u(2,5,:,:)=di*gfac*ux(1,:,:)
      fy_ux(2,1,:,:)=-gfac*Ve
      fy_ux(2,2,:,:)=di*gfac*u(2,:,:)
      fy_ux(2,6,:,:)=gfac*u(2,:,:)
      fy_uy(2,2,:,:)=-eta_local/x
      fy_uy(2,4,:,:)=-2*di*nu_t*Rinv/x*gfac
      fy_uy(2,5,:,:)=2*di**2*nu_t*Rinv/x*gfac
      
      fx_u(3,2,:,:)=hfac*uy(1,:,:)
      fx_u(3,3,:,:)=-gfac*uy(6,:,:)
      fx_u(3,4,:,:)=-hfac*uy(6,:,:) - nu_t*Rinv*x*kfac
      fx_u(3,5,:,:)=gfac*uy(1,:,:) + di*nu_t*Rinv*x*kfac
      fx_ux(3,3,:,:)=-mu_t*x*gfac
      fx_ux(3,4,:,:)=nu_t*x*hfac
      fx_ux(3,5,:,:)=-di*nu_t*x*hfac
      fx_uy(3,1,:,:)=gfac*(u(5,:,:) + 2*Rinv*gfac*u(2,:,:))
      fx_uy(3,6,:,:)=-gfac*(u(3,:,:) + 2*Rinv*gfac*u(4,:,:))
      
      fy_u(3,2,:,:)=-hfac*(x*Rinv*u(2,:,:) + ux(1,:,:))
      fy_u(3,3,:,:)=gfac*ux(6,:,:)
      fy_u(3,4,:,:)=hfac*(x*Rinv*u(4,:,:)+ux(6,:,:))
      fy_u(3,5,:,:)=-gfac*ux(1,:,:)
      fy_ux(3,1,:,:)=-gfac*u(5,:,:) - hfac*u(2,:,:)
      fy_ux(3,6,:,:)=gfac*u(3,:,:) + hfac*u(4,:,:)
      fy_uy(3,3,:,:)=-mu_t/x
      fy_uy(3,4,:,:)=2*nu_t*Rinv/x*gfac
      fy_uy(3,5,:,:)=-2*di*nu_t*Rinv/x*gfac

      s_u(3,2,:,:)=-2*Rinv*kfac*uy(1,:,:)
      s_u(3,4,:,:)=2*Rinv*kfac*uy(6,:,:)
      s_uy(3,1,:,:)=-2*Rinv*kfac*u(2,:,:)
      s_uy(3,6,:,:)=2*Rinv*kfac*u(4,:,:)

      fx_ux(4,4,:,:)=-(mu_t+nu_t)*x*gfac
      fx_ux(4,5,:,:)=di*nu_t*x*gfac
      fy_uy(4,4,:,:)=-(mu_t+nu_t)/x
      fy_uy(4,5,:,:)=di*nu_t/x

      s_u(4,3,:,:)=mu_t*x*hfac
      s_u(4,4,:,:)=-nu_t*kfac
      s_u(4,5,:,:)=di*nu_t*kfac
      s_ux(4,1,:,:)=gfac*uy(2,:,:)
      s_ux(4,2,:,:)=-gfac*uy(1,:,:)
      s_ux(4,4,:,:)=gfac*uy(6,:,:)
      s_ux(4,6,:,:)=-gfac*uy(4,:,:)
      s_uy(4,1,:,:)=-gfac*ux(2,:,:)
      s_uy(4,2,:,:)=gfac*ux(1,:,:)
      s_uy(4,4,:,:)=-gfac*ux(6,:,:)
      s_uy(4,6,:,:)=gfac*ux(4,:,:)
      
      fx_ux(5,1,:,:)=x*gfac
      fy_uy(5,1,:,:)=1/x
      s_u(5,2,:,:)=x*hfac
      s_u(5,5,:,:)=x*gfac

      fx_ux(6,6,:,:)=x*gfac
      fy_uy(6,6,:,:)=1/x
      s_u(6,3,:,:)=x*gfac
      s_u(6,4,:,:)=x*hfac
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
      USE fourfieldhlx_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      gfac=1/(1+(x*Rinv)**2)

      mass(1,1,:,:)=x*gfac
      mass(2,2,:,:)=x*gfac
      mass(2,1,:,:)=-2*x*Rinv*gfac**2
      mass(3,3,:,:)=x*gfac
      mass(3,4,:,:)=4*Rinv*x*gfac**2
      mass(4,4,:,:)=x*gfac
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
      USE fourfieldhlx_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: ksi,phi
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: x,y

      REAL(r8) :: rs
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      rs=1/(1+((1-r_s)/r_s)**third)
      rs=rs+1.5*gr_curve*(r_s-rs)
      x = ((ksi-rs)**3+5*gr_curve*ksi+rs**3)/(1+5*gr_curve+3*rs*(rs-1))
      y = twopi*phi
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
      USE fourfieldhlx_mod
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
      USE fourfieldhlx_mod
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
