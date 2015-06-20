c-----------------------------------------------------------------------
c     file heat.f.
c     contains specifications for anisotropic heat conduction.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. heat_mod.
c     1. heat_Dtensor
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
c     subprogram 0. heat_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE heat_mod
      USE local_mod
      IMPLICIT NONE

      CHARACTER(16) :: init_type="."
      REAL(r8) :: lx=0,phi=0,kperp=0,tau=0,lambda=0,ly=0,kpar=one
      REAL(r8), PARAMETER :: x1=0.252,x2=3.748,y1=0.7755,y2=1.122,rs=0.1

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. heat_Dtensor.
c     Define the 4 components of the anisotropic conductivity tensor, d.
c     The tensor is defined by kpar in the x' direction and kperp in the
c     y' direction, where x' and y' c.s. is aligned with the principle
c     axis of D.  In the x' and y' c.s., Dp is defined.  A rotation
c     matrix R defines the rotation to x' and y'.  To map Dp onto the
c     x and y c.s., note that Dp*Rp*x = Rt*b, where Rt is the transpose
c     of the rotation matrix.  So, (R*Dp*Rt)*x = b, and it is clear that
c     (Rt*Dp*R) is the desired d.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE heat_Dtensor(x,y,d)
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y      
      REAL(r8), DIMENSION(0:SIZE(x,1)-1,0:SIZE(x,2)-1,2,2), 
     $     INTENT(OUT) :: d

      REAL(r8), DIMENSION(0:SIZE(x,1)-1,0:SIZE(x,2)-1) :: phir
c-----------------------------------------------------------------------
c     define phi based on the geometry.
c     - for "source_sink", phi is parallel to the walls.
c     - for "polar", phi is in the azimuthal direction.
c     - for "nimrod", phi is perpendicular to contours of 
c     cos(pi*x)*cos(pi*y).
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("source_sink")
         phir = ATAN( 0.5*0.347*10.0/COSH(10.0*(x-lx/2.0))**2 )
      CASE("polar")
         phir = ATAN2(y,x) + pi*half
         WHERE(x==0 .AND. y==0) phir=0
      CASE("nimrod")
         phir = ATAN2(SIN(pi*y)*COS(pi*x),SIN(pi*x)*COS(pi*y)) + pi*half
      CASE DEFAULT
         phir = phi
      END SELECT
c-----------------------------------------------------------------------
c     calculate the heat conduction tensor.
c-----------------------------------------------------------------------
      d(:,:,1,1)=kpar*COS(phir)**2+kperp*SIN(phir)**2
      d(:,:,1,2)=(kpar-kperp)*SIN(phir)*COS(phir)
      d(:,:,2,1)=d(:,:,1,2)
      d(:,:,2,2)=kpar*SIN(phir)**2+kperp*COS(phir)**2
      IF(init_type=="polar" .OR. init_type=="nimrod")THEN
         WHERE (x==0 .AND. y==0)
            d(:,:,1,1)=kpar
            d(:,:,1,2)=0.
            d(:,:,2,2)=kpar
            d(:,:,2,1)=0.
        ENDWHERE
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE heat_Dtensor
      END MODULE heat_mod
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
      USE heat_mod
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
      NAMELIST/heat_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $  tmax,nstep,lx,init_type,phi,kperp,kpar,tau,lambda,ly
c-----------------------------------------------------------------------
c     definitions.
c-----------------------------------------------------------------------
c     phi: angle in degrees from x-axis defining direction of parallel
c     heat conduction.
c     kperp: constant defining perpendicular heat conduction.
c     kpar: constant defining parallel heat conduction.
c     tau: time constant for time dependence of the "left" BC.
c     lambda: parameter in the Gaussian in the "left" BC.
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=heat_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="heat"
      
      phi=phi*pi/180._r8
      SELECT CASE(init_type)
      CASE("simple","poly4")
         lambda=lambda/COS(phi)
      END SELECT

      nqty=1
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
      USE heat_mod
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
      CALL MPI_Bcast(kperp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kpar,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(phi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(tau,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
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
      USE heat_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      REAL(r8), DIMENSION(0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) :: ypp
c-----------------------------------------------------------------------
c     initial condition.
c-----------------------------------------------------------------------
      ui=0
c-----------------------------------------------------------------------
c     simple  : time-dependent gaussian source.
c     poly4   : time-dependent 4th order polynomial source.
c     simpler : trans-domain gaussian initial condition.
c     polar : radius-dependent gaussian.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("simpler")
         ypp=-(xpi-lx*half)*SIN(phi)+(ypi-ly*half)*COS(phi)
         ui(1,:,:)=EXP(-(ypp/lambda)**2)
      CASE("polar")
         ui(1,:,:) = EXP(-(SQRT(xpi**2+ypi**2)-tau)**2/lambda**2)
      CASE("square_circle")
         ui(1,:,:) = EXP(-(xpi**2+ypi**2)/lambda**2)
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
      USE heat_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      top%bc_type="zeroflux"
      bottom%bc_type="zeroflux"
      left%bc_type="zeroflux"
      right%bc_type="zeroflux"

      left%static=.TRUE.
      right%static=.TRUE.
      top%static=.TRUE.
      bottom%static=.TRUE.

      SELECT CASE(init_type)
      CASE("simple","poly4")
         left%bc_type="robin"
      CASE("cubit_heat")
         right%bc_type="robin"
      CASE("nimrod")
         left%bc_type="robin"
         right%bc_type="robin"
         top%bc_type="robin"
         bottom%bc_type="robin"
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
      USE heat_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      INTEGER :: i,j
      REAL(r8) :: area_left
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
      SELECT CASE(init_type)
      CASE("simple")
         ! adjust gaussian so that width perpendicular to the field is
         ! always the same.
         SELECT CASE(lrtb)
         CASE("left")
            c(1,:,:)=u(1,:,:)
     $           -EXP(-(y-.5*ly)**2/lambda**2)*(1-EXP(-t/tau))
         END SELECT
      CASE("poly4")
         ! adjust polynomial so that width perpendicular to the field is
         ! always the same.
         SELECT CASE(lrtb)
         CASE("left")
            c(1,:,:)=u(1,:,:)-1/lambda**4*((y-1)**2-lambda**2)**2
     $             *(1-EXP(-t/tau))
            DO i=1,SIZE(y,1)
               DO j=1,SIZE(y,2)
                  IF(ABS(y(i,j)-1)>lambda)THEN
                     c(1,i,j)=u(1,i,j)
                  ENDIF
               ENDDO
            ENDDO
         END SELECT
      CASE("cubit_heat")
         SELECT CASE(lrtb)
         CASE("right")
            area_left = (1.398 + 0.347*TANH(-4.0))
            c(1,:,:)=u(1,:,:)-EXP(-y**2/lambda**2)
     $           *(1-EXP(-t/tau))
         END SELECT
      CASE("nimrod")
         c(1,:,:)=u(1,:,:)
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
      USE heat_mod
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
      CASE("simple","poly4")
         SELECT CASE(lrtb)
         CASE("left")
            c_u(1,1,:,:)=1._r8
         END SELECT
      CASE("cubit_heat")
         SELECT CASE(lrtb)
         CASE("right")
            c_u(1,1,:,:) = 1._r8
         END SELECT
      CASE("nimrod")
         c_u(1,1,:,:) = 1._r8
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
      USE heat_mod
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
      USE heat_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8), DIMENSION(0:SIZE(x,1)-1,0:SIZE(x,2)-1,2,2) :: d
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: sshape,ramp
      REAL(r8), PARAMETER :: alpha = 5.0
c-----------------------------------------------------------------------
c     calculate the heat conduction tensor.
c-----------------------------------------------------------------------
      CALL heat_Dtensor(x,y,d)
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     heat equation fluxes
c-----------------------------------------------------------------------
      fx(1,:,:)=-d(:,:,1,1)*ux(1,:,:)-d(:,:,1,2)*uy(1,:,:)
      fy(1,:,:)=-d(:,:,2,1)*ux(1,:,:)-d(:,:,2,2)*uy(1,:,:)

      SELECT CASE(init_type)
      CASE("source_sink")
         IF (t <= alpha) THEN
            ramp=0.5*(1-COS(2*pi*t/alpha))
         ELSE
            ramp=0.0
         ENDIF
         sshape = 0.10*ramp*(EXP(-((x-x1)**2 + (y-y1)**2)/rs**2)
     $        - EXP(-((x-x2)**2 + (y-y2)**2)/rs**2))
         s(1,:,:) = sshape
      CASE("nimrod")
         s(1,:,:) = two*pi**2*COS(pi*x)*COS(pi*y)
      END SELECT
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
      USE heat_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8), DIMENSION(0:SIZE(x,1)-1,0:SIZE(x,2)-1,2,2) :: d
c-----------------------------------------------------------------------
c     calculate the heat conduction tensor.
c-----------------------------------------------------------------------
      CALL heat_Dtensor(x,y,d)
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
c     flux derivative.
c-----------------------------------------------------------------------
      fx_ux(1,1,:,:)=-d(:,:,1,1)
      fx_uy(1,1,:,:)=-d(:,:,1,2)

      fy_ux(1,1,:,:)=-d(:,:,2,1)
      fy_uy(1,1,:,:)=-d(:,:,2,2)
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
      USE heat_mod
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
      SUBROUTINE physics_grid(x,y,ksi,eta)
      USE heat_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,eta

      INTEGER :: ix,iy
      REAL(r8) :: p1,q1,f1,f2,d1,d2
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("simple","poly4","simpler")
         ksi=lx*x
         eta=ly*y
      CASE("polar")
         ksi=lx*(x-.5)
         eta=ly*(y-.5)
      CASE("square_circle")
         DO ix=1,SIZE(x,1)
            DO iy=1,SIZE(x,2)

               p1=two*ABS(x(ix,iy)-half)
               q1=two*ABS(y(ix,iy)-half)

               d1=(SQRT((one+one/TAN(half*pi*p1))**2 
     $              - (COS(.25*pi*p1))**2)-SIN(.25*pi*p1))
               d2=(SQRT((one+one/TAN(half*pi*q1))**2 
     $              - (COS(.25*pi*q1))**2)-SIN(.25*pi*q1))

               f1=-(COS(half*pi*p1) 
     $              + SQRT((one+SIN(pi*p1))/(COS(.25*pi*p1))**2 
     $              - (SIN(half*pi*p1))**2))
               f2=-(COS(half*pi*q1) 
     $              + SQRT((one+SIN(pi*q1))/(COS(.25*pi*q1))**2 
     $              - (SIN(half*pi*q1))**2))

               IF(p1==one)THEN
                  ksi(ix,iy)=COS(.25*pi*q1)
                  eta(ix,iy)=SIN(.25*pi*q1)
               ELSEIF(q1==one)THEN
                  ksi(ix,iy)=SIN(.25*pi*p1)
                  eta(ix,iy)=COS(.25*pi*p1)
               ELSEIF(p1==zero .AND. q1==zero)THEN
                  ksi(ix,iy)=zero
                  eta(ix,iy)=zero
               ELSEIF(p1==zero)THEN
                  ksi(ix,iy)=zero
                  eta(ix,iy)=(one+one/TAN(half*pi*q1))-d2
               ELSEIF(q1==zero)THEN
                  ksi(ix,iy)=(one+one/TAN(half*pi*p1))-d1
                  eta(ix,iy)=zero
               ELSE
                  ksi(ix,iy)=(SQRT(4.*d1**2*(one-f1/d1**2-f2/d2**2)
     $                 - (f1-f2)**2/d2**2)-(two*d1+d1/d2**2*(f1-f2)))
     $                 /(two*(one+d1**2/d2**2))
                  eta(ix,iy)=ksi(ix,iy)*d1/d2 + half*(f1-f2)/d2
               ENDIF
               ksi(ix,iy)=ksi(ix,iy)*SIGN(one,x(ix,iy)-half)
               eta(ix,iy)=eta(ix,iy)*SIGN(one,y(ix,iy)-half)
            ENDDO
         ENDDO
      CASE("source_sink")
         ksi = lx*x
         eta=0.5*(y+(1.3980 + 0.347*TANH(10.0*(ksi - (lx/2.0-0*y/4)) )))
      CASE("nimrod")
         ksi = x-half
         eta = y-half
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
      USE heat_mod
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
      USE heat_mod
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
