c-----------------------------------------------------------------------
c     file wave.f.
c     contains specifications for sound waves.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. wave_mod.
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
c     subprogram 0. wave_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE wave_mod
      USE local_mod
      IMPLICIT NONE

      INTEGER :: wave_knx=1,wave_kny=1
      REAL(r8) :: wave_t0=0,omega,kx,ky
      REAL(r8), DIMENSION(4) :: u0fac

      END MODULE wave_mod
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
      USE wave_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: physics_type
      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty,nqty_schur
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep,exit_flag
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax

      INTEGER :: myios
      REAL(r8) :: tfac
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/wave_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,du_diagnose,wave_knx,wave_kny,wave_t0
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=wave_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      nqty=3
      nqty_schur=1
      physics_type="wave"
c-----------------------------------------------------------------------
c     rescale dt-related input.
c-----------------------------------------------------------------------
      kx=twopi*wave_knx
      ky=twopi*wave_kny
      omega=SQRT(kx**2+ky**2)
      tfac=twopi/omega
      dt=dt*tfac
      dtmax=dtmax*tfac
      tmax=tmax*tfac
      wave_t0=wave_t0*tfac
c-----------------------------------------------------------------------
c     compute u0fac.
c-----------------------------------------------------------------------
      IF(xperiodic)THEN
         u0fac(1)=1
         u0fac(2)=0
      ELSE
         u0fac(1)=half
         u0fac(2)=half
      ENDIF
      IF(yperiodic)THEN
         u0fac(3)=1
         u0fac(4)=0
      ELSE
         u0fac(3)=half
         u0fac(4)=half
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_input
c-----------------------------------------------------------------------
c     subprogram b. physics_init_parameters
c     sets up special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_init_parameters(static,ground,adapt_qty)
      USE wave_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     broadcast integer input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(wave_knx,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(wave_kny,1,MPI_INTEGER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast local module variables.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(kx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ky,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(omega,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(wave_t0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(u0fac,SIZE(u0fac),MPI_DOUBLE_PRECISION,0,comm,ierr)
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
      USE wave_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: iqty,i
      COMPLEX(r8), PARAMETER :: ifac=(0,1)
      COMPLEX(r8) :: tfac
      COMPLEX(r8), DIMENSION(SIZE(xpi,1),SIZE(xpi,2)) :: expfacx,expfacy
c-----------------------------------------------------------------------
c     compute expfac.
c-----------------------------------------------------------------------
      expfacx(:,:)=EXP(ifac*(kx*xpi))
      expfacy(:,:)=EXP(ifac*(ky*ypi))
      tfac=EXP(-ifac*omega*wave_t0)
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      ui(1,:,:)=kx*(u0fac(1)*expfacx-u0fac(2)*CONJG(expfacx))
     $        *(u0fac(3)*expfacy+u0fac(4)*CONJG(expfacy))*tfac
      ui(2,:,:)=ky*(u0fac(1)*expfacx+u0fac(2)*CONJG(expfacx))
     $        *(u0fac(3)*expfacy-u0fac(4)*CONJG(expfacy))*tfac
      ui(3,:,:)=omega*(u0fac(1)*expfacx+u0fac(2)*CONJG(expfacx))
     $        *(u0fac(3)*expfacy+u0fac(4)*CONJG(expfacy))*tfac
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
      USE wave_mod
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
      USE wave_mod
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
      USE wave_mod
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
      USE wave_mod
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
c-----------------------------------------------------------------------
c     set nonzero values.
c-----------------------------------------------------------------------
      SELECT CASE(lrtb)
      CASE("left","right")
         mass(1,1,:,:) = one
         mass_x(2,2,:,:) = one
         mass_x(3,3,:,:) = one
      CASE("top","bottom")
         mass_y(1,1,:,:) = one
         mass(2,2,:,:) = one
         mass_y(3,3,:,:) = one
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
      SUBROUTINE physics_rhs(t,x,y,u,ux,uy,fx,fy,s,first)
      USE wave_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     wave equation.
c-----------------------------------------------------------------------
      fx(1,:,:)=u(3,:,:)
      fy(2,:,:)=u(3,:,:)
      fx(3,:,:)=u(1,:,:)
      fy(3,:,:)=u(2,:,:)
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
      USE wave_mod
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
c     ideal MHD formulation.
c-----------------------------------------------------------------------
      fx_u(1,3,:,:)=1
      fy_u(2,3,:,:)=1
      fx_u(3,1,:,:)=1
      fy_u(3,2,:,:)=1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_drdu
c-----------------------------------------------------------------------
c     subprogram l. physics_mass.
c     computes contributions to the mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_mass(x,y,mass,mass_x,mass_y)
      USE wave_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify mass matrices.
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
      USE wave_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta
c-----------------------------------------------------------------------
c     y positions.
c-----------------------------------------------------------------------
      ksi=x
      eta=y
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
      USE wave_mod
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
c     compute nonzero terms.
c-----------------------------------------------------------------------
      fx_ux(1,1,:,:)=-hfac**2
      fy_uy(1,1,:,:)=-hfac**2
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
      USE wave_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_dealloc
c-----------------------------------------------------------------------
c     subprogram q. physics_main.
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
