c-----------------------------------------------------------------------
c     file MAST-init10.f.
c     Initial condition setup for MAST with 10 variables
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. MASTinit_mod.
c     1. MASTinit_equil.
c     2. Btor.
c     3. ff_Bz.
c     4. ff_jz.
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
c     n. physics_schur
c     o. physics_dealloc
c     p. physics_main.
c-----------------------------------------------------------------------
c     subprogram 0. MASTinit_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE MASTinit_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL :: cylinder=.FALSE.

      CHARACTER(16) :: init_type="."
      INTEGER :: cyl_fac=0
      REAL(r8), PARAMETER :: rpost=0.
      REAL(r8) :: lambda=1.,beta0=0.,beta_e=0.,Bguide=0.,lx=0.,ly=0.,
     $     epsilon=0.,gr_curve=0.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. MASTinit_equil.
c     computes equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE MASTinit_equil(x,y,u)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), PARAMETER :: Rcen=0.0_r8,Zcen=1.0_r8
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: Bz1i,Bz2i,jz1i,jz2i,
     $     rminor
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      SELECT CASE(init_type)
c      CASE("cartFF")
      CASE("cylFF", "cartFF")
         u(1,:,:)=LOG(one)
         u(2,:,:)=0.0_r8

         rminor = SQRT((x-Zcen)**2 + (y-Rcen)**2)
         CALL ff_Bz(rminor,y,lambda,epsilon,Bguide,Bz1i)
         CALL ff_jz(rminor,lambda,epsilon,jz1i)

         rminor = SQRT((x+Zcen)**2 + (y-Rcen)**2)
         CALL ff_Bz(rminor,y,lambda,epsilon,Bguide,Bz2i)
         CALL ff_jz(rminor,lambda,epsilon,jz2i)
      
         u(3,:,:)=-(Bz1i + Bz2i)         
         u(7,:,:)=jz1i+jz2i
         u(8,:,:)=beta0*0.25
         u(9,:,:)=beta0*0.25
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE MASTinit_equil
c-----------------------------------------------------------------------
c     subprogram 2. Btor.
c     prescribes toroidal guide field.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION Btor(Bguide,r) RESULT(Bz)

      REAL(r8), INTENT(IN) :: Bguide
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: r
      REAL(r8), DIMENSION(SIZE(r,1),SIZE(r,2)) :: Bz
c-----------------------------------------------------------------------
c     compute Bz.
c-----------------------------------------------------------------------           
      SELECT CASE(init_type)
      CASE("cylFF")
         Bz = Bguide*0.85_r8/r
      CASE("cartFF") 
         Bz = Bguide
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION Btor
c-----------------------------------------------------------------------
c     subprogram 3. ff_Bz.
c     computes Bz for one force free plasmoid.
c     The background field is removed from each plasmoid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ff_Bz(rminor,rmajor,w,epsilon,Bguide,Bz)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rminor,rmajor
      REAL(r8), INTENT(IN) :: w,epsilon,Bguide
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: Bz
c-----------------------------------------------------------------------
c     compute Bz.
c-----------------------------------------------------------------------
      WHERE(rminor <= w)
         Bz(:,:)=SQRT(Btor(Bguide,rmajor)**2 
     $        + epsilon**2*ABS(47._r8*w**2/360._r8
     $        - rminor**2/2._r8 + 3._r8*rminor**4/(4._r8*w**2)
     $        - 5._r8*rminor**6/(9._r8*w**4) 
     $        + 5._r8*rminor**8/(24._r8*w**6)
     $        - rminor**10/(30._r8*w**8))) - Btor(Bguide,rmajor)
      ELSEWHERE
         Bz(:,:)=0.
      END WHERE
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ff_Bz
c-----------------------------------------------------------------------
c     subprogram 4. ff_jz.
c     computes jz for one force free plasmoid.
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ff_jz(r,w,epsilon,jz)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: r
      REAL(r8), INTENT(IN) :: w,epsilon
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: jz
c-----------------------------------------------------------------------
c     compute jz
c-----------------------------------------------------------------------
      WHERE(r <= w)
         jz(:,:)=epsilon*(1._r8 - r**2/w**2)**2
      ELSEWHERE
         jz(:,:)=0._r8
      END WHERE
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ff_jz
      END MODULE MASTinit_mod
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
      USE MASTinit_mod
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
      NAMELIST/MASTinit_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,init_type,cylinder,beta0,
     $     beta_e,lx,ly,lambda,Bguide,epsilon
c-----------------------------------------------------------------------
c     Sample namelist.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=MASTinit_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="MASTinit"

      nqty=9
      nqty_schur=0
      nstep=1

      SELECT CASE(init_type)
      CASE("cylFF")
         cylinder=.TRUE.
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
      USE MASTinit_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(2)=.TRUE.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(cylinder,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta_e,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Bguide,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      IF(cylinder)cyl_fac=1

      SELECT CASE(init_type)
      CASE("cartFF","cylFF")
c         beta_e=0
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
      USE MASTinit_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u
c-----------------------------------------------------------------------
c     MAST initial conditions.
c-----------------------------------------------------------------------
      CALL MASTinit_equil(x,y,u)

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
      USE MASTinit_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
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
      CASE("cartFF","cylFF")
         top%static(2)=.TRUE.
         bottom%bc_type(2)="zeroflux"
         left%bc_type(2)="zeroflux"
         right%static(2)=.TRUE.
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
      USE MASTinit_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      REAL(r8) :: Bv
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
      Bv=0.06_r8 !Specify vertical field.
 
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("cartFF","cylFF")
         SELECT CASE(lrtb)
         CASE("top","right")
c            c(2,:,:)=u(2,:,:)-0.5_r8*Bv*y
             c(2,:,:)=u(2,:,:)
c         CASE("left","bottom")
c            c(2,:,:)=u(2,:,:)-0.5_r8*Bv*y
c             c(2,:,:)=u(2,:,:)
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
      USE MASTinit_mod
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
c-----------------------------------------------------------------------
c     boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("cartFF","cylFF")
         SELECT CASE(lrtb)
         CASE("top","right")
            c_u(2,2,:,:)=one
c         CASE("left","right")
c            c_u(2,2,:,:)=one
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
      USE MASTinit_mod
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
      CASE("cartFF","cylFF")
         SELECT CASE(lrtb)
         CASE("top","bottom","left","right")
            mass(1,1,:,:)=one
            mass(3,3,:,:)=one
            mass(4,4,:,:)=one
            mass(5,5,:,:)=one
            mass(6,6,:,:)=one
            mass(7,7,:,:)=one
            mass(8,8,:,:)=one
            mass(9,9,:,:)=one
c            mass(10,10,:,:)=one
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
      USE MASTinit_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_faci
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: b1,b2,j3
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c----------Cyl to cartesian relationships:------------------------------
c     1: z --> x
c     2: r --> y
c     3: phi --> z
c-----------------------------------------------------------------------
      r_faci=y**(-cyl_fac)
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      j3 = u(7,:,:)
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     poloidal magnetic flux equation.
c-----------------------------------------------------------------------
      s(2,:,:)=j3
      fx(2,:,:)=b2
      fy(2,:,:)=-b1
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     out-of-plane electron momentum (current) equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     total pressure equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     in plane currents for hyperresistive term in 3.
c-----------------------------------------------------------------------
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
      USE MASTinit_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_faci
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
c     cylindrical to cartesian relationships:
c     1: z --> x
c     2: r --> y
c     3: phi --> z
c-----------------------------------------------------------------------
      r_faci=y**(-cyl_fac)
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     poloidal magnetic flux equation.
c-----------------------------------------------------------------------
      s_u(2,7,:,:)=one
      fx_ux(2,2,:,:)=one
      fy_u(2,2,:,:)=cyl_fac*r_faci
      fy_uy(2,2,:,:)=one
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------  
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------  
c-----------------------------------------------------------------------
c     out-of-plane electron momentum (current) equation.
c-----------------------------------------------------------------------  
c-----------------------------------------------------------------------
c     total pressure equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     in plane currents for hyperresistive term in 3.
c-----------------------------------------------------------------------
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
      USE MASTinit_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      mass(1,1,:,:)=one
      mass(3,3,:,:)=one
      mass(4,4,:,:)=one
      mass(5,5,:,:)=one
      mass(6,6,:,:)=one
      mass(7,7,:,:)=one
      mass(8,8,:,:)=one
      mass(9,9,:,:)=one
c      mass(10,10,:,:)=one
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
      USE MASTinit_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,phi
      REAL(r8) :: rcurve
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("cartFF")
         rcurve=0.25_r8
c          ksi = lx*x/two
          ksi=lx*(10.*x**3 + x)/(10. + one)
c         ksi=SIGN(half*lx*((two*x - one)**2 + rcurve*ABS(two*x - one))
c     $        /(one + rcurve),two*x - one) + rpost
c          phi = ly*y/two
          phi = ly*(y**2 + y)/(one + one)
c         phi=SIGN(half*ly*((two*y - one)**2 + gr_curve*ABS(two*y - one))
c     $        /(one + gr_curve),two*y - one) 
      CASE("cylFF")
         rcurve=0.25_r8         ! Was 0.1 for R=1.1
         ksi=SIGN(half*ly*((two*x - one)**2 + gr_curve*ABS(two*x - one))
     $        /(one + gr_curve),two*x - one)
         phi=lx*y + rpost
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
      USE MASTinit_mod
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
c     deallocates.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_dealloc
      USE MASTinit_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     deallocate appropriate arrays.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      END SELECT
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
