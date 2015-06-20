c-----------------------------------------------------------------------
c     file cmhd.f.
c     contains specifications for compressible MHD model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. cmhd_mod.
c     1. cmhd_equil.
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
c     o. physics_dealloc
c     p. physics_main.
c-----------------------------------------------------------------------
c     subprogram 0. cmhd_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE cmhd_mod
      USE extra_mod
      USE transport_mod
      IMPLICIT NONE

      LOGICAL :: source=.FALSE.,cylinder=.FALSE.,if_kinvisc=.FALSE.,
     $     flux_inflow=.FALSE.

      CHARACTER(16) :: init_type=".",equil_type=".",equilfile=".",
     $     coilfile=".",interp="bicube",eta_case=".",kappa_case="."
      INTEGER :: cyl_fac=0
      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,qe=1.602e-19,
     $     me=9.109e-31,mi=3.345e-27,mp=1.673e-27,ep0=8.854e-12,
     $     mu0=4.e-7*pi,chod_const=0.1,kappa_min=0.,kappa_max=1.e8

      REAL(r8) :: eta=0.,r_eta=1.e10,etavac=1.,v_chod_norm=1.,
     $     etac_norm=1.,etas_norm=1.,mu=0.,kappa_par=0.,kappa_perp=0.,
     $     ddiff=0.,beta0=0.,rhomin=0.,pmin=0.,Tmin=1.,n0=1.e20,
     $     b0=1.,L0=1.,bx=0.,lx=0.,ly=0.,lambda_psi=0.,zmin=0.,
     $     kx=0.,ky=0.,ksq=0.,epsilon=0.,gamma_fac=1.,ke_norm=1.,
     $     xe_norm=1.,ki_norm=1.,xi_norm=1.,gr_curve=1.,nu=0.,
     $     eps1=0.,eps2=0.,alpha=0.,b0_mt=1., n0_mt=0.3,n_bgr=0.1 

      TYPE(coil_type) :: coils
      TYPE(bicube_type) :: equil_bc
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: equil,equilxy

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. cmhd_equil.
c     sets equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cmhd_equil(x,y,u,ux,uy,deriv)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: deriv

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: coshy_psi,
     $     sinhy_psi,coshx_psi,sinhx_psi, bfunc
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: eq
c-----------------------------------------------------------------------
c     compute cosh.
c-----------------------------------------------------------------------
      coshy_psi=COSH(y/lambda_psi)
      coshx_psi=COSH(x/lambda_psi)
c-----------------------------------------------------------------------
c        compute sinh.
c-----------------------------------------------------------------------
      sinhx_psi=SINH(x/lambda_psi)
      sinhy_psi=SINH(y/lambda_psi)
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      
      SELECT CASE(init_type)
c Elena asymptotic magnetotail equilibria
c Elena in sel.in add b0 - asymptotic mf outside current sheet        
      CASE("MAGNETOTAIL")
         bfunc = EXP(eps1*(x/lambda_psi
     $        + alpha/eps2*(one + TANH(eps2/lambda_psi*(x - 0.7*lx)))))
         u(1,:,:) = n0_mt/(bfunc*COSH(y/(lambda_psi*bfunc)))**2 + n_bgr
         u(2,:,:) = lambda_psi*b0_mt*LOG(bfunc*COSH(y/lambda_psi/bfunc))
         u(3,:,:) = b0_mt**2/(two*n0_mt)*u(1,:,:)
c Elena               
      CASE("GEM","GEM_open")
         u(1,:,:)=one/coshy_psi**2 + .2_r8
         u(2,:,:)=-lambda_psi*LOG(coshy_psi)
         u(3,:,:)=beta0*u(1,:,:)
         u(6,:,:)=-one/(lambda_psi*coshy_psi**2)
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=-two*TANH(y/lambda_psi)/(lambda_psi*coshy_psi**2)
         uy(2,:,:)=-TANH(y/lambda_psi)
         uy(3,:,:)=beta0*uy(1,:,:)
         uy(6,:,:)=-u(6,:,:)*two*TANH(y/lambda_psi)/lambda_psi
      CASE("frc","frc_rm")
c-----------------------------------------------------------------------
c     eq(1) is flux, eq(2) is pressure, and eq(3) is jphi
c-----------------------------------------------------------------------
         CALL extra_interp(interp,x,y,equil_bc,equil,equilxy,eq)
c-----------------------------------------------------------------------
c     apply cosine smoothing to floor pressure (pmin).
c-----------------------------------------------------------------------
         WHERE(eq(2,:,:) < pmin)
            u(3,:,:)=pmin
         ELSEWHERE(eq(2,:,:) >= pmin .AND. eq(2,:,:) < two*pmin)
            u(3,:,:)=pmin 
     $           + half*pmin*(one - COS(pi*(eq(2,:,:)-pmin)/pmin))
         ELSEWHERE
            u(3,:,:)=eq(2,:,:)
         ENDWHERE
c-----------------------------------------------------------------------
c     apply cosine smoothing to floor density (rhomin).
c-----------------------------------------------------------------------
         u(1,:,:) = SQRT(eq(2,:,:))
         WHERE(eq(2,:,:) <= 0)u(1,:,:)=0
         WHERE(u(1,:,:) < rhomin)
            u(1,:,:) = rhomin
         ELSEWHERE(u(1,:,:) >= rhomin .AND. u(1,:,:) < two*rhomin)
            u(1,:,:) = rhomin
     $           + half*rhomin*(one - COS(pi*(u(1,:,:)-rhomin)/rhomin))
         ENDWHERE

         u(2,:,:) = -eq(1,:,:)/y
         WHERE(y == 0)u(2,:,:) = zero
         u(6,:,:) = eq(3,:,:)
      CASE("alfven")
         u(1,:,:) = one
         u(3,:,:) = one
      CASE("islands")
         u(1,:,:) = rhomin + b0**2*COS(pi*x)**2/beta0 
     $        *(0.625_r8*SIN(twopi*y)**2 + half*COS(twopi*y)**2
     $        - half*(TANH(y/lambda_psi)*COS(twopi*y)
     $        + SIN(twopi*y)/(twopi*lambda_psi
     $        *COSH(y/lambda_psi)**2))**2)
         u(2,:,:) = b0*TANH(y/lambda_psi)*COS(pi*x)*SIN(twopi*y)/twopi
         u(3,:,:) = beta0*u(1,:,:)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cmhd_equil
      END MODULE cmhd_mod
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
      USE cmhd_mod
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
      NAMELIST/cmhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,init_type,source,cylinder,eta,eta_case,
     $     r_eta,etavac,mu,nu,if_kinvisc,kappa_par,kappa_perp,
     $     ddiff,beta0,rhomin,pmin,b0,n0,L0,bx,lx,ly,lambda_psi,
     $     zmin,epsilon,equil_type,equilfile,interp,coilfile,
     $     flux_inflow,kappa_case,gr_curve,
     $     eps1, eps2, alpha, b0_mt, n0_mt, n_bgr
c-----------------------------------------------------------------------
c     Sample namelist.
c-----------------------------------------------------------------------
ccc&cmhd_list
ccc
ccc	cylinder=f
ccc
ccc	init_type="GEM"
ccc	source=f
ccc	equil_type="."
ccc	equilfile="."
ccc     interp="."
ccc	coilfile="."
ccc
ccc	lx=25.6
ccc	ly=12.8
ccc     lambda_psi=0.5
ccc	zmin=0.
ccc
ccc	beta0=0.5
ccc
ccc	epsilon=0.1
ccc
ccc	eta=2.e-3
ccc	eta_case="."
ccc	r_eta=0.
ccc	etavac=0.
ccc
ccc	mu=1.e-3
ccc	nu=0.
ccc	if_kinvisc=f
ccc     kappa_case="."
ccc	kappa_par=2.e-3
ccc	kappa_perp=2.e-3
ccc	ddiff=0.
ccc
ccc	rhomin=0.
ccc	pmin=0.
ccc     bx=0.
ccc
ccc	b0=1.
ccc	n0=1.
ccc	L0=1.
ccc/
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=cmhd_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="cmhd"
      nqty=6
      nqty_schur=0
c-----------------------------------------------------------------------
c     set cylinder to true.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("frc","frc_rm")
         cylinder=.TRUE.
      CASE("islands")
         cylinder=.FALSE.
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
      USE cmhd_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty

      REAL(r8) :: tnorm
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(6)=.TRUE.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(equil_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(equilfile,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(interp,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(coilfile,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(eta_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(kappa_case,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(cylinder,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(if_kinvisc,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(flux_inflow,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(r_eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(etavac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_par,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_perp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ddiff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rhomin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(L0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(zmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eps1,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eps2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0_mt,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0_mt,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n_bgr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      IF(cylinder)cyl_fac=1
c-----------------------------------------------------------------------
c     compute normalizations in MKS units.
c     v_chod_norm normalizes the ve/vs term in Chodura resistivity.
c     etac_norm absorbs the constants in front of Chodura resistivity.
c     etas_norm absorbs the constants in front of Spitzer resistivity,
c     with electron temperature being half of total temperature.
c
c     for Braginskii thermal conduction:
c     - ke/ki_norm absorb constants in front of Braginskii
c     perpendicular and parallel thermal conduction.
c     - xe/xi_norm absorb constants in front of omega*tau. 
c-----------------------------------------------------------------------
      tnorm=L0*SQRT(mu0*n0*mi)/b0
      v_chod_norm=SQRT(mi/(mu0*n0))/(qe*L0)
      etac_norm=me/(qe*L0*b0*SQRT(ep0*mu0))
      etas_norm=5.e-5*17.*tnorm*(two*n0*mu0*qe)**1.5/(mu0*L0**2*b0**3)
     $     /half**1.5
      ke_norm = 3.56e21*(b0**2/(2*n0*mu0*qe))**2.5*tnorm/(2*n0*L0**2)
      xe_norm = 3.56e21*(b0**2/(2*n0*mu0*qe))**1.5*b0/n0
      ki_norm = 1.97e-7/SQRT(mi*mp)*(b0**2/(2*n0*mu0*qe))**2.5
     $     *tnorm/(2*n0*L0**2)
      xi_norm = 1.97e-7/SQRT(mi*mp)*(b0**2/(2*n0*mu0*qe))**1.5*b0/n0

      SELECT CASE(init_type)
      CASE("GEM","GEM_open")
         lambda_psi=.5_r8
         lx=25.6_r8
         ly=12.8_r8
         beta0=half
         kx=twopi/lx
         ky=pi/ly
      CASE("alfven")
         kx=twopi/lx
      CASE("frc","frc_rm")
         SELECT CASE(equil_type)
         CASE("gm")
            CALL extra_read_marklin(equilfile,interp,1._r8,
     $           equil_bc,equil,equilxy)
         CASE("db")
            CALL extra_read_barnes(equilfile,interp,
     $           equil_bc,equil,equilxy)
         END SELECT
         CALL extra_coilalloc(coilfile,coils)
c-----------------------------------------------------------------------
c     normalize the coil times, voltages, and positions.
c     e.g. tnorm*t0=t ==> tnorm=t/t0
c-----------------------------------------------------------------------
         coils%tfire=coils%tfire/tnorm
         coils%tqtr=coils%tqtr/tnorm
         coils%tcrow=coils%tcrow/tnorm
         coils%vc0=coils%vc0/(L0**2*b0/tnorm)
         coils%zb=coils%zb/L0
         coils%ze=coils%ze/L0
c-----------------------------------------------------------------------
c     set rhomin and pmin to default values if not set by user.
c-----------------------------------------------------------------------
         IF(pmin==0)pmin=5.e-3
         IF(rhomin==0)rhomin=SQRT(pmin)
         Tmin=pmin/rhomin
      END SELECT
      gamma_fac=gamma/(gamma-one)
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
      SUBROUTINE physics_init(x,y,u)
      USE cmhd_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: cosx,cosy
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     cmhd initial conditions.
c-----------------------------------------------------------------------
      u=zero
      cosx=COS(kx*x)
      cosy=COS(ky*y)

      CALL cmhd_equil(x,y,u,ux,uy,.FALSE.)

      SELECT CASE(init_type)
      CASE("GEM")
         u(2,:,:) = u(2,:,:) + epsilon*cosx*cosy
      CASE("GEM_open")
         WHERE(x <= lx/4 .AND. y <= ly/2)
     $        u(2,:,:) = u(2,:,:) + epsilon*cosx*cosy
      CASE("alfven")
c-----------------------------------------------------------------------
c     this flux gives a perturbed By.
c     use "bx" input variable to add constant bx.
c-----------------------------------------------------------------------
         u(2,:,:) = u(2,:,:) + b0*epsilon/kx*cosx
         u(5,:,:) = u(5,:,:) + epsilon*SIN(kx*x)
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
      USE cmhd_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     set boundary conditions.
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
      CASE ("MAGNETOTAIL")
         edge_order=(/4,1,3,2/)
         top%bc_type(1)="natural"
         top%bc_type(6)="natural"
         top%static(3:5)=.TRUE.  
c symetric boundary         
         bottom%bc_type(1:4)="zeroflux"
         bottom%bc_type(6)="zeroflux"
         bottom%static(5)=.TRUE.
         
         left%bc_type(1)="natural"
         left%bc_type(6)="natural"
         left%static(3:5)=.TRUE.
         
         right%bc_type(1)="natural"
         right%bc_type(6)="natural"
         right%static(4:5)=.TRUE.
         
      CASE("GEM")
         top%bc_type(1)="natural"
         top%static(3:6)=.TRUE.

         bottom%static=.TRUE.
         left%static=.TRUE.
         right%static=.TRUE.
      CASE("GEM_open")
         top%static=.TRUE.
         right%static=.TRUE.
         left%static=.TRUE.
         bottom%static=.TRUE.
c-----------------------------------------------------------------------
c     "frc_rm" is an frc problem with b.c.'s like those used by Richard
c     Milroy in NIMROD simulations.
c-----------------------------------------------------------------------
      CASE("frc_rm")
         top%static=.TRUE.
         top%static(1:3)=.FALSE.

         bottom%static=.TRUE.

         left%static=.TRUE.
         left%static(1:3)=.FALSE.

         right%static=.TRUE.
         right%static(1:3)=.FALSE.
      CASE("frc")
         top%static=.TRUE.
         top%static(2)=.FALSE.

         bottom%static=.TRUE.

         left%static=.TRUE.
         left%static(2)=.FALSE.

         right%static=.TRUE.
         right%static(2)=.FALSE.
      CASE("islands")
         edge_order=(/2,3,1,4/)
         top%bc_type="robin"
         top%bc_type(1)="zeroflux"
         top%bc_type(4)="zeroflux"
         top%static=.TRUE.
         top%static(2)=.FALSE.

         bottom%bc_type="zeroflux"
         bottom%bc_type(5)="robin"
         bottom%static(5)=.TRUE.

         left%bc_type="zeroflux"
         left%bc_type(4)="robin"
         left%static(4)=.TRUE.

         right%bc_type="robin"
         right%bc_type(1)="zeroflux"
         right%bc_type(5)="zeroflux"
         right%static=.TRUE.
         right%static(2)=.FALSE.

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
      USE cmhd_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: volt,eta_local
      REAL(r8), DIMENSION(2,SIZE(x,1),SIZE(x,2)) :: vi
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
c-----------------------------------------------------------------------
c     boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("MAGNETOTAIL")
         SELECT CASE(lrtb)
         CASE("top")
            c(3,:,:)=uy(3,:,:)*u(1,:,:)-u(3,:,:)*uy(1,:,:) !grad_y(T)=0
            c(4,:,:)=u(4,:,:) !nVx=0
ccc            c(4,:,:)=uy(4,:,:)*u(1,:,:)-u(4,:,:)*uy(1,:,:) !grad_y(Vx)=0
            c(5,:,:)=u(5,:,:) ! nVy=0
ccc            c(6,:,:)=u(6,:,:) !jz=0
         CASE("bottom")
            c(5,:,:) = u(5,:,:) !nVy=0 
         CASE("left")
            c(3,:,:)=ux(3,:,:)*u(1,:,:)-u(3,:,:)*ux(1,:,:) !grad_x(T)=0
            c(4,:,:)=u(4,:,:) !nVx=0
            c(5,:,:)=ux(5,:,:)*u(1,:,:)-u(5,:,:)*ux(1,:,:) !grad_x(Vy)=0
         CASE("right")
            c(3,:,:)=ux(3,:,:)*u(1,:,:)-u(3,:,:)*ux(1,:,:) !grad_x(T)=0
            c(4,:,:)=u(4,:,:) !nVx=0
            c(5,:,:)=ux(5,:,:)*u(1,:,:)-u(5,:,:)*ux(1,:,:) !grad_x(Vy)=0
         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top")
            c(3,:,:)=uy(3,:,:)*u(1,:,:)-u(3,:,:)*uy(1,:,:)            
            c(4,:,:)=uy(4,:,:)*u(1,:,:)-u(4,:,:)*uy(1,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=uy(6,:,:)
         CASE("left","right")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=ux(2,:,:)
            c(3,:,:)=ux(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=ux(5,:,:)
            c(6,:,:)=ux(6,:,:)
         END SELECT
      CASE("GEM_open")
         SELECT CASE(lrtb)
         CASE("top")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=uyy(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=uy(5,:,:)*u(1,:,:)-uy(1,:,:)*u(5,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=uy(4,:,:)*u(1,:,:)-u(4,:,:)*uy(1,:,:)
            c(5,:,:)=u(5,:,:)
         CASE("left")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=ux(2,:,:)
            c(3,:,:)=ux(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=ux(5,:,:)*u(1,:,:)-u(5,:,:)*ux(1,:,:)
         CASE("right")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=uxx(2,:,:)
            c(3,:,:)=ux(3,:,:)
            c(4,:,:)=ux(4,:,:)*u(1,:,:)-ux(1,:,:)*u(4,:,:)
            c(5,:,:)=u(5,:,:)
         END SELECT
      CASE("frc_rm")
         CALL extra_coileval(t,x,coils,volt)

         DO i=1,2
            vi(i,:,:)=u(i+3,:,:)/u(1,:,:)
         ENDDO

         SELECT CASE(lrtb)
         CASE("top")
            c(2,:,:)=-volt/(twopi*y)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=-vi(2,:,:)*(u(2,:,:)/y + uy(2,:,:)) 
     $           + volt/(twopi*y)
            c(6,:,:)=u(6,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
         CASE("left","right")
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
         END SELECT
      CASE("frc")
         CALL extra_coileval(t,x,coils,volt)
         CALL transport_seteta(eta_case,ly-y,u(1,:,:),u(3,:,:),u(6,:,:),
     $        chod_const,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,
     $        etavac,eta_local)

         DO i=1,2
            vi(i,:,:)=u(i+3,:,:)/u(1,:,:)
         ENDDO

         SELECT CASE(lrtb)
         CASE("top")
            c(1,:,:)=u(1,:,:)-rhomin
            c(2,:,:)=-volt/(twopi*y)
            c(3,:,:)=u(3,:,:)-Tmin*u(1,:,:)
            c(4,:,:)=u(4,:,:)
c-----------------------------------------------------------------------
c           flux convection b.c.
c-----------------------------------------------------------------------
            IF(flux_inflow)THEN
               c(5,:,:)=-vi(2,:,:)*(u(2,:,:)/y + uy(2,:,:)) 
     $              + volt/(twopi*y)
               c(6,:,:)=u(6,:,:)
c-----------------------------------------------------------------------
c           flux diffusion b.c.
c-----------------------------------------------------------------------
            ELSE
               c(5,:,:)=u(5,:,:)
               c(6,:,:)=eta_local*u(6,:,:) + volt/(twopi*y)
            ENDIF
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
         CASE("left","right")
            c(1,:,:)=u(1,:,:)-rhomin
            c(3,:,:)=u(3,:,:)-Tmin*u(1,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
         END SELECT
      CASE("islands")
         c(3,:,:) = u(3,:,:) - beta0*u(1,:,:)
         SELECT CASE(lrtb)
         CASE("top")
            c(5,:,:) = u(5,:,:)
            c(6,:,:) = uy(6,:,:)
         CASE("bottom")
            c(5,:,:) = u(5,:,:)
         CASE("left")
            c(4,:,:) = u(4,:,:)
         CASE("right")
            c(4,:,:) = u(4,:,:)
            c(6,:,:) = ux(6,:,:)
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
      USE cmhd_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: eta_local,
     $     eta_u1,eta_u3,eta_u6
      REAL(r8), DIMENSION(2,SIZE(x,1),SIZE(x,2)) :: vi
      REAL(r8), DIMENSION(2,2,SIZE(x,1),SIZE(x,2)) :: vi_u
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
      CASE("MAGNETOTAIL")
         SELECT CASE(lrtb)
         CASE("top")
            c_u(3,1,:,:)=uy(3,:,:)
            c_u(3,3,:,:)=-uy(1,:,:)
            c_uy(3,1,:,:)=-u(3,:,:)
            c_uy(3,3,:,:)=u(1,:,:)

c            c_u(4,1,:,:)=uy(4,:,:)
c            c_u(4,4,:,:)=-uy(1,:,:)
c            c_uy(4,1,:,:)=-u(4,:,:)
c            c_uy(4,4,:,:)=u(1,:,:)
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
         CASE("bottom")
            c_u(5,5,:,:)=one
         CASE("left")
            c_u(3,1,:,:)=ux(3,:,:)
            c_u(3,3,:,:)=-ux(1,:,:)
            c_ux(3,1,:,:)=-u(3,:,:)
            c_ux(3,3,:,:)=u(1,:,:)

            c_u(4,4,:,:)=one

            c_u(5,1,:,:)=ux(5,:,:)
            c_u(5,5,:,:)=-ux(1,:,:)
            c_ux(5,1,:,:)=-u(5,:,:)
            c_ux(5,5,:,:)=u(1,:,:)
         CASE("right")
            c_u(3,1,:,:)=ux(3,:,:)
            c_u(3,3,:,:)=-ux(1,:,:)
            c_ux(3,1,:,:)=-u(3,:,:)
            c_ux(3,3,:,:)=u(1,:,:)            

            c_u(4,4,:,:)=one

            c_u(5,1,:,:)=ux(5,:,:)
            c_u(5,5,:,:)=-ux(1,:,:)
            c_ux(5,1,:,:)=-u(5,:,:)
            c_ux(5,5,:,:)=u(1,:,:)
         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top")
            c_u(3,1,:,:)=uy(3,:,:)
            c_u(3,3,:,:)=-uy(1,:,:)
            c_uy(3,1,:,:)=-u(3,:,:)
            c_uy(3,3,:,:)=u(1,:,:)

            c_u(4,1,:,:)=uy(4,:,:)
            c_u(4,4,:,:)=-uy(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)
            c_uy(4,4,:,:)=u(1,:,:)

            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_uy(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_uy(6,6,:,:)=one
         CASE("left","right")
            c_ux(1,1,:,:)=one
            c_ux(2,2,:,:)=one
            c_ux(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_ux(5,5,:,:)=one
            c_ux(6,6,:,:)=one
         END SELECT
      CASE("GEM_open")
         SELECT CASE(lrtb)
         CASE("top")
            c_uy(1,1,:,:)=one
            c_uyy(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_u(4,4,:,:)=one

            c_u(5,1,:,:)=uy(5,:,:)
            c_u(5,5,:,:)=-uy(1,:,:)
            c_uy(5,1,:,:)=-u(5,:,:)
            c_uy(5,5,:,:)=u(1,:,:)
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_uy(2,2,:,:)=one
            c_uy(3,3,:,:)=one

            c_u(4,1,:,:)=uy(4,:,:)
            c_u(4,4,:,:)=-uy(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)
            c_uy(4,4,:,:)=u(1,:,:)

            c_u(5,5,:,:)=one
         CASE("left")
            c_ux(1,1,:,:)=one
            c_ux(2,2,:,:)=one
            c_ux(3,3,:,:)=one
            c_u(4,4,:,:)=one

            c_u(5,1,:,:)=ux(5,:,:)
            c_u(5,5,:,:)=-ux(1,:,:)
            c_ux(5,1,:,:)=-u(5,:,:)
            c_ux(5,5,:,:)=u(1,:,:)
         CASE("right")
            c_ux(1,1,:,:)=one
            c_uxx(2,2,:,:)=one
            c_ux(3,3,:,:)=one

            c_u(4,1,:,:)=ux(4,:,:)
            c_u(4,4,:,:)=-ux(1,:,:)
            c_ux(4,1,:,:)=-u(4,:,:)
            c_ux(4,4,:,:)=u(1,:,:)

            c_u(5,5,:,:)=one
        END SELECT
      CASE("frc_rm")
         DO i=1,2
            vi(i,:,:)=u(i+3,:,:)/u(1,:,:)
            vi_u(i,1,:,:)=-u(i+3,:,:)/u(1,:,:)**2
            vi_u(i,2,:,:)=one/u(1,:,:)
         ENDDO

         SELECT CASE(lrtb)
         CASE("top")
            c_u(4,4,:,:)=one
            c_u(5,1,:,:)=-vi_u(2,1,:,:)*(u(2,:,:)/y + uy(2,:,:))
            c_u(5,2,:,:)=-vi(2,:,:)/y
            c_u(5,5,:,:)=-vi_u(2,2,:,:)*(u(2,:,:)/y + uy(2,:,:))
            c_uy(5,2,:,:)=-vi(2,:,:)
            c_u(6,6,:,:)=one
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
         CASE("left","right")
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
         END SELECT
      CASE("frc")
         CALL transport_seteta_u(eta_case,ly-y,u(1,:,:),u(3,:,:),
     $        u(6,:,:),chod_const,etas_norm,etac_norm,v_chod_norm,
     $        ly-r_eta,eta,etavac,eta_local,eta_u1,eta_u3,eta_u6)

         DO i=1,2
            vi(i,:,:)=u(i+3,:,:)/u(1,:,:)
            vi_u(i,1,:,:)=-u(i+3,:,:)/u(1,:,:)**2
            vi_u(i,2,:,:)=one/u(1,:,:)
         ENDDO

         SELECT CASE(lrtb)
         CASE("top")
            c_u(1,1,:,:)=one
            c_u(3,1,:,:)=-Tmin
            c_u(3,3,:,:)=one
            c_u(4,4,:,:)=one
c-----------------------------------------------------------------------
c           flux convection b.c.
c-----------------------------------------------------------------------
            IF(flux_inflow)THEN
               c_u(5,1,:,:)=-vi_u(2,1,:,:)*(u(2,:,:)/y + uy(2,:,:))
               c_u(5,2,:,:)=-vi(2,:,:)/y
               c_u(5,5,:,:)=-vi_u(2,2,:,:)*(u(2,:,:)/y + uy(2,:,:))
               c_uy(5,2,:,:)=-vi(2,:,:)
               c_u(6,6,:,:)=one
c-----------------------------------------------------------------------
c           flux diffusion b.c.
c-----------------------------------------------------------------------
            ELSE
               c_u(5,5,:,:)=one
               c_u(6,1,:,:)=u(6,:,:)*eta_u1
               c_u(6,3,:,:)=u(6,:,:)*eta_u3
               c_u(6,6,:,:)=(u(6,:,:)*eta_u6 + eta_local)
            ENDIF
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
         CASE("left","right")
            c_u(1,1,:,:)=one

            c_u(3,1,:,:)=-Tmin
            c_u(3,3,:,:)=one

            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one

         END SELECT
      CASE("islands")
         c_u(3,1,:,:) = -beta0
         c_u(3,3,:,:) = one
         SELECT CASE(lrtb)
         CASE("top")
            c_u(5,5,:,:) = one
            c_uy(6,6,:,:) = one
         CASE("bottom")
            c_u(5,5,:,:) = one
         CASE("left")
            c_u(4,4,:,:) = one
         CASE("right")
            c_u(4,4,:,:) = one
            c_ux(6,6,:,:) = one
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
      USE cmhd_mod
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
      CASE("MAGNETOTAIL")
         SELECT CASE(lrtb)
         CASE("top")
ccc            mass(1,1,:,:)=one
            mass(2,2,:,:)=one
ccc            mass(3,3,:,:)=one
         CASE("right")
ccc            mass(1,1,:,:)=one
            mass(2,2,:,:)=one
ccc            mass(3,3,:,:)=one   
         CASE("left")
ccc            mass(1,1,:,:)=one
            mass_x(2,2,:,:)=one
ccc            mass(3,3,:,:)=one   
         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top")
            DO iqty=1,SIZE(mass,1)
               mass(iqty,iqty,:,:)=one
            ENDDO
         END SELECT
      CASE("frc_rm")
         SELECT CASE(lrtb)
         CASE("top")
            DO iqty=1,3
               mass(iqty,iqty,:,:)=one
            ENDDO
         CASE("left","right")
            DO iqty=1,3
               mass(iqty,iqty,:,:)=one
            ENDDO
         END SELECT
      CASE("frc")
         SELECT CASE(lrtb)
         CASE("top","left","right")
            mass(2,2,:,:)=one
         END SELECT
      CASE("islands")
         SELECT CASE(lrtb)
         CASE("top","right")
            mass(2,2,:,:)=one
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
      USE cmhd_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Tix,Tiy,kperp,kfac,
     $     eta_local,r_fac,r_faci,visc,nil,b1,b2,Bsq,n_inv,nx_inv,ny_inv
      REAL(r8), DIMENSION(2,SIZE(u,2),SIZE(u,3)) :: vi,vix,viy
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: BdotT
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     cylindrical to cartesian relationships:
c     1: z --> x
c     2: r --> y
c     3: phi --> z
c-----------------------------------------------------------------------
      r_fac=y**cyl_fac
      r_faci=y**(-cyl_fac)
c-----------------------------------------------------------------------
c     inverse density and derivatives.
c-----------------------------------------------------------------------
      n_inv=one/u(1,:,:)
      nx_inv=-ux(1,:,:)*n_inv**2
      ny_inv=-uy(1,:,:)*n_inv**2
c-----------------------------------------------------------------------
c     viscosity.
c-----------------------------------------------------------------------
      IF(if_kinvisc)THEN
         visc=mu*u(1,:,:)
      ELSE
         visc=mu
      ENDIF
c-----------------------------------------------------------------------
c     resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta(eta_case,ly-y,u(1,:,:),u(3,:,:),u(6,:,:),
     $     chod_const,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,
     $     etavac,eta_local)
c-----------------------------------------------------------------------
c     temperature gradients.
c-----------------------------------------------------------------------
      Tix = ux(3,:,:)*n_inv + u(3,:,:)*nx_inv
      Tiy = uy(3,:,:)*n_inv + u(3,:,:)*ny_inv
c-----------------------------------------------------------------------
c     magnetic fields.
c-----------------------------------------------------------------------
      nil=0
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT(b1,b2,nil,Tix,Tiy,BdotT)
         CALL transport_kbrag(u(1,:,:),u(3,:,:),half,Bsq,ke_norm,
     $        ki_norm,xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac)
      CASE("anisotropic")
         CALL transport_BdotT(b1,b2,nil,Tix,Tiy,BdotT)
         CALL transport_setkaniso(kappa_par,kappa_perp,Bsq,kperp,kfac)
      CASE DEFAULT
         BdotT=0
         kperp=kappa_par
         kfac=0
      END SELECT
c-----------------------------------------------------------------------
c     velocities and their gradients.
c-----------------------------------------------------------------------
      DO i=1,2
         vi(i,:,:)=u(i+3,:,:)*n_inv
         vix(i,:,:)=ux(i+3,:,:)*n_inv + u(i+3,:,:)*nx_inv
         viy(i,:,:)=uy(i+3,:,:)*n_inv + u(i+3,:,:)*ny_inv
      ENDDO
c-----------------------------------------------------------------------
c     equations for density and -Az (cartesian) or -Aphi (cylindrical).
c-----------------------------------------------------------------------
      fx(1,:,:)=r_fac*(u(4,:,:) - ddiff*ux(1,:,:))
      fy(1,:,:)=r_fac*(u(5,:,:) - ddiff*uy(1,:,:))

      fx(2,:,:) = r_fac*nu*ux(6,:,:)
      fy(2,:,:) = r_fac*nu*uy(6,:,:)
      s(2,:,:) = r_fac*(vi(2,:,:)*(b1+bx) - vi(1,:,:)*b2
     $     + eta_local*u(6,:,:)) + cyl_fac*r_faci*nu*u(6,:,:)
c-----------------------------------------------------------------------
c     pressure equation.
c-----------------------------------------------------------------------
      fx(3,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi(1,:,:)
     $     - kfac*BdotT(1,:,:) - kperp*Tix)
      fy(3,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi(2,:,:)
     $     - kfac*BdotT(2,:,:) - kperp*Tiy)
      s(3,:,:)=r_fac*(vi(1,:,:)*ux(3,:,:) + vi(2,:,:)*uy(3,:,:)
     $     + eta_local*u(6,:,:)**2
     $     + visc*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 + viy(1,:,:)**2
     $     + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:))
     $     + nu*(ux(6,:,:)**2 + uy(6,:,:)**2))
     $     + cyl_fac*r_faci*(two*visc*vi(2,:,:)**2+nu*u(6,:,:)**2)
c-----------------------------------------------------------------------
c     momentum equations.
c-----------------------------------------------------------------------
      fx(4,:,:)=r_fac*(u(4,:,:)*vi(1,:,:)
     $     - two*visc*vix(1,:,:) + u(3,:,:))

      fy(4,:,:)=r_fac*(u(4,:,:)*vi(2,:,:) 
     $     - visc*(viy(1,:,:) + vix(2,:,:)))
      
      s(4,:,:)=-r_fac*u(6,:,:)*b2

      fx(5,:,:)=r_fac*(u(5,:,:)*vi(1,:,:)
     $     - visc*(vix(2,:,:) + viy(1,:,:)))
      
      fy(5,:,:)=r_fac*(u(5,:,:)*vi(2,:,:)
     $     - two*visc*viy(2,:,:))

      s(5,:,:)=r_fac*(u(6,:,:)*(b1+bx) - uy(3,:,:))
     $     - cyl_fac*two*visc*vi(2,:,:)*r_faci
c-----------------------------------------------------------------------
c     current.
c-----------------------------------------------------------------------
      fx(6,:,:)=b2
      fy(6,:,:)=-b1
      s(6,:,:)=u(6,:,:)
c-----------------------------------------------------------------------
c     initial equilibrium source term.
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("GEM","GEM_open")
            CALL cmhd_equil(x,y,u0,u0x,u0y,.TRUE.)
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
      USE cmhd_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: r_fac,r_faci,Tix,Tiy,
     $     Tix_un,Tiy_un,Ti_un,kperp,kfac,kpar_u1,kperp_bsq,kperp_u1,
     $     kperp_u3,kpar_u3,eta_local,visc,visc_rho,eta_u1,eta_u3,
     $     eta_u6,nil,b1,b2,Bsq,n_inv,nx_inv,ny_inv
      REAL(r8), DIMENSION(2,SIZE(u,2),SIZE(u,3)) :: vi,vix,viy,
     $     vi_un,vix_un,viy_un
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: nil3,BdotT,BdotT_Tx,
     $     BdotT_Ty,BdotT_b1,BdotT_b2
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
      r_fac=y**cyl_fac
      r_faci=y**(-cyl_fac)
c-----------------------------------------------------------------------
c     inverse density and derivatives.
c-----------------------------------------------------------------------
      n_inv=one/u(1,:,:)
      nx_inv=-ux(1,:,:)*n_inv**2
      ny_inv=-uy(1,:,:)*n_inv**2
c-----------------------------------------------------------------------
c     viscosity.
c-----------------------------------------------------------------------
      IF(if_kinvisc)THEN
         visc=mu*u(1,:,:)
         visc_rho=mu
      ELSE
         visc=mu
         visc_rho=zero
      ENDIF
c-----------------------------------------------------------------------
c     resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta_u(eta_case,ly-y,u(1,:,:),u(3,:,:),u(6,:,:),
     $     chod_const,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,
     $     etavac,eta_local,eta_u1,eta_u3,eta_u6) 
c-----------------------------------------------------------------------
c     temperature gradients and derivatives.
c-----------------------------------------------------------------------
      Tix = ux(3,:,:)*n_inv + u(3,:,:)*nx_inv
      Tiy = uy(3,:,:)*n_inv + u(3,:,:)*ny_inv
      Ti_un=-u(3,:,:)*n_inv**2
      Tix_un=-(two*u(3,:,:)*nx_inv + ux(3,:,:)*n_inv)*n_inv
      Tiy_un=-(two*u(3,:,:)*ny_inv + uy(3,:,:)*n_inv)*n_inv
c-----------------------------------------------------------------------
c     magnetic fields.
c-----------------------------------------------------------------------
      nil=0
      nil3=0
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT_u(b1,b2,nil,Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,nil3,BdotT_Tx,BdotT_Ty)
         CALL transport_kbrag_u(u(1,:,:),u(3,:,:),half,Bsq,ke_norm,
     $        ki_norm,xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac,
     $        kpar_u1,kpar_u3,kperp_u1,kperp_u3,kperp_bsq)
      CASE("anisotropic")
         kpar_u1=0
         kpar_u3=0
         kperp_u1=0
         kperp_u3=0
         CALL transport_BdotT_u(b1,b2,nil,Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,nil3,BdotT_Tx,BdotT_Ty)
         CALL transport_setkaniso_u(kappa_par,kappa_perp,Bsq,
     $        kperp,kfac,kperp_bsq)
      CASE DEFAULT
         BdotT=0
         BdotT_b1=0
         BdotT_b2=0
         BdotT_Tx=0
         BdotT_Ty=0
         kperp=kappa_par
         kfac=0
         kpar_u1=0
         kpar_u3=0
         kperp_u1=0
         kperp_u3=0
         kperp_bsq=0
      END SELECT
c-----------------------------------------------------------------------
c     velocities and their derivatives.
c-----------------------------------------------------------------------
      DO i=1,2
         vi(i,:,:)=u(i+3,:,:)*n_inv
         vix(i,:,:)=ux(i+3,:,:)*n_inv + u(i+3,:,:)*nx_inv
         viy(i,:,:)=uy(i+3,:,:)*n_inv + u(i+3,:,:)*ny_inv

         vi_un(i,:,:)=-vi(i,:,:)*n_inv
         vix_un(i,:,:) = -(vi(i,:,:)*nx_inv + vix(i,:,:)*n_inv)
         viy_un(i,:,:) = -(vi(i,:,:)*ny_inv + viy(i,:,:)*n_inv)
      ENDDO
c-----------------------------------------------------------------------
c     equations for density and -Az (cartesian) or -Aphi (cylindrical).
c-----------------------------------------------------------------------
      fx_u(1,4,:,:)=r_fac
      fy_u(1,5,:,:)=r_fac
      fx_ux(1,1,:,:)=-r_fac*ddiff
      fy_uy(1,1,:,:)=-r_fac*ddiff

      fx_ux(2,6,:,:) = r_fac*nu
      fy_uy(2,6,:,:) = r_fac*nu
      s_u(2,1,:,:)=r_fac*(vi_un(2,:,:)*(b1+bx) - vi_un(1,:,:)*b2
     $     + eta_u1*u(6,:,:))
      s_u(2,2,:,:)=-cyl_fac*vi(2,:,:)
      s_u(2,3,:,:)=r_fac*eta_u3*u(6,:,:)
      s_u(2,4,:,:)=-r_fac*n_inv*b2
      s_u(2,5,:,:)=r_fac*n_inv*(b1+bx)
      s_u(2,6,:,:)=r_fac*(eta_local + u(6,:,:)*eta_u6)
     $     + cyl_fac*r_faci*nu
      s_ux(2,2,:,:)=-r_fac*vi(1,:,:)
      s_uy(2,2,:,:)=-r_fac*vi(2,:,:)
c-----------------------------------------------------------------------
c     pressure equation.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi_un(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*Tix_un 
     $     - kfac*BdotT_Ty(1,:,:)*Tiy_un
     $     - (kpar_u1 - kperp_u1)*BdotT(1,:,:) - kperp_u1*Tix)
      fx_u(3,2,:,:)=cyl_fac*(kfac*BdotT_b1(1,:,:)
     $     - two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_u(3,3,:,:)=r_fac*(gamma_fac*vi(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*nx_inv
     $     - kfac*BdotT_Ty(1,:,:)*ny_inv
     $     - (kpar_u3 - kperp_u3)*BdotT(1,:,:) - kperp_u3*Tix)
      fx_u(3,4,:,:)=r_fac*gamma_fac*u(3,:,:)*n_inv

      fx_ux(3,1,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*Ti_un
      fx_ux(3,2,:,:)=r_fac*(-kfac*BdotT_b2(1,:,:)
     $     + two*b2*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_ux(3,3,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*n_inv

      fx_uy(3,1,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*Ti_un
      fx_uy(3,2,:,:)=r_fac*(kfac*BdotT_b1(1,:,:)
     $     - two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_uy(3,3,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*n_inv

      fy_u(3,1,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi_un(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*Tiy_un
     $     - kfac*BdotT_Tx(2,:,:)*Tix_un
     $     - (kpar_u1 - kperp_u1)*BdotT(2,:,:) - kperp_u1*Tiy)
      fy_u(3,2,:,:)=cyl_fac*(kfac*BdotT_b1(2,:,:)
     $     - two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(3,3,:,:)=r_fac*(gamma_fac*vi(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*ny_inv
     $     - kfac*BdotT_Tx(2,:,:)*nx_inv
     $     - (kpar_u3 - kperp_u3)*BdotT(2,:,:) - kperp_u3*Tiy)
      fy_u(3,5,:,:)=r_fac*gamma_fac*u(3,:,:)*n_inv

      fy_ux(3,1,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*Ti_un
      fy_ux(3,2,:,:)=r_fac*(-kfac*BdotT_b2(2,:,:)
     $     + two*b2*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_ux(3,3,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*n_inv

      fy_uy(3,1,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*Ti_un
      fy_uy(3,2,:,:)=r_fac*(kfac*BdotT_b1(2,:,:)
     $     - two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_uy(3,3,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*n_inv

      s_u(3,1,:,:)=r_fac*(vi_un(1,:,:)*ux(3,:,:) 
     $     + vi_un(2,:,:)*uy(3,:,:)
     $     + eta_u1*u(6,:,:)**2
     $     + two*visc*(two*vix(1,:,:)*vix_un(1,:,:) 
     $     + two*viy(2,:,:)*viy_un(2,:,:) 
     $     + viy(1,:,:)*viy_un(1,:,:) + vix(2,:,:)*vix_un(2,:,:)
     $     + viy(1,:,:)*vix_un(2,:,:) + vix(2,:,:)*viy_un(1,:,:))
     $     + visc_rho*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:))
     $     + cyl_fac*two*r_faci**2*(visc*two*vi(2,:,:)*vi_un(2,:,:)
     $     + visc_rho*vi(2,:,:)**2))
      s_u(3,3,:,:)=r_fac*eta_u3*u(6,:,:)**2
      s_u(3,4,:,:)=r_fac*(ux(3,:,:)*n_inv + two*visc*
     $     (two*vix(1,:,:)*nx_inv + (viy(1,:,:) + vix(2,:,:))*ny_inv))
      s_u(3,5,:,:)=r_fac*(uy(3,:,:)*n_inv + two*visc*
     $     (two*viy(2,:,:)*ny_inv + (vix(2,:,:) + viy(1,:,:))*nx_inv)
     $     + cyl_fac*4._r8*visc*r_faci**2*vi(2,:,:)*n_inv)
      s_u(3,6,:,:)=r_fac*(two*eta_local*u(6,:,:) + eta_u6*u(6,:,:)**2)
     $     + cyl_fac*r_faci*two*nu*u(6,:,:)

      s_ux(3,1,:,:)=r_fac*two*visc*(two*vix(1,:,:)*vi_un(1,:,:)  
     $     + vi_un(2,:,:)*(vix(2,:,:) + viy(1,:,:)))
      s_ux(3,3,:,:)=r_fac*vi(1,:,:)
      s_ux(3,4,:,:)=r_fac*visc*4._r8*vix(1,:,:)*n_inv
      s_ux(3,5,:,:)=r_fac*visc*two*n_inv*(vix(2,:,:) + viy(1,:,:))
      s_ux(3,6,:,:)=r_fac*two*nu*ux(6,:,:)

      s_uy(3,1,:,:)=r_fac*two*visc*(two*viy(2,:,:)*vi_un(2,:,:)  
     $     + vi_un(1,:,:)*(viy(1,:,:) + vix(2,:,:)))
      s_uy(3,3,:,:)=r_fac*vi(2,:,:)
      s_uy(3,4,:,:)=r_fac*visc*two*n_inv*(viy(1,:,:) + vix(2,:,:))
      s_uy(3,5,:,:)=r_fac*visc*4._r8*viy(2,:,:)*n_inv
      s_uy(3,6,:,:)=r_fac*two*nu*uy(6,:,:)
c-----------------------------------------------------------------------
c     momentum equations.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(1,:,:) 
     $     - two*visc*vix_un(1,:,:) - two*visc_rho*vix(1,:,:))
      fx_u(4,3,:,:)=r_fac
      fx_u(4,4,:,:)=r_fac*two*(vi(1,:,:) - visc*nx_inv)
      fx_ux(4,1,:,:)=-r_fac*visc*two*vi_un(1,:,:)
      fx_ux(4,4,:,:)=-r_fac*visc*two*n_inv

      fy_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(2,:,:) 
     $     - visc*(viy_un(1,:,:) + vix_un(2,:,:))
     $     - visc_rho*(viy(1,:,:) + vix(2,:,:)))
      fy_u(4,4,:,:)=r_fac*(vi(2,:,:) - visc*ny_inv)
      fy_u(4,5,:,:)=r_fac*(vi(1,:,:) - visc*nx_inv)
      fy_ux(4,1,:,:)=-r_fac*visc*vi_un(2,:,:)
      fy_ux(4,5,:,:)=-r_fac*visc*n_inv
      fy_uy(4,1,:,:)=-r_fac*visc*vi_un(1,:,:)
      fy_uy(4,4,:,:)=-r_fac*visc*n_inv

      s_u(4,6,:,:)=-r_fac*b2
      s_ux(4,2,:,:)=-r_fac*u(6,:,:)

      fx_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(1,:,:) 
     $     - visc*(vix_un(2,:,:) + viy_un(1,:,:))
     $     - visc_rho*(vix(2,:,:) + viy(1,:,:)))
      fx_u(5,4,:,:)=r_fac*(vi(2,:,:) - visc*ny_inv)
      fx_u(5,5,:,:)=r_fac*(vi(1,:,:) - visc*nx_inv)
      fx_ux(5,1,:,:)=-r_fac*visc*vi_un(2,:,:)
      fx_ux(5,5,:,:)=-r_fac*visc*n_inv
      fx_uy(5,1,:,:)=-r_fac*visc*vi_un(1,:,:)
      fx_uy(5,4,:,:)=-r_fac*visc*n_inv
      
      fy_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(2,:,:) 
     $     - two*visc*viy_un(2,:,:) - two*visc_rho*viy(2,:,:))
      fy_u(5,5,:,:)=r_fac*two*(vi(2,:,:) - visc*ny_inv)
      fy_uy(5,1,:,:)=-r_fac*visc*two*vi_un(2,:,:)
      fy_uy(5,5,:,:)=-r_fac*visc*two*n_inv

      s_u(5,1,:,:)=-cyl_fac*two*r_faci*(visc*vi_un(2,:,:)
     $     + visc_rho*vi(2,:,:))
      s_u(5,2,:,:)=-cyl_fac*u(6,:,:)
      s_u(5,5,:,:)=-cyl_fac*two*visc*r_faci*n_inv
      s_u(5,6,:,:)=r_fac*(b1+bx)
      s_uy(5,2,:,:)=-r_fac*u(6,:,:)
      s_uy(5,3,:,:)=-r_fac
c-----------------------------------------------------------------------
c     current.
c-----------------------------------------------------------------------
      fx_ux(6,2,:,:)=one
      fy_u(6,2,:,:)=cyl_fac*r_faci
      fy_uy(6,2,:,:)=one
      s_u(6,6,:,:)=one
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
      USE cmhd_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c        cylindrical to cartesian relationships:
c        1: z --> x
c        2: r --> y
c        3: phi --> z
c-----------------------------------------------------------------------
      r_fac=y**cyl_fac

      mass(1,1,:,:)=r_fac
      mass(2,2,:,:)=r_fac
      mass(3,3,:,:)=r_fac/(gamma-one)
      mass(4,4,:,:)=r_fac
      mass(5,5,:,:)=r_fac
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
      SUBROUTINE physics_grid(x,y,ksi,etag)
      USE cmhd_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,etag

      REAL(r8) :: fac
c-----------------------------------------------------------------------
c     set grid according to init_type.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("frc","alfven","frc_rm")
         IF(gr_curve == 0.)THEN
            ksi=lx*x+zmin
            etag=ly*y
         ELSE
            fac=(2.1_r8 - zmin)/lx
            ksi=((x - fac)**3 + gr_curve*x + fac**3)*lx/
     $           ((one - fac)**3 + gr_curve + fac**3) + zmin
            etag=ly*(y**2 + gr_curve*y)/(one + gr_curve)
         ENDIF
      CASE("MAGNETOTAIL")
         ksi=lx*half*(3._r8*x-one)
         etag=ly*half*(gr_curve*y**3 + y)/(one+gr_curve)
      CASE("GEM")
         ksi=lx*half*x
         etag=ly*half*(y**4+gr_curve*y)/(one+gr_curve)
      CASE("GEM_open")
         ksi=two*lx*x
         etag=two*ly*(y**2+gr_curve*y)/(one+gr_curve)
      CASE("islands")
         ksi=half*x
         etag=half*(gr_curve*y**4 + y)/(one+gr_curve)
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
      USE cmhd_mod
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
      USE cmhd_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     deallocate appropriate arrays.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("frc","frc_rm")
         CALL extra_equil_dealloc(interp,equil_bc,equil,equilxy)
      END SELECT

      SELECT CASE(init_type)
      CASE("frc","frc_rm")
         CALL extra_coildealloc(coils)
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