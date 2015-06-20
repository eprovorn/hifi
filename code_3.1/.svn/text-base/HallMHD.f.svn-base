c-----------------------------------------------------------------------
c     file HallMHD.f.
c     contains specifications for Hall MHD model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. HallMHD_mod.
c     1. HallMHD_equil.
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
c     subprogram 0. HallMHD_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE HallMHD_mod
      USE extra_mod
      USE transport_mod
      IMPLICIT NONE

      LOGICAL :: source=.FALSE.,cylinder=.FALSE.,
     $     flux_inflow=.FALSE.
      CHARACTER(16) :: init_type=".",equil_type=".",
     $     equilfile=".",coilfile=".",interp="bicube",eta_case="."
      INTEGER :: cyl_fac=0
      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,qe=1.602e-19,
     $     me=9.109e-31,mi=3.345e-27,ep0=8.854e-12,mu0=4.e-7*pi,
     $     chod_const=0.1
      REAL(r8) :: di=0.,nu=0.,eta=0.,r_eta=1.e10,etavac=1.,
     $     v_chod_norm=1.,etac_norm=1.,etas_norm=1.,mu=0.,kappa_par=0.,
     $     kappa_perp=0.,Dn=0.,beta0=0.,beta_e=0.,rhomin=0.,pmin=0.,
     $     pmax=1.,Tmin=1.,n0=1.,T0=1.,b0=1.,Bguide=0.,L0=1.,lx=0.,
     $     ly=0.,lambda_psi=0.,lambda_phi=0.,zmin=-1.,kx=0.,ky=0.,
     $     ksq=0.,epsilon=0.,alfven=0.,mach=0.,gamma_fac=1.,gr_curve=0.

      TYPE(coil_type) :: coils
      TYPE(bicube_type) :: equil_bc
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: equil,equilxy

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. HallMHD_equil.
c     computes equilibrium for magnetic reconnection.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHD_equil(x,y,u,ux,uy,deriv)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: deriv

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: coshy_psi,coshy_phi,
     $     sinhy_psi,sinhy_phi,coshx_psi,coshx_phi,sinhx_psi,sinhx_phi,
     $     r0
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: eq
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
      CASE("Harris")
         u(1,:,:)=1
         u(2,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(4,:,:)=-mach*sinhy_phi/coshy_phi
         u(7,:,:)=alfven/(lambda_psi*coshy_psi**2)
         u(8,:,:)=beta0
      CASE("Huba")
         u(1,:,:)=one/coshy_psi**2
         u(3,:,:)=-sinhy_psi/coshy_psi
         u(8,:,:)=half*u(1,:,:)
      CASE("GEM","GEM_open")
         u(1,:,:)=one/coshy_psi**2 + .2_r8
         u(2,:,:)=-lambda_psi*LOG(coshy_psi)
         u(3,:,:)=Bguide
         u(6,:,:)=-(beta0-beta_e)/beta0*di/(lambda_psi*coshy_psi**2)
         u(7,:,:)=-one/(lambda_psi*coshy_psi**2)
         u(8,:,:)=beta0*u(1,:,:)
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=-two*TANH(y/lambda_psi)/(lambda_psi*coshy_psi**2)
         uy(2,:,:)=-TANH(y/lambda_psi)
         uy(6,:,:)=-u(6,:,:)*two*TANH(y/lambda_psi)/lambda_psi
         uy(7,:,:)=-u(7,:,:)*two*TANH(y/lambda_psi)/lambda_psi
         uy(8,:,:)=beta0*uy(1,:,:)
      CASE("frc")
c-----------------------------------------------------------------------
c     eq(1) is flux, eq(2) is pressure, and eq(3) is jphi
c-----------------------------------------------------------------------
         CALL extra_interp(interp,x,y,equil_bc,equil,equilxy,eq)
c-----------------------------------------------------------------------
c     apply cosine smoothing to floor pressure (pmin).
c-----------------------------------------------------------------------
         WHERE(eq(2,:,:) < pmin)
            eq(2,:,:)=pmin
         ELSEWHERE(eq(2,:,:) >= pmin .AND. eq(2,:,:) < two*pmin)
            eq(2,:,:)=pmin 
     $           + half*pmin*(one - COS(pi*(eq(2,:,:)-pmin)/pmin))
         END WHERE
c-----------------------------------------------------------------------
c     apply cosine smoothing to floor density (rhomin).
c-----------------------------------------------------------------------
         u(1,:,:) = SQRT(eq(2,:,:))*b0**2/(two*mu0*n0*qe*T0)
         WHERE(u(1,:,:) < rhomin)
            u(1,:,:) = rhomin
         ELSEWHERE(u(1,:,:) >= rhomin .AND. u(1,:,:) < two*rhomin)
            u(1,:,:) = rhomin
     $           + half*rhomin*(one - COS(pi*(u(1,:,:)-rhomin)/rhomin))
         END WHERE

         u(2,:,:) = -eq(1,:,:)/y
         WHERE(y == 0)u(2,:,:) = zero
         u(6,:,:) = di*eq(3,:,:)*(one-beta_e/beta0)
         u(7,:,:) = eq(3,:,:)
         u(8,:,:) = eq(2,:,:)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHD_equil
      END MODULE HallMHD_mod
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
      USE HallMHD_mod
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
      NAMELIST/HallMHD_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,init_type,source,cylinder,di,eta,
     $     eta_case,r_eta,etavac,mu,nu,kappa_par,kappa_perp,Dn,beta0,
     $     beta_e,rhomin,pmin,pmax,n0,T0,b0,Bguide,L0,lx,ly,lambda_psi,
     $     lambda_phi,zmin,epsilon,mach,alfven,equil_type,equilfile,
     $     interp,coilfile,flux_inflow
c-----------------------------------------------------------------------
c     Sample namelist.
c-----------------------------------------------------------------------
c$$$&HallMHD_list
c$$$
c$$$	cylinder=t
c$$$
c$$$	init_type="frc"
c$$$	source=f
c$$$	equil_type="gm"
c$$$	equilfile="rm_equil.data"
c$$$	coilfile="frc_coils.csv"
c$$$	flux_inflow=f
c$$$
c$$$	lx=21.43
c$$$	ly=1.
c$$$    lambda_psi=0.5
c$$$    lambda_phi=0.5
c$$$	zmin=-1.
c$$$
c$$$	beta0=1.
c$$$	beta_e=0.5
c$$$	epsilon=0.
c$$$
c$$$	di=.1
c$$$	nu=2.e-5
c$$$
c$$$	eta=2.01e-3
c$$$	eta_case="spitzer-chodura"
c$$$	r_eta=0.85
c$$$	etavac=100.
c$$$
c$$$	mu=1.00e-1
c$$$	kappa_par=20.0
c$$$	kappa_perp=4.02e-4
c$$$	Dn=2.00e-3
c$$$
c$$$	rhomin=5.e-2
c$$$	pmin=5.e-3
c$$$	pmax=1.57
c$$$
c$$$	b0=0.0199
c$$$	n0=7.44e19
c$$$	T0=16.56
c$$$
c$$$	Bguide=0.
c$$$	alfven=0.
c$$$	mach=0.
c$$$/
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=HallMHD_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="HallMHD"

      nqty=8
      nqty_schur=0

      SELECT CASE(init_type)
      CASE("frc")
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
      USE HallMHD_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty

      REAL(r8) :: tnorm
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(7)=.TRUE.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(equil_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(equilfile,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(interp,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(coilfile,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(eta_case,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(cylinder,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(flux_inflow,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(r_eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(etavac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_par,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_perp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Dn,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta_e,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rhomin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(T0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Bguide,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(L0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_phi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(zmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alfven,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      IF(cylinder)cyl_fac=1
c-----------------------------------------------------------------------
c     compute normalizations in MKS units.
c     v_chod_norm normalizes the ve/vs term in Chodura resistivity.
c     etac_norm absorbs the constants in front of Chodura resistivity.
c     etas_norm absorbs the constants in front of Spitzer resistivity.
c-----------------------------------------------------------------------
      tnorm=L0*SQRT(mu0*n0*mi)/b0
      v_chod_norm=one/(qe*L0)*SQRT(mi/(mu0*n0))
      etac_norm=me/(qe*L0*b0*SQRT(ep0*mu0))
      etas_norm=5.e-5*17.*tnorm*(two*n0*mu0*qe)**1.5/(mu0*L0**2*b0**3)
     $     /(beta_e/beta0)**1.5
      SELECT CASE(init_type)
      CASE("Harris")
         kx=twopi/lx
         ky=pi
      CASE("GEM","GEM_open")
         lambda_psi=.5_r8
         lx=25.6_r8
         ly=12.8_r8
         di=one
         beta0=half
         beta_e=one/12.
         kx=twopi/lx
         ky=pi/ly
      CASE("frc")
         beta0=one
         beta_e=half
         SELECT CASE(equil_type)
         CASE("gm")
            CALL extra_read_marklin(equilfile,interp,pmax,
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
      USE HallMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: cosx,cosy,
     $     sinx,siny,coshy_psi,tanhx_phi
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     HallMHD initial conditions.
c-----------------------------------------------------------------------
      CALL HallMHD_equil(x,y,u,ux,uy,.FALSE.)

      sinx=SIN(kx*x)
      siny=SIN(ky*y)
      cosx=COS(kx*x)
      cosy=COS(ky*y)
      coshy_psi=COSH(y/lambda_psi)
      tanhx_phi=TANH(x/lambda_phi)

      SELECT CASE(init_type)
      CASE("Harris")
         u(2,:,:)=u(2,:,:)-epsilon*cosx*cosy*alfven/ksq
         u(4,:,:)=u(4,:,:)-epsilon*cosx*siny*mach*ky/ksq
         u(5,:,:)=u(5,:,:)+epsilon*sinx*cosy*mach*kx/ksq
         u(7,:,:)=u(7,:,:)+epsilon*cosx*cosy*alfven
      CASE("Huba")
         u(1,:,:)=u(1,:,:)-half*(one-epsilon)*(one+tanhx_phi)
     $        /coshy_psi**2
         u(8,:,:)=u(1,:,:)/(one+epsilon-(one-epsilon)*tanhx_phi)
      CASE("GEM")
         u(2,:,:)=u(2,:,:)+epsilon*cosx*cosy
      CASE("GEM_open")
         WHERE(x <= lx/4 .AND. y <= ly/2)
     $        u(2,:,:)=u(2,:,:)+epsilon*cosx*cosy
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
      USE HallMHD_mod
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
      CASE("Harris")
         top%bc_type(1)="natural"
         top%static(4:6)=.TRUE.
         top%bc_type(7:8)="natural"

         bottom%bc_type(1)="natural"
         bottom%static(4:6)=.TRUE.
         bottom%bc_type(7:8)="natural"
      CASE("Huba")
         left%bc_type(7)="natural"
         right%bc_type(7)="natural"
         top%bc_type(7)="natural"
         bottom%bc_type(7)="natural"
      CASE("GEM")
         top%static(2:8)=.TRUE.
         bottom%static(2:8)=.TRUE.
         top%bc_type(1)="natural"
         bottom%bc_type(1)="natural"
      CASE("GEM_open")
         top%static=.TRUE.
         right%static=.TRUE.
         left%static=.TRUE.
         bottom%static=.TRUE.
      CASE("frc")
         top%static=.TRUE.
         top%static(2)=.FALSE.

         bottom%static=.TRUE.

         left%static=.TRUE.
         left%static(2)=.FALSE.

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
      USE HallMHD_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: volt,eta_local,jtot
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0

      jtot = SQRT((uy(3,:,:) + cyl_fac/y*u(3,:,:))**2 + ux(3,:,:)**2
     $     + u(7,:,:)**2)
      WHERE(y == 0)jtot = SQRT(uy(3,:,:)**2 + ux(3,:,:)**2
     $     + u(7,:,:)**2)

      CALL transport_seteta(eta_case,ly-y,u(1,:,:),u(8,:,:),jtot,
     $     chod_const,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,
     $     etavac,eta_local)
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Harris")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(4,:,:)=uy(4,:,:)*u(1,:,:)-u(4,:,:)*uy(1,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=uy(6,:,:)*u(1,:,:)-u(6,:,:)*uy(1,:,:)
         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(2,:,:)=u(2,:,:)+lambda_psi*LOG(COSH(y/lambda_psi))
            c(3,:,:)=-di*(u(3,:,:)*ux(3,:,:) + beta_e/beta0*ux(8,:,:)) 
     $           + eta_local*uy(3,:,:)*u(1,:,:)
            c(4,:,:)=uy(4,:,:) - uy(1,:,:)*u(4,:,:)/u(1,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=uy(6,:,:) - uy(1,:,:)*u(6,:,:)/u(1,:,:)
            c(7,:,:)=uy(7,:,:) - uy(1,:,:)*u(7,:,:)/u(1,:,:)
            c(8,:,:)=uy(8,:,:) - uy(1,:,:)*u(8,:,:)/u(1,:,:)
         END SELECT
      CASE("GEM_open")
         SELECT CASE(lrtb)
         CASE("top")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=uyy(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=uy(5,:,:)*u(1,:,:)-uy(1,:,:)*u(5,:,:)
            c(6,:,:)=uy(6,:,:)*u(1,:,:)-uy(1,:,:)*u(6,:,:)
            c(7,:,:)=uy(7,:,:)*u(1,:,:)-uy(1,:,:)*u(7,:,:)
            c(8,:,:)=uy(8,:,:)
         CASE("bottom")
c-----------------------------------------------------------------------
c     even boundary conditions
c-----------------------------------------------------------------------
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=uy(2,:,:)
            c(4,:,:)=uy(4,:,:)*u(1,:,:)-u(4,:,:)*uy(1,:,:)
            c(6,:,:)=uy(6,:,:)*u(1,:,:)-u(6,:,:)*uy(1,:,:)
            c(7,:,:)=uy(7,:,:)*u(1,:,:)-u(7,:,:)*uy(1,:,:)
            c(8,:,:)=uy(8,:,:)
c-----------------------------------------------------------------------
c     odd boundary conditions
c-----------------------------------------------------------------------
            c(3,:,:)=u(3,:,:)
            c(5,:,:)=u(5,:,:)
         CASE("left")
c-----------------------------------------------------------------------
c     even boundary conditions
c-----------------------------------------------------------------------
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=ux(2,:,:)
            c(5,:,:)=ux(5,:,:)*u(1,:,:)-u(5,:,:)*ux(1,:,:)
            c(6,:,:)=ux(6,:,:)*u(1,:,:)-u(6,:,:)*ux(1,:,:)
            c(7,:,:)=ux(7,:,:)*u(1,:,:)-u(7,:,:)*ux(1,:,:)
            c(8,:,:)=ux(8,:,:)
c-----------------------------------------------------------------------
c     odd boundary conditions
c-----------------------------------------------------------------------
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=u(4,:,:)
         CASE("right")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=uxx(2,:,:)
            c(3,:,:)=ux(3,:,:)
            c(4,:,:)=ux(4,:,:)*u(1,:,:)-ux(1,:,:)*u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=ux(6,:,:)*u(1,:,:)-ux(1,:,:)*u(6,:,:)
            c(7,:,:)=ux(7,:,:)*u(1,:,:)-ux(1,:,:)*u(7,:,:)
            c(8,:,:)=ux(8,:,:)
         END SELECT
      CASE("frc")
         CALL extra_coileval(t,x,coils,volt)

         SELECT CASE(lrtb)
         CASE("top")
            c(1,:,:)=u(1,:,:)-rhomin
            c(2,:,:)=-volt/(twopi*y)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(6,:,:)=u(6,:,:)
c-----------------------------------------------------------------------
c           flux convection b.c.
c-----------------------------------------------------------------------
            IF(flux_inflow)THEN
               c(5,:,:)=-u(5,:,:)*(u(2,:,:)/y + uy(2,:,:)) 
     $              + volt/(twopi*y)*u(1,:,:)
               c(7,:,:)=u(7,:,:)
c-----------------------------------------------------------------------
c           flux diffusion b.c.
c-----------------------------------------------------------------------
            ELSE
               c(5,:,:)=u(5,:,:)
               c(7,:,:)=eta_local*u(7,:,:) + volt/(twopi*y)
            ENDIF            
            c(8,:,:)=u(8,:,:)-Tmin*u(1,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=u(7,:,:)
            c(8,:,:)=uy(8,:,:)
         CASE("left","right")
            c(1,:,:)=u(1,:,:)-rhomin
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=u(7,:,:)
            c(8,:,:)=u(8,:,:)-Tmin*u(1,:,:)
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
      USE HallMHD_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: volt,eta_local,jtot,
     $     eta_rho,eta_p,eta_j,j_u3,j_ux3,j_uy3,j_u7
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
      jtot = SQRT((uy(3,:,:) + cyl_fac/y*u(3,:,:))**2 + ux(3,:,:)**2
     $     + u(7,:,:)**2)
      WHERE(y == 0)jtot = SQRT(uy(3,:,:)**2 + ux(3,:,:)**2
     $     + u(7,:,:)**2)
      j_u3 = cyl_fac/y*(uy(3,:,:) + cyl_fac/y*u(3,:,:))/jtot
      WHERE(y == 0)j_u3 = 0._r8
      j_ux3 = ux(3,:,:)/jtot
      j_uy3 = (uy(3,:,:) + cyl_fac/y*u(3,:,:))/jtot
      WHERE(y == 0)j_uy3 = uy(3,:,:)/jtot
      j_u7 = u(7,:,:)/jtot
      WHERE(jtot == 0)
         j_u3 = 0
         j_ux3 = 0
         j_uy3 = 0
         j_u7 = 0
      END WHERE

      CALL transport_seteta_u(eta_case,ly-y,u(1,:,:),u(8,:,:),jtot,
     $     chod_const,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,
     $     etavac,eta_local,eta_rho,eta_p,eta_j)

      SELECT CASE(init_type)
      CASE("Harris")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_u(4,1,:,:)=uy(4,:,:)
            c_u(4,4,:,:)=-uy(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)
            c_uy(4,4,:,:)=u(1,:,:)
            
            c_u(5,5,:,:)=one
            
            c_u(6,1,:,:)=uy(6,:,:)
            c_u(6,6,:,:)=-uy(1,:,:)
            c_uy(6,1,:,:)=-u(6,:,:)
            c_uy(6,6,:,:)=u(1,:,:)
         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_u(2,2,:,:)=one
            c_u(3,1,:,:)=uy(3,:,:)*(eta_local + eta_rho*u(1,:,:))
            c_u(3,3,:,:)=-di*ux(3,:,:) + eta_j*j_u3*uy(3,:,:)*u(1,:,:)
            c_u(3,7,:,:)=eta_j*j_u7*uy(3,:,:)*u(1,:,:)
            c_u(3,8,:,:)=eta_p*uy(3,:,:)*u(1,:,:)
            c_ux(3,3,:,:)=-di*u(3,:,:) + eta_j*j_ux3*uy(3,:,:)*u(1,:,:)
            c_ux(3,8,:,:)=-di*beta_e/beta0
            c_uy(3,3,:,:)=u(1,:,:)*(eta_local + eta_j*j_uy3*uy(3,:,:))
            c_u(4,1,:,:)=uy(1,:,:)*u(4,:,:)/u(1,:,:)**2
            c_u(4,4,:,:)=-uy(1,:,:)/u(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)/u(1,:,:)
            c_uy(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,1,:,:)=uy(1,:,:)*u(6,:,:)/u(1,:,:)**2
            c_u(6,6,:,:)=-uy(1,:,:)/u(1,:,:)
            c_uy(6,1,:,:)=-u(6,:,:)/u(1,:,:)
            c_uy(6,6,:,:)=one
            c_u(7,1,:,:)=uy(1,:,:)*u(7,:,:)/u(1,:,:)**2
            c_u(7,7,:,:)=-uy(1,:,:)/u(1,:,:)
            c_uy(7,1,:,:)=-u(7,:,:)/u(1,:,:)
            c_uy(7,7,:,:)=one
            c_u(8,1,:,:)=uy(1,:,:)*u(8,:,:)/u(1,:,:)**2
            c_u(8,8,:,:)=-uy(1,:,:)/u(1,:,:)
            c_uy(8,1,:,:)=-u(8,:,:)/u(1,:,:)
            c_uy(8,8,:,:)=one
        END SELECT
      CASE("GEM_open")
         SELECT CASE(lrtb)
         CASE("top")
            c_uy(1,1,:,:)=one
            c_uyy(2,2,:,:)=one

            c_u(3,1,:,:)=uy(3,:,:)
            c_u(3,3,:,:)=-uy(1,:,:)
            c_uy(3,1,:,:)=-u(3,:,:)
            c_uy(3,3,:,:)=u(1,:,:)

            c_u(4,4,:,:)=one

            c_u(5,1,:,:)=uy(5,:,:)
            c_u(5,5,:,:)=-uy(1,:,:)
            c_uy(5,1,:,:)=-u(5,:,:)
            c_uy(5,5,:,:)=u(1,:,:)

            c_u(6,1,:,:)=uy(6,:,:)
            c_u(6,6,:,:)=-uy(1,:,:)
            c_uy(6,1,:,:)=-u(6,:,:)
            c_uy(6,6,:,:)=u(1,:,:)

            c_u(7,1,:,:)=uy(7,:,:)
            c_u(7,7,:,:)=-uy(1,:,:)
            c_uy(7,1,:,:)=-u(7,:,:)
            c_uy(7,7,:,:)=u(1,:,:)

            c_uy(8,8,:,:)=one
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_uy(2,2,:,:)=one
            c_u(3,3,:,:)=one

            c_u(4,1,:,:)=uy(4,:,:)
            c_u(4,4,:,:)=-uy(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)
            c_uy(4,4,:,:)=u(1,:,:)

            c_u(5,5,:,:)=one

            c_u(6,1,:,:)=uy(6,:,:)
            c_u(6,6,:,:)=-uy(1,:,:)
            c_uy(6,1,:,:)=-u(6,:,:)
            c_uy(6,6,:,:)=u(1,:,:)

            c_u(7,1,:,:)=uy(7,:,:)
            c_u(7,7,:,:)=-uy(1,:,:)
            c_uy(7,1,:,:)=-u(7,:,:)
            c_uy(7,7,:,:)=u(1,:,:)

            c_uy(8,8,:,:)=one
         CASE("left")
            c_ux(1,1,:,:)=one
            c_ux(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_u(4,4,:,:)=one

            c_u(5,1,:,:)=ux(5,:,:)
            c_u(5,5,:,:)=-ux(1,:,:)
            c_ux(5,1,:,:)=-u(5,:,:)
            c_ux(5,5,:,:)=u(1,:,:)

            c_u(6,1,:,:)=ux(6,:,:)
            c_u(6,6,:,:)=-ux(1,:,:)
            c_ux(6,1,:,:)=-u(6,:,:)
            c_ux(6,6,:,:)=u(1,:,:)

            c_u(7,1,:,:)=ux(7,:,:)
            c_u(7,7,:,:)=-ux(1,:,:)
            c_ux(7,1,:,:)=-u(7,:,:)
            c_ux(7,7,:,:)=u(1,:,:)

            c_ux(8,8,:,:)=one
         CASE("right")
            c_ux(1,1,:,:)=one
            c_uxx(2,2,:,:)=one
            c_ux(3,3,:,:)=one

            c_u(4,1,:,:)=ux(4,:,:)
            c_u(4,4,:,:)=-ux(1,:,:)
            c_ux(4,1,:,:)=-u(4,:,:)
            c_ux(4,4,:,:)=u(1,:,:)

            c_u(5,5,:,:)=one

            c_u(6,1,:,:)=ux(6,:,:)
            c_u(6,6,:,:)=-ux(1,:,:)
            c_ux(6,1,:,:)=-u(6,:,:)
            c_ux(6,6,:,:)=u(1,:,:)

            c_u(7,1,:,:)=ux(7,:,:)
            c_u(7,7,:,:)=-ux(1,:,:)
            c_ux(7,1,:,:)=-u(7,:,:)
            c_ux(7,7,:,:)=u(1,:,:)

            c_ux(8,8,:,:)=one
        END SELECT
      CASE("frc")
         CALL extra_coileval(t,x,coils,volt)

         SELECT CASE(lrtb)
         CASE("top")
            c_u(1,1,:,:)=one
            c_u(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_u(6,6,:,:)=one
c-----------------------------------------------------------------------
c           flux convection b.c.
c-----------------------------------------------------------------------
            IF(flux_inflow)THEN
               c_u(5,1,:,:)=volt/(twopi*y)
               c_u(5,2,:,:)=-u(5,:,:)/y
               c_u(5,5,:,:)=-(u(2,:,:)/y + uy(2,:,:)) 
               c_uy(5,2,:,:)=-u(5,:,:)
               c_u(7,7,:,:)=one
c-----------------------------------------------------------------------
c           flux diffusion b.c.
c-----------------------------------------------------------------------
            ELSE
               c_u(5,5,:,:)=one
               c_u(7,1,:,:)=eta_rho*u(7,:,:)
               c_u(7,3,:,:)=eta_j*j_u3*u(7,:,:) 
               c_u(7,7,:,:)=eta_j*j_u7*u(7,:,:) + eta_local
               c_u(7,8,:,:)=eta_p*u(7,:,:)
               c_ux(7,3,:,:)=eta_j*j_ux3*u(7,:,:) 
               c_uy(7,3,:,:)=eta_j*j_uy3*u(7,:,:)
            ENDIF
            c_u(8,1,:,:)=-Tmin
            c_u(8,8,:,:)=one
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
            c_u(7,7,:,:)=one
            c_uy(8,8,:,:)=one
         CASE("left","right")
            c_u(1,1,:,:)=one
            c_u(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
            c_u(7,7,:,:)=one
            c_u(8,1,:,:)=-Tmin
            c_u(8,8,:,:)=one
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
      USE HallMHD_mod
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
         DO iqty=1,SIZE(mass,1)
            mass(iqty,iqty,:,:)=one
         ENDDO
      CASE("Huba")
         DO iqty=1,SIZE(mass,1)
            mass_x(iqty,iqty,:,:)=nhat(1,:,:)
            mass_y(iqty,iqty,:,:)=nhat(2,:,:)
         ENDDO
      CASE("frc")
         SELECT CASE(lrtb)
         CASE("top","left","right")
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
      USE HallMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Tix,Tiy,kperp,kfac,
     $     ve,vex,vey,r_fac,r_faci,j1,j2,j3,jtot,eta_local,b1,b2,Bsq,
     $     n_inv,nx_inv,ny_inv
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: vi,vix,viy,BdotT
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
c     temperature gradients.
c-----------------------------------------------------------------------
      Tix = ux(8,:,:)*n_inv + u(8,:,:)*nx_inv
      Tiy = uy(8,:,:)*n_inv + u(8,:,:)*ny_inv
c-----------------------------------------------------------------------
c     magnetic fields and currents.
c-----------------------------------------------------------------------
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2 + u(3,:,:)**2
      j1 = uy(3,:,:) + cyl_fac*r_faci*u(3,:,:)
      j2 = -ux(3,:,:)
      j3 = u(7,:,:)
      jtot = SQRT(j1**2 + j2**2 + j3**2)
c-----------------------------------------------------------------------
c     resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta(eta_case,ly-y,u(1,:,:),u(8,:,:),jtot,
     $     chod_const,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,
     $     etavac,eta_local)
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      IF(kappa_par /= kappa_perp)THEN
         CALL transport_BdotT(b1,b2,u(3,:,:),Tix,Tiy,BdotT)
         CALL transport_setkaniso(kappa_par,kappa_perp,Bsq,kperp,kfac)
      ELSE
         BdotT=0
         kperp=kappa_par
         kfac=0
      ENDIF
c-----------------------------------------------------------------------
c     velocities and their gradients.
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:)=u(i+3,:,:)*n_inv
         vix(i,:,:)=ux(i+3,:,:)*n_inv + u(i+3,:,:)*nx_inv
         viy(i,:,:)=uy(i+3,:,:)*n_inv + u(i+3,:,:)*ny_inv
      ENDDO

      ve=(u(6,:,:)-di*u(7,:,:))*n_inv
      vex=(ux(6,:,:)-di*ux(7,:,:))*n_inv + (u(6,:,:)-di*u(7,:,:))*nx_inv
      vey=(uy(6,:,:)-di*uy(7,:,:))*n_inv + (u(6,:,:)-di*u(7,:,:))*ny_inv
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
      fx(1,:,:) = r_fac*(u(4,:,:) - Dn*ux(1,:,:))
      fy(1,:,:) = r_fac*(u(5,:,:) - Dn*uy(1,:,:))
c-----------------------------------------------------------------------
c     poloidal magnetic flux equation.
c-----------------------------------------------------------------------
      fx(2,:,:) = -r_fac*di*nu*vex*n_inv
      fy(2,:,:) = -r_fac*di*nu*vey*n_inv
         
      s(2,:,:) = r_fac*(vi(2,:,:)*b1 - vi(1,:,:)*b2
     $     - di*(j2*b1 - j1*b2)*n_inv + eta_local*j3
     $     - di*nu*(vex*nx_inv + vey*ny_inv))
     $     - cyl_fac*r_faci*di*nu*ve*n_inv
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
      fx(3,:,:) = u(3,:,:)*vi(1,:,:) - ve*b1
     $     - di*(beta_e/beta0)*uy(8,:,:)*n_inv + eta_local*j2

      fy(3,:,:) = u(3,:,:)*vi(2,:,:) - ve*b2
     $     + di*(beta_e/beta0)*ux(8,:,:)*n_inv - eta_local*j1

      s(3,:,:) = di*u(3,:,:)*(j2*ny_inv + j1*nx_inv) 
     $     + cyl_fac*two*r_faci*di*u(3,:,:)*ux(3,:,:)*n_inv
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(4,:,:) = r_fac*(u(4,:,:)*vi(1,:,:) - two*mu*vix(1,:,:) 
     $     + half*u(3,:,:)**2 + u(8,:,:))
         
      fy(4,:,:) = r_fac*(u(4,:,:)*vi(2,:,:) 
     $     - mu*(viy(1,:,:) + vix(2,:,:)))

      s(4,:,:) = -r_fac*j3*b2
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(5,:,:) = r_fac*(u(5,:,:)*vi(1,:,:) 
     $     - mu*(vix(2,:,:) + viy(1,:,:)))
      
      fy(5,:,:) = r_fac*(u(5,:,:)*vi(2,:,:) + half*u(3,:,:)**2
     $     - two*mu*viy(2,:,:))
      
      s(5,:,:) = r_fac*(j3*b1 - uy(8,:,:)) + cyl_fac*(u(6,:,:)*vi(3,:,:)
     $     - half*u(3,:,:)**2 - two*mu*vi(2,:,:)*r_faci)
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(6,:,:) = r_fac*(u(6,:,:)*vi(1,:,:) - nu*vex - mu*vix(3,:,:))
      
      fy(6,:,:) = r_fac*(u(6,:,:)*vi(2,:,:) - nu*vey - mu*viy(3,:,:))

      s(6,:,:) = r_fac*(j1*b2 - j2*b1)
     $     - cyl_fac*(u(6,:,:)*vi(2,:,:) + r_faci*(mu*vi(3,:,:)+nu*ve))
c-----------------------------------------------------------------------
c     out-of-plane electron momentum (current) equation.
c-----------------------------------------------------------------------
      fx(7,:,:)=b2
      fy(7,:,:)=-b1
      s(7,:,:)=j3
c-----------------------------------------------------------------------
c     total pressure equation.
c-----------------------------------------------------------------------
      fx(8,:,:)=r_fac*(gamma_fac*u(8,:,:)*vi(1,:,:)
     $     - kfac*BdotT(1,:,:) - kperp*Tix)
      fy(8,:,:)=r_fac*(gamma_fac*u(8,:,:)*vi(2,:,:)
     $     - kfac*BdotT(2,:,:) - kperp*Tiy)
      s(8,:,:)=r_fac*(vi(1,:,:)*ux(8,:,:) + vi(2,:,:)*uy(8,:,:)
     $     + eta_local*jtot**2 
     $     + mu*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:)
     $     + vix(3,:,:)**2 + viy(3,:,:)**2) + nu*(vex**2 + vey**2)
     $     + cyl_fac*r_faci**2*(mu*(two*vi(2,:,:)**2 + vi(3,:,:)**2)
     $     + nu*ve**2))
c-----------------------------------------------------------------------
c     initial equilibrium source term.
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("GEM","GEM_open")
            CALL HallMHD_equil(x,y,u0,u0x,u0y,.TRUE.)
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
      USE HallMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci,Tix,Tiy,
     $     Tix_un,Tiy_un,Ti_un,kperp,kfac,kappa_bsq,ve,vex,vey,j1,j2,j3,
     $     jtot,eta_local,b1,b2,Bsq,eta_rho,eta_p,eta_j,
     $     j_u3,j_ux3,j_uy3,j_u7,n_inv,nx_inv,ny_inv,ve_un,vex_un,vey_un
      REAL(r8), DIMENSION(3,SIZE(x,1),SIZE(x,2)) :: vi,vix,viy,
     $     vi_un,vix_un,viy_un,BdotT,BdotT_b1,BdotT_b2,BdotT_b3,
     $     BdotT_Tx,BdotT_Ty
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
c     temperature gradients and derivatives.
c-----------------------------------------------------------------------
      Tix = ux(8,:,:)*n_inv + u(8,:,:)*nx_inv
      Tiy = uy(8,:,:)*n_inv + u(8,:,:)*ny_inv
      Ti_un=-u(8,:,:)*n_inv**2
      Tix_un=-(two*u(8,:,:)*nx_inv + ux(8,:,:)*n_inv)*n_inv
      Tiy_un=-(two*u(8,:,:)*ny_inv + uy(8,:,:)*n_inv)*n_inv
c-----------------------------------------------------------------------
c     magnetic fields and currents.
c-----------------------------------------------------------------------
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2 + u(3,:,:)**2
      j1 = uy(3,:,:) + cyl_fac*r_faci*u(3,:,:)
      j2 = -ux(3,:,:)
      j3 = u(7,:,:)
      jtot = SQRT(j1**2 + j2**2 + j3**2)

      j_u3 = cyl_fac*r_faci*(uy(3,:,:) + cyl_fac*r_faci*u(3,:,:))/jtot
      j_ux3 = ux(3,:,:)/jtot
      j_uy3 = (uy(3,:,:) + cyl_fac*r_faci*u(3,:,:))/jtot
      j_u7 = u(7,:,:)/jtot
      WHERE(jtot == 0)
         j_u3 = 0
         j_ux3 = 0
         j_uy3 = 0
         j_u7 = 0
      END WHERE
c-----------------------------------------------------------------------
c     resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta_u(eta_case,ly-y,u(1,:,:),u(8,:,:),jtot,
     $     chod_const,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,
     $     etavac,eta_local,eta_rho,eta_p,eta_j)
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      IF(kappa_par /= kappa_perp)THEN
         CALL transport_BdotT_u(b1,b2,u(3,:,:),Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty)
         CALL transport_setkaniso_u(kappa_par,kappa_perp,Bsq,kperp,kfac,
     $        kappa_bsq)
      ELSE
         BdotT=0
         BdotT_b1=0
         BdotT_b2=0
         BdotT_b3=0
         BdotT_Tx=0
         BdotT_Ty=0
         kperp=kappa_par
         kfac=0
         kappa_bsq=0
      ENDIF
c-----------------------------------------------------------------------
c     velocities and their derivatives.
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:)=u(i+3,:,:)*n_inv
         vix(i,:,:)=ux(i+3,:,:)*n_inv + u(i+3,:,:)*nx_inv
         viy(i,:,:)=uy(i+3,:,:)*n_inv + u(i+3,:,:)*ny_inv

         vi_un(i,:,:)=-vi(i,:,:)*n_inv
         vix_un(i,:,:) = -(vi(i,:,:)*nx_inv + vix(i,:,:)*n_inv)
         viy_un(i,:,:) = -(vi(i,:,:)*ny_inv + viy(i,:,:)*n_inv)
      ENDDO

      ve=(u(6,:,:)-di*u(7,:,:))*n_inv
      vex=(ux(6,:,:)-di*ux(7,:,:))*n_inv + (u(6,:,:)-di*u(7,:,:))*nx_inv
      vey=(uy(6,:,:)-di*uy(7,:,:))*n_inv + (u(6,:,:)-di*u(7,:,:))*ny_inv
      
      ve_un=-ve*n_inv
      vex_un = -(ve*nx_inv + vex*n_inv)
      vey_un = -(ve*ny_inv + vey*n_inv)
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
      fx_u(1,4,:,:) = r_fac
      fx_ux(1,1,:,:) = -r_fac*Dn
      fy_u(1,5,:,:) = r_fac
      fy_uy(1,1,:,:) = -r_fac*Dn
c-----------------------------------------------------------------------
c     poloidal magnetic flux equation.
c-----------------------------------------------------------------------
      fx_u(2,1,:,:)=-r_fac*di*nu*(vex_un - vex*n_inv)*n_inv
      fx_u(2,6,:,:)=-r_fac*di*nu*nx_inv*n_inv
      fx_u(2,7,:,:)= r_fac*di**2*nu*nx_inv*n_inv
      fx_ux(2,1,:,:)=-r_fac*di*nu*ve_un*n_inv
      fx_ux(2,6,:,:)=-r_fac*di*nu*n_inv**2
      fx_ux(2,7,:,:)= r_fac*di**2*nu*n_inv**2

      fy_u(2,1,:,:)=-r_fac*di*nu*(vey_un - vey*n_inv)*n_inv
      fy_u(2,6,:,:)=-r_fac*di*nu*ny_inv*n_inv
      fy_u(2,7,:,:)= r_fac*di**2*nu*ny_inv*n_inv
      fy_uy(2,1,:,:)=-r_fac*di*nu*ve_un*n_inv
      fy_uy(2,6,:,:)=-r_fac*di*nu*n_inv**2
      fy_uy(2,7,:,:)=r_fac*di**2*nu*n_inv**2

      s_u(2,1,:,:)=r_fac*(vi_un(2,:,:)*b1 - vi_un(1,:,:)*b2
     $     + di*(j2*b1 - j1*b2)*n_inv**2 + eta_rho*j3
     $     - di*nu*((vex_un*nx_inv + vey_un*ny_inv)
     $     - two*(vex*nx_inv + vey*ny_inv)*n_inv))
     $     + cyl_fac*r_faci*two*di*nu*ve*n_inv**2
      s_u(2,2,:,:)=-cyl_fac*(vi(2,:,:)-di*j2*n_inv)
      s_u(2,3,:,:)=di*cyl_fac*b2*n_inv + r_fac*eta_j*j_u3*j3
      s_u(2,4,:,:)=-r_fac*b2*n_inv
      s_u(2,5,:,:)=r_fac*b1*n_inv
      s_u(2,6,:,:)= -r_fac*di*nu*(nx_inv**2 + ny_inv**2)
     $     - cyl_fac*r_faci*di*nu*n_inv**2
      s_u(2,7,:,:)=r_fac*(di**2*nu*(nx_inv**2 + ny_inv**2)
     $     + eta_j*j_u7*j3 + eta_local)
     $     + cyl_fac*r_faci*di**2*nu*n_inv**2
      s_u(2,8,:,:)=r_fac*eta_p*j3
      s_ux(2,1,:,:)=r_fac*di*nu*(vex + ve_un*ux(1,:,:))*n_inv**2
      s_ux(2,2,:,:)=r_fac*(di*j1*n_inv - vi(1,:,:))
      s_ux(2,3,:,:)=r_fac*(di*b1*n_inv + eta_j*j_ux3*j3)
      s_ux(2,6,:,:)=-r_fac*di*nu*nx_inv*n_inv
      s_ux(2,7,:,:)=r_fac*di**2*nu*nx_inv*n_inv
      s_uy(2,1,:,:)=r_fac*di*nu*(vey + ve_un*uy(1,:,:))*n_inv**2
      s_uy(2,2,:,:)=-r_fac*(vi(2,:,:)-di*j2*n_inv)
      s_uy(2,3,:,:)=r_fac*(di*b2*n_inv + eta_j*j_uy3*j3)
      s_uy(2,6,:,:)=-r_fac*di*nu*ny_inv*n_inv
      s_uy(2,7,:,:)=r_fac*di**2*nu*ny_inv*n_inv
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=u(3,:,:)*vi_un(1,:,:) - ve_un*b1 
     $     + di*(beta_e/beta0)*uy(8,:,:)*n_inv**2 +  eta_rho*j2
      fx_u(3,2,:,:)=cyl_fac*r_faci*ve 
      fx_u(3,3,:,:)=vi(1,:,:) + eta_j*j_u3*j2
      fx_u(3,4,:,:)=u(3,:,:)*n_inv
      fx_u(3,6,:,:)=-n_inv*b1
      fx_u(3,7,:,:)=di*n_inv*b1 + eta_j*j_u7*j2
      fx_u(3,8,:,:)= eta_p*j2
      fx_ux(3,3,:,:)=-eta_local + eta_j*j_ux3*j2
      fx_uy(3,2,:,:)=ve
      fx_uy(3,3,:,:)=eta_j*j_uy3*j2
      fx_uy(3,8,:,:)=-di*(beta_e/beta0)*n_inv
         
      fy_u(3,1,:,:)=u(3,:,:)*vi_un(2,:,:) - ve_un*b2 
     $     - di*(beta_e/beta0)*ux(8,:,:)*n_inv**2 - eta_rho*j1
      fy_u(3,3,:,:)=vi(2,:,:) - eta_local*cyl_fac*r_faci - eta_j*j_u3*j1
      fy_u(3,5,:,:)=u(3,:,:)*n_inv
      fy_u(3,6,:,:)=-n_inv*b2
      fy_u(3,7,:,:)=di*n_inv*b2 - eta_j*j_u7*j1
      fy_u(3,8,:,:)= -eta_p*j1
      fy_ux(3,2,:,:)=-ve
      fy_ux(3,3,:,:)=-eta_j*j_ux3*j1
      fy_ux(3,8,:,:)=di*(beta_e/beta0)*n_inv
      fy_uy(3,3,:,:)=-eta_local - eta_j*j_uy3*j1
         
      s_u(3,1,:,:)=-two*di*u(3,:,:)*(j2*ny_inv + j1*nx_inv)*n_inv
     $     - cyl_fac*two*r_faci*di*u(3,:,:)*ux(3,:,:)*n_inv**2
      s_u(3,3,:,:)=di*(j2*ny_inv + j1*nx_inv)
     $     + cyl_fac*r_faci*di*(u(3,:,:)*nx_inv + two*ux(3,:,:)*n_inv)
      s_ux(3,1,:,:)=-di*u(3,:,:)*j1*n_inv**2
      s_ux(3,3,:,:)=-di*u(3,:,:)*ny_inv
     $     + cyl_fac*two*r_faci*di*u(3,:,:)*n_inv
      s_uy(3,1,:,:)=-di*u(3,:,:)*j2*n_inv**2
      s_uy(3,3,:,:)=di*u(3,:,:)*nx_inv
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(1,:,:) - mu*two*vix_un(1,:,:))
      fx_u(4,3,:,:)=r_fac*u(3,:,:)
      fx_u(4,4,:,:)=r_fac*two*(vi(1,:,:) - mu*nx_inv)
      fx_u(4,8,:,:)=r_fac
      fx_ux(4,1,:,:)=-r_fac*mu*two*vi_un(1,:,:)
      fx_ux(4,4,:,:)=-r_fac*mu*two*n_inv

      fy_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(2,:,:)
     $     - mu*(viy_un(1,:,:) + vix_un(2,:,:)))
      fy_u(4,4,:,:)=r_fac*(vi(2,:,:) - mu*ny_inv)
      fy_u(4,5,:,:)=r_fac*(vi(1,:,:) - mu*nx_inv)
      fy_ux(4,1,:,:)=-r_fac*mu*vi_un(2,:,:)
      fy_ux(4,5,:,:)=-r_fac*mu*n_inv
      fy_uy(4,1,:,:)=-r_fac*mu*vi_un(1,:,:)
      fy_uy(4,4,:,:)=-r_fac*mu*n_inv

      s_u(4,7,:,:)=-r_fac*b2
      s_ux(4,2,:,:)=-r_fac*j3
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------      
      fx_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(1,:,:)
     $     - mu*(vix_un(2,:,:) + viy_un(1,:,:)))
      fx_u(5,4,:,:)=r_fac*(u(5,:,:)*n_inv - mu*ny_inv)
      fx_u(5,5,:,:)=r_fac*(vi(1,:,:)-mu*nx_inv)
      fx_ux(5,1,:,:)=-r_fac*mu*vi_un(2,:,:)
      fx_ux(5,5,:,:)=-r_fac*mu*n_inv
      fx_uy(5,1,:,:)=-r_fac*mu*vi_un(1,:,:)
      fx_uy(5,4,:,:)=-r_fac*mu*n_inv
      
      fy_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(2,:,:) - mu*two*viy_un(2,:,:))
      fy_u(5,3,:,:)=r_fac*u(3,:,:)
      fy_u(5,5,:,:)=r_fac*two*(vi(2,:,:) - mu*ny_inv)
      fy_uy(5,1,:,:)=-r_fac*mu*two*vi_un(2,:,:)
      fy_uy(5,5,:,:)=-r_fac*mu*two*n_inv

      s_u(5,1,:,:)=cyl_fac*(u(6,:,:)*vi_un(3,:,:) 
     $     - r_faci*two*mu*vi_un(2,:,:))
      s_u(5,2,:,:)=-cyl_fac*j3
      s_u(5,3,:,:)=-cyl_fac*u(3,:,:)
      s_u(5,5,:,:)=-cyl_fac*r_faci*two*mu*n_inv
      s_u(5,6,:,:)=cyl_fac*two*vi(3,:,:)
      s_u(5,7,:,:)=r_fac*b1
      s_uy(5,2,:,:)=-r_fac*j3
      s_uy(5,8,:,:)=-r_fac
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------      
      fx_u(6,1,:,:)=r_fac*(u(6,:,:)*vi_un(1,:,:) - mu*vix_un(3,:,:)
     $     - nu*vex_un)
      fx_u(6,4,:,:)=r_fac*u(6,:,:)*n_inv
      fx_u(6,6,:,:)=r_fac*(vi(1,:,:) - (mu+nu)*nx_inv)
      fx_u(6,7,:,:)= r_fac*di*nu*nx_inv
      fx_ux(6,1,:,:)=-r_fac*(mu*vi_un(3,:,:) + nu*ve_un)
      fx_ux(6,6,:,:)=-r_fac*(mu+nu)*n_inv
      fx_ux(6,7,:,:)= r_fac*di*nu*n_inv
      
      fy_u(6,1,:,:)=r_fac*(u(6,:,:)*vi_un(2,:,:) - mu*viy_un(3,:,:)
     $     - nu*vey_un)
      fy_u(6,5,:,:)=r_fac*u(6,:,:)*n_inv
      fy_u(6,6,:,:)=r_fac*(vi(2,:,:) - (mu+nu)*ny_inv)
      fy_u(6,7,:,:)= r_fac*di*nu*ny_inv
      fy_uy(6,1,:,:)=-r_fac*(mu*vi_un(3,:,:) + nu*ve_un)
      fy_uy(6,6,:,:)=-r_fac*(mu+nu)*n_inv
      fy_uy(6,7,:,:)= r_fac*di*nu*n_inv
      
      s_u(6,1,:,:)=cyl_fac*(vi(2,:,:)*vi(3,:,:) 
     $     - r_faci*(mu*vi_un(3,:,:) + nu*ve_un))
      s_u(6,2,:,:)=cyl_fac*j2
      s_u(6,3,:,:)=cyl_fac*b2
      s_u(6,5,:,:)=-cyl_fac*vi(3,:,:)
      s_u(6,6,:,:)=-cyl_fac*(vi(2,:,:) + r_faci*(mu+nu)*n_inv)
      s_u(6,7,:,:)= cyl_fac*r_faci*di*nu*n_inv
      s_ux(6,2,:,:)=r_fac*j1
      s_ux(6,3,:,:)=r_fac*b1
      s_uy(6,2,:,:)=r_fac*j2
      s_uy(6,3,:,:)=r_fac*b2
c-----------------------------------------------------------------------
c     out-of-plane electron momentum (current) equation.
c-----------------------------------------------------------------------      
      fx_ux(7,2,:,:)=one
      fy_u(7,2,:,:)=cyl_fac*r_faci
      fy_uy(7,2,:,:)=one
      s_u(7,7,:,:)=one
c-----------------------------------------------------------------------
c     total pressure equation.
c-----------------------------------------------------------------------
      fx_u(8,1,:,:)=r_fac*(gamma_fac*u(8,:,:)*vi_un(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*Tix_un 
     $     - kfac*BdotT_Ty(1,:,:)*Tiy_un)
      fx_u(8,2,:,:)=-cyl_fac*(-kfac*BdotT_b1(1,:,:)
     $     + two*b1*kappa_bsq*(BdotT(1,:,:) - Tix))
      fx_u(8,3,:,:)=r_fac*(-kfac*BdotT_b3(1,:,:)
     $     + two*u(3,:,:)*kappa_bsq*(BdotT(1,:,:) - Tix))
      fx_u(8,4,:,:)=r_fac*gamma_fac*u(8,:,:)*n_inv
      fx_u(8,8,:,:)=r_fac*(gamma_fac*vi(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*nx_inv
     $     - kfac*BdotT_Ty(1,:,:)*ny_inv)

      fx_ux(8,1,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*Ti_un
      fx_ux(8,2,:,:)=r_fac*(-kfac*BdotT_b2(1,:,:)
     $     + two*b2*kappa_bsq*(BdotT(1,:,:) - Tix))
      fx_ux(8,8,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*n_inv

      fx_uy(8,1,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*Ti_un
      fx_uy(8,2,:,:)=-r_fac*(-kfac*BdotT_b1(1,:,:)
     $     + two*b1*kappa_bsq*(BdotT(1,:,:) - Tix))
      fx_uy(8,8,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*n_inv

      fy_u(8,1,:,:)=r_fac*(gamma_fac*u(8,:,:)*vi_un(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*Tiy_un
     $     - kfac*BdotT_Tx(2,:,:)*Tix_un)
      fy_u(8,2,:,:)=-cyl_fac*(-kfac*BdotT_b1(2,:,:)
     $     + two*b1*kappa_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(8,3,:,:)=r_fac*(-kfac*BdotT_b3(2,:,:)
     $     + two*u(3,:,:)*kappa_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(8,5,:,:)=r_fac*gamma_fac*u(8,:,:)*n_inv
      fy_u(8,8,:,:)=r_fac*(gamma_fac*vi(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*ny_inv
     $     - kfac*BdotT_Tx(2,:,:)*nx_inv)

      fy_ux(8,1,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*Ti_un
      fy_ux(8,2,:,:)=r_fac*(-kfac*BdotT_b2(2,:,:)
     $     + two*b2*kappa_bsq*(BdotT(2,:,:) - Tiy))
      fy_ux(8,8,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*n_inv

      fy_uy(8,1,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*Ti_un
      fy_uy(8,2,:,:)=-r_fac*(-kfac*BdotT_b1(2,:,:)
     $     + two*b1*kappa_bsq*(BdotT(2,:,:) - Tiy))
      fy_uy(8,8,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*n_inv

      s_u(8,1,:,:)=r_fac*(vi_un(1,:,:)*ux(8,:,:)
     $     + vi_un(2,:,:)*uy(8,:,:) + eta_rho*jtot**2
     $     + mu*two*(two*vix(1,:,:)*vix_un(1,:,:) 
     $     + two*viy(2,:,:)*viy_un(2,:,:) 
     $     + viy(1,:,:)*viy_un(1,:,:) + vix(2,:,:)*vix_un(2,:,:)
     $     + viy(1,:,:)*vix_un(2,:,:) + vix(2,:,:)*viy_un(1,:,:)
     $     + vix(3,:,:)*vix_un(3,:,:) + viy(3,:,:)*viy_un(3,:,:))
     $     + nu*two*(vex*vex_un + vey*vey_un)
     $     + cyl_fac*two*r_faci**2
     $     *(mu*(two*vi(2,:,:)*vi_un(2,:,:) + vi(3,:,:)*vi_un(3,:,:))
     $     + nu*ve*ve_un))
      s_u(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_u3
      s_u(8,4,:,:)=r_fac*(ux(8,:,:)*n_inv + mu*(4*vix(1,:,:)*nx_inv
     $     + two*viy(1,:,:)*ny_inv + two*vix(2,:,:)*ny_inv))
      s_u(8,5,:,:)=r_fac*(uy(8,:,:)*n_inv + mu*(4*viy(2,:,:)*ny_inv
     $     + two*vix(2,:,:)*nx_inv + two*viy(1,:,:)*nx_inv)
     $     + cyl_fac*4._r8*mu*r_faci**2*vi(2,:,:)*n_inv)
      s_u(8,6,:,:)=r_fac*two*(mu*(vix(3,:,:)*nx_inv + viy(3,:,:)*ny_inv)
     $     + nu*(vex*nx_inv + vey*ny_inv)
     $     + cyl_fac*r_faci**2*(mu*vi(3,:,:)+nu*ve)*n_inv)
      s_u(8,7,:,:)=r_fac*((eta_local*two + eta_j*jtot)*jtot*j_u7
     $     - di*nu*two*(vex*nx_inv + vey*ny_inv
     $     + cyl_fac*r_faci**2*ve*n_inv))
      s_u(8,8,:,:)=r_fac*eta_p*jtot**2

      s_ux(8,1,:,:)=r_fac*(mu*two*(two*vix(1,:,:)*vi_un(1,:,:)  
     $     + vi_un(2,:,:)*(vix(2,:,:) + viy(1,:,:))
     $     + vix(3,:,:)*vi_un(3,:,:)) + two*nu*vex*ve_un)
      s_ux(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_ux3
      s_ux(8,4,:,:)=r_fac*mu*4._r8*vix(1,:,:)*n_inv
      s_ux(8,5,:,:)=r_fac*mu*two*n_inv*(vix(2,:,:) + viy(1,:,:))
      s_ux(8,6,:,:)=r_fac*two*(mu*vix(3,:,:) + nu*vex)*n_inv
      s_ux(8,7,:,:)=-r_fac*di*nu*two*vex*n_inv
      s_ux(8,8,:,:)=r_fac*vi(1,:,:)

      s_uy(8,1,:,:)=r_fac*(mu*two*(two*viy(2,:,:)*vi_un(2,:,:)  
     $     + vi_un(1,:,:)*(viy(1,:,:) + vix(2,:,:))
     $     + viy(3,:,:)*vi_un(3,:,:)) + nu*two*vey*ve_un)
      s_uy(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_uy3
      s_uy(8,4,:,:)=r_fac*mu*two*n_inv*(viy(1,:,:) + vix(2,:,:))
      s_uy(8,5,:,:)=r_fac*mu*4._r8*viy(2,:,:)*n_inv
      s_uy(8,6,:,:)=r_fac*two*(mu*viy(3,:,:) + nu*vey)*n_inv
      s_uy(8,7,:,:)=-r_fac*di*nu*two*vey*n_inv
      s_uy(8,8,:,:)=r_fac*vi(2,:,:)
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
      USE HallMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y

      REAL(r8), DIMENSION(SIZE(y,1),SIZE(y,2)) :: r_fac
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
c        cylindrical to cartesian relationships:
c        1: z --> x
c        2: r --> y
c        3: phi --> z
c-----------------------------------------------------------------------
      r_fac=y**cyl_fac

      mass(1,1,:,:)=r_fac
      mass(2,2,:,:)=r_fac
      mass(4,4,:,:)=r_fac
      mass(5,5,:,:)=r_fac
      mass(6,6,:,:)=r_fac
      mass(8,8,:,:)=r_fac/(gamma-one)
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
      USE HallMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,phi
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Harris")
         ksi=lx*(x-half)
         phi=y-half
      CASE("Huba")
         ksi=lx*(x-half)
         phi=ly*(y-half)
      CASE("GEM")
         ksi=SIGN(half*lx*((two*x - one)**2 + gr_curve*ABS(two*x - one))
     $        /(one + gr_curve),two*x - one)
         phi=SIGN(half*ly*((two*y - one)**2 + gr_curve*ABS(two*y - one))
     $        /(one + gr_curve),two*y - one)
      CASE("GEM_open")
         ksi=two*lx*x
         phi=two*ly*(y**2+gr_curve*y)/(one+gr_curve)
      CASE("frc")
         ksi=lx*x + zmin/L0
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
      USE HallMHD_mod
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
      USE HallMHD_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     deallocate appropriate arrays.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("frc")
         CALL extra_equil_dealloc(interp,equil_bc,equil,equilxy)
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
