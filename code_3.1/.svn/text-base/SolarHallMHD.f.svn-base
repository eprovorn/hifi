c-----------------------------------------------------------------------
c     file SolarHallMHD.f.
c     contains specifications for Hall MHD model with gravity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. SolarHallMHD_mod.
c     1. SolarHallMHD_equil.
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
c     subprogram 0. SolarHallMHD_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE SolarHallMHD_mod
      USE transport_mod
      USE extra_mod
      IMPLICIT NONE

      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,qe=1.602e-19,
     $     me=9.109e-31,mp=1.673e-27,ep0=8.854e-12,mu0=4.e-7*pi,
     $     R0=6.96e8,k_B=1.381e-23,g0=2.74e2,Zeff=1.,chod_const=0.1,
     $     kappa_min=0.,kappa_max=1.e8

      LOGICAL :: source=.FALSE.,cylinder=.FALSE.
      CHARACTER(16) :: init_type=".",eta_case=".",kappa_case="."
      INTEGER :: cyl_fac=0
      REAL(r8) :: di=0.,nu=0.,eta=0.,etavac=1.,
     $     v_chod_norm=1.,etac_norm=1.,etas_norm=1.,mu=0.,kappa_par=0.,
     $     kappa_perp=0.,Dn=0.,beta0=1.,Te_frac=0.,n0=1.,T0=1.,b0=1.,
     $     Rmax=0.,gamma_fac=1.,gravity=1.,v_peak=0.,v_period=1.,
     $     v_angle=0.,ke_norm=1.,ki_norm=1.,xe_norm=1.,xi_norm=1.,
     $     mu_fac=0.,alpha=1.,lambda=1.,lx=1.,ly=1.,rscale=1.,
     $     phiscale=1.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. SolarHallMHD_equil.
c     computes equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE SolarHallMHD_equil(x,y,u,ux,uy,deriv)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: deriv

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: rsq,coshy_psi
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0

      rsq = x**2 + y**2
      SELECT CASE(init_type)      
      CASE("breakout","half_breakout")
         CALL RANDOM_NUMBER(u)
         u = 1.e-20*u

         u(1,:,:) = rsq**(half*(one-mu_fac))
         u(2,:,:) = -y/rsq**1.5*(one + alpha*(4.*x**2 - y**2)/rsq**2)
         u(8,:,:) = (T0*mu0*n0*k_B/b0**2)*u(1,:,:)/SQRT(rsq)
         IF(.NOT. deriv)RETURN
         ux(1,:,:) = (one-mu_fac)*x*rsq**(-half*(one+mu_fac))
         uy(1,:,:) = (one-mu_fac)*y*rsq**(-half*(one+mu_fac))
         ux(2,:,:) = x*y/rsq**2.5*(3. 
     $        + alpha*(20.*x**2 - 15.*y**2)/rsq**2)
         uy(2,:,:) = (2.*y**2 - x**2)/rsq**2.5
     $        + alpha*(35.*x**2*y**2 - 4.*rsq**2)/rsq**4.5
         ux(8,:,:) = -(T0*mu0*n0*k_B/b0**2)*mu_fac
     $        *x*rsq**(-half*mu_fac-one)
         uy(8,:,:) = -(T0*mu0*n0*k_B/b0**2)*mu_fac
     $        *y*rsq**(-half*mu_fac-one)
      CASE("GEM","GEM9")
         CALL RANDOM_NUMBER(u)
         u = 1.e-20*u
         coshy_psi=COSH(y/lambda)

         u(1,:,:)=one/coshy_psi**2 + .2_r8
         u(2,:,:)=-lambda*LOG(coshy_psi)
         u(7,:,:)=-one/(lambda*coshy_psi**2)
         u(8,:,:)=beta0*u(1,:,:)
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=-two*TANH(y/lambda)/(lambda*coshy_psi**2)
         uy(2,:,:)=-TANH(y/lambda)
         uy(7,:,:)=-u(7,:,:)*two*TANH(y/lambda)/lambda
         uy(8,:,:)=beta0*uy(1,:,:)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SolarHallMHD_equil
      END MODULE SolarHallMHD_mod
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
      USE SolarHallMHD_mod
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
      NAMELIST/SolarHallMHD_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,rscale,phiscale,init_type,source,cylinder,
     $     di,eta,eta_case,etavac,mu,nu,kappa_case,kappa_par,
     $     kappa_perp,Dn,Te_frac,n0,T0,b0,Rmax,v_peak,v_period,v_angle,
     $     alpha,lambda,beta0,lx,ly
c-----------------------------------------------------------------------
c     Sample namelist.
c-----------------------------------------------------------------------
c$$$&SolarHallMHD_list
c$$$
c$$$	cylinder=t
c$$$    rscale=0.
c$$$    phiscale=0.
c$$$
c$$$	init_type="breakout"
c$$$	source=f
c$$$
c$$$    alpha=1.
c$$$	Rmax=30.0
c$$$
c$$$    v_peak=8682.52
c$$$    v_period=1.e5
c$$$    v_angle=18
c$$$
c$$$	Te_frac=0.5
c$$$
c$$$    lambda=1.
c$$$    beta0=1.
c$$$    lx=1.
c$$$    ly=1.
c$$$
c$$$	di=.1
c$$$	nu=2.e-5
c$$$
c$$$	mu=1.00e-1
c$$$
c$$$	eta=2.01e-3
c$$$	eta_case="spitzer-chodura"
c$$$	etavac=100.
c$$$
c$$$    kappa_case="anisotropic"
c$$$	kappa_par=20.0
c$$$	kappa_perp=4.02e-4
c$$$	Dn=2.00e-3
c$$$
c$$$	b0=1.e-4
c$$$	n0=2.e14
c$$$	T0=2.e6
c$$$/
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=SolarHallMHD_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="SolarHallMHD"

      nqty=8
      nqty_schur=0
      IF(init_type=="GEM9")nqty=9
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
      USE SolarHallMHD_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty

      REAL(r8) :: tnorm
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(eta_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(kappa_case,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(cylinder,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(rscale,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(phiscale,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(etavac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_par,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_perp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Dn,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Te_frac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(T0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Rmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(v_peak,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(v_period,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(v_angle,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(7)=.TRUE.
      adapt_qty(7)=.FALSE.
      IF(init_type=="GEM9")THEN
         static(9)=.TRUE.
         adapt_qty(9)=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     compute normalizations in MKS units.
c     v_chod_norm normalizes the ve/vs term in Chodura resistivity.
c     etac_norm absorbs the constants in front of Chodura resistivity.
c     etas_norm absorbs the constants in front of Spitzer resistivity.
c-----------------------------------------------------------------------
      tnorm=R0*SQRT(mu0*n0*mp*Zeff)/b0
      gravity=g0*tnorm**2/R0
      v_chod_norm=one/(qe*R0)*SQRT(mp*Zeff/(mu0*n0))
      etac_norm=me/(qe*R0*b0*SQRT(ep0*mu0))
      etas_norm=5.e-5*17.*tnorm*(two*n0*mu0*qe)**1.5/(mu0*R0**2*b0**3)
      ke_norm = 3.56e21*(b0**2/(two*n0*mu0*qe))**2.5
     $     *tnorm/(two*n0*R0**2)
      xe_norm = 3.56e21*(b0**2/(two*n0*mu0*qe))**1.5*b0/n0
      ki_norm = 1.97e-7*(b0**2/(two*n0*mu0*qe))**2.5
     $     *SQRT(mu0/(n0*mp))/(two*R0*b0)
      xi_norm = 1.97e-7*(b0**2/(two*n0*mu0*qe))**1.5
     $     *b0/(n0*mp*SQRT(Zeff))
      gamma_fac=gamma/(gamma-one)

      SELECT CASE(init_type)
      CASE("breakout","half_breakout")
         cylinder=.TRUE.
         di = zero
         T0 = 4.e6
         v_period = 1.e5/tnorm
         v_peak = v_peak*(pi**2/4._r8)/v_period
         v_angle = pi/15._r8
         mu_fac=g0*R0*mp*Zeff/(k_B*T0)
      CASE("GEM","GEM9")
         cylinder=.FALSE.
         gravity=0.
         lambda=.5_r8
         lx=25.6_r8
         ly=12.8_r8
         di=one
         beta0=half
         Te_frac=one/6._r8
      END SELECT

      IF(cylinder)cyl_fac=1
      IF(Te_frac > 1.)Te_frac=1.
      IF(Te_frac < 0.)Te_frac=0.
      IF(Te_frac > 0)THEN
         etas_norm = etas_norm/Te_frac**1.5
      ELSE
         etas_norm=0.
      ENDIF
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
      USE SolarHallMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     SolarHallMHD initial conditions.
c-----------------------------------------------------------------------
      CALL SolarHallMHD_equil(x,y,u,ux,uy,.FALSE.)

      SELECT CASE(init_type)
      CASE("GEM","GEM9")
         u(2,:,:)=u(2,:,:) - alpha*COS(twopi/lx*x)*COS(pi/ly*y)
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
      USE SolarHallMHD_mod
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
      CASE("breakout","half_breakout")
         top%static=.TRUE.
         top%bc_type(2)="robin+"
         top%bc_type(7)="robin"

         bottom%static=.TRUE.
         bottom%bc_type(2)="robin+"
         bottom%bc_type(7)="robin+"

         left%static(1)=.TRUE.
         left%bc_type(2)="robin+"
         left%bc_type(3)="normflux"
         left%static(4:8)=.TRUE.
         left%bc_type(7)="natural"

         right%static(1)=.TRUE.
         right%bc_type(2)="robin+"
         right%static(3:8)=.TRUE.
         right%bc_type(7)="robin+"
      CASE("GEM")
         top%static(2:8)=.TRUE.
         top%bc_type(1)="natural"
         top%bc_type(4)="zeroflux"
         top%bc_type(6)="zeroflux"
         top%bc_type(8)="zeroflux"
         bottom%static=.TRUE.
         left%static=.TRUE.
         right%static=.TRUE.
      CASE("GEM9")
         top%static(2:9)=.TRUE.
         top%bc_type(1)="natural"
         top%bc_type(4)="zeroflux"
         top%bc_type(6)="zeroflux"
         top%bc_type(8)="zeroflux"
         bottom%static=.TRUE.
         left%static=.TRUE.
         right%static=.TRUE.
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
      USE SolarHallMHD_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci
c-----------------------------------------------------------------------
c     cylindrical to cartesian relationships:
c     1: z --> x
c     2: r --> y
c     3: phi --> z
c-----------------------------------------------------------------------
      r_fac=y**cyl_fac
      r_faci=y**(-cyl_fac)
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
c-----------------------------------------------------------------------
c     boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("breakout")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=u(7,:,:)
            c(8,:,:)=uy(8,:,:)
         CASE("left")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:)
     $           + (one - mu_fac)*u(1,:,:)
            c(3,:,:)=nhat(2,:,:)*(cyl_fac*u(3,:,:)+r_fac*uy(3,:,:))
     $           + nhat(1,:,:)*r_fac*ux(3,:,:)
            c(4,:,:)=(nhat(1,:,:)*nhat(2,:,:)*(ux(4,:,:) - uy(5,:,:))
     $           + nhat(2,:,:)**2*uy(4,:,:) - nhat(1,:,:)**2*ux(5,:,:))
     $           *u(1,:,:)
     $           - (nhat(1,:,:)*nhat(2,:,:)
     $           *(ux(1,:,:)*u(4,:,:) - uy(1,:,:)*u(5,:,:))
     $           + nhat(2,:,:)**2*uy(1,:,:)*u(4,:,:) 
     $           - nhat(1,:,:)**2*ux(1,:,:)*u(5,:,:))
            c(5,:,:)=nhat(1,:,:)*u(4,:,:) + nhat(2,:,:)*u(5,:,:)
            c(6,:,:)=u(6,:,:)/u(1,:,:)
            IF(t <= v_period)THEN
               WHERE(y > 0. .AND. ABS(ATAN(x/y)) < v_angle)
                  c(6,:,:)=u(6,:,:)/u(1,:,:)
     $                 - v_peak*(ATAN(x/y)**2 - v_angle**2)**2
     $                 *x/SQRT(x**2+y**2)
     $                 *half*(one-COS(twopi/v_period*t))
               END WHERE
            ENDIF
            c(7,:,:)=ux(7,:,:)*nhat(1,:,:) + uy(7,:,:)*nhat(2,:,:)
            c(8,:,:)=nhat(1,:,:)*ux(8,:,:) + nhat(2,:,:)*uy(8,:,:)
     $           - gravity*u(1,:,:)
         CASE("right")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:)
     $           - (one - mu_fac)*u(1,:,:)/Rmax
            c(3,:,:)=nhat(2,:,:)*(cyl_fac*u(3,:,:)+r_fac*uy(3,:,:))
     $           + nhat(1,:,:)*r_fac*ux(3,:,:)
            c(4,:,:)=nhat(1,:,:)*(ux(4,:,:)*u(1,:,:)-ux(1,:,:)*u(4,:,:))
     $           + nhat(2,:,:)*(uy(4,:,:)*u(1,:,:)-uy(1,:,:)*u(4,:,:))
            c(5,:,:)=nhat(1,:,:)*(ux(5,:,:)*u(1,:,:)-ux(1,:,:)*u(5,:,:))
     $           + nhat(2,:,:)*(uy(5,:,:)*u(1,:,:)-uy(1,:,:)*u(5,:,:))
            c(6,:,:)=nhat(1,:,:)*(ux(6,:,:)*u(1,:,:)-ux(1,:,:)*u(6,:,:))
     $           + nhat(2,:,:)*(uy(6,:,:)*u(1,:,:)-uy(1,:,:)*u(6,:,:))
            c(7,:,:)=u(7,:,:)
            c(8,:,:)=nhat(1,:,:)*ux(8,:,:) + nhat(2,:,:)*uy(8,:,:)
     $           + gravity*u(1,:,:)/Rmax**2
         END SELECT
      CASE("half_breakout")
         SELECT CASE(lrtb)
         CASE("bottom")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=ux(2,:,:)
            c(3,:,:)=ux(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=ux(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=ux(7,:,:)
            c(8,:,:)=ux(8,:,:)
         CASE("top")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=u(7,:,:)
            c(8,:,:)=uy(8,:,:)
         CASE("left")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:)
     $           + (one - mu_fac)*u(1,:,:)
            IF(t <= v_period)THEN
               WHERE(y > 0. .AND. ABS(ATAN(x/y)) < v_angle)
                  c(3,:,:) = v_peak*x**2*(v_angle**2 - ATAN(x/y)**2)**2
     $                 *(one - COS(twopi/v_period*t))
     $                 *(one + two*alpha*(one + x**2 - 4._r8*y**2))
               END WHERE
            ENDIF
            c(4,:,:)=(nhat(1,:,:)*nhat(2,:,:)*(ux(4,:,:) - uy(5,:,:))
     $           + nhat(2,:,:)**2*uy(4,:,:) - nhat(1,:,:)**2*ux(5,:,:))
     $           *u(1,:,:)
     $           - (nhat(1,:,:)*nhat(2,:,:)
     $           *(ux(1,:,:)*u(4,:,:) - uy(1,:,:)*u(5,:,:))
     $           + nhat(2,:,:)**2*uy(1,:,:)*u(4,:,:) 
     $           - nhat(1,:,:)**2*ux(1,:,:)*u(5,:,:))
            c(5,:,:)=nhat(1,:,:)*u(4,:,:) + nhat(2,:,:)*u(5,:,:)
            c(6,:,:)=nhat(1,:,:)*(ux(6,:,:)*u(1,:,:)-ux(1,:,:)*u(6,:,:))
     $           + nhat(2,:,:)*(uy(6,:,:)*u(1,:,:)-uy(1,:,:)*u(6,:,:))
            c(7,:,:)=ux(7,:,:)*nhat(1,:,:) + uy(7,:,:)*nhat(2,:,:)
            c(8,:,:)=nhat(1,:,:)*ux(8,:,:) + nhat(2,:,:)*uy(8,:,:)
     $           - gravity*u(1,:,:)
     $           + u(3,:,:)*(nhat(1,:,:)*ux(3,:,:)
     $           + nhat(2,:,:)*(cyl_fac*r_faci*u(3,:,:)+uy(3,:,:)))
     $           + u(7,:,:)*(nhat(1,:,:)*ux(2,:,:)
     $           + nhat(2,:,:)*(cyl_fac*r_faci*u(2,:,:)+uy(2,:,:)))

         CASE("right")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:)
     $           - (one - mu_fac)*u(1,:,:)/Rmax
            c(3,:,:)=nhat(2,:,:)*(cyl_fac*u(3,:,:)+r_fac*uy(3,:,:))
     $           + nhat(1,:,:)*r_fac*ux(3,:,:)
            c(4,:,:)=(nhat(1,:,:)*nhat(2,:,:)*(ux(4,:,:) - uy(5,:,:))
     $           + nhat(2,:,:)**2*uy(4,:,:) - nhat(1,:,:)**2*ux(5,:,:))
     $           *u(1,:,:)
     $           - (nhat(1,:,:)*nhat(2,:,:)
     $           *(ux(1,:,:)*u(4,:,:) - uy(1,:,:)*u(5,:,:))
     $           + nhat(2,:,:)**2*uy(1,:,:)*u(4,:,:) 
     $           - nhat(1,:,:)**2*ux(1,:,:)*u(5,:,:))
            c(5,:,:)=nhat(1,:,:)*u(4,:,:) + nhat(2,:,:)*u(5,:,:)
            c(6,:,:)=nhat(1,:,:)*(ux(6,:,:)*u(1,:,:)-ux(1,:,:)*u(6,:,:))
     $           + nhat(2,:,:)*(uy(6,:,:)*u(1,:,:)-uy(1,:,:)*u(6,:,:))
            c(7,:,:)=u(7,:,:)
            c(8,:,:)=u(8,:,:)-(T0*mu0*n0*k_B/b0**2)*u(1,:,:)/Rmax

         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top")
            c(2,:,:)=u(2,:,:)+lambda*LOG(COSH(y/lambda))
            c(3,:,:)=-di*(u(3,:,:)*ux(3,:,:) + Te_frac*ux(8,:,:)) 
     $           + eta*uy(3,:,:)*u(1,:,:)
            c(4,:,:)=uy(4,:,:)*u(1,:,:) - uy(1,:,:)*u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=uy(6,:,:)*u(1,:,:) - uy(1,:,:)*u(6,:,:)
            c(7,:,:)=uy(7,:,:)*u(1,:,:) - uy(1,:,:)*u(7,:,:)
            c(8,:,:)=uy(8,:,:)*u(1,:,:) - uy(1,:,:)*u(8,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=uy(6,:,:)
            c(7,:,:)=uy(7,:,:)
            c(8,:,:)=uy(8,:,:)
         CASE("left","right")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=ux(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=ux(5,:,:)
            c(6,:,:)=ux(6,:,:)
            c(7,:,:)=ux(7,:,:)
            c(8,:,:)=ux(8,:,:)
         END SELECT
      CASE("GEM9")
         SELECT CASE(lrtb)
         CASE("top")
            c(2,:,:)=u(2,:,:)+lambda*LOG(COSH(y/lambda))
            c(3,:,:)=-di*(u(3,:,:)*ux(3,:,:) + Te_frac*ux(8,:,:)) 
     $           + eta*uy(3,:,:)*u(1,:,:) - di**2*nu*uy(9,:,:)
            c(4,:,:)=uy(4,:,:)*u(1,:,:) - uy(1,:,:)*u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=uy(6,:,:)*u(1,:,:) - uy(1,:,:)*u(6,:,:)
            c(7,:,:)=uy(7,:,:)*u(1,:,:) - uy(1,:,:)*u(7,:,:)
            c(8,:,:)=uy(8,:,:)*u(1,:,:) - uy(1,:,:)*u(8,:,:)
            c(9,:,:)=u(9,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=uy(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=uy(6,:,:)
            c(7,:,:)=uy(7,:,:)
            c(8,:,:)=uy(8,:,:)
            c(9,:,:)=u(9,:,:)
         CASE("left","right")
            c(1,:,:)=ux(1,:,:)
            c(2,:,:)=ux(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=ux(5,:,:)
            c(6,:,:)=ux(6,:,:)
            c(7,:,:)=ux(7,:,:)
            c(8,:,:)=ux(8,:,:)
            c(9,:,:)=u(9,:,:)
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
      USE SolarHallMHD_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci
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
c     cylindrical to cartesian relationships:
c     1: z --> x
c     2: r --> y
c     3: phi --> z
c-----------------------------------------------------------------------
      r_fac=y**cyl_fac
      r_faci=y**(-cyl_fac)
c-----------------------------------------------------------------------
c     boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("breakout")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_uy(1,1,:,:) = one
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_uy(4,4,:,:) = one
            c_u(5,5,:,:) = one
            c_u(6,6,:,:) = one
            c_u(7,7,:,:) = one
            c_uy(8,8,:,:) = one
         CASE("left")
            c_u(1,1,:,:) = (one - mu_fac)
            c_ux(1,1,:,:) = nhat(1,:,:)
            c_uy(1,1,:,:) = nhat(2,:,:)
            c_u(3,3,:,:) = nhat(2,:,:)*cyl_fac
            c_ux(3,3,:,:) = nhat(1,:,:)*r_fac
            c_uy(3,3,:,:) = nhat(2,:,:)*r_fac
            c_u(4,1,:,:) = nhat(1,:,:)*nhat(2,:,:)*(ux(4,:,:)-uy(5,:,:))
     $           + nhat(2,:,:)**2*uy(4,:,:) - nhat(1,:,:)**2*ux(5,:,:)
            c_u(4,4,:,:) = -nhat(1,:,:)*nhat(2,:,:)*ux(1,:,:)
     $           - nhat(2,:,:)**2*uy(1,:,:)
            c_u(4,5,:,:) = nhat(1,:,:)*nhat(2,:,:)*uy(1,:,:)
     $           + nhat(1,:,:)**2*ux(1,:,:)
            c_ux(4,1,:,:) = -nhat(1,:,:)*nhat(2,:,:)*u(4,:,:)
     $           + nhat(1,:,:)**2*u(5,:,:)
            c_ux(4,4,:,:) = nhat(1,:,:)*nhat(2,:,:)*u(1,:,:)
            c_ux(4,5,:,:) = -nhat(1,:,:)**2*u(1,:,:)
            c_uy(4,1,:,:) = nhat(1,:,:)*nhat(2,:,:)*u(5,:,:)
     $           - nhat(2,:,:)**2*u(4,:,:)
            c_uy(4,4,:,:) = nhat(2,:,:)**2*u(1,:,:)
            c_uy(4,5,:,:) = -nhat(1,:,:)*nhat(2,:,:)*u(1,:,:)
            c_u(5,4,:,:) = nhat(1,:,:)
            c_u(5,5,:,:) = nhat(2,:,:)
            c_u(6,1,:,:) = -u(6,:,:)/u(1,:,:)**2 
            c_u(6,6,:,:) = one/u(1,:,:)
            c_ux(7,7,:,:) = nhat(1,:,:)
            c_uy(7,7,:,:) = nhat(2,:,:)
            c_u(8,1,:,:) = -gravity
            c_ux(8,8,:,:) = nhat(1,:,:)
            c_uy(8,8,:,:) = nhat(2,:,:)
         CASE("right")
            c_u(1,1,:,:) = -(one - mu_fac)/Rmax
            c_ux(1,1,:,:) = nhat(1,:,:)
            c_uy(1,1,:,:) = nhat(2,:,:)
            c_u(3,3,:,:) = nhat(2,:,:)*cyl_fac
            c_ux(3,3,:,:) = nhat(1,:,:)*r_fac
            c_uy(3,3,:,:) = nhat(2,:,:)*r_fac
            c_u(4,1,:,:) = nhat(1,:,:)*ux(4,:,:)+nhat(2,:,:)*uy(4,:,:)
            c_u(4,4,:,:) =-nhat(1,:,:)*ux(1,:,:)-nhat(2,:,:)*uy(1,:,:)
            c_ux(4,1,:,:) = -nhat(1,:,:)*u(4,:,:)
            c_uy(4,1,:,:) = -nhat(2,:,:)*u(4,:,:)
            c_ux(4,4,:,:) = nhat(1,:,:)*u(1,:,:)
            c_uy(4,4,:,:) = nhat(2,:,:)*u(1,:,:)
            c_u(5,1,:,:) = nhat(1,:,:)*ux(5,:,:)+nhat(2,:,:)*uy(5,:,:)
            c_u(5,5,:,:) =-nhat(1,:,:)*ux(1,:,:)-nhat(2,:,:)*uy(1,:,:)
            c_ux(5,1,:,:) = -nhat(1,:,:)*u(5,:,:)
            c_uy(5,1,:,:) = -nhat(2,:,:)*u(5,:,:)
            c_ux(5,5,:,:) = nhat(1,:,:)*u(1,:,:)
            c_uy(5,5,:,:) = nhat(2,:,:)*u(1,:,:)
            c_u(6,1,:,:) = nhat(1,:,:)*ux(6,:,:)+nhat(2,:,:)*uy(6,:,:)
            c_u(6,6,:,:) =-nhat(1,:,:)*ux(1,:,:)-nhat(2,:,:)*uy(1,:,:)
            c_ux(6,1,:,:) = -nhat(1,:,:)*u(6,:,:)
            c_uy(6,1,:,:) = -nhat(2,:,:)*u(6,:,:)
            c_ux(6,6,:,:) = nhat(1,:,:)*u(1,:,:)
            c_uy(6,6,:,:) = nhat(2,:,:)*u(1,:,:)
            c_u(7,7,:,:) = one
            c_u(8,1,:,:) = gravity/Rmax**2
            c_ux(8,8,:,:) = nhat(1,:,:)
            c_uy(8,8,:,:) = nhat(2,:,:)
         END SELECT
      CASE("half_breakout")
         SELECT CASE(lrtb)
         CASE("bottom")
            c_ux(1,1,:,:) = one
            c_ux(2,2,:,:) = one
            c_ux(3,3,:,:) = one
            c_u(4,4,:,:) = one
            c_ux(5,5,:,:) = one
            c_u(6,6,:,:) = one
            c_ux(7,7,:,:) = one
            c_ux(8,8,:,:) = one
         CASE("top")
            c_uy(1,1,:,:) = one
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_uy(4,4,:,:) = one
            c_u(5,5,:,:) = one
            c_u(6,6,:,:) = one
            c_u(7,7,:,:) = one
            c_uy(8,8,:,:) = one
         CASE("left")
            c_u(1,1,:,:) = (one - mu_fac)
            c_ux(1,1,:,:) = nhat(1,:,:)
            c_uy(1,1,:,:) = nhat(2,:,:)
            c_u(4,1,:,:) = nhat(1,:,:)*nhat(2,:,:)*(ux(4,:,:)-uy(5,:,:))
     $           + nhat(2,:,:)**2*uy(4,:,:) - nhat(1,:,:)**2*ux(5,:,:)
            c_u(4,4,:,:) = -nhat(1,:,:)*nhat(2,:,:)*ux(1,:,:)
     $           - nhat(2,:,:)**2*uy(1,:,:)
            c_u(4,5,:,:) = nhat(1,:,:)*nhat(2,:,:)*uy(1,:,:)
     $           + nhat(1,:,:)**2*ux(1,:,:)
            c_ux(4,1,:,:) = -nhat(1,:,:)*nhat(2,:,:)*u(4,:,:)
     $           + nhat(1,:,:)**2*u(5,:,:)
            c_ux(4,4,:,:) = nhat(1,:,:)*nhat(2,:,:)*u(1,:,:)
            c_ux(4,5,:,:) = -nhat(1,:,:)**2*u(1,:,:)
            c_uy(4,1,:,:) = nhat(1,:,:)*nhat(2,:,:)*u(5,:,:)
     $           - nhat(2,:,:)**2*u(4,:,:)
            c_uy(4,4,:,:) = nhat(2,:,:)**2*u(1,:,:)
            c_uy(4,5,:,:) = -nhat(1,:,:)*nhat(2,:,:)*u(1,:,:)
            c_u(5,4,:,:) = nhat(1,:,:)
            c_u(5,5,:,:) = nhat(2,:,:)
            c_u(6,1,:,:) = nhat(1,:,:)*ux(6,:,:)+nhat(2,:,:)*uy(6,:,:)
            c_u(6,6,:,:) =-nhat(1,:,:)*ux(1,:,:)-nhat(2,:,:)*uy(1,:,:)
            c_ux(6,1,:,:) = -nhat(1,:,:)*u(6,:,:)
            c_uy(6,1,:,:) = -nhat(2,:,:)*u(6,:,:)
            c_ux(6,6,:,:) = nhat(1,:,:)*u(1,:,:)
            c_uy(6,6,:,:) = nhat(2,:,:)*u(1,:,:)
            c_ux(7,7,:,:) = nhat(1,:,:)
            c_uy(7,7,:,:) = nhat(2,:,:)
            c_u(8,1,:,:) = -gravity
            c_u(8,2,:,:) = u(7,:,:)*nhat(2,:,:)*cyl_fac*r_faci
            c_u(8,3,:,:) = nhat(1,:,:)*ux(3,:,:)
     $           + nhat(2,:,:)*(two*cyl_fac*r_faci*u(3,:,:)+uy(3,:,:))
            c_u(8,7,:,:) = (nhat(1,:,:)*ux(2,:,:)
     $           + nhat(2,:,:)*(cyl_fac*r_faci*u(2,:,:)+uy(2,:,:)))
            c_ux(8,2,:,:) = u(7,:,:)*nhat(1,:,:)
            c_ux(8,3,:,:) = u(3,:,:)*nhat(1,:,:)
            c_ux(8,8,:,:) = nhat(1,:,:)
            c_uy(8,2,:,:) = u(7,:,:)*nhat(2,:,:)
            c_uy(8,3,:,:) = u(3,:,:)*nhat(2,:,:)
            c_uy(8,8,:,:) = nhat(2,:,:)

         CASE("right")
            c_u(1,1,:,:) = -(one - mu_fac)/Rmax
            c_ux(1,1,:,:) = nhat(1,:,:)
            c_uy(1,1,:,:) = nhat(2,:,:)

            c_u(3,3,:,:) = nhat(2,:,:)*cyl_fac
            c_ux(3,3,:,:) = nhat(1,:,:)*r_fac
            c_uy(3,3,:,:) = nhat(2,:,:)*r_fac

            c_u(4,1,:,:) = nhat(1,:,:)*nhat(2,:,:)*(ux(4,:,:)-uy(5,:,:))
     $           + nhat(2,:,:)**2*uy(4,:,:) - nhat(1,:,:)**2*ux(5,:,:)
            c_u(4,4,:,:) = -nhat(1,:,:)*nhat(2,:,:)*ux(1,:,:)
     $           - nhat(2,:,:)**2*uy(1,:,:)
            c_u(4,5,:,:) = nhat(1,:,:)*nhat(2,:,:)*uy(1,:,:)
     $           + nhat(1,:,:)**2*ux(1,:,:)
            c_ux(4,1,:,:) = -nhat(1,:,:)*nhat(2,:,:)*u(4,:,:)
     $           + nhat(1,:,:)**2*u(5,:,:)
            c_ux(4,4,:,:) = nhat(1,:,:)*nhat(2,:,:)*u(1,:,:)
            c_ux(4,5,:,:) = -nhat(1,:,:)**2*u(1,:,:)
            c_uy(4,1,:,:) = nhat(1,:,:)*nhat(2,:,:)*u(5,:,:)
     $           - nhat(2,:,:)**2*u(4,:,:)
            c_uy(4,4,:,:) = nhat(2,:,:)**2*u(1,:,:)
            c_uy(4,5,:,:) = -nhat(1,:,:)*nhat(2,:,:)*u(1,:,:)
            c_u(5,4,:,:) = nhat(1,:,:)
            c_u(5,5,:,:) = nhat(2,:,:)

            c_u(6,1,:,:) = nhat(1,:,:)*ux(6,:,:)+nhat(2,:,:)*uy(6,:,:)
            c_u(6,6,:,:) =-nhat(1,:,:)*ux(1,:,:)-nhat(2,:,:)*uy(1,:,:)
            c_ux(6,1,:,:) = -nhat(1,:,:)*u(6,:,:)
            c_uy(6,1,:,:) = -nhat(2,:,:)*u(6,:,:)
            c_ux(6,6,:,:) = nhat(1,:,:)*u(1,:,:)
            c_uy(6,6,:,:) = nhat(2,:,:)*u(1,:,:)

            c_u(7,7,:,:) = one

            c_u(8,1,:,:) = -(T0*mu0*n0*k_B/b0**2)/Rmax
            c_u(8,8,:,:) = one
         END SELECT
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top")
            c_u(2,2,:,:)=one
            c_u(3,1,:,:)=eta*uy(3,:,:)
            c_u(3,3,:,:)=-di*ux(3,:,:)
            c_ux(3,3,:,:)=-di*u(3,:,:)
            c_ux(3,8,:,:)=-di*Te_frac
            c_uy(3,3,:,:)=eta*u(1,:,:)
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
            c_u(8,1,:,:)=uy(8,:,:)
            c_u(8,8,:,:)=-uy(1,:,:)
            c_uy(8,1,:,:)=-u(8,:,:)
            c_uy(8,8,:,:)=u(1,:,:)
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_uy(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_uy(6,6,:,:)=one
            c_uy(7,7,:,:)=one
            c_uy(8,8,:,:)=one
         CASE("left","right")
            c_ux(1,1,:,:)=one
            c_ux(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_ux(5,5,:,:)=one
            c_ux(6,6,:,:)=one
            c_ux(7,7,:,:)=one
            c_ux(8,8,:,:)=one
        END SELECT
      CASE("GEM9")
         SELECT CASE(lrtb)
         CASE("top")
            c_u(2,2,:,:)=one
            c_u(3,1,:,:)=eta*uy(3,:,:)
            c_u(3,3,:,:)=-di*ux(3,:,:)
            c_ux(3,3,:,:)=-di*u(3,:,:)
            c_ux(3,8,:,:)=-di*Te_frac
            c_uy(3,3,:,:)=eta*u(1,:,:)
            c_uy(3,9,:,:)=-di**2*nu
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
            c_u(8,1,:,:)=uy(8,:,:)
            c_u(8,8,:,:)=-uy(1,:,:)
            c_uy(8,1,:,:)=-u(8,:,:)
            c_uy(8,8,:,:)=u(1,:,:)
            c_u(9,9,:,:)=one
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_uy(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_uy(6,6,:,:)=one
            c_uy(7,7,:,:)=one
            c_uy(8,8,:,:)=one
            c_u(9,9,:,:)=one
         CASE("left","right")
            c_ux(1,1,:,:)=one
            c_ux(2,2,:,:)=one
            c_u(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_ux(5,5,:,:)=one
            c_ux(6,6,:,:)=one
            c_ux(7,7,:,:)=one
            c_ux(8,8,:,:)=one
            c_u(9,9,:,:)=one
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
      USE SolarHallMHD_mod
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
      CASE("breakout","half_breakout")
         SELECT CASE(lrtb)
         CASE("left","right")
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
      USE SolarHallMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Tix,Tiy,kperp,kfac,
     $     ve,vex,vey,r_fac,r_faci,j1,j2,j3,jtot,eta_local,b1,b2,Bsq,
     $     n_inv,nx_inv,ny_inv,gx_vec,gy_vec,r_sphr
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
c     compute components of gravitational vector.
c-----------------------------------------------------------------------
      gx_vec=gravity*x/(x**2+y**2)**1.5
      gy_vec=gravity*y/(x**2+y**2)**1.5
      WHERE((x**2 + y**2) == 0)
         gx_vec = zero
         gy_vec = zero
      END WHERE
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
      r_sphr=SQRT(x**2 + y**2)-one
      SELECT CASE(eta_case)
      CASE("spitzer")
         CALL transport_seteta("spitzer-chodura",r_sphr,u(1,:,:),
     $        u(8,:,:),jtot,0._r8,etas_norm,0._r8,0._r8,
     $        0._r8,eta,etavac,eta_local)
      CASE("spitzer-chodura")
         CALL transport_seteta(eta_case,r_sphr,u(1,:,:),
     $        u(8,:,:),jtot,chod_const,etas_norm,etac_norm,v_chod_norm,
     $        0._r8,eta,etavac,eta_local)
      CASE DEFAULT
         eta_local=eta
      END SELECT
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT(b1,b2,u(3,:,:),Tix,Tiy,BdotT)
         CALL transport_kbrag(u(1,:,:),u(8,:,:),Te_frac,Bsq,
     $        ke_norm,ki_norm,xe_norm,xi_norm,kappa_min,kappa_max,
     $        kperp,kfac)
      CASE("anisotropic")
         CALL transport_BdotT(b1,b2,u(3,:,:),Tix,Tiy,BdotT)
         CALL transport_setkaniso(kappa_par,kappa_perp,Bsq,kperp,kfac)
         kperp=kperp*u(1,:,:)
         kfac=kfac*u(1,:,:)
      CASE("scalar")
         BdotT=0
         kperp=kappa_par
         kfac=0
      CASE DEFAULT
         BdotT=0
         kperp=kappa_par*u(1,:,:)
         kfac=0
      END SELECT
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
      fx(2,:,:) = -r_fac*di*nu*vex
      fy(2,:,:) = -r_fac*di*nu*vey
      
      s(2,:,:) = r_fac*(vi(2,:,:)*b1 - vi(1,:,:)*b2
     $     - di*(j2*b1 - j1*b2)*n_inv + eta_local*j3
     $     + di*nu*n_inv*(vex*ux(1,:,:) + vey*uy(1,:,:)))
     $     - cyl_fac*r_faci*di*nu*ve
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
      fx(3,:,:) = u(3,:,:)*vi(1,:,:) - ve*b1
     $     - di*Te_frac*uy(8,:,:)*n_inv + eta_local*j2

      fy(3,:,:) = u(3,:,:)*vi(2,:,:) - ve*b2
     $     + di*Te_frac*ux(8,:,:)*n_inv - eta_local*j1

      s(3,:,:) = di*u(3,:,:)*(j2*ny_inv + j1*nx_inv) 
     $     + cyl_fac*two*r_faci*di*u(3,:,:)*ux(3,:,:)*n_inv

      IF(init_type=="GEM9")THEN
         fx(3,:,:)=fx(3,:,:) + di**2*nu*ux(9,:,:)*n_inv
         fy(3,:,:)=fy(3,:,:) 
     $        + di**2*nu*(cyl_fac*r_faci*u(9,:,:) + uy(9,:,:))*n_inv
      ENDIF
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(4,:,:) = r_fac*(u(4,:,:)*vi(1,:,:) - two*mu*u(1,:,:)*vix(1,:,:) 
     $     + half*u(3,:,:)**2 + u(8,:,:))
         
      fy(4,:,:) = r_fac*(u(4,:,:)*vi(2,:,:) 
     $     - mu*u(1,:,:)*(viy(1,:,:) + vix(2,:,:)))

      s(4,:,:) = -r_fac*(j3*b2 + gx_vec*u(1,:,:))
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(5,:,:) = r_fac*(u(5,:,:)*vi(1,:,:) 
     $     - mu*u(1,:,:)*(vix(2,:,:) + viy(1,:,:)))
      
      fy(5,:,:) = r_fac*(u(5,:,:)*vi(2,:,:) + half*u(3,:,:)**2
     $     - two*mu*u(1,:,:)*viy(2,:,:))
      
      s(5,:,:) = r_fac*(j3*b1 - uy(8,:,:) - gy_vec*u(1,:,:)) 
     $     + cyl_fac*(u(6,:,:)*vi(3,:,:) - half*u(3,:,:)**2 
     $     - two*mu*u(5,:,:)*r_faci)
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(6,:,:) = r_fac*(u(6,:,:)*vi(1,:,:) 
     $     - u(1,:,:)*(nu*vex + mu*vix(3,:,:)))
      
      fy(6,:,:) = r_fac*(u(6,:,:)*vi(2,:,:) 
     $     - u(1,:,:)*(nu*vey + mu*viy(3,:,:)))

      s(6,:,:) = r_fac*(j1*b2 - j2*b1) - cyl_fac*(u(6,:,:)*vi(2,:,:) 
     $     + r_faci*((mu+nu)*u(6,:,:) - di*nu*u(7,:,:)))
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
     $     + u(1,:,:)*(mu*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:)
     $     + vix(3,:,:)**2 + viy(3,:,:)**2) + nu*(vex**2 + vey**2)))
     $     + cyl_fac*r_faci*(mu*(two*vi(2,:,:)*u(5,:,:) 
     $     + vi(3,:,:)*u(6,:,:)) + nu*ve*(u(6,:,:)-di*u(7,:,:)))
c-----------------------------------------------------------------------
c     del^2(Bz).
c-----------------------------------------------------------------------
      IF(init_type=="GEM9")THEN
         fx(9,:,:)=ux(3,:,:)
         fy(9,:,:)=cyl_fac*r_faci*u(3,:,:) + uy(3,:,:)
         s(9,:,:)=u(9,:,:)
      ENDIF
c-----------------------------------------------------------------------
c     initial equilibrium source term.
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("breakout","half_breakout","GEM")
            CALL SolarHallMHD_equil(x,y,u0,u0x,u0y,.TRUE.)
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
      USE SolarHallMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci,Tix,Tiy,
     $     Tix_un,Tiy_un,Ti_un,kperp,kfac,ve,vex,vey,j1,j2,j3,
     $     jtot,eta_local,b1,b2,Bsq,eta_rho,eta_p,eta_j,j_u3,j_ux3,
     $     j_uy3,j_u7,n_inv,nx_inv,ny_inv,ve_un,vex_un,vey_un,
     $     gx_vec,gy_vec,r_sphr,kpar_un,kpar_p,kperp_un,kperp_p,
     $     kperp_bsq
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
c     compute components of gravitational vector.
c-----------------------------------------------------------------------
      gx_vec=gravity*x/(x**2+y**2)**1.5
      gy_vec=gravity*y/(x**2+y**2)**1.5
      WHERE((x**2 + y**2) == 0)
         gx_vec = zero
         gy_vec = zero
      END WHERE
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
      r_sphr=SQRT(x**2 + y**2)-one
      SELECT CASE(eta_case)
      CASE("spitzer")
         CALL transport_seteta_u("spitzer-chodura",r_sphr,u(1,:,:),
     $        u(8,:,:),jtot,0._r8,etas_norm,0._r8,0._r8,0._r8,
     $        eta,etavac,eta_local,eta_rho,eta_p,eta_j)
      CASE("spitzer-chodura")
         CALL transport_seteta_u(eta_case,r_sphr,u(1,:,:),
     $        u(8,:,:),jtot,chod_const,etas_norm,etac_norm,v_chod_norm,
     $        0._r8,eta,etavac,eta_local,eta_rho,eta_p,eta_j)
         WHERE(eta_local == etavac)
            eta_rho=0.
            eta_p=0.
            eta_j=0.
         END WHERE
      CASE DEFAULT
         eta_local=eta
         eta_rho=0.
         eta_p=0.
         eta_j=0.
      END SELECT
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT_u(b1,b2,u(3,:,:),Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty)
         CALL transport_kbrag_u(u(1,:,:),u(8,:,:),Te_frac,Bsq,
     $        ke_norm,ki_norm,xe_norm,xi_norm,kappa_min,kappa_max,
     $        kperp,kfac,kpar_un,kpar_p,kperp_un,kperp_p,kperp_bsq)
      CASE("anisotropic")
         CALL transport_BdotT_u(b1,b2,u(3,:,:),Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty)
         CALL transport_setkaniso_u(kappa_par,kappa_perp,Bsq,kperp,kfac,
     $        kperp_bsq)
         kperp_un=kperp
         kperp_p=0
         kperp=kperp*u(1,:,:)
         kperp_bsq=kperp_bsq*u(1,:,:)
         kfac=kfac*u(1,:,:)
         kpar_un=kappa_par
         kpar_p=0
      CASE("scalar")
         BdotT=0
         BdotT_b1=0
         BdotT_b2=0
         BdotT_b3=0
         BdotT_Tx=0
         BdotT_Ty=0
         kperp=kappa_par
         kfac=0
         kpar_un=0
         kpar_p=0
         kperp_un=0
         kperp_p=0
         kperp_bsq=0
      CASE DEFAULT
         BdotT=0
         BdotT_b1=0
         BdotT_b2=0
         BdotT_b3=0
         BdotT_Tx=0
         BdotT_Ty=0
         kperp=kappa_par*u(1,:,:)
         kfac=0
         kpar_un=0
         kpar_p=0
         kperp_un=kappa_par
         kperp_p=0
         kperp_bsq=0
      END SELECT
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
      fx_u(2,1,:,:)=-r_fac*di*nu*vex_un
      fx_u(2,6,:,:)=-r_fac*di*nu*nx_inv
      fx_u(2,7,:,:)= r_fac*di**2*nu*nx_inv
      fx_ux(2,1,:,:)=-r_fac*di*nu*ve_un
      fx_ux(2,6,:,:)=-r_fac*di*nu*n_inv
      fx_ux(2,7,:,:)= r_fac*di**2*nu*n_inv
      
      fy_u(2,1,:,:)=-r_fac*di*nu*vey_un
      fy_u(2,6,:,:)=-r_fac*di*nu*ny_inv
      fy_u(2,7,:,:)= r_fac*di**2*nu*ny_inv
      fy_uy(2,1,:,:)=-r_fac*di*nu*ve_un
      fy_uy(2,6,:,:)=-r_fac*di*nu*n_inv
      fy_uy(2,7,:,:)=r_fac*di**2*nu*n_inv
      
      s_u(2,1,:,:)=r_fac*(vi_un(2,:,:)*b1 - vi_un(1,:,:)*b2
     $     + di*(j2*b1 - j1*b2)*n_inv**2 + eta_rho*j3
     $     + di*nu*(n_inv*(vex_un*ux(1,:,:) + vey_un*uy(1,:,:))
     $     + (vex*nx_inv + vey*ny_inv)))
     $     - cyl_fac*r_faci*di*nu*ve_un
      s_u(2,2,:,:)=-cyl_fac*(vi(2,:,:)-di*j2*n_inv)
      s_u(2,3,:,:)=di*cyl_fac*b2*n_inv + r_fac*eta_j*j_u3*j3
      s_u(2,4,:,:)=-r_fac*b2*n_inv
      s_u(2,5,:,:)=r_fac*b1*n_inv
      s_u(2,6,:,:)= -r_fac*di*nu*u(1,:,:)*(nx_inv**2 + ny_inv**2)
     $     - cyl_fac*r_faci*di*nu*n_inv
      s_u(2,7,:,:)=r_fac*(di**2*nu*u(1,:,:)*(nx_inv**2 + ny_inv**2)
     $     + eta_j*j_u7*j3 + eta_local)
     $     + cyl_fac*r_faci*di**2*nu*n_inv
      s_u(2,8,:,:)=r_fac*eta_p*j3
      s_ux(2,1,:,:)=r_fac*di*nu*(vex + ve_un*ux(1,:,:))*n_inv
      s_ux(2,2,:,:)=r_fac*(di*j1*n_inv - vi(1,:,:))
      s_ux(2,3,:,:)=r_fac*(di*b1*n_inv + eta_j*j_ux3*j3)
      s_ux(2,6,:,:)=-r_fac*di*nu*nx_inv
      s_ux(2,7,:,:)=r_fac*di**2*nu*nx_inv
      s_uy(2,1,:,:)=r_fac*di*nu*(vey + ve_un*uy(1,:,:))*n_inv
      s_uy(2,2,:,:)=-r_fac*(vi(2,:,:)-di*j2*n_inv)
      s_uy(2,3,:,:)=r_fac*(di*b2*n_inv + eta_j*j_uy3*j3)
      s_uy(2,6,:,:)=-r_fac*di*nu*ny_inv
      s_uy(2,7,:,:)=r_fac*di**2*nu*ny_inv
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=u(3,:,:)*vi_un(1,:,:) - ve_un*b1 
     $     + di*Te_frac*uy(8,:,:)*n_inv**2 +  eta_rho*j2
      fx_u(3,2,:,:)=cyl_fac*r_faci*ve 
      fx_u(3,3,:,:)=vi(1,:,:) + eta_j*j_u3*j2
      fx_u(3,4,:,:)=u(3,:,:)*n_inv
      fx_u(3,6,:,:)=-n_inv*b1
      fx_u(3,7,:,:)=di*n_inv*b1 + eta_j*j_u7*j2
      fx_u(3,8,:,:)= eta_p*j2
      fx_ux(3,3,:,:)=-eta_local + eta_j*j_ux3*j2
      fx_uy(3,2,:,:)=ve
      fx_uy(3,3,:,:)=eta_j*j_uy3*j2
      fx_uy(3,8,:,:)=-di*Te_frac*n_inv
         
      fy_u(3,1,:,:)=u(3,:,:)*vi_un(2,:,:) - ve_un*b2 
     $     - di*Te_frac*ux(8,:,:)*n_inv**2 - eta_rho*j1
      fy_u(3,3,:,:)=vi(2,:,:) - eta_local*cyl_fac*r_faci - eta_j*j_u3*j1
      fy_u(3,5,:,:)=u(3,:,:)*n_inv
      fy_u(3,6,:,:)=-n_inv*b2
      fy_u(3,7,:,:)=di*n_inv*b2 - eta_j*j_u7*j1
      fy_u(3,8,:,:)= -eta_p*j1
      fy_ux(3,2,:,:)=-ve
      fy_ux(3,3,:,:)=-eta_j*j_ux3*j1
      fy_ux(3,8,:,:)=di*Te_frac*n_inv
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

      IF(init_type=="GEM9")THEN
         fx_u(3,1,:,:)=fx_u(3,1,:,:) - di**2*nu*ux(9,:,:)*n_inv**2
         fx_ux(3,9,:,:)=di**2*nu*n_inv
         fy_u(3,1,:,:)=fy_u(3,1,:,:) 
     $        - di**2*nu*(cyl_fac*r_faci*u(9,:,:) + uy(9,:,:))*n_inv**2
         fy_u(3,9,:,:)=di**2*nu*cyl_fac*r_faci*n_inv
         fy_uy(3,9,:,:)=di**2*nu*n_inv
      ENDIF
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(1,:,:) 
     $     - mu*two*vi(1,:,:)*ux(1,:,:)*n_inv)
      fx_u(4,3,:,:)=r_fac*u(3,:,:)
      fx_u(4,4,:,:)=r_fac*two*(vi(1,:,:) + mu*ux(1,:,:)*n_inv)
      fx_u(4,8,:,:)=r_fac
      fx_ux(4,1,:,:)= r_fac*mu*two*vi(1,:,:)
      fx_ux(4,4,:,:)=-r_fac*mu*two

      fy_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(2,:,:)
     $     - mu*(vi(1,:,:)*uy(1,:,:) + vi(2,:,:)*ux(1,:,:))*n_inv)
      fy_u(4,4,:,:)=r_fac*(vi(2,:,:) + mu*uy(1,:,:)*n_inv)
      fy_u(4,5,:,:)=r_fac*(vi(1,:,:) + mu*ux(1,:,:)*n_inv)
      fy_ux(4,1,:,:)= r_fac*mu*vi(2,:,:)
      fy_ux(4,5,:,:)=-r_fac*mu
      fy_uy(4,1,:,:)= r_fac*mu*vi(1,:,:)
      fy_uy(4,4,:,:)=-r_fac*mu

      s_u(4,1,:,:)=-r_fac*gx_vec
      s_u(4,7,:,:)=-r_fac*b2
      s_ux(4,2,:,:)=-r_fac*j3
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(1,:,:)
     $     - mu*(vi(2,:,:)*ux(1,:,:) + vi(1,:,:)*uy(1,:,:))*n_inv)
      fx_u(5,4,:,:)=r_fac*(vi(2,:,:) + mu*uy(1,:,:)*n_inv)
      fx_u(5,5,:,:)=r_fac*(vi(1,:,:) + mu*ux(1,:,:)*n_inv)
      fx_ux(5,1,:,:)= r_fac*mu*vi(2,:,:)
      fx_ux(5,5,:,:)=-r_fac*mu
      fx_uy(5,1,:,:)= r_fac*mu*vi(1,:,:)
      fx_uy(5,4,:,:)=-r_fac*mu
      
      fy_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(2,:,:) 
     $     - mu*two*vi(2,:,:)*uy(1,:,:)*n_inv)
      fy_u(5,3,:,:)=r_fac*u(3,:,:)
      fy_u(5,5,:,:)=r_fac*two*(vi(2,:,:) + mu*uy(1,:,:)*n_inv)
      fy_uy(5,1,:,:)= r_fac*mu*two*vi(2,:,:)
      fy_uy(5,5,:,:)=-r_fac*mu*two

      s_u(5,1,:,:)=-r_fac*gy_vec + cyl_fac*u(6,:,:)*vi_un(3,:,:)
      s_u(5,2,:,:)=-cyl_fac*j3
      s_u(5,3,:,:)=-cyl_fac*u(3,:,:)
      s_u(5,5,:,:)=-cyl_fac*r_faci*two*mu
      s_u(5,6,:,:)=cyl_fac*two*vi(3,:,:)
      s_u(5,7,:,:)=r_fac*b1
      s_uy(5,2,:,:)=-r_fac*j3
      s_uy(5,8,:,:)=-r_fac
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx_u(6,1,:,:)=r_fac*(u(6,:,:)*vi_un(1,:,:) 
     $     - (mu*vi(3,:,:) + nu*ve)*ux(1,:,:)*n_inv)
      fx_u(6,4,:,:)=r_fac*u(6,:,:)*n_inv
      fx_u(6,6,:,:)=r_fac*(vi(1,:,:) + (mu+nu)*ux(1,:,:)*n_inv)
      fx_u(6,7,:,:)= -r_fac*di*nu*ux(1,:,:)*n_inv
      fx_ux(6,1,:,:)= r_fac*(mu*vi(3,:,:) + nu*ve)
      fx_ux(6,6,:,:)=-r_fac*(mu+nu)
      fx_ux(6,7,:,:)= r_fac*di*nu
      
      fy_u(6,1,:,:)=r_fac*(u(6,:,:)*vi_un(2,:,:)
     $     - (mu*vi(3,:,:) + nu*ve)*uy(1,:,:)*n_inv)
      fy_u(6,5,:,:)=r_fac*u(6,:,:)*n_inv
      fy_u(6,6,:,:)=r_fac*(vi(2,:,:) + (mu+nu)*uy(1,:,:)*n_inv)
      fy_u(6,7,:,:)= -r_fac*di*nu*uy(1,:,:)*n_inv
      fy_uy(6,1,:,:)= r_fac*(mu*vi(3,:,:) + nu*ve)
      fy_uy(6,6,:,:)=-r_fac*(mu+nu)
      fy_uy(6,7,:,:)= r_fac*di*nu
      
      s_u(6,1,:,:)=cyl_fac*vi(2,:,:)*vi(3,:,:)
      s_u(6,2,:,:)=cyl_fac*j2
      s_u(6,3,:,:)=cyl_fac*b2
      s_u(6,5,:,:)=-cyl_fac*vi(3,:,:)
      s_u(6,6,:,:)=-cyl_fac*(vi(2,:,:) + r_faci*(mu+nu))
      s_u(6,7,:,:)= cyl_fac*r_faci*di*nu
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
     $     - kfac*BdotT_Ty(1,:,:)*Tiy_un
     $     - (kpar_un - kperp_un)*BdotT(1,:,:) - kperp_un*Tix)
      fx_u(8,2,:,:)=-cyl_fac*(-kfac*BdotT_b1(1,:,:)
     $     + two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_u(8,3,:,:)=r_fac*(-kfac*BdotT_b3(1,:,:)
     $     + two*u(3,:,:)*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_u(8,4,:,:)=r_fac*gamma_fac*u(8,:,:)*n_inv
      fx_u(8,8,:,:)=r_fac*(gamma_fac*vi(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*nx_inv
     $     - kfac*BdotT_Ty(1,:,:)*ny_inv
     $     - (kpar_p - kperp_p)*BdotT(1,:,:) - kperp_p*Tix)

      fx_ux(8,1,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*Ti_un
      fx_ux(8,2,:,:)=r_fac*(-kfac*BdotT_b2(1,:,:)
     $     + two*b2*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_ux(8,8,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*n_inv

      fx_uy(8,1,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*Ti_un
      fx_uy(8,2,:,:)=-r_fac*(-kfac*BdotT_b1(1,:,:)
     $     + two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_uy(8,8,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*n_inv

      fy_u(8,1,:,:)=r_fac*(gamma_fac*u(8,:,:)*vi_un(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*Tiy_un
     $     - kfac*BdotT_Tx(2,:,:)*Tix_un
     $     - (kpar_un - kperp_un)*BdotT(2,:,:) - kperp_un*Tiy)
      fy_u(8,2,:,:)=-cyl_fac*(-kfac*BdotT_b1(2,:,:)
     $     + two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(8,3,:,:)=r_fac*(-kfac*BdotT_b3(2,:,:)
     $     + two*u(3,:,:)*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(8,5,:,:)=r_fac*gamma_fac*u(8,:,:)*n_inv
      fy_u(8,8,:,:)=r_fac*(gamma_fac*vi(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*ny_inv
     $     - kfac*BdotT_Tx(2,:,:)*nx_inv
     $     - (kpar_p - kperp_p)*BdotT(2,:,:) - kperp_p*Tiy)

      fy_ux(8,1,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*Ti_un
      fy_ux(8,2,:,:)=r_fac*(-kfac*BdotT_b2(2,:,:)
     $     + two*b2*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_ux(8,8,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*n_inv

      fy_uy(8,1,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*Ti_un
      fy_uy(8,2,:,:)=-r_fac*(-kfac*BdotT_b1(2,:,:)
     $     + two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_uy(8,8,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*n_inv

      s_u(8,1,:,:)=r_fac*(vi_un(1,:,:)*ux(8,:,:)
     $     + vi_un(2,:,:)*uy(8,:,:) + eta_rho*jtot**2
     $     + mu*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:)
     $     + vix(3,:,:)**2 + viy(3,:,:)**2) + nu*(vex**2 + vey**2)
     $     + u(1,:,:)*(mu*two*(two*vix(1,:,:)*vix_un(1,:,:) 
     $     + two*viy(2,:,:)*viy_un(2,:,:) 
     $     + viy(1,:,:)*viy_un(1,:,:) + vix(2,:,:)*vix_un(2,:,:)
     $     + viy(1,:,:)*vix_un(2,:,:) + vix(2,:,:)*viy_un(1,:,:)
     $     + vix(3,:,:)*vix_un(3,:,:) + viy(3,:,:)*viy_un(3,:,:))
     $     + nu*two*(vex*vex_un + vey*vey_un)))
     $     - cyl_fac*r_faci*(mu*(two*vi(2,:,:)**2 + vi(3,:,:)**2)
     $     + nu*ve**2)
      s_u(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_u3
      s_u(8,4,:,:)=r_fac*(ux(8,:,:)*n_inv 
     $     - mu*n_inv*(4._r8*vix(1,:,:)*ux(1,:,:)
     $     + two*(vix(2,:,:) + viy(1,:,:))*uy(1,:,:)))
      s_u(8,5,:,:)=r_fac*(uy(8,:,:)*n_inv 
     $     - mu*n_inv*(4._r8*viy(2,:,:)*uy(1,:,:)
     $     + two*(vix(2,:,:) + viy(1,:,:))*ux(1,:,:)))
     $     + cyl_fac*4._r8*mu*r_faci*vi(2,:,:)
      s_u(8,6,:,:)=-r_fac*two*n_inv
     $     *(mu*(vix(3,:,:)*ux(1,:,:) + viy(3,:,:)*uy(1,:,:))
     $     + nu*(vex*ux(1,:,:) + vey*uy(1,:,:)))
     $     + cyl_fac*two*r_faci*(mu*vi(3,:,:)+nu*ve)
      s_u(8,7,:,:)=r_fac*((eta_local*two + eta_j*jtot)*jtot*j_u7
     $     + di*nu*two*n_inv*(vex*ux(1,:,:) + vey*uy(1,:,:)))
     $     - cyl_fac*two*r_faci*di*nu*ve
      s_u(8,8,:,:)=r_fac*eta_p*jtot**2

      s_ux(8,1,:,:)= -r_fac*two*(mu*(two*vix(1,:,:)*vi(1,:,:)  
     $     + vi(2,:,:)*(vix(2,:,:) + viy(1,:,:))
     $     + vix(3,:,:)*vi(3,:,:)) + nu*vex*ve)
      s_ux(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_ux3
      s_ux(8,4,:,:)=r_fac*mu*4._r8*vix(1,:,:)
      s_ux(8,5,:,:)=r_fac*mu*two*(vix(2,:,:) + viy(1,:,:))
      s_ux(8,6,:,:)=r_fac*two*(mu*vix(3,:,:) + nu*vex)
      s_ux(8,7,:,:)=-r_fac*di*nu*two*vex
      s_ux(8,8,:,:)=r_fac*vi(1,:,:)

      s_uy(8,1,:,:)= -r_fac*two*(mu*(two*viy(2,:,:)*vi(2,:,:)  
     $     + vi(1,:,:)*(viy(1,:,:) + vix(2,:,:))
     $     + viy(3,:,:)*vi(3,:,:)) + nu*vey*ve)
      s_uy(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_uy3
      s_uy(8,4,:,:)=r_fac*mu*two*(viy(1,:,:) + vix(2,:,:))
      s_uy(8,5,:,:)=r_fac*mu*4._r8*viy(2,:,:)
      s_uy(8,6,:,:)=r_fac*two*(mu*viy(3,:,:) + nu*vey)
      s_uy(8,7,:,:)=-r_fac*di*nu*two*vey
      s_uy(8,8,:,:)=r_fac*vi(2,:,:)
c-----------------------------------------------------------------------
c     del^2(Bz).
c-----------------------------------------------------------------------
      IF(init_type=="GEM9")THEN
         fx_ux(9,3,:,:)=one
         fy_u(9,3,:,:)=cyl_fac*r_faci
         fy_uy(9,3,:,:)=one
         s_u(9,9,:,:)=one
      ENDIF
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
      USE SolarHallMHD_mod
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
      SUBROUTINE physics_grid(ksi,phi,x,y)
      USE SolarHallMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: ksi,phi
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: x,y
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("breakout")
         x = SIN(pi*(phiscale*(phi-half)**3 + (phi-half))/
     $        (0.25_r8*phiscale+one))
     $        *((Rmax-one)*(rscale*ksi**4 + ksi)/(rscale + one) + one)
         y = COS(pi*(phiscale*(phi-half)**3 + (phi-half))/
     $        (0.25_r8*phiscale+one))
     $        *((Rmax-one)*(rscale*ksi**4 + ksi)/(rscale + one) + one)
      CASE("half_breakout")
         x = SIN(half*pi*(phiscale*phi**2 + phi)/
     $        (phiscale+one))
     $        *((Rmax-one)*(rscale*ksi**2 + ksi)/(rscale + one) + one)
         y = COS(half*pi*(phiscale*phi**2 + phi)/
     $        (phiscale+one))
     $        *((Rmax-one)*(rscale*ksi**2 + ksi)/(rscale + one) + one)
      CASE("GEM","GEM9")
         x = half*lx*(ksi**2 + rscale*ksi)/(one + rscale)
         y = half*ly*(phi**2 + rscale*phi)/(one + rscale)
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
      USE SolarHallMHD_mod
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
      USE SolarHallMHD_mod
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
