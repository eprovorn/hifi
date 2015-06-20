c-----------------------------------------------------------------------
c     file SolarLnHMHD.f.
c     contains specifications for Hall MHD model with gravity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. SolarLnHMHD_mod.
c     1. SolarLnHMHD_equil.
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
c     subprogram 0. SolarLnHMHD_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE SolarLnHMHD_mod
      USE transport_mod
      USE extra_mod
      IMPLICIT NONE

      LOGICAL :: source=.FALSE.,cylinder=.FALSE.
      CHARACTER(16) :: init_type=".",eta_case=".",kappa_case=".",
     $     mu_case="."
      INTEGER :: cyl_fac=0
      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,qe=1.602e-19,
     $     me=9.109e-31,mp=1.673e-27,ep0=8.854e-12,mu0=4.e-7*pi,
     $     R0=6.96e8,k_B=1.381e-23,g0=2.74e2,m_eff=1.,chod_const=0.1,
     $     kappa_min=0.,kappa_max=1.e8
      REAL(r8) :: di=0.,nu=0.,eta=0.,etavac=1.,v_chod_norm=1.,
     $     etac_norm=1.,etas_norm=1.,mu=0.,mu_min=0.,kappa_par=0.,
     $     kappa_perp=0.,Dn=0.,Te_frac=0.,n0=1.,T0=1.,b0=1.,
     $     Rmax=0.,gamma_fac=1.,gravity=1.,v_peak=0.,v_period=1.,
     $     v_angle=0.,ke_norm=1.,ki_norm=1.,xe_norm=1.,xi_norm=1.,
     $     mu_fac=0.,alpha=1.,rscale=1.,phiscale=1.,hyper_eta=0.,
     $     j_crit=1.e20

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. SolarLnHMHD_equil.
c     computes equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE SolarLnHMHD_equil(x,y,u,ux,uy,deriv)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: deriv

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: rsq
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0

      rsq = x**2 + y**2
      SELECT CASE(init_type)      
      CASE("half_breakout")
         CALL RANDOM_NUMBER(u)
         u = 1.e-20*u

         u(1,:,:) = half*(one-mu_fac)*LOG(rsq)
         u(2,:,:) = -y/rsq**1.5*(one + alpha*(4.*x**2 - y**2)/rsq**2)
         u(8,:,:) = (T0*mu0*n0*k_B/b0**2)*EXP(u(1,:,:))/SQRT(rsq)
         IF(.NOT. deriv)RETURN
         ux(1,:,:) = (one-mu_fac)*x/rsq
         uy(1,:,:) = (one-mu_fac)*y/rsq
         ux(2,:,:) = x*y/rsq**2.5*(3. 
     $        + alpha*(20.*x**2 - 15.*y**2)/rsq**2)
         uy(2,:,:) = (2.*y**2 - x**2)/rsq**2.5
     $        + alpha*(35.*x**2*y**2 - 4.*rsq**2)/rsq**4.5
         ux(8,:,:) = -(T0*mu0*n0*k_B/b0**2)*mu_fac
     $        *x*rsq**(-half*mu_fac-one)
         uy(8,:,:) = -(T0*mu0*n0*k_B/b0**2)*mu_fac
     $        *y*rsq**(-half*mu_fac-one)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SolarLnHMHD_equil
      END MODULE SolarLnHMHD_mod
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
      USE SolarLnHMHD_mod
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
     $     di,eta,eta_case,etavac,nu,kappa_case,kappa_par,kappa_perp,Dn,
     $     mu_case,mu,mu_min,Te_frac,n0,T0,b0,Rmax,v_peak,v_period,
     $     v_angle,alpha,hyper_eta,j_crit
c-----------------------------------------------------------------------
c     Sample namelist.
c-----------------------------------------------------------------------
c$$$&SolarLnHMHD_list
c$$$
c$$$	cylinder=t
c$$$    rscale=0.
c$$$    phiscale=0.
c$$$
c$$$	init_type="half_breakout"
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
c$$$	di=.1
c$$$	nu=2.e-5
c$$$
c$$$    mu_case="scalar"
c$$$	mu=1.00e-1
c$$$    mu_min=1.e-4
c$$$
c$$$	eta_case="spitzer-chodura"
c$$$	eta=2.e-3
c$$$	etavac=100.
c$$$    j_crit=1.e4
c$$$
c$$$    kappa_case="anisotropic"
c$$$	kappa_par=20.0
c$$$	kappa_perp=4.02e-4
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
      physics_type="SolarLnHMHD"

      nqty=8
      nqty_schur=0
      IF(hyper_eta > 0.)nqty=9

      SELECT CASE(init_type)
      CASE("half_breakout")
         cylinder=.TRUE.
      END SELECT

      IF(Te_frac > 1.)Te_frac=1.
      IF(Te_frac < 0.)Te_frac=0.
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
      USE SolarLnHMHD_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty

      REAL(r8) :: tnorm
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(eta_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(kappa_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(mu_case,16,MPI_CHARACTER,0,comm,ierr)
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
      CALL MPI_Bcast(mu_min,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_par,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_perp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Te_frac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(T0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Rmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(v_peak,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(v_period,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(v_angle,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(hyper_eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(j_crit,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(7)=.TRUE.
      adapt_qty(7)=.FALSE.
      IF(hyper_eta > 0.)THEN
         static(9)=.TRUE.
         adapt_qty(9)=.FALSE.
      ENDIF
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
      tnorm=R0*SQRT(mu0*n0*mp*m_eff)/b0
      gravity=g0*tnorm**2/R0
      v_chod_norm=one/(qe*R0)*SQRT(mp*m_eff/(mu0*n0))
      etac_norm=me/(qe*R0*b0*SQRT(ep0*mu0))
      etas_norm=5.e-5*17.*tnorm*(two*n0*mu0*qe)**1.5/(mu0*R0**2*b0**3)
      IF(Te_frac > 0)THEN
         etas_norm = etas_norm/Te_frac**1.5
      ELSE
         etas_norm=0.
      ENDIF
      ke_norm = 3.56e21*(b0**2/(two*n0*mu0*qe))**2.5
     $     *tnorm/(two*n0*R0**2)
      xe_norm = 3.56e21*(b0**2/(two*n0*mu0*qe))**1.5*b0/n0
      ki_norm = 1.97e-7*(b0**2/(two*n0*mu0*qe))**2.5
     $     *SQRT(mu0/(n0*mp))/(two*R0*b0)
      xi_norm = 1.97e-7*(b0**2/(two*n0*mu0*qe))**1.5
     $     *b0/(n0*mp*SQRT(m_eff))
      gamma_fac=gamma/(gamma-one)

      SELECT CASE(init_type)
      CASE("half_breakout")
         di = zero
         T0 = 4.e6
         v_period = v_period/tnorm
         v_peak = v_peak*(pi**2/4._r8)/v_period
         v_angle = v_angle*pi/180._r8
         mu_fac=g0*R0*mp*m_eff/(k_B*T0)
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
      USE SolarLnHMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     SolarLnHMHD initial conditions.
c-----------------------------------------------------------------------
      CALL SolarLnHMHD_equil(x,y,u,ux,uy,.FALSE.)

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
      USE SolarLnHMHD_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      edge_order=(/1,3,2,4/)
      top%bc_type="robin"
      top%static=.FALSE.
      bottom%bc_type="robin"
      bottom%static=.FALSE.
      left%bc_type="robin"
      left%static=.FALSE.
      right%bc_type="robin"
      right%static=.FALSE.

      SELECT CASE(init_type)
      CASE("half_breakout")
         top%static=.TRUE.

         bottom%static=.TRUE.
         bottom%bc_type(1:3)="zeroflux"
         bottom%bc_type(5)="zeroflux"
         bottom%bc_type(7:8)="zeroflux"

         left%bc_type(1)="normflux"
         left%bc_type(2)="zeroflux"
         left%bc_type(3)="normflux"
         left%static(4:6)=.TRUE.
         left%bc_type(8)="zeroflux"

         right%static(1)=.TRUE.
         right%static(3:6)=.TRUE.
         right%bc_type(7)="natural"
         right%static(8)=.TRUE.

         IF(hyper_eta > 0.)THEN
            bottom%bc_type(9)="zeroflux"
            left%bc_type(9)="zeroflux"
            right%static(9)=.TRUE.
         ENDIF
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
      USE SolarLnHMHD_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci
c----------------------------------------------------------------------
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
      CASE("half_breakout")
         SELECT CASE(lrtb)
         CASE("bottom")
            c(4,:,:)=u(4,:,:)
            c(6,:,:)=u(6,:,:)
         CASE("top")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=u(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=u(7,:,:)
            c(8,:,:)=uy(8,:,:)
            IF(hyper_eta > 0.)THEN
               c(9,:,:)=u(9,:,:)
            ENDIF
         CASE("left")
            WHERE((nhat(1,:,:)*ux(1,:,:) + nhat(2,:,:)*uy(1,:,:)) < 0)
               c(1,:,:) = -1.e-2*r_fac
     $              *(nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:))**2
            END WHERE
            IF(t <= v_period)THEN
               WHERE(y > 0. .AND. ABS(ATAN(x/y)) < v_angle)
                  c(3,:,:) = v_peak*x**2*(v_angle**2 - ATAN(x/y)**2)**2
     $                 *(one - COS(twopi/v_period*t))
     $                 *(one + two*alpha*(one + x**2 - 4._r8*y**2))
               END WHERE
            ENDIF
            c(4,:,:)=(nhat(1,:,:)*nhat(2,:,:)
     $           *(ux(4,:,:) - uy(5,:,:) - cyl_fac*r_faci*u(5,:,:))
     $           + nhat(2,:,:)**2*uy(4,:,:) - nhat(1,:,:)**2*ux(5,:,:))
     $           - (nhat(1,:,:)*nhat(2,:,:)
     $           *(ux(1,:,:)*u(4,:,:) - uy(1,:,:)*u(5,:,:))
     $           + nhat(2,:,:)**2*uy(1,:,:)*u(4,:,:) 
     $           - nhat(1,:,:)**2*ux(1,:,:)*u(5,:,:))
            c(5,:,:)=nhat(1,:,:)*u(4,:,:) + nhat(2,:,:)*u(5,:,:)
            c(6,:,:)=nhat(1,:,:)*(ux(6,:,:) - ux(1,:,:)*u(6,:,:)) 
     $           + nhat(2,:,:)*(uy(6,:,:) - uy(1,:,:)*u(6,:,:)
     $           + cyl_fac*r_faci*u(6,:,:))
ccc            c(7,:,:)=nhat(1,:,:)*ux(7,:,:)
ccc     $           + nhat(2,:,:)*(uy(7,:,:) + cyl_fac*r_faci*u(7,:,:))
ccc            c(8,:,:)=u(8,:,:)-(T0*mu0*n0*k_B/b0**2)*EXP(u(1,:,:))
         CASE("right")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:)
     $           - (one - mu_fac)/Rmax
            c(3,:,:)=nhat(2,:,:)*(cyl_fac*r_faci*u(3,:,:) + uy(3,:,:))
     $           + nhat(1,:,:)*ux(3,:,:)
            c(4,:,:)=nhat(1,:,:)*(ux(4,:,:) - ux(1,:,:)*u(4,:,:)) 
     $           + nhat(2,:,:)*(uy(4,:,:) - uy(1,:,:)*u(4,:,:))
            c(5,:,:)=nhat(1,:,:)*(ux(5,:,:) - ux(1,:,:)*u(5,:,:))
     $           + nhat(2,:,:)*(uy(5,:,:) - uy(1,:,:)*u(5,:,:)
     $           + cyl_fac*r_faci*u(5,:,:))
            c(6,:,:)=nhat(1,:,:)*(ux(6,:,:) - ux(1,:,:)*u(6,:,:))
     $           + nhat(2,:,:)*(uy(6,:,:) - uy(1,:,:)*u(6,:,:)
     $           + cyl_fac*r_faci*u(6,:,:))
            c(8,:,:)=u(8,:,:)-(T0*mu0*n0*k_B/b0**2)*EXP(u(1,:,:))/Rmax
            IF(hyper_eta > 0.)THEN
               c(9,:,:)=u(9,:,:)
            ENDIF
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
      USE SolarLnHMHD_mod
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
      CASE("half_breakout")
         SELECT CASE(lrtb)
         CASE("bottom")
            c_u(4,4,:,:) = one
            c_u(6,6,:,:) = one
         CASE("top")
            c_uy(1,1,:,:) = one
            c_u(2,2,:,:) = one
            c_u(3,3,:,:) = one
            c_uy(4,4,:,:) = one
            c_u(5,5,:,:) = one
            c_u(6,6,:,:) = one
            c_u(7,7,:,:) = one
            c_uy(8,8,:,:) = one
            IF(hyper_eta > 0.)THEN
               c_u(9,9,:,:) = one
            ENDIF
         CASE("left")
            WHERE((nhat(1,:,:)*ux(1,:,:) + nhat(2,:,:)*uy(1,:,:)) < 0)
               c_ux(1,1,:,:) = -1.e-2*r_fac*two*nhat(1,:,:)
     $              *(nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:))
               c_uy(1,1,:,:) = -1.e-2*r_fac*two*nhat(2,:,:)
     $              *(nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:))
            END WHERE

            c_u(4,4,:,:) = -nhat(1,:,:)*nhat(2,:,:)*ux(1,:,:)
     $           - nhat(2,:,:)**2*uy(1,:,:)
            c_u(4,5,:,:) = nhat(1,:,:)*nhat(2,:,:)
     $           *(uy(1,:,:) - cyl_fac*r_faci)
     $           + nhat(1,:,:)**2*ux(1,:,:)
            c_ux(4,1,:,:) = -nhat(1,:,:)*nhat(2,:,:)*u(4,:,:)
     $           + nhat(1,:,:)**2*u(5,:,:)
            c_ux(4,4,:,:) = nhat(1,:,:)*nhat(2,:,:)
            c_ux(4,5,:,:) = -nhat(1,:,:)**2
            c_uy(4,1,:,:) = nhat(1,:,:)*nhat(2,:,:)*u(5,:,:)
     $           - nhat(2,:,:)**2*u(4,:,:)
            c_uy(4,4,:,:) = nhat(2,:,:)**2
            c_uy(4,5,:,:) = -nhat(1,:,:)*nhat(2,:,:)

            c_u(5,4,:,:) = nhat(1,:,:)
            c_u(5,5,:,:) = nhat(2,:,:)

            c_u(6,6,:,:) =-nhat(1,:,:)*ux(1,:,:)
     $           - nhat(2,:,:)*(uy(1,:,:) - cyl_fac*r_faci)
            c_ux(6,1,:,:) = -nhat(1,:,:)*u(6,:,:)
            c_uy(6,1,:,:) = -nhat(2,:,:)*u(6,:,:)
            c_ux(6,6,:,:) = nhat(1,:,:)
            c_uy(6,6,:,:) = nhat(2,:,:)

ccc            c_u(7,7,:,:) = nhat(2,:,:)*r_faci*cyl_fac
ccc            c_ux(7,7,:,:) = nhat(1,:,:)
ccc            c_uy(7,7,:,:) = nhat(2,:,:)

ccc            c_u(8,1,:,:) = -(T0*mu0*n0*k_B/b0**2)*EXP(u(1,:,:))
ccc            c_u(8,8,:,:) = one
         CASE("right")
            c_ux(1,1,:,:) = nhat(1,:,:)
            c_uy(1,1,:,:) = nhat(2,:,:)

            c_u(3,3,:,:) = nhat(2,:,:)*cyl_fac*r_faci
            c_ux(3,3,:,:) = nhat(1,:,:)
            c_uy(3,3,:,:) = nhat(2,:,:)

            c_u(4,4,:,:) =-nhat(1,:,:)*ux(1,:,:)
     $           - nhat(2,:,:)*uy(1,:,:)
            c_ux(4,1,:,:) = -nhat(1,:,:)*u(4,:,:)
            c_uy(4,1,:,:) = -nhat(2,:,:)*u(4,:,:)
            c_ux(4,4,:,:) = nhat(1,:,:)
            c_uy(4,4,:,:) = nhat(2,:,:)

            c_u(5,5,:,:) =-nhat(1,:,:)*ux(1,:,:)
     $           - nhat(2,:,:)*(uy(1,:,:) - cyl_fac*r_faci)
            c_ux(5,1,:,:) = -nhat(1,:,:)*u(5,:,:)
            c_uy(5,1,:,:) = -nhat(2,:,:)*u(5,:,:)
            c_ux(5,5,:,:) = nhat(1,:,:)
            c_uy(5,5,:,:) = nhat(2,:,:)

            c_u(6,6,:,:) =-nhat(1,:,:)*ux(1,:,:)
     $           - nhat(2,:,:)*(uy(1,:,:) - cyl_fac*r_faci)
            c_ux(6,1,:,:) = -nhat(1,:,:)*u(6,:,:)
            c_uy(6,1,:,:) = -nhat(2,:,:)*u(6,:,:)
            c_ux(6,6,:,:) = nhat(1,:,:)
            c_uy(6,6,:,:) = nhat(2,:,:)

            c_u(8,1,:,:) = -(T0*mu0*n0*k_B/b0**2)*EXP(u(1,:,:))/Rmax
            c_u(8,8,:,:) = one
            IF(hyper_eta > 0.)THEN
               c_u(9,9,:,:) = one
            ENDIF
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
      USE SolarLnHMHD_mod
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
      CASE("half_breakout")
         SELECT CASE(lrtb)
         CASE("left")
            mass(7,2,:,:)=one
         CASE("right")
            mass(2,2,:,:)=cyl_fac*y**(-cyl_fac)*nhat(2,:,:)
            mass_x(2,2,:,:)=nhat(1,:,:)
            mass_y(2,2,:,:)=nhat(2,:,:)
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
      USE SolarLnHMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Tix,Tiy,kperp,kfac,
     $     ve,vex,vey,r_fac,r_faci,j1,j2,j3,jtot,eta_local,mu_local,
     $     b1,b2,Bsq,rho,rho_inv,rhox_inv,rhoy_inv,gx_vec,gy_vec,
     $     r_sphr
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
      rho=EXP(u(1,:,:))
      rho_inv=EXP(-u(1,:,:))
      rhox_inv=-ux(1,:,:)*rho_inv
      rhoy_inv=-uy(1,:,:)*rho_inv
c-----------------------------------------------------------------------
c     temperature gradients.
c-----------------------------------------------------------------------
      Tix = rho_inv*(ux(8,:,:) - u(8,:,:)*ux(1,:,:))
      Tiy = rho_inv*(uy(8,:,:) - u(8,:,:)*uy(1,:,:))
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
         CALL transport_seteta("spitzer-chodura",r_sphr,rho,
     $        u(8,:,:),jtot,0._r8,etas_norm,0._r8,0._r8,0._r8,
     $        eta,etavac,eta_local)
      CASE("spitzer-chodura")
         CALL transport_seteta(eta_case,r_sphr,rho,u(8,:,:),jtot,
     $        chod_const,etas_norm,etac_norm,v_chod_norm,0._r8,
     $        eta,etavac,eta_local)
      CASE("enhanced")
         eta_local = eta
         WHERE(jtot > .95*j_crit .AND. jtot < 1.05*j_crit)
            eta_local = eta + (etavac-eta)*half
     $           *(one - COS(10.*pi*(jtot/j_crit-.95)))
         ELSEWHERE(jtot >= 1.05*j_crit)
            eta_local = etavac
         END WHERE
      CASE DEFAULT
         eta_local = eta
      END SELECT
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT(b1,b2,u(3,:,:),Tix,Tiy,BdotT)
         CALL transport_kbrag(rho,u(8,:,:),Te_frac,Bsq,ke_norm,ki_norm,
     $        xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac)
         kperp = kperp + kappa_par
      CASE("anisotropic")
         CALL transport_BdotT(b1,b2,u(3,:,:),Tix,Tiy,BdotT)
         CALL transport_setkaniso(kappa_par,kappa_perp,Bsq,kperp,kfac)
         kperp=kperp*rho
         kfac=kfac*rho
      CASE("scalar")
         BdotT=0
         kperp=kappa_par
         kfac=0
      CASE DEFAULT
         BdotT=0
         kperp=kappa_par*rho
         kfac=0
      END SELECT
c-----------------------------------------------------------------------
c     viscosity.
c-----------------------------------------------------------------------
      SELECT CASE(mu_case)
      CASE("scalar")
         mu_local = mu
      CASE DEFAULT
         mu_local = mu*rho + mu_min
      END SELECT
c-----------------------------------------------------------------------
c     velocities and their gradients.
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:)=u(i+3,:,:)*rho_inv
         vix(i,:,:)=ux(i+3,:,:)*rho_inv - vi(i,:,:)*ux(1,:,:)
         viy(i,:,:)=uy(i+3,:,:)*rho_inv - vi(i,:,:)*uy(1,:,:)
      ENDDO

      ve=(u(6,:,:)-di*u(7,:,:))*rho_inv
      vex=(ux(6,:,:)-di*ux(7,:,:))*rho_inv - ve*ux(1,:,:)
      vey=(uy(6,:,:)-di*uy(7,:,:))*rho_inv - ve*uy(1,:,:)
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
      fx(1,:,:) = r_fac*vi(1,:,:)
      fy(1,:,:) = r_fac*vi(2,:,:)
      s(1,:,:) = -r_fac*(vi(1,:,:)*ux(1,:,:) + vi(2,:,:)*uy(1,:,:))
c-----------------------------------------------------------------------
c     poloidal magnetic flux equation.
c-----------------------------------------------------------------------
      fx(2,:,:) = -r_fac*(di*nu*vex - hyper_eta*ux(7,:,:))
      fy(2,:,:) = -r_fac*(di*nu*vey - hyper_eta*uy(7,:,:))
      
      s(2,:,:) = r_fac*(vi(2,:,:)*b1 - vi(1,:,:)*b2
     $     - di*(j2*b1 - j1*b2)*rho_inv + eta_local*j3
     $     + di*nu*(vex*ux(1,:,:) + vey*uy(1,:,:)))
     $     - cyl_fac*r_faci*(di*nu*ve - hyper_eta*u(7,:,:))
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
      fx(3,:,:) = u(3,:,:)*vi(1,:,:) - ve*b1
     $     - di*Te_frac*uy(8,:,:)*rho_inv + eta_local*j2

      fy(3,:,:) = u(3,:,:)*vi(2,:,:) - ve*b2
     $     + di*Te_frac*ux(8,:,:)*rho_inv - eta_local*j1

      s(3,:,:) = di*u(3,:,:)*(j2*rhoy_inv + j1*rhox_inv) 
     $     + cyl_fac*two*r_faci*di*u(3,:,:)*ux(3,:,:)*rho_inv

      IF(hyper_eta > 0.)THEN
         fx(3,:,:)=fx(3,:,:) + hyper_eta*ux(9,:,:)
         fy(3,:,:)=fy(3,:,:) 
     $        + hyper_eta*(cyl_fac*r_faci*u(9,:,:) + uy(9,:,:))
      ENDIF
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(4,:,:) = r_fac*(u(4,:,:)*vi(1,:,:) - two*mu_local*vix(1,:,:) 
     $     + half*u(3,:,:)**2 + u(8,:,:))
         
      fy(4,:,:) = r_fac*(u(4,:,:)*vi(2,:,:) 
     $     - mu_local*(viy(1,:,:) + vix(2,:,:)))

      s(4,:,:) = -r_fac*(j3*b2 + gx_vec*rho)
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(5,:,:) = r_fac*(u(5,:,:)*vi(1,:,:) 
     $     - mu_local*(vix(2,:,:) + viy(1,:,:)))
      
      fy(5,:,:) = r_fac*(u(5,:,:)*vi(2,:,:) + half*u(3,:,:)**2
     $     - two*mu_local*viy(2,:,:))
      
      s(5,:,:) = r_fac*(j3*b1 - uy(8,:,:) - gy_vec*rho) 
     $     + cyl_fac*(u(6,:,:)*vi(3,:,:) - half*u(3,:,:)**2 
     $     - two*mu_local*vi(2,:,:)*r_faci)
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(6,:,:) = r_fac*(u(6,:,:)*vi(1,:,:) 
     $     - rho*nu*vex - mu_local*vix(3,:,:))
      
      fy(6,:,:) = r_fac*(u(6,:,:)*vi(2,:,:) 
     $     - rho*nu*vey - mu_local*viy(3,:,:))

      s(6,:,:) = r_fac*(j1*b2 - j2*b1) - cyl_fac*(u(6,:,:)*vi(2,:,:) 
     $     + r_faci*(mu_local*vi(3,:,:) + nu*rho*ve))
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
     $     + mu_local*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + (viy(1,:,:) + vix(2,:,:))**2
     $     + vix(3,:,:)**2 + viy(3,:,:)**2) 
     $     + nu*rho*(vex**2 + vey**2))
     $     + cyl_fac*r_faci*(mu_local*(two*vi(2,:,:)**2 + vi(3,:,:)**2) 
     $     + nu*rho*ve**2)
c-----------------------------------------------------------------------
c     del^2(out-of-plane B).
c-----------------------------------------------------------------------
      IF(hyper_eta > 0.)THEN
         fx(9,:,:)=ux(3,:,:)
         fy(9,:,:)=cyl_fac*r_faci*u(3,:,:) + uy(3,:,:)
         s(9,:,:)=u(9,:,:)
      ENDIF
c-----------------------------------------------------------------------
c     initial equilibrium source term.
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("half_breakout")
            CALL SolarLnHMHD_equil(x,y,u0,u0x,u0y,.TRUE.)
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
      USE SolarLnHMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i
      REAL(r8) :: mu_diff
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci,Tix,Tiy,
     $     Ti_un,kperp,kfac,ve,vex,vey,j1,j2,j3,jtot,eta_local,b1,b2,
     $     Bsq,eta_rho,eta_p,eta_j,j_u3,j_ux3,j_uy3,j_u7,rho,rho_inv,
     $     rhox_inv,rhoy_inv,gx_vec,gy_vec,r_sphr,kpar_un,kpar_p,
     $     kperp_un,kperp_p,kperp_bsq,mu_local
      REAL(r8), DIMENSION(3,SIZE(x,1),SIZE(x,2)) :: vi,vix,viy,
     $     BdotT,BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty
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
      rho=EXP(u(1,:,:))
      rho_inv=EXP(-u(1,:,:))
      rhox_inv=-ux(1,:,:)*rho_inv
      rhoy_inv=-uy(1,:,:)*rho_inv
c-----------------------------------------------------------------------
c     temperature gradients and derivatives.
c-----------------------------------------------------------------------
      Tix = ux(8,:,:)*rho_inv + u(8,:,:)*rhox_inv
      Tiy = uy(8,:,:)*rho_inv + u(8,:,:)*rhoy_inv
      Ti_un=-u(8,:,:)*rho_inv
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
         CALL transport_seteta_u("spitzer-chodura",r_sphr,rho,u(8,:,:),
     $        jtot,0._r8,etas_norm,0._r8,0._r8,0._r8,eta,etavac,
     $        eta_local,eta_rho,eta_p,eta_j)
         eta_rho = eta_rho*rho
      CASE("spitzer-chodura")
         CALL transport_seteta_u(eta_case,r_sphr,rho,u(8,:,:),jtot,
     $        chod_const,etas_norm,etac_norm,v_chod_norm,0._r8,
     $        eta,etavac,eta_local,eta_rho,eta_p,eta_j)
         eta_rho = eta_rho*rho
         WHERE(eta_local == etavac)
            eta_rho=0.
            eta_p=0.
            eta_j=0.
         END WHERE
      CASE("enhanced")
         eta_local = eta
         WHERE(jtot > .95*j_crit .AND. jtot < 1.05*j_crit)
            eta_local = eta + (etavac-eta)*half
     $           *(one - COS(10.*pi*(jtot/j_crit-.95)))
         ELSEWHERE(jtot >= 1.05*j_crit)
            eta_local = etavac
         END WHERE
         eta_rho=0.
         eta_p=0.
         eta_j=0.
         WHERE(jtot > .95*j_crit .AND. jtot < 1.05*j_crit)
            eta_j = (etavac-eta)*half*10.*pi/j_crit
     $           *SIN(10.*pi*(jtot/j_crit-.95))
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
         CALL transport_kbrag_u(rho,u(8,:,:),Te_frac,Bsq,ke_norm,
     $        ki_norm,xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac,
     $        kpar_un,kpar_p,kperp_un,kperp_p,kperp_bsq)
         kperp = kperp + kappa_par
         kpar_un = kpar_un*rho
         kperp_un = kperp_un*rho
      CASE("anisotropic")
         CALL transport_BdotT_u(b1,b2,u(3,:,:),Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty)
         CALL transport_setkaniso_u(kappa_par,kappa_perp,Bsq,kperp,kfac,
     $        kperp_bsq)
         kperp_un=kperp*rho
         kperp_p=0
         kperp=kperp*rho
         kperp_bsq=kperp_bsq*rho
         kfac=kfac*rho
         kpar_un=kappa_par*rho
         kpar_p=0
      CASE("scalar")
         kperp=kappa_par
         BdotT=0
         BdotT_b1=0
         BdotT_b2=0
         BdotT_b3=0
         BdotT_Tx=0
         BdotT_Ty=0
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
         kperp=kappa_par*rho
         kfac=0
         kpar_un=0
         kpar_p=0
         kperp_un=kappa_par*rho
         kperp_p=0
         kperp_bsq=0
      END SELECT
c-----------------------------------------------------------------------
c     viscosity.
c-----------------------------------------------------------------------
      SELECT CASE(mu_case)
      CASE("scalar")
         mu_local = mu
         mu_diff = -mu
      CASE DEFAULT
         mu_local = mu*rho + mu_min
         mu_diff = -mu_min
      END SELECT
c-----------------------------------------------------------------------
c     velocities and their derivatives.
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:) = u(i+3,:,:)*rho_inv
         vix(i,:,:) = ux(i+3,:,:)*rho_inv - vi(i,:,:)*ux(1,:,:)
         viy(i,:,:) = uy(i+3,:,:)*rho_inv - vi(i,:,:)*uy(1,:,:)
      ENDDO

      ve = (u(6,:,:)-di*u(7,:,:))*rho_inv
      vex = (ux(6,:,:)-di*ux(7,:,:))*rho_inv - ve*ux(1,:,:)
      vey = (uy(6,:,:)-di*uy(7,:,:))*rho_inv - ve*uy(1,:,:)
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
      fx_u(1,1,:,:) = -r_fac*vi(1,:,:)
      fx_u(1,4,:,:) = r_fac*rho_inv
      fy_u(1,1,:,:) = -r_fac*vi(2,:,:)
      fy_u(1,5,:,:) = r_fac*rho_inv
      s_u(1,1,:,:) = r_fac*(vi(1,:,:)*ux(1,:,:) + vi(2,:,:)*uy(1,:,:))
      s_u(1,4,:,:) = r_fac*rhox_inv
      s_u(1,5,:,:) = r_fac*rhoy_inv
      s_ux(1,1,:,:) = -r_fac*vi(1,:,:) 
      s_uy(1,1,:,:) = -r_fac*vi(2,:,:)
c-----------------------------------------------------------------------
c     poloidal magnetic flux equation.
c-----------------------------------------------------------------------
      fx_u(2,1,:,:) = r_fac*di*nu*vex
      fx_u(2,6,:,:) =-r_fac*di*nu*rhox_inv
      fx_u(2,7,:,:) = r_fac*di**2*nu*rhox_inv
      fx_ux(2,1,:,:) = r_fac*di*nu*ve
      fx_ux(2,6,:,:) =-r_fac*di*nu*rho_inv
      fx_ux(2,7,:,:) = r_fac*(di**2*nu*rho_inv + hyper_eta)
      
      fy_u(2,1,:,:) = r_fac*di*nu*vey
      fy_u(2,6,:,:) =-r_fac*di*nu*rhoy_inv
      fy_u(2,7,:,:) = r_fac*di**2*nu*rhoy_inv
      fy_uy(2,1,:,:) = r_fac*di*nu*ve
      fy_uy(2,6,:,:) =-r_fac*di*nu*rho_inv
      fy_uy(2,7,:,:) = r_fac*(di**2*nu*rho_inv + hyper_eta)
      
      s_u(2,1,:,:)=r_fac*(-vi(2,:,:)*b1 + vi(1,:,:)*b2
     $     + di*(j2*b1 - j1*b2)*rho_inv + eta_rho*j3
     $     - di*nu*(vex*ux(1,:,:) + vey*uy(1,:,:)))
     $     + cyl_fac*r_faci*di*nu*ve
      s_u(2,2,:,:)=-cyl_fac*(vi(2,:,:)-di*j2*rho_inv)
      s_u(2,3,:,:)=di*cyl_fac*b2*rho_inv + r_fac*eta_j*j_u3*j3
      s_u(2,4,:,:)=-r_fac*b2*rho_inv
      s_u(2,5,:,:)=r_fac*b1*rho_inv
      s_u(2,6,:,:)=r_fac*di*nu*(rhox_inv*ux(1,:,:) + rhoy_inv*uy(1,:,:))
     $     - cyl_fac*r_faci*di*nu*rho_inv
      s_u(2,7,:,:)=r_fac*(-di**2*nu*(rhox_inv*ux(1,:,:) 
     $     + rhoy_inv*uy(1,:,:)) + eta_j*j_u7*j3 + eta_local)
     $     + cyl_fac*r_faci*(di**2*nu*rho_inv + hyper_eta)
      s_u(2,8,:,:)=r_fac*eta_p*j3
      s_ux(2,1,:,:)=r_fac*di*nu*(vex - ve*ux(1,:,:))
      s_ux(2,2,:,:)=r_fac*(di*j1*rho_inv - vi(1,:,:))
      s_ux(2,3,:,:)=r_fac*(di*b1*rho_inv + eta_j*j_ux3*j3)
      s_ux(2,6,:,:)=-r_fac*di*nu*rhox_inv
      s_ux(2,7,:,:)=r_fac*di**2*nu*rhox_inv
      s_uy(2,1,:,:)=r_fac*di*nu*(vey - ve*uy(1,:,:))
      s_uy(2,2,:,:)=-r_fac*(vi(2,:,:)-di*j2*rho_inv)
      s_uy(2,3,:,:)=r_fac*(di*b2*rho_inv + eta_j*j_uy3*j3)
      s_uy(2,6,:,:)=-r_fac*di*nu*rhoy_inv
      s_uy(2,7,:,:)=r_fac*di**2*nu*rhoy_inv
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=-u(3,:,:)*vi(1,:,:) + ve*b1 
     $     + di*Te_frac*uy(8,:,:)*rho_inv +  eta_rho*j2
      fx_u(3,2,:,:)=cyl_fac*r_faci*ve 
      fx_u(3,3,:,:)=vi(1,:,:) + eta_j*j_u3*j2
      fx_u(3,4,:,:)=u(3,:,:)*rho_inv
      fx_u(3,6,:,:)=-rho_inv*b1
      fx_u(3,7,:,:)=di*rho_inv*b1 + eta_j*j_u7*j2
      fx_u(3,8,:,:)= eta_p*j2
      fx_ux(3,3,:,:)=-eta_local + eta_j*j_ux3*j2
      fx_uy(3,2,:,:)=ve
      fx_uy(3,3,:,:)=eta_j*j_uy3*j2
      fx_uy(3,8,:,:)=-di*Te_frac*rho_inv
         
      fy_u(3,1,:,:)=-u(3,:,:)*vi(2,:,:) + ve*b2 
     $     - di*Te_frac*ux(8,:,:)*rho_inv - eta_rho*j1
      fy_u(3,3,:,:)=vi(2,:,:) - eta_local*cyl_fac*r_faci - eta_j*j_u3*j1
      fy_u(3,5,:,:)=u(3,:,:)*rho_inv
      fy_u(3,6,:,:)=-rho_inv*b2
      fy_u(3,7,:,:)=di*rho_inv*b2 - eta_j*j_u7*j1
      fy_u(3,8,:,:)= -eta_p*j1
      fy_ux(3,2,:,:)=-ve
      fy_ux(3,3,:,:)=-eta_j*j_ux3*j1
      fy_ux(3,8,:,:)=di*Te_frac*rho_inv
      fy_uy(3,3,:,:)=-eta_local - eta_j*j_uy3*j1
         
      s_u(3,1,:,:)=-di*u(3,:,:)*(j2*rhoy_inv + j1*rhox_inv)
     $     - cyl_fac*two*r_faci*di*u(3,:,:)*ux(3,:,:)*rho_inv
      s_u(3,3,:,:)=di*(j2*rhoy_inv + j1*rhox_inv)
     $     + cyl_fac*r_faci*di*(u(3,:,:)*rhox_inv 
     $     + two*ux(3,:,:)*rho_inv)
      s_ux(3,1,:,:)=-di*u(3,:,:)*j1*rho_inv
      s_ux(3,3,:,:)=-di*u(3,:,:)*rhoy_inv
     $     + cyl_fac*two*r_faci*di*u(3,:,:)*rho_inv
      s_uy(3,1,:,:)=-di*u(3,:,:)*j2*rho_inv
      s_uy(3,3,:,:)=di*u(3,:,:)*rhox_inv

      IF(hyper_eta > 0.)THEN
         fx_ux(3,9,:,:)=hyper_eta
         fy_u(3,9,:,:)=hyper_eta*cyl_fac*r_faci
         fy_uy(3,9,:,:)=hyper_eta
      ENDIF
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:)=-r_fac*(u(4,:,:)*vi(1,:,:)
     $     + two*mu_diff*vix(1,:,:)) 
      fx_u(4,3,:,:)=r_fac*u(3,:,:)
      fx_u(4,4,:,:)=r_fac*two*(vi(1,:,:) - mu_local*rhox_inv)
      fx_u(4,8,:,:)=r_fac
      fx_ux(4,1,:,:)= r_fac*two*mu_local*vi(1,:,:)
      fx_ux(4,4,:,:)=-r_fac*two*mu_local*rho_inv

      fy_u(4,1,:,:)=-r_fac*(u(4,:,:)*vi(2,:,:)
     $     + mu_diff*(viy(1,:,:) + vix(2,:,:))) 
      fy_u(4,4,:,:)=r_fac*(vi(2,:,:) - mu_local*rhoy_inv)
      fy_u(4,5,:,:)=r_fac*(vi(1,:,:) - mu_local*rhox_inv)
      fy_ux(4,1,:,:)= r_fac*mu_local*vi(2,:,:)
      fy_ux(4,5,:,:)=-r_fac*mu_local*rho_inv
      fy_uy(4,1,:,:)= r_fac*mu_local*vi(1,:,:)
      fy_uy(4,4,:,:)=-r_fac*mu_local*rho_inv

      s_u(4,1,:,:)=-r_fac*gx_vec*rho
      s_u(4,7,:,:)=-r_fac*b2
      s_ux(4,2,:,:)=-r_fac*j3
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------      
      fx_u(5,1,:,:)=-r_fac*(u(5,:,:)*vi(1,:,:)
     $     + mu_diff*(vix(2,:,:) + viy(1,:,:))) 
      fx_u(5,4,:,:)=r_fac*(vi(2,:,:) - mu_local*rhoy_inv)
      fx_u(5,5,:,:)=r_fac*(vi(1,:,:) - mu_local*rhox_inv)
      fx_ux(5,1,:,:)= r_fac*mu_local*vi(2,:,:)
      fx_ux(5,5,:,:)=-r_fac*mu_local*rho_inv
      fx_uy(5,1,:,:)= r_fac*mu_local*vi(1,:,:)
      fx_uy(5,4,:,:)=-r_fac*mu_local*rho_inv
      
      fy_u(5,1,:,:)=-r_fac*(u(5,:,:)*vi(2,:,:)
     $     + two*mu_diff*viy(2,:,:)) 
      fy_u(5,3,:,:)= r_fac*u(3,:,:)
      fy_u(5,5,:,:)= r_fac*two*(vi(2,:,:) - mu_local*rhoy_inv)
      fy_uy(5,1,:,:)= r_fac*two*mu_local*vi(2,:,:)
      fy_uy(5,5,:,:)=-r_fac*two*mu_local*rho_inv

      s_u(5,1,:,:)=-r_fac*gy_vec*rho - cyl_fac*(u(6,:,:)*vi(3,:,:) 
     $     + two*r_faci*mu_diff*vi(2,:,:))
      s_u(5,2,:,:)=-cyl_fac*j3
      s_u(5,3,:,:)=-cyl_fac*u(3,:,:)
      s_u(5,5,:,:)=-cyl_fac*r_faci*two*mu_local*rho_inv
      s_u(5,6,:,:)=cyl_fac*two*vi(3,:,:)
      s_u(5,7,:,:)=r_fac*b1
      s_uy(5,2,:,:)=-r_fac*j3
      s_uy(5,8,:,:)=-r_fac
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------      
      fx_u(6,1,:,:)=-r_fac*(u(6,:,:)*vi(1,:,:) 
     $     + mu_diff*vix(3,:,:)) 
      fx_u(6,4,:,:)=r_fac*u(6,:,:)*rho_inv
      fx_u(6,6,:,:)=r_fac*(vi(1,:,:) - mu_local*rhox_inv + nu*ux(1,:,:))
      fx_u(6,7,:,:)= -r_fac*di*nu*ux(1,:,:)
      fx_ux(6,1,:,:)= r_fac*(mu_local*vi(3,:,:) + nu*ve*rho)
      fx_ux(6,6,:,:)=-r_fac*(mu_local*rho_inv + nu)
      fx_ux(6,7,:,:)= r_fac*di*nu
      
      fy_u(6,1,:,:)=-r_fac*(u(6,:,:)*vi(2,:,:)
     $     + mu_diff*viy(3,:,:)) 
      fy_u(6,5,:,:)=r_fac*u(6,:,:)*rho_inv
      fy_u(6,6,:,:)=r_fac*(vi(2,:,:) - mu_local*rhoy_inv + nu*uy(1,:,:))
      fy_u(6,7,:,:)= -r_fac*di*nu*uy(1,:,:)
      fy_uy(6,1,:,:)= r_fac*(mu_local*vi(3,:,:) + nu*ve*rho)
      fy_uy(6,6,:,:)=-r_fac*(mu_local*rho_inv + nu)
      fy_uy(6,7,:,:)= r_fac*di*nu
      
      s_u(6,1,:,:)=cyl_fac*(u(6,:,:)*vi(2,:,:)
     $     - r_faci*mu_diff*vi(3,:,:))
      s_u(6,2,:,:)=cyl_fac*j2
      s_u(6,3,:,:)=cyl_fac*b2
      s_u(6,5,:,:)=-cyl_fac*vi(3,:,:)
      s_u(6,6,:,:)=-cyl_fac*(vi(2,:,:) + r_faci*(mu_local*rho_inv + nu))
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
      fx_u(8,1,:,:)=r_fac*(-gamma_fac*u(8,:,:)*vi(1,:,:)
     $     + (kfac*BdotT_Tx(1,:,:) + kperp)*Tix
     $     + kfac*BdotT_Ty(1,:,:)*Tiy
     $     - (kpar_un - kperp_un)*BdotT(1,:,:) - kperp_un*Tix)
      fx_u(8,2,:,:)=-cyl_fac*(-kfac*BdotT_b1(1,:,:)
     $     + two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_u(8,3,:,:)=r_fac*(-kfac*BdotT_b3(1,:,:)
     $     + two*u(3,:,:)*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_u(8,4,:,:)=r_fac*gamma_fac*u(8,:,:)*rho_inv
      fx_u(8,8,:,:)=r_fac*(gamma_fac*vi(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*rhox_inv
     $     - kfac*BdotT_Ty(1,:,:)*rhoy_inv
     $     - (kpar_p - kperp_p)*BdotT(1,:,:) - kperp_p*Tix)

      fx_ux(8,1,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*Ti_un
      fx_ux(8,2,:,:)=r_fac*(-kfac*BdotT_b2(1,:,:)
     $     + two*b2*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_ux(8,8,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*rho_inv

      fx_uy(8,1,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*Ti_un
      fx_uy(8,2,:,:)=-r_fac*(-kfac*BdotT_b1(1,:,:)
     $     + two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_uy(8,8,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*rho_inv

      fy_u(8,1,:,:)=r_fac*(-gamma_fac*u(8,:,:)*vi(2,:,:)
     $     + (kfac*BdotT_Ty(2,:,:) + kperp)*Tiy
     $     + kfac*BdotT_Tx(2,:,:)*Tix
     $     - (kpar_un - kperp_un)*BdotT(2,:,:) - kperp_un*Tiy)
      fy_u(8,2,:,:)=-cyl_fac*(-kfac*BdotT_b1(2,:,:)
     $     + two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(8,3,:,:)=r_fac*(-kfac*BdotT_b3(2,:,:)
     $     + two*u(3,:,:)*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(8,5,:,:)=r_fac*gamma_fac*u(8,:,:)*rho_inv
      fy_u(8,8,:,:)=r_fac*(gamma_fac*vi(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*rhoy_inv
     $     - kfac*BdotT_Tx(2,:,:)*rhox_inv
     $     - (kpar_p - kperp_p)*BdotT(2,:,:) - kperp_p*Tiy)

      fy_ux(8,1,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*Ti_un
      fy_ux(8,2,:,:)=r_fac*(-kfac*BdotT_b2(2,:,:)
     $     + two*b2*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_ux(8,8,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*rho_inv

      fy_uy(8,1,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*Ti_un
      fy_uy(8,2,:,:)=-r_fac*(-kfac*BdotT_b1(2,:,:)
     $     + two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_uy(8,8,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*rho_inv

      s_u(8,1,:,:)=r_fac*(-vi(1,:,:)*ux(8,:,:)
     $     - vi(2,:,:)*uy(8,:,:) + eta_rho*jtot**2
     $     + (mu_diff - mu_local)
     $     *(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:)
     $     + vix(3,:,:)**2 + viy(3,:,:)**2) 
     $     - nu*rho*(vex**2 + vey**2))
     $     + cyl_fac*r_faci
     $     *((mu_diff - mu_local)*(two*vi(2,:,:)**2 + vi(3,:,:)**2) 
     $     - nu*rho*ve**2)
      s_u(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_u3
      s_u(8,4,:,:)=r_fac*(ux(8,:,:)*rho_inv 
     $     + mu_local*(4._r8*vix(1,:,:)*rhox_inv
     $     + two*(vix(2,:,:) + viy(1,:,:))*rhoy_inv))
      s_u(8,5,:,:)=r_fac*(uy(8,:,:)*rho_inv 
     $     + mu_local*(4._r8*viy(2,:,:)*rhoy_inv
     $     + two*(vix(2,:,:) + viy(1,:,:))*rhox_inv))
     $     + cyl_fac*r_faci*4._r8*mu_local*vi(2,:,:)*rho_inv
      s_u(8,6,:,:)=r_fac*two
     $     *(mu_local*(vix(3,:,:)*rhox_inv + viy(3,:,:)*rhoy_inv)
     $     - nu*(vex*ux(1,:,:) + vey*uy(1,:,:)))
     $     + cyl_fac*r_faci*two*(mu_local*vi(3,:,:)*rho_inv + nu*ve)
      s_u(8,7,:,:)=r_fac*((eta_local*two + eta_j*jtot)*jtot*j_u7
     $     + di*nu*two*(vex*ux(1,:,:) + vey*uy(1,:,:)))
     $     - cyl_fac*two*r_faci*di*nu*ve
      s_u(8,8,:,:)=r_fac*eta_p*jtot**2

      s_ux(8,1,:,:)= -r_fac*two*(mu_local*(two*vix(1,:,:)*vi(1,:,:)  
     $     + (vix(2,:,:) + viy(1,:,:))*vi(2,:,:)
     $     + vix(3,:,:)*vi(3,:,:)) + nu*rho*vex*ve)
      s_ux(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_ux3
      s_ux(8,4,:,:)=r_fac*4._r8*mu_local*vix(1,:,:)*rho_inv
      s_ux(8,5,:,:)=r_fac*two*mu_local*(vix(2,:,:) + viy(1,:,:))*rho_inv
      s_ux(8,6,:,:)=r_fac*two*(mu_local*vix(3,:,:)*rho_inv + nu*vex)
      s_ux(8,7,:,:)=-r_fac*di*nu*two*vex
      s_ux(8,8,:,:)=r_fac*vi(1,:,:)

      s_uy(8,1,:,:)= -r_fac*two*(mu_local*(two*viy(2,:,:)*vi(2,:,:)  
     $     + (viy(1,:,:) + vix(2,:,:))*vi(1,:,:)
     $     + viy(3,:,:)*vi(3,:,:)) + nu*rho*vey*ve)
      s_uy(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_uy3
      s_uy(8,4,:,:)=r_fac*two*mu_local*(viy(1,:,:) + vix(2,:,:))*rho_inv
      s_uy(8,5,:,:)=r_fac*4._r8*mu_local*viy(2,:,:)*rho_inv
      s_uy(8,6,:,:)=r_fac*two*(mu_local*viy(3,:,:)*rho_inv + nu*vey)
      s_uy(8,7,:,:)=-r_fac*di*nu*two*vey
      s_uy(8,8,:,:)=r_fac*vi(2,:,:)
c-----------------------------------------------------------------------
c     del^2(out-of-plane B).
c-----------------------------------------------------------------------
      IF(hyper_eta > 0.)THEN
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
      USE SolarLnHMHD_mod
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
      USE SolarLnHMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: ksi,phi
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: x,y
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("half_breakout")
         x = SIN(half*pi*(phiscale*phi**2 + phi)/
     $        (phiscale+one))
     $        *((Rmax-one)*(rscale*ksi**2 + ksi)/(rscale + one) + one)
         y = COS(half*pi*(phiscale*phi**2 + phi)/
     $        (phiscale+one))
     $        *((Rmax-one)*(rscale*ksi**2 + ksi)/(rscale + one) + one)
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
      USE SolarLnHMHD_mod
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
      USE SolarLnHMHD_mod
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
