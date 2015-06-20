c-----------------------------------------------------------------------
c     file XMHD.f.
c     contains specifications for extended MHD model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. XMHD_mod.
c     1. XMHD_equil.
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
c     subprogram 0. physics_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE XMHD_mod
      USE transport_mod
      USE extra_mod
      IMPLICIT NONE

      REAL(r8), PARAMETER :: mass_r=5.44617e-4_r8,
     $     gamma=5._r8/3._r8,gfac=gamma/(gamma-one),Bmin=1.e-8_r8,
     $     kappa_min=0.,kappa_max=1.e8

      LOGICAL :: source=.FALSE.
      CHARACTER(16) :: init_type="."
      REAL(r8) :: T_i=1.,T_e=1.,L0=1.,B0=1.,di=1.,
     $     lambda_psi=1.,lx=1.,ly=1.,epsilon=0.,kappa_i=0.,kappa_e=0.,
     $     xi_norm=0.,xe_norm=0.,d_inv=1.,B_eq=0.,kx=1.,ky=1.,rhomin=0.,
     $     resist=0.51_r8,eta_i0=-.96_r8,eta_e0=-3.9757e-4_r8,
     $     fheat=0.71_r8,ieheat=3._r8

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. XMHD_equil.
c     computes equilibrium for magnetic reconnection.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHD_equil(x,y,u,ux,uy,deriv)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: deriv

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: coshy_psi,tanhy_psi
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      coshy_psi=COSH(y/lambda_psi)
      tanhy_psi=TANH(y/lambda_psi)

      SELECT CASE(init_type)
      CASE("GEM")
         u(1,:,:)=one/coshy_psi**2 + rhomin
         u(4,:,:)=-two*T_i/(d_inv*lambda_psi*u(1,:,:)*coshy_psi**2)
         u(7,:,:)=lambda_psi*LOG(coshy_psi)
         u(8,:,:)=T_i*u(1,:,:)
         u(9,:,:)=T_e*u(1,:,:)
         u(12,:,:)=two*T_e/(d_inv*lambda_psi*u(1,:,:)*coshy_psi**2)
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=-two*tanhy_psi/(lambda_psi*coshy_psi**2)
         uy(4,:,:)=-u(4,:,:)*(two*tanhy_psi/lambda_psi
     $        +uy(1,:,:)/u(1,:,:))
         uy(7,:,:)=tanhy_psi
         uy(8,:,:)=T_i*uy(1,:,:)
         uy(9,:,:)=T_e*uy(1,:,:)
         uy(12,:,:)=-u(12,:,:)*(two*tanhy_psi/lambda_psi
     $        +uy(1,:,:)/u(1,:,:))
      CASE("uniform")
         u(1,:,:) = one
         u(5,:,:) = -two*lambda_psi*ATAN(EXP(y/lambda_psi))
         u(7,:,:) = lambda_psi*LOG(coshy_psi)
         u(8,:,:) = T_i*u(1,:,:)
         u(9,:,:) = T_e*u(1,:,:)
         u(10,:,:) = tanhy_psi/(d_inv*lambda_psi*u(1,:,:)*coshy_psi)
         u(12,:,:) = one/(d_inv*lambda_psi*u(1,:,:)*coshy_psi**2)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHD_equil
      END MODULE XMHD_mod
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
      USE XMHD_mod
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
      NAMELIST/XMHD_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,du_diagnose,source,init_type,
     $     T_i,T_e,L0,B0,di,lambda_psi,epsilon,lx,ly,
     $     resist,eta_i0,eta_e0,fheat,ieheat
c-----------------------------------------------------------------------
c     Sample namelist.
c-----------------------------------------------------------------------
c$$$&XMHD_list
c$$$
c$$$	init_type="GEM"
c$$$	source=f
c$$$
c$$$    di=1.                  \\ (c/omega_pi)/L0
c$$$
c$$$	T_i=.25
c$$$	T_e=.25
c$$$	
c$$$	L0=1.                  \\ unit length in cm
c$$$    B0=1.e3                \\ unit B-field in gauss
c$$$	lx=25.6
c$$$	ly=12.8
c$$$	lambda_psi=.5
c$$$
c$$$	epsilon=1.e-2
c$$$	
c$$$	fheat=0.
c$$$/
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=XMHD_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="XMHD"
c-----------------------------------------------------------------------
c     set number of dependent variables.
c-----------------------------------------------------------------------
      nqty=12
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
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_init_parameters(static,ground,adapt_qty)
      USE XMHD_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(10:12)=.TRUE.
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
      CALL MPI_Bcast(T_i,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(T_e,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(L0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(B0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(resist,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta_i0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta_e0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(fheat,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ieheat,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     special assignments.
c-----------------------------------------------------------------------
      B_eq=0

      SELECT CASE(init_type)
      CASE("DispTest1a","DispTest1b","DispTest1c")
         B_eq=one
         kx=.1_r8
         ky=1._r8
      CASE("DispTest2a","DispTest2b","DispTest2c")
         B_eq=one
         kx=.1_r8
         ky=.1_r8
      CASE("DispTest3a","DispTest3b","DispTest3c")
         B_eq=100._r8
         kx=.1_r8
         ky=.1_r8
      CASE("DispTest4")
         B_eq=one
         kx=one
         ky=one
      CASE("GEM","uniform")
         rhomin=0.2_r8
         lambda_psi=.5_r8
         lx=25.6_r8
         ly=12.8_r8
         kx=twopi/lx
         ky=pi/ly
      END SELECT
c-----------------------------------------------------------------------
c     initialize physical constants.
c     kappa_i = tau_i/tau_A
c     kappa_e = (m_i/m_e)*(tau_e/tau_A)
c-----------------------------------------------------------------------
      d_inv=one/di
      kappa_i = 2.4076e-11_r8*(B0**4*L0**5*di**6)
      kappa_e = 7.2950e-10_r8*(B0**4*L0**5*di**6)
      xi_norm = kappa_i*d_inv
      xe_norm = kappa_e*d_inv
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
      USE XMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: ui

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(xpi,1),SIZE(xpi,2)) :: cosx,cosy
      REAL(r8), DIMENSION(SIZE(ui,1),SIZE(ui,2),SIZE(ui,3)) :: uix,uiy
c-----------------------------------------------------------------------
c     initial conditions.
c-----------------------------------------------------------------------
      CALL XMHD_equil(xpi,ypi,ui,uix,uiy,.FALSE.)
      cosx=COS(kx*xpi)
      cosy=COS(ky*ypi)
c-----------------------------------------------------------------------
c     linear dispersion relation tests, initial conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
c-----------------------------------------------------------------------
c     parameter set 1
c-----------------------------------------------------------------------
      CASE("DispTest1a")
         ui(1,:,:)=one+epsilon*9.98961e-1*COS(kx*xpi+ky*ypi) 
         ui(2,:,:)=epsilon*1.86355e-2*COS(kx*xpi+ky*ypi-3.56409e-1) 
         ui(3,:,:)=epsilon*6.72296e-7*COS(kx*xpi+ky*ypi+2.83207) 
         ui(4,:,:)=epsilon*1.24116e-4*COS(kx*xpi+ky*ypi+1.4641)
         ui(5,:,:)=epsilon*3.69746e-3*COS(kx*xpi+ky*ypi-2.90787) 
         ui(6,:,:)=epsilon*3.6997e-2*COS(kx*xpi+ky*ypi-2.90792) 
         ui(7,:,:)=epsilon*3.34868e-4*COS(kx*xpi+ky*ypi+1.48559)
         ui(8,:,:)=T_i+epsilon*1.71437e-4*COS(kx*xpi+ky*ypi-8.27389e-2) 
         ui(9,:,:)=T_e+epsilon*1.66967e-4*COS(kx*xpi+ky*ypi-8.37158e-2) 
         ui(10,:,:)=epsilon*1.86349e-2*COS(kx*xpi+ky*ypi-3.56429e-1) 
         ui(11,:,:)=epsilon*6.08562e-7*COS(kx*xpi+ky*ypi+2.77842) 
         ui(12,:,:)=epsilon*1.5414e-5*COS(kx*xpi+ky*ypi+1.31168) 
      CASE("DispTest1b")
         ui(1,:,:)=one+epsilon*3.21581e-2*COS(kx*xpi+ky*ypi+4.67001e-3) 
         ui(2,:,:)=epsilon*7.5308e-6*COS(kx*xpi+ky*ypi-2.32219) 
         ui(3,:,:)=epsilon*3.21378e-3*COS(kx*xpi+ky*ypi-1.62165e-3) 
         ui(4,:,:)=epsilon*9.89368e-2*COS(kx*xpi+ky*ypi-1.57723)
         ui(5,:,:)=epsilon*1.50668e-4*COS(kx*xpi+ky*ypi+1.55887) 
         ui(6,:,:)=epsilon*9.89133e-1*COS(kx*xpi+ky*ypi) 
         ui(7,:,:)=epsilon*3.05845e-4*COS(kx*xpi+ky*ypi-1.56371)
         ui(8,:,:)=T_i+epsilon*5.37099e-6*COS(kx*xpi+ky*ypi+2.03487e-2) 
         ui(9,:,:)=T_e+epsilon*6.85804e-6*COS(kx*xpi+ky*ypi-1.16239e-2) 
         ui(10,:,:)=epsilon*3.18444e-2*COS(kx*xpi+ky*ypi-1.69621e-3) 
         ui(11,:,:)=epsilon*2.88277e-5*COS(kx*xpi+ky*ypi-1.25048e-2) 
         ui(12,:,:)=epsilon*9.88374e-2*COS(kx*xpi+ky*ypi-1.57724) 
      CASE("DispTest1c")
         ui(1,:,:)=one+epsilon*4.85615e-1*COS(kx*xpi+ky*ypi+8.22472e-4)
         ui(2,:,:)=epsilon*2.61695e-5*COS(kx*xpi+ky*ypi-7.65792e-1) 
         ui(3,:,:)=epsilon*4.88224e-1*COS(kx*xpi+ky*ypi-1.38757e-6) 
         ui(4,:,:)=epsilon*1.49084e-3*COS(kx*xpi+ky*ypi+1.57042)
         ui(5,:,:)=epsilon*4.81868e-6*COS(kx*xpi+ky*ypi-2.63533) 
         ui(6,:,:)=epsilon*1.58665e-1*COS(kx*xpi+ky*ypi-3.14012) 
         ui(7,:,:)=epsilon*4.86093e-1*COS(kx*xpi+ky*ypi-1.57012)
         ui(8,:,:)=T_i+epsilon*8.0936e-5*COS(kx*xpi+ky*ypi+8.08789e-4) 
         ui(9,:,:)=T_e+epsilon*8.09116e-5*COS(kx*xpi+ky*ypi+5.25993e-4) 
         ui(10,:,:)=epsilon*5.08885e-3*COS(kx*xpi+ky*ypi-3.1367) 
         ui(11,:,:)=epsilon*4.88734e-1*COS(kx*xpi+ky*ypi) 
         ui(12,:,:)=epsilon*1.59578e-1*COS(kx*xpi+ky*ypi+1.57146) 
c-----------------------------------------------------------------------
c     parameter set 2
c-----------------------------------------------------------------------
      CASE("DispTest2a")
         ui(1,:,:)=one+epsilon*9.99655e-1*COS(kx*xpi+ky*ypi) 
         ui(2,:,:)=epsilon*1.81611e-2*COS(kx*xpi+ky*ypi-4.83217e-2) 
         ui(3,:,:)=epsilon*3.00348e-6*COS(kx*xpi+ky*ypi+3.01067) 
         ui(4,:,:)=epsilon*1.20718e-5*COS(kx*xpi+ky*ypi+1.47097)
         ui(5,:,:)=epsilon*3.6836e-3*COS(kx*xpi+ky*ypi+3.07893) 
         ui(6,:,:)=epsilon*3.68579e-3*COS(kx*xpi+ky*ypi+3.07889) 
         ui(7,:,:)=epsilon*1.64974e-3*COS(kx*xpi+ky*ypi+1.48654)
         ui(8,:,:)=T_i+epsilon*1.66932e-4*COS(kx*xpi+ky*ypi-8.18364e-2) 
         ui(9,:,:)=T_e+epsilon*1.62938e-4*COS(kx*xpi+ky*ypi-8.2536e-2) 
         ui(10,:,:)=epsilon*1.81611e-2*COS(kx*xpi+ky*ypi-4.83216e-2) 
         ui(11,:,:)=epsilon*2.99642e-6*COS(kx*xpi+ky*ypi+3.01068) 
         ui(12,:,:)=epsilon*1.45817e-6*COS(kx*xpi+ky*ypi+1.35734) 
      CASE("DispTest2b")
         ui(1,:,:)=one+epsilon*6.36133e-3*COS(kx*xpi+ky*ypi+1.06904e-4) 
         ui(2,:,:)=epsilon*5.94464e-7*COS(kx*xpi+ky*ypi-4.36866e-1) 
         ui(3,:,:)=epsilon*6.3525e-3*COS(kx*xpi+ky*ypi+1.84498e-5) 
         ui(4,:,:)=epsilon*9.90358e-2*COS(kx*xpi+ky*ypi-1.57091)
         ui(5,:,:)=epsilon*1.49447e-5*COS(kx*xpi+ky*ypi+1.56841) 
         ui(6,:,:)=epsilon*9.89605e-1*COS(kx*xpi+ky*ypi) 
         ui(7,:,:)=epsilon*3.17184e-2*COS(kx*xpi+ky*ypi-1.57064)
         ui(8,:,:)=T_i+epsilon*1.06146e-6*COS(kx*xpi+ky*ypi+7.96938e-3) 
         ui(9,:,:)=T_e+epsilon*1.21001e-6*COS(kx*xpi+ky*ypi-1.04513e-2) 
         ui(10,:,:)=epsilon*3.18707e-3*COS(kx*xpi+ky*ypi-9.40177e-5) 
         ui(11,:,:)=epsilon*3.16597e-3*COS(kx*xpi+ky*ypi+5.22191e-5) 
         ui(12,:,:)=epsilon*9.88315e-2*COS(kx*xpi+ky*ypi-1.57091) 
      CASE("DispTest2c")
         ui(1,:,:)=one+epsilon*9.70188e-2*COS(kx*xpi+ky*ypi+1.57079) 
         ui(2,:,:)=epsilon*2.32734e-5*COS(kx*xpi+ky*ypi+1.41265) 
         ui(3,:,:)=epsilon*1.37298e-1*COS(kx*xpi+ky*ypi+1.57071) 
         ui(4,:,:)=epsilon*6.22171e-3*COS(kx*xpi+ky*ypi-3.14155)
         ui(5,:,:)=epsilon*4.85115e-6*COS(kx*xpi+ky*ypi-1.3816) 
         ui(6,:,:)=epsilon*8.81634e-2*COS(kx*xpi+ky*ypi-1.57068) 
         ui(7,:,:)=epsilon*9.72029e-1*COS(kx*xpi+ky*ypi)
         ui(8,:,:)=T_i+epsilon*1.61711e-5*COS(kx*xpi+ky*ypi+1.57068) 
         ui(9,:,:)=T_e+epsilon*1.61589e-5*COS(kx*xpi+ky*ypi+1.56871) 
         ui(10,:,:)=epsilon*2.60914e-4*COS(kx*xpi+ky*ypi-1.55663) 
         ui(11,:,:)=epsilon*1.37582e-1*COS(kx*xpi+ky*ypi+1.57071) 
         ui(12,:,:)=epsilon*1.24816e-2*COS(kx*xpi+ky*ypi-3.14157)
c-----------------------------------------------------------------------
c     parameter set 3
c-----------------------------------------------------------------------
      CASE("DispTest3a")
         ui(1,:,:)=one+epsilon*9.99657e-1*COS(kx*xpi+ky*ypi) 
         ui(2,:,:)=epsilon*1.81597e-2*COS(kx*xpi+ky*ypi-4.83112e-2) 
         ui(3,:,:)=epsilon*3.00325e-10*COS(kx*xpi+ky*ypi+3.01069) 
         ui(4,:,:)=epsilon*1.20661e-7*COS(kx*xpi+ky*ypi+1.47102)
         ui(5,:,:)=epsilon*3.68331e-3*COS(kx*xpi+ky*ypi+3.07894) 
         ui(6,:,:)=epsilon*3.68331e-3*COS(kx*xpi+ky*ypi+3.07894) 
         ui(7,:,:)=epsilon*1.64948e-5*COS(kx*xpi+ky*ypi+1.48656)
         ui(8,:,:)=T_i+epsilon*1.66933e-4*COS(kx*xpi+ky*ypi-8.18296e-2) 
         ui(9,:,:)=T_e+epsilon*1.62939e-4*COS(kx*xpi+ky*ypi-8.25298e-2) 
         ui(10,:,:)=epsilon*1.81597e-2*COS(kx*xpi+ky*ypi-4.83112e-2) 
         ui(11,:,:)=epsilon*2.99619e-10*COS(kx*xpi+ky*ypi+3.0107) 
         ui(12,:,:)=epsilon*1.4542e-8*COS(kx*xpi+ky*ypi+1.35727) 
      CASE("DispTest3b")
         ui(1,:,:)=one+epsilon*6.36143e-5*COS(kx*xpi+ky*ypi+1.06993e-6) 
         ui(2,:,:)=epsilon*1.73425e-6*COS(kx*xpi+ky*ypi-3.14017) 
         ui(3,:,:)=epsilon*6.35487e-3*COS(kx*xpi+ky*ypi+1.77638e-7) 
         ui(4,:,:)=epsilon*9.90377e-2*COS(kx*xpi+ky*ypi-1.5708)
         ui(5,:,:)=epsilon*5.78081e-7*COS(kx*xpi+ky*ypi+2.61374e-1) 
         ui(6,:,:)=epsilon*9.89624e-1*COS(kx*xpi+ky*ypi) 
         ui(7,:,:)=epsilon*3.17415e-2*COS(kx*xpi+ky*ypi-1.57079)
         ui(8,:,:)=T_i+epsilon*1.06024e-8*COS(kx*xpi+ky*ypi+8.11356e-5) 
         ui(9,:,:)=T_e+epsilon*1.21127e-8*COS(kx*xpi+ky*ypi-1.05662e-4) 
         ui(10,:,:)=epsilon*3.18485e-3*COS(kx*xpi+ky*ypi-9.26927e-7) 
         ui(11,:,:)=epsilon*3.16828e-3*COS(kx*xpi+ky*ypi+5.08122e-7) 
         ui(12,:,:)=epsilon*9.88333e-2*COS(kx*xpi+ky*ypi-1.5708) 
      CASE("DispTest3c")
         ui(1,:,:)=one+epsilon*9.74625e-4*COS(kx*xpi+ky*ypi+1.5708) 
         ui(2,:,:)=epsilon*1.61764e-7*COS(kx*xpi+ky*ypi-1.34285) 
         ui(3,:,:)=epsilon*1.37937e-1*COS(kx*xpi+ky*ypi-1.5708) 
         ui(4,:,:)=epsilon*6.25342e-3*COS(kx*xpi+ky*ypi+3.14159)
         ui(5,:,:)=epsilon*5.13251e-8*COS(kx*xpi+ky*ypi+1.38865) 
         ui(6,:,:)=epsilon*8.86006e-2*COS(kx*xpi+ky*ypi+1.5708)
         ui(7,:,:)=epsilon*9.76637e-1*COS(kx*xpi+ky*ypi)
         ui(8,:,:)=T_i+epsilon*1.62437e-7*COS(kx*xpi+ky*ypi+1.5708) 
         ui(9,:,:)=T_e+epsilon*1.62342e-7*COS(kx*xpi+ky*ypi+1.57082) 
         ui(10,:,:)=epsilon*2.85136e-4*COS(kx*xpi+ky*ypi+1.57067) 
         ui(11,:,:)=epsilon*1.38223e-1*COS(kx*xpi+ky*ypi-1.5708) 
         ui(12,:,:)=epsilon*1.2543e-2*COS(kx*xpi+ky*ypi+3.14159)
c-----------------------------------------------------------------------
c     parameter set 4
c-----------------------------------------------------------------------
      CASE("DispTest4")
         ui(1,:,:)=one+epsilon*9.99722e-1*COS(kx*xpi+ky*ypi)
         ui(2,:,:)=epsilon*1.65314e-2*COS(kx*xpi+ky*ypi-1.66695e-1)
         ui(3,:,:)=epsilon*2.6895e-6*COS(kx*xpi+ky*ypi+2.80983)
         ui(4,:,:)=epsilon*8.87881e-5*COS(kx*xpi+ky*ypi+1.33658)
         ui(5,:,:)=epsilon*2.17906e-3*COS(kx*xpi+ky*ypi+3.06992)
         ui(6,:,:)=epsilon*2.18048e-3*COS(kx*xpi+ky*ypi+3.06979)
         ui(7,:,:)=epsilon*1.34923e-4*COS(kx*xpi+ky*ypi+1.3739)
         ui(8,:,:)=T_i+epsilon*1.63458e-4*COS(kx*xpi+ky*ypi-1.93861e-1)
         ui(9,:,:)=T_e+epsilon*1.06311e-4*COS(kx*xpi+ky*ypi-1.43143e-1)
         ui(10,:,:)=epsilon*1.6531e-2*COS(kx*xpi+ky*ypi-1.66692e-1)
         ui(11,:,:)=epsilon*2.22431e-6*COS(kx*xpi+ky*ypi+2.79605)
         ui(12,:,:)=epsilon*3.78728e-6*COS(kx*xpi+ky*ypi+3.09222e-1)
c-----------------------------------------------------------------------
c     magnetic reconnection, initial conditions.
c-----------------------------------------------------------------------
      CASE("GEM","uniform")
         ui(7,:,:)=ui(7,:,:)+epsilon*cosx*cosy
c-----------------------------------------------------------------------
c     abort if not recognized.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT

      DO i=2,7
         ui(i,:,:)=ui(1,:,:)*ui(i,:,:)
      ENDDO
      DO i=10,12
         ui(i,:,:)=ui(1,:,:)*ui(i,:,:)
      ENDDO

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
      USE XMHD_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     GEM magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("GEM","uniform")
c-----------------------------------------------------------------------
c     density
c-----------------------------------------------------------------------
         top%bc_type(1)="natural"
         bottom%bc_type(1)="natural"
c-----------------------------------------------------------------------
c     ion momentum
c-----------------------------------------------------------------------
         top%bc_type(2:4)="robin"
         bottom%bc_type(2:4)="robin"
         top%static(2:4)=.TRUE.
         bottom%static(2:4)=.TRUE.
c-----------------------------------------------------------------------
c     vector potential
c-----------------------------------------------------------------------
         top%bc_type(5:7)="robin"
         bottom%bc_type(5:7)="robin"
         top%static(5)=.FALSE.
         bottom%static(5)=.FALSE.
         top%static(6)=.TRUE.
         bottom%static(6)=.TRUE.
         top%static(7)=.FALSE.
         bottom%static(7)=.FALSE.
c-----------------------------------------------------------------------
c     pressure
c-----------------------------------------------------------------------
         top%bc_type(8)="zeroflux"
         bottom%bc_type(8)="zeroflux"
         top%bc_type(9)="normflux"
         bottom%bc_type(9)="normflux"
c-----------------------------------------------------------------------
c     electron momentum
c-----------------------------------------------------------------------
         top%bc_type(10:12)="robin"
         bottom%bc_type(10:12)="robin"
         top%static(10:12)=.TRUE.
         bottom%static(10:12)=.TRUE.
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
      USE XMHD_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     zero out output.
c-----------------------------------------------------------------------
      c=0
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("GEM","uniform")
         SELECT CASE(lrtb)
         CASE("bottom","top")
c-----------------------------------------------------------------------
c     ion momentum
c-----------------------------------------------------------------------
            c(2,:,:) = uy(2,:,:)*u(1,:,:)-uy(1,:,:)*u(2,:,:)
            c(3,:,:) = u(3,:,:)
            c(4,:,:) = uy(4,:,:)*u(1,:,:)-uy(1,:,:)*u(4,:,:)
c-----------------------------------------------------------------------
c     vector potential
c-----------------------------------------------------------------------
            c(5,:,:) = -u(5,:,:)*(ux(2,:,:)+uy(3,:,:))/u(1,:,:)
            c(6,:,:) = uy(6,:,:)*u(1,:,:)-uy(1,:,:)*u(6,:,:)
            c(7,:,:) = -u(7,:,:)*(ux(2,:,:)+uy(3,:,:))/u(1,:,:)
c-----------------------------------------------------------------------
c     pressure
c-----------------------------------------------------------------------
            WHERE(nhat(2,:,:)*u(11,:,:) > 0)
               c(9,:,:) = nhat(2,:,:)*gfac*u(9,:,:)*u(11,:,:)/u(1,:,:)
            ELSEWHERE
               c(9,:,:) = nhat(2,:,:)*gfac*T_e*u(11,:,:)
            END WHERE
c-----------------------------------------------------------------------
c     electron momentum
c-----------------------------------------------------------------------
            c(10,:,:) = u(2,:,:)-u(10,:,:)
            c(11,:,:) = uy(3,:,:)-uy(11,:,:)
            c(12,:,:) = u(4,:,:)-u(12,:,:)
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
      USE XMHD_mod
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
c     give values to c_u,c_ux,c_uy.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("GEM","uniform")
         SELECT CASE(lrtb)
         CASE("bottom","top")
c-----------------------------------------------------------------------
c     ion momentum
c-----------------------------------------------------------------------
            c_u(2,1,:,:)=uy(2,:,:)
            c_u(2,2,:,:)=-uy(1,:,:)
            c_uy(2,1,:,:)=-u(2,:,:)
            c_uy(2,2,:,:)=u(1,:,:)
            c_u(3,3,:,:)=one
            c_u(4,1,:,:)=uy(4,:,:)
            c_u(4,4,:,:)=-uy(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)
            c_uy(4,4,:,:)=u(1,:,:)
c-----------------------------------------------------------------------
c     vector potential
c-----------------------------------------------------------------------
            c_u(5,1,:,:)=u(5,:,:)*(ux(2,:,:)+uy(3,:,:))/u(1,:,:)**2
            c_u(5,5,:,:)=-(ux(2,:,:)+uy(3,:,:))/u(1,:,:)
            c_ux(5,2,:,:)=-u(5,:,:)/u(1,:,:)
            c_uy(5,3,:,:)=-u(5,:,:)/u(1,:,:)
            c_u(6,1,:,:)=uy(6,:,:)
            c_u(6,6,:,:)=-uy(1,:,:)
            c_uy(6,1,:,:)=-u(6,:,:)
            c_uy(6,6,:,:)=u(1,:,:)
            c_u(7,1,:,:)=u(7,:,:)*(ux(2,:,:)+uy(3,:,:))/u(1,:,:)**2
            c_u(7,7,:,:)=-(ux(2,:,:)+uy(3,:,:))/u(1,:,:)
            c_ux(7,2,:,:)=-u(7,:,:)/u(1,:,:)
            c_uy(7,3,:,:)=-u(7,:,:)/u(1,:,:)
c-----------------------------------------------------------------------
c     pressure
c-----------------------------------------------------------------------
            WHERE(nhat(2,:,:)*u(11,:,:) > 0)
               c_u(9,1,:,:) = -nhat(2,:,:)*gfac*u(9,:,:)*u(11,:,:)
     $              /u(1,:,:)**2
               c_u(9,9,:,:) = nhat(2,:,:)*gfac*u(11,:,:)/u(1,:,:)
               c_u(9,11,:,:) = nhat(2,:,:)*gfac*u(9,:,:)/u(1,:,:)
            ELSEWHERE
               c_u(9,11,:,:) = nhat(2,:,:)*gfac*T_e
            END WHERE
c-----------------------------------------------------------------------
c     electron momentum
c-----------------------------------------------------------------------
            c_u(10,2,:,:) = one
            c_u(10,10,:,:) = -one
            c_uy(11,3,:,:) = one
            c_uy(11,11,:,:) = -one
            c_u(12,4,:,:) = one
            c_u(12,12,:,:) = -one
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
      USE XMHD_mod
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

      SELECT CASE(init_type)
      CASE("GEM","uniform")
         SELECT CASE(lrtb)
         CASE("top","bottom")
c-----------------------------------------------------------------------
c     vector potential
c-----------------------------------------------------------------------
            mass(5,5,:,:)=one
            mass(7,7,:,:)=one
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
      USE XMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: invBsq,bx,by,bz,
     $     Ti,Tix,Tiy,Te,Tex,Tey,n_inv,nx_inv,ny_inv,
     $     kfaci,kperpi,kface,kperpe
      REAL(r8), DIMENSION(3,SIZE(x,1),SIZE(x,2)) :: vi,vix,viy,
     $     ve,vex,vey,A,Ax,Ay,BdotTi,BdotTe

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     inverse density.
c-----------------------------------------------------------------------
      n_inv=one/u(1,:,:)
      nx_inv=-ux(1,:,:)*n_inv**2
      ny_inv=-uy(1,:,:)*n_inv**2
c-----------------------------------------------------------------------
c     velocities and vector potential
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:) = u(i+1,:,:)*n_inv
         vix(i,:,:) = ux(i+1,:,:)*n_inv + u(i+1,:,:)*nx_inv
         viy(i,:,:) = uy(i+1,:,:)*n_inv + u(i+1,:,:)*ny_inv
         A(i,:,:) = u(i+4,:,:)*n_inv
         Ax(i,:,:) = ux(i+4,:,:)*n_inv + u(i+4,:,:)*nx_inv
         Ay(i,:,:) = uy(i+4,:,:)*n_inv + u(i+4,:,:)*ny_inv
         ve(i,:,:) = u(i+9,:,:)*n_inv
         vex(i,:,:) = ux(i+9,:,:)*n_inv + u(i+9,:,:)*nx_inv
         vey(i,:,:) = uy(i+9,:,:)*n_inv + u(i+9,:,:)*ny_inv
      ENDDO
c-----------------------------------------------------------------------
c     B-field.
c-----------------------------------------------------------------------
      Ay(3,:,:) = Ay(3,:,:)+B_eq
      bx = Ay(3,:,:)
      by = -Ax(3,:,:)
      bz = Ax(2,:,:) - Ay(1,:,:)
      invBsq = one/(bx**2 + by**2 + bz**2)

      WHERE((bx**2 + by**2 + bz**2) < Bmin**2)invBsq = zero
c-----------------------------------------------------------------------
c     Temperatures.
c-----------------------------------------------------------------------
      Ti = u(8,:,:)*n_inv
      Tix = ux(8,:,:)*n_inv + u(8,:,:)*nx_inv
      Tiy = uy(8,:,:)*n_inv + u(8,:,:)*ny_inv
      Te = u(9,:,:)*n_inv
      Tex = ux(9,:,:)*n_inv + u(9,:,:)*nx_inv
      Tey = uy(9,:,:)*n_inv + u(9,:,:)*ny_inv

      CALL transport_BdotT(bx,by,bz,Tix,Tiy,BdotTi)
      CALL transport_BdotT(bx,by,bz,Tex,Tey,BdotTe)
      CALL transport_kbrag(u(1,:,:),u(8,:,:),zero,
     $     (bx**2 + by**2 + bz**2),zero,kappa_i,zero,xi_norm,kappa_min,
     $     kappa_max,kperpi,kfaci)
      CALL transport_kbrag(u(1,:,:),u(9,:,:),one,
     $     (bx**2 + by**2 + bz**2),kappa_e,zero,xe_norm,zero,kappa_min,
     $     kappa_max,kperpe,kface)
c-----------------------------------------------------------------------
c     density
c-----------------------------------------------------------------------
      fx(1,:,:)=u(2,:,:)
      fy(1,:,:)=u(3,:,:)
c-----------------------------------------------------------------------
c     total x-momentum
c-----------------------------------------------------------------------
      fx(2,:,:)=u(8,:,:)+u(9,:,:)+u(2,:,:)*vi(1,:,:)
     $     +mass_r*u(10,:,:)*ve(1,:,:)
     $     +eta_i0*kappa_i*Ti**2.5*(one+third)*vix(1,:,:)
     $     +eta_e0*kappa_e*Te**2.5*(one+third)*vex(1,:,:)

      fy(2,:,:)=u(3,:,:)*vi(1,:,:)+mass_r*u(11,:,:)*ve(1,:,:)
     $     +eta_i0*kappa_i*Ti**2.5*(viy(1,:,:)+third*vix(2,:,:))
     $     +eta_e0*kappa_e*Te**2.5*(vey(1,:,:)+third*vex(2,:,:))

      s(2,:,:)=((u(3,:,:)-u(11,:,:))*(Ax(2,:,:)-Ay(1,:,:))
     $     +(u(4,:,:)-u(12,:,:))*Ax(3,:,:))*d_inv
c-----------------------------------------------------------------------
c     total y-momentum
c-----------------------------------------------------------------------
      fx(3,:,:)=u(2,:,:)*vi(2,:,:)+mass_r*u(10,:,:)*ve(2,:,:)
     $     +eta_i0*kappa_i*Ti**2.5*(vix(2,:,:)+third*viy(1,:,:))
     $     +eta_e0*kappa_e*Te**2.5*(vex(2,:,:)+third*vey(1,:,:))

      fy(3,:,:)=u(8,:,:)+u(9,:,:)+u(3,:,:)*vi(2,:,:)
     $     +mass_r*u(11,:,:)*ve(2,:,:)
     $     +eta_i0*kappa_i*Ti**2.5*(one+third)*viy(2,:,:)
     $     +eta_e0*kappa_e*Te**2.5*(one+third)*vey(2,:,:)

      s(3,:,:)=((u(2,:,:)-u(10,:,:))*(Ay(1,:,:)-Ax(2,:,:))
     $     +(u(4,:,:)-u(12,:,:))*Ay(3,:,:))*d_inv
c-----------------------------------------------------------------------
c     total z-momentum
c-----------------------------------------------------------------------
      fx(4,:,:)=u(2,:,:)*vi(3,:,:)+mass_r*u(10,:,:)*ve(3,:,:)
     $     + eta_i0*kappa_i*Ti**2.5*vix(3,:,:)
     $     + eta_e0*kappa_e*Te**2.5*vex(3,:,:)

      fy(4,:,:)=u(3,:,:)*vi(3,:,:)+mass_r*u(11,:,:)*ve(3,:,:)
     $     + eta_i0*kappa_i*Ti**2.5*viy(3,:,:)
     $     + eta_e0*kappa_e*Te**2.5*vey(3,:,:)

      s(4,:,:)=-((u(2,:,:)-u(10,:,:))*Ax(3,:,:)
     $     + (u(3,:,:)-u(11,:,:))*Ay(3,:,:))*d_inv
c-----------------------------------------------------------------------
c     electron x-momentum
c-----------------------------------------------------------------------
      fx(5,:,:)=-u(9,:,:)-u(10,:,:)*(mass_r*ve(1,:,:)-d_inv*A(1,:,:))
     $     - eta_e0*kappa_e*Te**2.5*(one+third)*vex(1,:,:)

      fy(5,:,:)=-u(11,:,:)*(mass_r*ve(1,:,:)-d_inv*A(1,:,:))
     $     - eta_e0*kappa_e*Te**2.5*(vey(1,:,:)+third*vex(2,:,:))

      s(5,:,:)=(u(10,:,:)*Ax(1,:,:)+u(11,:,:)*Ax(2,:,:)
     $     + u(12,:,:)*Ax(3,:,:))*d_inv
     $     - resist*(u(2,:,:)-u(10,:,:))*u(1,:,:)/(kappa_e*Te**1.5)
     $     + fheat*u(1,:,:)*BdotTe(1,:,:)
c-----------------------------------------------------------------------
c     electron y-momentum
c-----------------------------------------------------------------------
      fx(6,:,:)=-u(10,:,:)*(mass_r*ve(2,:,:)-d_inv*A(2,:,:))
     $     - eta_e0*kappa_e*Te**2.5*(vex(2,:,:)+third*vey(1,:,:))
      fy(6,:,:)=-u(9,:,:)-u(11,:,:)*(mass_r*ve(2,:,:)-d_inv*A(2,:,:))
     $     - eta_e0*kappa_e*Te**2.5*(one+third)*vey(2,:,:)

      s(6,:,:)=(u(10,:,:)*Ay(1,:,:)+u(11,:,:)*Ay(2,:,:)
     $     + u(12,:,:)*Ay(3,:,:))*d_inv
     $     - resist*(u(3,:,:)-u(11,:,:))*u(1,:,:)/(kappa_e*Te**1.5)
     $     + fheat*u(1,:,:)*BdotTe(2,:,:)
c-----------------------------------------------------------------------
c     electron z-momentum
c-----------------------------------------------------------------------
      fx(7,:,:)=-u(10,:,:)*(mass_r*ve(3,:,:)-d_inv*A(3,:,:))
     $     - eta_e0*kappa_e*Te**2.5*vex(3,:,:)
      fy(7,:,:)=-u(11,:,:)*(mass_r*ve(3,:,:)-d_inv*A(3,:,:))
     $     - eta_e0*kappa_e*Te**2.5*vey(3,:,:)

      s(7,:,:)=-d_inv*B_eq*u(11,:,:)
     $     - resist*(u(4,:,:)-u(12,:,:))*u(1,:,:)/(kappa_e*Te**1.5)
     $     + fheat*u(1,:,:)*BdotTe(3,:,:)
c-----------------------------------------------------------------------
c     ion pressure
c-----------------------------------------------------------------------
      fx(8,:,:) = gfac*u(8,:,:)*vi(1,:,:)
     $     - kfaci*BdotTi(1,:,:) - kperpi*Tix
      fy(8,:,:) = gfac*u(8,:,:)*vi(2,:,:)
     $     - kfaci*BdotTi(2,:,:) - kperpi*Tiy
      
      s(8,:,:)= vi(1,:,:)*ux(8,:,:)+vi(2,:,:)*uy(8,:,:)
     $     - eta_i0*kappa_i*Ti**2.5
     $     *((one+third)*(vix(1,:,:)**2+viy(2,:,:)**2)
     $     +two*third*vix(2,:,:)*viy(1,:,:)
     $     +viy(1,:,:)**2+vix(2,:,:)**2+vix(3,:,:)**2+viy(3,:,:)**2)
     $     + ieheat*(u(9,:,:)-u(8,:,:))*u(1,:,:)/(kappa_e*Te**1.5)
c-----------------------------------------------------------------------
c     electron pressure
c-----------------------------------------------------------------------
      fx(9,:,:)=gfac*u(9,:,:)*ve(1,:,:)
     $     - kface*BdotTe(1,:,:) - kperpe*Tex
     $     + fheat*Te*invBsq*((u(10,:,:)-u(2,:,:))*bx
     $     + (u(11,:,:)-u(3,:,:))*by+(u(12,:,:)-u(4,:,:))*bz)*bx
      fy(9,:,:)=gfac*u(9,:,:)*ve(2,:,:)
     $     - kface*BdotTe(2,:,:) - kperpe*Tey
     $     + fheat*Te*invBsq*((u(10,:,:)-u(2,:,:))*bx
     $     + (u(11,:,:)-u(3,:,:))*by+(u(12,:,:)-u(4,:,:))*bz)*by
      
      s(9,:,:) = ve(1,:,:)*ux(9,:,:)+ve(2,:,:)*uy(9,:,:)
     $     - eta_e0*kappa_e*Te**2.5
     $     *((one+third)*(vex(1,:,:)**2+vey(2,:,:)**2)
     $     + two*third*vex(2,:,:)*vey(1,:,:)
     $     + vey(1,:,:)**2+vex(2,:,:)**2+vex(3,:,:)**2+vey(3,:,:)**2)
     $     + (ieheat*(u(8,:,:)-u(9,:,:))*u(1,:,:)
     $     + resist*((u(10,:,:)-u(2,:,:))**2 + (u(11,:,:)-u(3,:,:))**2 
     $     + (u(12,:,:)-u(4,:,:))**2))/(kappa_e*Te**1.5)
     $     + fheat*((u(10,:,:)-u(2,:,:))*BdotTe(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe(2,:,:) 
     $     + (u(12,:,:)-u(4,:,:))*BdotTe(3,:,:))
c-----------------------------------------------------------------------
c     x-component of Ampere's Law
c-----------------------------------------------------------------------
      fx(10,:,:)=-Ay(2,:,:)
      fy(10,:,:)=Ay(1,:,:)
      s(10,:,:)=d_inv*(u(10,:,:)-u(2,:,:))
c-----------------------------------------------------------------------
c     y-component of Ampere's Law
c-----------------------------------------------------------------------
      fx(11,:,:)=Ax(2,:,:)
      fy(11,:,:)=-Ax(1,:,:)
      s(11,:,:)=d_inv*(u(11,:,:)-u(3,:,:))      
c-----------------------------------------------------------------------
c     z-component of Ampere's Law
c-----------------------------------------------------------------------
      fx(12,:,:)=Ax(3,:,:)
      fy(12,:,:)=Ay(3,:,:)
      s(12,:,:)=d_inv*(u(12,:,:)-u(4,:,:))

      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("GEM")
            CALL XMHD_equil(x,y,u0,u0x,u0y,.TRUE.)
            DO i=2,7
               u0x(i,:,:)=u0x(1,:,:)*u0(i,:,:)+u0(1,:,:)*u0x(i,:,:)
               u0y(i,:,:)=u0y(1,:,:)*u0(i,:,:)+u0(1,:,:)*u0y(i,:,:)
               u0(i,:,:)=u0(1,:,:)*u0(i,:,:)
            ENDDO
            DO i=10,12
               u0x(i,:,:)=u0x(1,:,:)*u0(i,:,:)+u0(1,:,:)*u0x(i,:,:)
               u0y(i,:,:)=u0y(1,:,:)*u0(i,:,:)+u0(1,:,:)*u0y(i,:,:)
               u0(i,:,:)=u0(1,:,:)*u0(i,:,:)
            ENDDO
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
      USE XMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i
      
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: invBsq,bx,by,bz,
     $     PidelV,PedelV,Ti,Tix,Tiy,Te,Tex,Tey,Ti_un,Tix_un,Tiy_un,
     $     Te_un,Tex_un,Tey_un,n_inv,nx_inv,ny_inv,delVpar,delVpar_un,
     $     delVpar_uxn,delVpar_uyn,delVsq,kfaci,kperpi,kface,kperpe,
     $     kpari_un,kpari_pi,kperpi_un,kperpi_pi,kperpi_bsq,
     $     kpare_un,kpare_pe,kperpe_un,kperpe_pe,kperpe_bsq
      REAL(r8), DIMENSION(3,SIZE(x,1),SIZE(x,2)) :: vi,vix,viy,
     $     ve,vex,vey,A,Ax,Ay,vi_un,vix_un,viy_un,
     $     ve_un,vex_un,vey_un,A_un,Ax_un,Ay_un,
     $     BdotTi,BdotTi_bx,BdotTi_by,BdotTi_bz,BdotTi_Tix,BdotTi_Tiy,
     $     BdotTe,BdotTe_bx,BdotTe_by,BdotTe_bz,BdotTe_Tex,BdotTe_Tey,
     $     BdotTi_un,BdotTi_uxn,BdotTi_uyn,
     $     BdotTe_un,BdotTe_uxn,BdotTe_uyn
      REAL(r8), DIMENSION(4,SIZE(x,1),SIZE(x,2)) :: Bsq_u,Bsq_ux,Bsq_uy
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
c     inverse density.
c-----------------------------------------------------------------------
      n_inv=one/u(1,:,:)
      nx_inv=-ux(1,:,:)*n_inv**2
      ny_inv=-uy(1,:,:)*n_inv**2
c-----------------------------------------------------------------------
c     velocities and vector potential
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:) = u(i+1,:,:)*n_inv
         vix(i,:,:) = ux(i+1,:,:)*n_inv + u(i+1,:,:)*nx_inv
         viy(i,:,:) = uy(i+1,:,:)*n_inv + u(i+1,:,:)*ny_inv
         A(i,:,:) = u(i+4,:,:)*n_inv
         Ax(i,:,:) = ux(i+4,:,:)*n_inv + u(i+4,:,:)*nx_inv
         Ay(i,:,:) = uy(i+4,:,:)*n_inv + u(i+4,:,:)*ny_inv
         ve(i,:,:) = u(i+9,:,:)*n_inv
         vex(i,:,:) = ux(i+9,:,:)*n_inv + u(i+9,:,:)*nx_inv
         vey(i,:,:) = uy(i+9,:,:)*n_inv + u(i+9,:,:)*ny_inv

         vi_un(i,:,:) = -vi(i,:,:)*n_inv
         vix_un(i,:,:) = -(vi(i,:,:)*nx_inv + vix(i,:,:)*n_inv)
         viy_un(i,:,:) = -(vi(i,:,:)*ny_inv + viy(i,:,:)*n_inv)
         A_un(i,:,:) = -A(i,:,:)*n_inv
         Ax_un(i,:,:) = -(A(i,:,:)*nx_inv + Ax(i,:,:)*n_inv)
         Ay_un(i,:,:) = -(A(i,:,:)*ny_inv + Ay(i,:,:)*n_inv)
         ve_un(i,:,:) = -ve(i,:,:)*n_inv
         vex_un(i,:,:) = -(ve(i,:,:)*nx_inv + vex(i,:,:)*n_inv)
         vey_un(i,:,:) = -(ve(i,:,:)*ny_inv + vey(i,:,:)*n_inv)
      ENDDO
c-----------------------------------------------------------------------
c     B-field.
c-----------------------------------------------------------------------
      Ay(3,:,:) = Ay(3,:,:) + B_eq
      bx = Ay(3,:,:)
      by = -Ax(3,:,:)
      bz = Ax(2,:,:) - Ay(1,:,:)
      invBsq = one/(bx**2 + by**2 + bz**2)

      Bsq_u(1,:,:) =
     $     two*((Ay(1,:,:) - Ax(2,:,:))*(Ay_un(1,:,:) - Ax_un(2,:,:))
     $     + Ax(3,:,:)*Ax_un(3,:,:) + Ay(3,:,:)*Ay_un(3,:,:))
      Bsq_u(2,:,:) = two*(Ay(1,:,:) - Ax(2,:,:))*ny_inv
      Bsq_u(3,:,:) = -two*(Ay(1,:,:) - Ax(2,:,:))*nx_inv
      Bsq_u(4,:,:) = two*(Ax(3,:,:)*nx_inv + Ay(3,:,:)*ny_inv)

      Bsq_ux = 0
      Bsq_ux(1,:,:) = two*(-(Ay(1,:,:) - Ax(2,:,:))*A_un(2,:,:)
     $     + Ax(3,:,:)*A_un(3,:,:))
      Bsq_ux(3,:,:) = -two*(Ay(1,:,:) - Ax(2,:,:))*n_inv
      Bsq_ux(4,:,:) = two*Ax(3,:,:)*n_inv

      Bsq_uy = 0
      Bsq_uy(1,:,:) = two*((Ay(1,:,:) - Ax(2,:,:))*A_un(1,:,:)
     $     + Ay(3,:,:)*A_un(3,:,:))
      Bsq_uy(2,:,:) = two*(Ay(1,:,:) - Ax(2,:,:))*n_inv
      Bsq_uy(4,:,:) = two*Ay(3,:,:)*n_inv

      WHERE((bx**2 + by**2 + bz**2) < Bmin**2)invBsq = zero
      DO i=1,4
         WHERE((bx**2 + by**2 + bz**2) < Bmin**2)
            Bsq_u(i,:,:) = zero
            Bsq_ux(i,:,:) = zero
            Bsq_uy(i,:,:) = zero
         END WHERE
      ENDDO
c-----------------------------------------------------------------------
c     Temperatures.
c-----------------------------------------------------------------------
      Ti = u(8,:,:)*n_inv
      Tix = ux(8,:,:)*n_inv + u(8,:,:)*nx_inv
      Tiy = uy(8,:,:)*n_inv + u(8,:,:)*ny_inv
      Te = u(9,:,:)*n_inv
      Tex = ux(9,:,:)*n_inv + u(9,:,:)*nx_inv
      Tey = uy(9,:,:)*n_inv + u(9,:,:)*ny_inv

      Ti_un = -Ti*n_inv
      Tix_un = -(Ti*nx_inv + Tix*n_inv)
      Tiy_un = -(Ti*ny_inv + Tiy*n_inv)
      Te_un = -Te*n_inv
      Tex_un = -(Te*nx_inv + Tex*n_inv)
      Tey_un = -(Te*ny_inv + Tey*n_inv)

      CALL transport_kbrag_u(u(1,:,:),u(8,:,:),zero,
     $     (bx**2 + by**2 + bz**2),zero,kappa_i,zero,xi_norm,kappa_min,
     $     kappa_max,kperpi,kfaci,kpari_un,kpari_pi,kperpi_un,kperpi_pi,
     $     kperpi_bsq)
      CALL transport_kbrag_u(u(1,:,:),u(9,:,:),one,
     $     (bx**2 + by**2 + bz**2),kappa_e,zero,xe_norm,zero,kappa_min,
     $     kappa_max,kperpe,kface,kpare_un,kpare_pe,kperpe_un,kperpe_pe,
     $     kperpe_bsq)
c-----------------------------------------------------------------------
c     complex gradient/velocity related expressions.
c-----------------------------------------------------------------------
      delVsq=(u(10,:,:)-u(2,:,:))**2+(u(11,:,:)-u(3,:,:))**2
     $     +(u(12,:,:)-u(4,:,:))**2
      PidelV=(one+third)*(vix(1,:,:)**2+viy(2,:,:)**2)
     $     +two*third*vix(2,:,:)*viy(1,:,:)
     $     +viy(1,:,:)**2+vix(2,:,:)**2+vix(3,:,:)**2+viy(3,:,:)**2
      PedelV=(one+third)*(vex(1,:,:)**2+vey(2,:,:)**2)
     $     +two*third*vex(2,:,:)*vey(1,:,:)
     $     +vey(1,:,:)**2+vex(2,:,:)**2+vex(3,:,:)**2+vey(3,:,:)**2

      CALL transport_BdotT_u(bx,by,bz,Tix,Tiy,BdotTi,
     $     BdotTi_bx,BdotTi_by,BdotTi_bz,BdotTi_Tix,BdotTi_Tiy)
      CALL transport_BdotT_u(bx,by,bz,Tex,Tey,BdotTe,
     $     BdotTe_bx,BdotTe_by,BdotTe_bz,BdotTe_Tex,BdotTe_Tey)
      DO i=1,3
         BdotTi_un(i,:,:) = BdotTi_bx(i,:,:)*Ay_un(3,:,:)
     $        - BdotTi_by(i,:,:)*Ax_un(3,:,:)
     $        + BdotTi_bz(i,:,:)*(Ax_un(2,:,:) - Ay_un(1,:,:))
     $        + BdotTi_Tix(i,:,:)*Tix_un + BdotTi_Tiy(i,:,:)*Tiy_un
         BdotTi_uxn(i,:,:) = - BdotTi_by(i,:,:)*A_un(3,:,:)
     $        + BdotTi_bz(i,:,:)*A_un(2,:,:) + BdotTi_Tix(i,:,:)*Ti_un
         BdotTi_uyn(i,:,:) = BdotTi_bx(i,:,:)*A_un(3,:,:)
     $        - BdotTi_bz(i,:,:)*A_un(1,:,:) + BdotTi_Tiy(i,:,:)*Ti_un
         BdotTe_un(i,:,:) = BdotTe_bx(i,:,:)*Ay_un(3,:,:)
     $        - BdotTe_by(i,:,:)*Ax_un(3,:,:)
     $        + BdotTe_bz(i,:,:)*(Ax_un(2,:,:) - Ay_un(1,:,:))
     $        + BdotTe_Tex(i,:,:)*Tex_un + BdotTe_Tey(i,:,:)*Tey_un
         BdotTe_uxn(i,:,:) = - BdotTe_by(i,:,:)*A_un(3,:,:)
     $        + BdotTe_bz(i,:,:)*A_un(2,:,:) + BdotTe_Tex(i,:,:)*Te_un
         BdotTe_uyn(i,:,:) = BdotTe_bx(i,:,:)*A_un(3,:,:)
     $        - BdotTe_bz(i,:,:)*A_un(1,:,:) + BdotTe_Tey(i,:,:)*Te_un
      ENDDO
      delVpar=((u(10,:,:)-u(2,:,:))*bx+(u(11,:,:)-u(3,:,:))*by+
     $     (u(12,:,:)-u(4,:,:))*bz)*invBsq
      delVpar_un=((u(10,:,:)-u(2,:,:))*Ay_un(3,:,:)
     $     - (u(11,:,:)-u(3,:,:))*Ax_un(3,:,:)
     $     + (u(12,:,:)-u(4,:,:))*(Ax_un(2,:,:) - Ay_un(1,:,:))
     $     - delVpar*Bsq_u(1,:,:))*invBsq
      delVpar_uxn=(-(u(11,:,:)-u(3,:,:))*A_un(3,:,:)
     $     + (u(12,:,:)-u(4,:,:))*A_un(2,:,:) 
     $     - delVpar*Bsq_ux(1,:,:))*invBsq
      delVpar_uyn=((u(10,:,:)-u(2,:,:))*A_un(3,:,:)
     $     - (u(12,:,:)-u(4,:,:))*A_un(1,:,:)
     $     - delVpar*Bsq_uy(1,:,:))*invBsq
c-----------------------------------------------------------------------
c     density.
c-----------------------------------------------------------------------
      fx_u(1,2,:,:)=one
      fy_u(1,3,:,:)=one
c-----------------------------------------------------------------------
c     x-momentum.
c-----------------------------------------------------------------------
      fx_u(2,1,:,:)=u(2,:,:)*vi_un(1,:,:)
     $     + mass_r*u(10,:,:)*ve_un(1,:,:)
     $     + eta_i0*kappa_i*(one+third)*Ti**1.5
     $     *(2.5_r8*Ti_un*vix(1,:,:) + Ti*vix_un(1,:,:))
     $     + eta_e0*kappa_e*(one+third)*Te**1.5
     $     *(2.5_r8*Te_un*vex(1,:,:) + Te*vex_un(1,:,:))
      fx_u(2,2,:,:)=vi(1,:,:) + u(2,:,:)*n_inv
     $     + eta_i0*kappa_i*(one+third)*Ti**2.5*nx_inv
      fx_u(2,8,:,:)=one + 2.5_r8*eta_i0*kappa_i*(one+third)
     $     *Ti**1.5*vix(1,:,:)*n_inv
      fx_u(2,9,:,:)=one + 2.5_r8*eta_e0*kappa_e*(one+third)
     $     *Te**1.5*vex(1,:,:)*n_inv
      fx_u(2,10,:,:)=mass_r*ve(1,:,:) + mass_r*u(10,:,:)*n_inv
     $     + eta_e0*kappa_e*(one+third)*Te**2.5*nx_inv

      fx_ux(2,1,:,:)=eta_i0*kappa_i*Ti**2.5*(one+third)*vi_un(1,:,:)
     $     + eta_e0*kappa_e*Te**2.5*(one+third)*ve_un(1,:,:)
      fx_ux(2,2,:,:)=eta_i0*kappa_i*Ti**2.5*(one+third)*n_inv
      fx_ux(2,10,:,:)=eta_e0*kappa_e*Te**2.5*(one+third)*n_inv

      fy_u(2,1,:,:)=u(3,:,:)*vi_un(1,:,:) 
     $     + mass_r*u(11,:,:)*ve_un(1,:,:)
     $     + eta_i0*kappa_i*Ti**1.5
     $     *(2.5_r8*Ti_un*(viy(1,:,:)+third*vix(2,:,:))
     $     + Ti*(viy_un(1,:,:)+third*vix_un(2,:,:)))
     $     + eta_e0*kappa_e*Te**1.5
     $     *(2.5_r8*Te_un*(vey(1,:,:)+third*vex(2,:,:))
     $     + Te*(vey_un(1,:,:)+third*vex_un(2,:,:)))
      fy_u(2,2,:,:)=u(3,:,:)*n_inv + eta_i0*kappa_i*Ti**2.5*ny_inv
      fy_u(2,3,:,:)=vi(1,:,:) + eta_i0*kappa_i*Ti**2.5*third*nx_inv
      fy_u(2,8,:,:)=2.5_r8*eta_i0*kappa_i*Ti**1.5*n_inv
     $     *(viy(1,:,:)+third*vix(2,:,:))
      fy_u(2,9,:,:)=2.5_r8*eta_e0*kappa_e*Te**1.5*n_inv
     $     *(vey(1,:,:)+third*vex(2,:,:))
      fy_u(2,10,:,:)=mass_r*u(11,:,:)*n_inv
     $     + eta_e0*kappa_e*Te**2.5*ny_inv
      fy_u(2,11,:,:)=mass_r*ve(1,:,:)
     $     + eta_e0*kappa_e*Te**2.5*third*nx_inv

      fy_ux(2,1,:,:)=eta_i0*kappa_i*Ti**2.5*third*vi_un(2,:,:)
     $     + eta_e0*kappa_e*Te**2.5*third*ve_un(2,:,:)
      fy_ux(2,3,:,:)=eta_i0*kappa_i*Ti**2.5*third*n_inv
      fy_ux(2,11,:,:)=eta_e0*kappa_e*Te**2.5*third*n_inv

      fy_uy(2,1,:,:)=eta_i0*kappa_i*Ti**2.5*vi_un(1,:,:)
     $     + eta_e0*kappa_e*Te**2.5*ve_un(1,:,:)
      fy_uy(2,2,:,:)=eta_i0*kappa_i*Ti**2.5*n_inv
      fy_uy(2,10,:,:)=eta_e0*kappa_e*Te**2.5*n_inv

      s_u(2,1,:,:)=((u(3,:,:)-u(11,:,:))*(Ax_un(2,:,:)-Ay_un(1,:,:))
     $     + (u(4,:,:)-u(12,:,:))*Ax_un(3,:,:))*d_inv
      s_u(2,3,:,:)=(Ax(2,:,:)-Ay(1,:,:))*d_inv
      s_u(2,4,:,:)=Ax(3,:,:)*d_inv
      s_u(2,5,:,:)=-(u(3,:,:)-u(11,:,:))*ny_inv*d_inv
      s_u(2,6,:,:)=(u(3,:,:)-u(11,:,:))*nx_inv*d_inv
      s_u(2,7,:,:)=(u(4,:,:)-u(12,:,:))*nx_inv*d_inv
      s_u(2,11,:,:)=-(Ax(2,:,:)-Ay(1,:,:))*d_inv
      s_u(2,12,:,:)=-Ax(3,:,:)*d_inv

      s_ux(2,1,:,:)=((u(3,:,:)-u(11,:,:))*A_un(2,:,:)
     $     + (u(4,:,:)-u(12,:,:))*A_un(3,:,:))*d_inv
      s_ux(2,6,:,:)=(u(3,:,:)-u(11,:,:))*n_inv*d_inv
      s_ux(2,7,:,:)=(u(4,:,:)-u(12,:,:))*n_inv*d_inv

      s_uy(2,1,:,:)=-(u(3,:,:)-u(11,:,:))*A_un(1,:,:)*d_inv
      s_uy(2,5,:,:)=-(u(3,:,:)-u(11,:,:))*n_inv*d_inv
c-----------------------------------------------------------------------
c     y-momentum.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=u(2,:,:)*vi_un(2,:,:) 
     $     + mass_r*u(10,:,:)*ve_un(2,:,:)
     $     + eta_i0*kappa_i*Ti**1.5
     $     *(2.5_r8*Ti_un*(vix(2,:,:)+third*viy(1,:,:))
     $     + Ti*(vix_un(2,:,:)+third*viy_un(1,:,:)))
     $     + eta_e0*kappa_e*Te**1.5
     $     *(2.5_r8*Te_un*(vex(2,:,:)+third*vey(1,:,:))
     $     + Te*(vex_un(2,:,:)+third*vey_un(1,:,:)))
      fx_u(3,2,:,:)=vi(2,:,:)
     $     + eta_i0*kappa_i*Ti**2.5*third*ny_inv
      fx_u(3,3,:,:)=u(2,:,:)*n_inv
     $     + eta_i0*kappa_i*Ti**2.5*nx_inv
      fx_u(3,8,:,:)=2.5_r8*eta_i0*kappa_i*Ti**1.5*n_inv
     $     *(vix(2,:,:)+third*viy(1,:,:))
      fx_u(3,9,:,:)=2.5_r8*eta_e0*kappa_e*Te**1.5*n_inv
     $     *(vex(2,:,:)+third*vey(1,:,:))
      fx_u(3,10,:,:)=mass_r*ve(2,:,:)
     $     + eta_e0*kappa_e*Te**2.5*third*ny_inv
      fx_u(3,11,:,:)=mass_r*u(10,:,:)*n_inv
     $     + eta_e0*kappa_e*Te**2.5*nx_inv

      fx_ux(3,1,:,:)=eta_i0*kappa_i*Ti**2.5*vi_un(2,:,:)
     $     + eta_e0*kappa_e*Te**2.5*ve_un(2,:,:)
      fx_ux(3,3,:,:)=eta_i0*kappa_i*Ti**2.5*n_inv
      fx_ux(3,11,:,:)=eta_e0*kappa_e*Te**2.5*n_inv

      fx_uy(3,1,:,:)=eta_i0*kappa_i*Ti**2.5*third*vi_un(1,:,:)
     $     + eta_e0*kappa_e*Te**2.5*third*ve_un(1,:,:)
      fx_uy(3,2,:,:)=eta_i0*kappa_i*Ti**2.5*third*n_inv
      fx_uy(3,10,:,:)=eta_e0*kappa_e*Te**2.5*third*n_inv

      fy_u(3,1,:,:)=u(3,:,:)*vi_un(2,:,:)
     $     + mass_r*u(11,:,:)*ve_un(2,:,:)
     $     + eta_i0*kappa_i*(one+third)*Ti**1.5
     $     *(2.5_r8*Ti_un*viy(2,:,:) + Ti*viy_un(2,:,:))
     $     + eta_e0*kappa_e*(one+third)*Te**1.5
     $     *(2.5_r8*Te_un*vey(2,:,:) + Te*vey_un(2,:,:))
      fy_u(3,3,:,:)=vi(2,:,:) + u(3,:,:)*n_inv
     $     + eta_i0*kappa_i*(one+third)*Ti**2.5*ny_inv
      fy_u(3,8,:,:)=one+2.5_r8*eta_i0*kappa_i*(one+third)
     $     *Ti**1.5*viy(2,:,:)*n_inv
      fy_u(3,9,:,:)=one+2.5_r8*eta_e0*kappa_e*(one+third)
     $     *Te**1.5*vey(2,:,:)*n_inv
      fy_u(3,11,:,:)=mass_r*ve(2,:,:) + mass_r*u(11,:,:)*n_inv
     $     + eta_e0*kappa_e*(one+third)*Te**2.5*ny_inv

      fy_uy(3,1,:,:)=eta_i0*kappa_i*Ti**2.5*(one+third)*vi_un(2,:,:)
     $     + eta_e0*kappa_e*Te**2.5*(one+third)*ve_un(2,:,:)
      fy_uy(3,3,:,:)=eta_i0*kappa_i*Ti**2.5*(one+third)*n_inv
      fy_uy(3,11,:,:)=eta_e0*kappa_e*Te**2.5*(one+third)*n_inv

      s_u(3,1,:,:)=((u(2,:,:)-u(10,:,:))*(Ay_un(1,:,:)-Ax_un(2,:,:))
     $     + (u(4,:,:)-u(12,:,:))*Ay_un(3,:,:))*d_inv
      s_u(3,2,:,:)=(Ay(1,:,:)-Ax(2,:,:))*d_inv
      s_u(3,4,:,:)=Ay(3,:,:)*d_inv
      s_u(3,5,:,:)=(u(2,:,:)-u(10,:,:))*ny_inv*d_inv
      s_u(3,6,:,:)=-(u(2,:,:)-u(10,:,:))*nx_inv*d_inv
      s_u(3,7,:,:)=(u(4,:,:)-u(12,:,:))*ny_inv*d_inv
      s_u(3,10,:,:)=-(Ay(1,:,:)-Ax(2,:,:))*d_inv
      s_u(3,12,:,:)=-Ay(3,:,:)*d_inv

      s_ux(3,1,:,:)=-(u(2,:,:)-u(10,:,:))*A_un(2,:,:)*d_inv
      s_ux(3,6,:,:)=-(u(2,:,:)-u(10,:,:))*n_inv*d_inv

      s_uy(3,1,:,:)=((u(2,:,:)-u(10,:,:))*A_un(1,:,:)
     $     + (u(4,:,:)-u(12,:,:))*A_un(3,:,:))*d_inv
      s_uy(3,5,:,:)=(u(2,:,:)-u(10,:,:))*n_inv*d_inv
      s_uy(3,7,:,:)=(u(4,:,:)-u(12,:,:))*n_inv*d_inv
c-----------------------------------------------------------------------
c     z-momentum.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:)=u(2,:,:)*vi_un(3,:,:)
     $     + mass_r*u(10,:,:)*ve_un(3,:,:)
     $     + eta_i0*kappa_i*Ti**1.5
     $     *(2.5_r8*Ti_un*vix(3,:,:) + Ti*vix_un(3,:,:))
     $     + eta_e0*kappa_e*Te**1.5
     $     *(2.5_r8*Te_un*vex(3,:,:) + Te*vex_un(3,:,:))
      fx_u(4,2,:,:)=vi(3,:,:)
      fx_u(4,4,:,:)=u(2,:,:)*n_inv
     $     + eta_i0*kappa_i*Ti**2.5*nx_inv
      fx_u(4,8,:,:)=2.5_r8*eta_i0*kappa_i
     $     *Ti**1.5*vix(3,:,:)*n_inv
      fx_u(4,9,:,:)=2.5_r8*eta_e0*kappa_e
     $     *Te**1.5*vex(3,:,:)*n_inv
      fx_u(4,10,:,:)=mass_r*ve(3,:,:)
      fx_u(4,12,:,:)=mass_r*u(10,:,:)*n_inv 
     $     + eta_e0*kappa_e*Te**2.5*nx_inv

      fx_ux(4,1,:,:)=eta_i0*kappa_i*Ti**2.5*vi_un(3,:,:)
     $     + eta_e0*kappa_e*Te**2.5*ve_un(3,:,:)
      fx_ux(4,4,:,:)=eta_i0*kappa_i*Ti**2.5*n_inv
      fx_ux(4,12,:,:)=eta_e0*kappa_e*Te**2.5*n_inv

      fy_u(4,1,:,:)=u(3,:,:)*vi_un(3,:,:)
     $     + mass_r*u(11,:,:)*ve_un(3,:,:)
     $     + eta_i0*kappa_i*Ti**1.5
     $     *(2.5_r8*Ti_un*viy(3,:,:) + Ti*viy_un(3,:,:))
     $     + eta_e0*kappa_e*Te**1.5
     $     *(2.5_r8*Te_un*vey(3,:,:) + Te*vey_un(3,:,:))
      fy_u(4,3,:,:)=vi(3,:,:)
      fy_u(4,4,:,:)=u(3,:,:)*n_inv + eta_i0*kappa_i*Ti**2.5*ny_inv
      fy_u(4,8,:,:)=2.5_r8*eta_i0*kappa_i*Ti**1.5*viy(3,:,:)*n_inv
      fy_u(4,9,:,:)=2.5_r8*eta_e0*kappa_e*Te**1.5*vey(3,:,:)*n_inv
      fy_u(4,11,:,:)=mass_r*ve(3,:,:)
      fy_u(4,12,:,:)=mass_r*u(11,:,:)*n_inv
     $     + eta_e0*kappa_e*Te**2.5*ny_inv

      fy_uy(4,1,:,:)=eta_i0*kappa_i*Ti**2.5*vi_un(3,:,:)
     $     + eta_e0*kappa_e*Te**2.5*ve_un(3,:,:)
      fy_uy(4,4,:,:)=eta_i0*kappa_i*Ti**2.5*n_inv
      fy_uy(4,12,:,:)=eta_e0*kappa_e*Te**2.5*n_inv

      s_u(4,1,:,:)=-((u(2,:,:)-u(10,:,:))*Ax_un(3,:,:)
     $     + (u(3,:,:)-u(11,:,:))*Ay_un(3,:,:))*d_inv
      s_u(4,2,:,:)=-Ax(3,:,:)*d_inv
      s_u(4,3,:,:)=-Ay(3,:,:)*d_inv
      s_u(4,7,:,:)=-((u(2,:,:)-u(10,:,:))*nx_inv
     $     + (u(3,:,:)-u(11,:,:))*ny_inv)*d_inv
      s_u(4,10,:,:)=Ax(3,:,:)*d_inv
      s_u(4,11,:,:)=Ay(3,:,:)*d_inv

      s_ux(4,1,:,:)=-(u(2,:,:)-u(10,:,:))*A_un(3,:,:)*d_inv
      s_ux(4,7,:,:)=-(u(2,:,:)-u(10,:,:))*n_inv*d_inv

      s_uy(4,1,:,:)=-(u(3,:,:)-u(11,:,:))*A_un(3,:,:)*d_inv
      s_uy(4,7,:,:)=-(u(3,:,:)-u(11,:,:))*n_inv*d_inv
c-----------------------------------------------------------------------
c     x-Ohm's law.
c-----------------------------------------------------------------------
      fx_u(5,1,:,:)=-u(10,:,:)*(mass_r*ve_un(1,:,:)-d_inv*A_un(1,:,:))
     $     - eta_e0*kappa_e*(one+third)*Te**1.5
     $     *(2.5_r8*Te_un*vex(1,:,:) + Te*vex_un(1,:,:))
      fx_u(5,5,:,:)=d_inv*u(10,:,:)*n_inv
      fx_u(5,9,:,:)=-one-2.5_r8*eta_e0*kappa_e*(one+third)
     $     *Te**1.5*vex(1,:,:)*n_inv
      fx_u(5,10,:,:)=-mass_r*ve(1,:,:)+d_inv*A(1,:,:)
     $     - mass_r*u(10,:,:)*n_inv
     $     - eta_e0*kappa_e*(one+third)*Te**2.5*nx_inv

      fx_ux(5,1,:,:)=-eta_e0*kappa_e*Te**2.5*(one+third)*ve_un(1,:,:)
      fx_ux(5,10,:,:)=-eta_e0*kappa_e*Te**2.5*(one+third)*n_inv

      fy_u(5,1,:,:)=-u(11,:,:)*(mass_r*ve_un(1,:,:)-d_inv*A_un(1,:,:))
     $     - eta_e0*kappa_e*Te**1.5
     $     *(2.5_r8*Te_un*(vey(1,:,:)+third*vex(2,:,:))
     $     + Te*(vey_un(1,:,:)+third*vex_un(2,:,:)))
      fy_u(5,5,:,:)=d_inv*u(11,:,:)*n_inv
      fy_u(5,9,:,:)=-2.5_r8*eta_e0*kappa_e*Te**1.5*n_inv
     $     *(vey(1,:,:)+third*vex(2,:,:))
      fy_u(5,10,:,:)=-mass_r*u(11,:,:)*n_inv
     $     - eta_e0*kappa_e*Te**2.5*ny_inv
      fy_u(5,11,:,:)=-mass_r*ve(1,:,:)+d_inv*A(1,:,:)
     $     - eta_e0*kappa_e*Te**2.5*third*nx_inv

      fy_ux(5,1,:,:)=-eta_e0*kappa_e*Te**2.5*third*ve_un(2,:,:)
      fy_ux(5,11,:,:)=-eta_e0*kappa_e*Te**2.5*third*n_inv

      fy_uy(5,1,:,:)=-eta_e0*kappa_e*Te**2.5*ve_un(1,:,:)
      fy_uy(5,10,:,:)=-eta_e0*kappa_e*Te**2.5*n_inv

      s_u(5,1,:,:)=(u(10,:,:)*Ax_un(1,:,:)+u(11,:,:)*Ax_un(2,:,:)
     $     + u(12,:,:)*Ax_un(3,:,:))*d_inv
     $     - 2.5_r8*resist*(u(2,:,:)-u(10,:,:))/(kappa_e*Te**1.5)
     $     + fheat*(BdotTe(1,:,:) + u(1,:,:)*BdotTe_un(1,:,:))
      s_u(5,2,:,:)=-resist*u(1,:,:)/(kappa_e*Te**1.5)
      s_u(5,5,:,:)=u(10,:,:)*nx_inv*d_inv
     $     - fheat*u(1,:,:)*BdotTe_bz(1,:,:)*ny_inv
      s_u(5,6,:,:)=u(11,:,:)*nx_inv*d_inv
     $     + fheat*u(1,:,:)*BdotTe_bz(1,:,:)*nx_inv
      s_u(5,7,:,:)=u(12,:,:)*nx_inv*d_inv + fheat*u(1,:,:)
     $     *(BdotTe_bx(1,:,:)*ny_inv - BdotTe_by(1,:,:)*nx_inv)
      s_u(5,9,:,:)=1.5_r8*resist*(u(2,:,:)-u(10,:,:))/(kappa_e*Te**2.5)
     $     + fheat*u(1,:,:)
     $     *(BdotTe_Tex(1,:,:)*nx_inv + BdotTe_Tey(1,:,:)*ny_inv)
      s_u(5,10,:,:)=Ax(1,:,:)*d_inv + resist*u(1,:,:)/(kappa_e*Te**1.5)
      s_u(5,11,:,:)=Ax(2,:,:)*d_inv
      s_u(5,12,:,:)=Ax(3,:,:)*d_inv

      s_ux(5,1,:,:)=(u(10,:,:)*A_un(1,:,:)+u(11,:,:)*A_un(2,:,:)
     $     + u(12,:,:)*A_un(3,:,:))*d_inv
     $     + fheat*u(1,:,:)*BdotTe_uxn(1,:,:)
      s_ux(5,5,:,:)=u(10,:,:)*n_inv*d_inv
      s_ux(5,6,:,:)=u(11,:,:)*n_inv*d_inv + fheat*BdotTe_bz(1,:,:)
      s_ux(5,7,:,:)=u(12,:,:)*n_inv*d_inv - fheat*BdotTe_by(1,:,:)
      s_ux(5,9,:,:)=fheat*BdotTe_Tex(1,:,:)

      s_uy(5,1,:,:)=fheat*u(1,:,:)*BdotTe_uyn(1,:,:)
      s_uy(5,5,:,:)=-fheat*BdotTe_bz(1,:,:)
      s_uy(5,7,:,:)=fheat*BdotTe_bx(1,:,:)
      s_uy(5,9,:,:)=fheat*BdotTe_Tey(1,:,:)
c-----------------------------------------------------------------------
c     y-Ohm's law.
c-----------------------------------------------------------------------
      fx_u(6,1,:,:)=-u(10,:,:)*(mass_r*ve_un(2,:,:)-d_inv*A_un(2,:,:))
     $     - eta_e0*kappa_e*Te**1.5
     $     *(2.5_r8*Te_un*(vex(2,:,:)+third*vey(1,:,:))
     $     + Te*(vex_un(2,:,:)+third*vey_un(1,:,:)))
      fx_u(6,6,:,:)=d_inv*u(10,:,:)*n_inv
      fx_u(6,9,:,:)=-2.5_r8*eta_e0*kappa_e*Te**1.5*n_inv
     $     *(vex(2,:,:)+third*vey(1,:,:))
      fx_u(6,10,:,:)=-mass_r*ve(2,:,:)+d_inv*A(2,:,:)
     $     - eta_e0*kappa_e*Te**2.5*third*ny_inv
      fx_u(6,11,:,:)=-mass_r*u(10,:,:)*n_inv
     $     - eta_e0*kappa_e*Te**2.5*nx_inv

      fx_ux(6,1,:,:)=-eta_e0*kappa_e*Te**2.5*ve_un(2,:,:)
      fx_ux(6,11,:,:)=-eta_e0*kappa_e*Te**2.5*n_inv

      fx_uy(6,1,:,:)=-eta_e0*kappa_e*Te**2.5*third*ve_un(1,:,:)
      fx_uy(6,10,:,:)=-eta_e0*kappa_e*Te**2.5*third*n_inv

      fy_u(6,1,:,:)=-u(11,:,:)*(mass_r*ve_un(2,:,:)-d_inv*A_un(2,:,:))
     $     - eta_e0*kappa_e*(one+third)*Te**1.5
     $     *(2.5_r8*Te_un*vey(2,:,:) + Te*vey_un(2,:,:))
      fy_u(6,6,:,:)=d_inv*u(11,:,:)*n_inv
      fy_u(6,9,:,:)=-one-2.5_r8*eta_e0*kappa_e*(one+third)
     $     *Te**1.5*vey(2,:,:)*n_inv
      fy_u(6,11,:,:)=-mass_r*ve(2,:,:)+d_inv*A(2,:,:)
     $     - mass_r*u(11,:,:)*n_inv
     $     - eta_e0*kappa_e*(one+third)*Te**2.5*ny_inv

      fy_uy(6,1,:,:)=-eta_e0*kappa_e*Te**2.5*(one+third)*ve_un(2,:,:)
      fy_uy(6,11,:,:)=-eta_e0*kappa_e*Te**2.5*(one+third)*n_inv

      s_u(6,1,:,:)=(u(10,:,:)*Ay_un(1,:,:)+u(11,:,:)*Ay_un(2,:,:)
     $     + u(12,:,:)*Ay_un(3,:,:))*d_inv
     $     - 2.5_r8*resist*(u(3,:,:)-u(11,:,:))/(kappa_e*Te**1.5)
     $     + fheat*(BdotTe(2,:,:) + u(1,:,:)*BdotTe_un(2,:,:))
      s_u(6,3,:,:)=-resist*u(1,:,:)/(kappa_e*Te**1.5)
      s_u(6,5,:,:)=u(10,:,:)*ny_inv*d_inv
     $     - fheat*u(1,:,:)*BdotTe_bz(2,:,:)*ny_inv
      s_u(6,6,:,:)=u(11,:,:)*ny_inv*d_inv
     $     + fheat*u(1,:,:)*BdotTe_bz(2,:,:)*nx_inv
      s_u(6,7,:,:)=u(12,:,:)*ny_inv*d_inv + fheat*u(1,:,:)
     $     *(BdotTe_bx(2,:,:)*ny_inv - BdotTe_by(2,:,:)*nx_inv)
      s_u(6,9,:,:)=1.5_r8*resist*(u(3,:,:)-u(11,:,:))/(kappa_e*Te**2.5)
     $     + fheat*u(1,:,:)
     $     *(BdotTe_Tex(2,:,:)*nx_inv + BdotTe_Tey(2,:,:)*ny_inv)
      s_u(6,10,:,:)=Ay(1,:,:)*d_inv
      s_u(6,11,:,:)=Ay(2,:,:)*d_inv + resist*u(1,:,:)/(kappa_e*Te**1.5)
      s_u(6,12,:,:)=Ay(3,:,:)*d_inv

      s_ux(6,1,:,:)=fheat*u(1,:,:)*BdotTe_uxn(2,:,:)
      s_ux(6,6,:,:)=fheat*BdotTe_bz(2,:,:)
      s_ux(6,7,:,:)=-fheat*BdotTe_by(2,:,:)
      s_ux(6,9,:,:)=fheat*BdotTe_Tex(2,:,:)

      s_uy(6,1,:,:)=(u(10,:,:)*A_un(1,:,:)+u(11,:,:)*A_un(2,:,:)
     $     + u(12,:,:)*A_un(3,:,:))*d_inv
     $     + fheat*u(1,:,:)*BdotTe_uyn(2,:,:)
      s_uy(6,5,:,:)=u(10,:,:)*n_inv*d_inv - fheat*BdotTe_bz(2,:,:)
      s_uy(6,6,:,:)=u(11,:,:)*n_inv*d_inv
      s_uy(6,7,:,:)=u(12,:,:)*n_inv*d_inv + fheat*BdotTe_bx(2,:,:)
      s_uy(6,9,:,:)=fheat*BdotTe_Tey(2,:,:)
c-----------------------------------------------------------------------
c     z-Ohm's law.
c-----------------------------------------------------------------------
      fx_u(7,1,:,:)=-u(10,:,:)*(mass_r*ve_un(3,:,:)-d_inv*A_un(3,:,:))
     $     - eta_e0*kappa_e*Te**1.5
     $     *(2.5_r8*Te_un*vex(3,:,:)+Te*vex_un(3,:,:))
      fx_u(7,7,:,:)=d_inv*u(10,:,:)*n_inv
      fx_u(7,9,:,:)=-2.5_r8*eta_e0*kappa_e*Te**1.5*n_inv
     $     *vex(3,:,:)
      fx_u(7,10,:,:)=-mass_r*ve(3,:,:)+d_inv*A(3,:,:)
      fx_u(7,12,:,:)=-mass_r*u(10,:,:)*n_inv
     $     - eta_e0*kappa_e*Te**2.5*nx_inv

      fx_ux(7,1,:,:)=-eta_e0*kappa_e*Te**2.5*ve_un(3,:,:)
      fx_ux(7,12,:,:)=-eta_e0*kappa_e*Te**2.5*n_inv

      fy_u(7,1,:,:)=-u(11,:,:)*(mass_r*ve_un(3,:,:)-d_inv*A_un(3,:,:))
     $     - eta_e0*kappa_e*Te**1.5
     $     *(2.5_r8*Te_un*vey(3,:,:)+Te*vey_un(3,:,:))
      fy_u(7,7,:,:)=d_inv*u(11,:,:)*n_inv
      fy_u(7,9,:,:)=-2.5_r8*eta_e0*kappa_e*Te**1.5*n_inv
     $     *vey(3,:,:)
      fy_u(7,11,:,:)=-mass_r*ve(3,:,:)+d_inv*A(3,:,:)
      fy_u(7,12,:,:)=-mass_r*u(11,:,:)*n_inv
     $     - eta_e0*kappa_e*Te**2.5*ny_inv

      fy_uy(7,1,:,:)=-eta_e0*kappa_e*Te**2.5*ve_un(3,:,:)
      fy_uy(7,12,:,:)=-eta_e0*kappa_e*Te**2.5*n_inv

      s_u(7,1,:,:)=-2.5_r8*resist*(u(4,:,:)-u(12,:,:))/(kappa_e*Te**1.5)
     $     + fheat*(BdotTe(3,:,:) + u(1,:,:)*BdotTe_un(3,:,:))
      s_u(7,4,:,:)=-resist*u(1,:,:)/(kappa_e*Te**1.5)
      s_u(7,5,:,:)=-fheat*u(1,:,:)*BdotTe_bz(3,:,:)*ny_inv
      s_u(7,6,:,:)=fheat*u(1,:,:)*BdotTe_bz(3,:,:)*nx_inv
      s_u(7,7,:,:)=fheat*u(1,:,:)
     $     *(BdotTe_bx(3,:,:)*ny_inv - BdotTe_by(3,:,:)*nx_inv)
      s_u(7,9,:,:)=1.5_r8*resist*(u(4,:,:)-u(12,:,:))/(kappa_e*Te**2.5)
     $     + fheat*u(1,:,:)
     $     *(BdotTe_Tex(3,:,:)*nx_inv + BdotTe_Tey(3,:,:)*ny_inv)
      s_u(7,11,:,:)=-d_inv*B_eq
      s_u(7,12,:,:)=resist*u(1,:,:)/(kappa_e*Te**1.5)

      s_ux(7,1,:,:)=fheat*u(1,:,:)*BdotTe_uxn(3,:,:)
      s_ux(7,6,:,:)=fheat*BdotTe_bz(3,:,:)
      s_ux(7,7,:,:)=-fheat*BdotTe_by(3,:,:)
      s_ux(7,9,:,:)=fheat*BdotTe_Tex(3,:,:)

      s_uy(7,1,:,:)=fheat*u(1,:,:)*BdotTe_uyn(3,:,:)
      s_uy(7,5,:,:)=-fheat*BdotTe_bz(3,:,:)
      s_uy(7,7,:,:)=fheat*BdotTe_bx(3,:,:)
      s_uy(7,9,:,:)=fheat*BdotTe_Tey(3,:,:)
c-----------------------------------------------------------------------
c     ion pressure.
c-----------------------------------------------------------------------
      fx_u(8,1,:,:) = gfac*u(8,:,:)*vi_un(1,:,:)
     $     - kfaci*BdotTi_un(1,:,:) - kperpi*Tix_un
     $     - (kpari_un - kperpi_un)*BdotTi(1,:,:) - kperpi_un*Tix
     $     + kperpi_bsq*Bsq_u(1,:,:)*(BdotTi(1,:,:) - Tix)
      fx_u(8,2,:,:) = gfac*Ti
      fx_u(8,5,:,:) = kfaci*BdotTi_bz(1,:,:)*ny_inv
     $     + kperpi_bsq*Bsq_u(2,:,:)*(BdotTi(1,:,:) - Tix)
      fx_u(8,6,:,:) = -kfaci*BdotTi_bz(1,:,:)*nx_inv    
     $     + kperpi_bsq*Bsq_u(3,:,:)*(BdotTi(1,:,:) - Tix)
      fx_u(8,7,:,:) = -kfaci
     $     *(BdotTi_bx(1,:,:)*ny_inv - BdotTi_by(1,:,:)*nx_inv)
     $     + kperpi_bsq*Bsq_u(4,:,:)*(BdotTi(1,:,:) - Tix)
      fx_u(8,8,:,:) = gfac*vi(1,:,:) - kfaci
     $     *(BdotTi_Tix(1,:,:)*nx_inv + BdotTi_Tiy(1,:,:)*ny_inv)
     $     - kperpi*nx_inv
     $     - (kpari_pi - kperpi_pi)*BdotTi(1,:,:) - kperpi_pi*Tix

      fx_ux(8,1,:,:) = -kfaci*BdotTi_uxn(1,:,:) - kperpi*Ti_un
     $     + kperpi_bsq*Bsq_ux(1,:,:)*(BdotTi(1,:,:) - Tix)
      fx_ux(8,6,:,:) = -kfaci*BdotTi_bz(1,:,:)*n_inv
     $     + kperpi_bsq*Bsq_ux(3,:,:)*(BdotTi(1,:,:) - Tix)
      fx_ux(8,7,:,:) = kfaci*BdotTi_by(1,:,:)*n_inv
     $     + kperpi_bsq*Bsq_ux(4,:,:)*(BdotTi(1,:,:) - Tix)
      fx_ux(8,8,:,:) = -kfaci*BdotTi_Tix(1,:,:)*n_inv - kperpi*n_inv

      fx_uy(8,1,:,:) = -kfaci*BdotTi_uyn(1,:,:)
     $     + kperpi_bsq*Bsq_uy(1,:,:)*(BdotTi(1,:,:) - Tix)
      fx_uy(8,5,:,:) = kfaci*BdotTi_bz(1,:,:)*n_inv
     $     + kperpi_bsq*Bsq_uy(2,:,:)*(BdotTi(1,:,:) - Tix)
      fx_uy(8,7,:,:) = -kfaci*BdotTi_bx(1,:,:)*n_inv
     $     + kperpi_bsq*Bsq_uy(4,:,:)*(BdotTi(1,:,:) - Tix)
      fx_uy(8,8,:,:) = -kfaci*BdotTi_Tiy(1,:,:)*n_inv

      fy_u(8,1,:,:) = gfac*u(8,:,:)*vi_un(2,:,:)
     $     - kfaci*BdotTi_un(2,:,:) - kperpi*Tiy_un
     $     - (kpari_un - kperpi_un)*BdotTi(2,:,:) - kperpi_un*Tiy
     $     + kperpi_bsq*Bsq_u(1,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_u(8,3,:,:) = gfac*Ti
      fy_u(8,5,:,:) = kfaci*BdotTi_bz(2,:,:)*ny_inv
     $     + kperpi_bsq*Bsq_u(2,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_u(8,6,:,:) = -kfaci*BdotTi_bz(2,:,:)*nx_inv    
     $     + kperpi_bsq*Bsq_u(3,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_u(8,7,:,:) = -kfaci
     $     *(BdotTi_bx(2,:,:)*ny_inv - BdotTi_by(2,:,:)*nx_inv)
     $     + kperpi_bsq*Bsq_u(4,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_u(8,8,:,:) = gfac*vi(2,:,:) - kfaci
     $     *(BdotTi_Tix(2,:,:)*nx_inv + BdotTi_Tiy(2,:,:)*ny_inv)
     $     - kperpi*ny_inv
     $     - (kpari_pi - kperpi_pi)*BdotTi(2,:,:) - kperpi_pi*Tiy

      fy_ux(8,1,:,:) = -kfaci*BdotTi_uxn(2,:,:)
     $     + kperpi_bsq*Bsq_ux(1,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_ux(8,6,:,:) = -kfaci*BdotTi_bz(2,:,:)*n_inv
     $     + kperpi_bsq*Bsq_ux(3,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_ux(8,7,:,:) = kfaci*BdotTi_by(2,:,:)*n_inv
     $     + kperpi_bsq*Bsq_ux(4,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_ux(8,8,:,:) = -kfaci*BdotTi_Tix(2,:,:)*n_inv

      fy_uy(8,1,:,:) = -kfaci*BdotTi_uyn(2,:,:) - kperpi*Ti_un
     $     + kperpi_bsq*Bsq_uy(1,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_uy(8,5,:,:) = kfaci*BdotTi_bz(2,:,:)*n_inv
     $     + kperpi_bsq*Bsq_uy(2,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_uy(8,7,:,:) = -kfaci*BdotTi_bx(2,:,:)*n_inv
     $     + kperpi_bsq*Bsq_uy(4,:,:)*(BdotTi(2,:,:) - Tiy)
      fy_uy(8,8,:,:) = -kfaci*BdotTi_Tiy(2,:,:)*n_inv - kperpi*n_inv

      s_u(8,1,:,:)=vi_un(1,:,:)*ux(8,:,:) + vi_un(2,:,:)*uy(8,:,:)
     $     - eta_i0*kappa_i*Ti**1.5*(2.5_r8*Ti_un*PidelV
     $     + two*Ti*((one+third)*(vix(1,:,:)*vix_un(1,:,:)
     $     + viy(2,:,:)*viy_un(2,:,:))
     $     + third*vix(2,:,:)*viy_un(1,:,:)
     $     + third*viy(1,:,:)*vix_un(2,:,:)
     $     + viy(1,:,:)*viy_un(1,:,:) + vix(2,:,:)*vix_un(2,:,:)
     $     + vix(3,:,:)*vix_un(3,:,:) + viy(3,:,:)*viy_un(3,:,:)))
     $     + 2.5_r8*ieheat*(u(9,:,:)-u(8,:,:))/(kappa_e*Te**1.5)
      s_u(8,2,:,:)=n_inv*ux(8,:,:) - two*eta_i0*kappa_i*Ti**2.5
     $     *((one+third)*vix(1,:,:)*nx_inv
     $     + third*vix(2,:,:)*ny_inv + viy(1,:,:)*ny_inv)
      s_u(8,3,:,:)=n_inv*uy(8,:,:) - two*eta_i0*kappa_i*Ti**2.5
     $     *((one+third)*viy(2,:,:)*ny_inv
     $     + third*viy(1,:,:)*nx_inv + vix(2,:,:)*nx_inv)
      s_u(8,4,:,:)=-two*eta_i0*kappa_i*Ti**2.5
     $     *(vix(3,:,:)*nx_inv + viy(3,:,:)*ny_inv)
      s_u(8,8,:,:)=-eta_i0*kappa_i*2.5_r8*Ti**1.5*n_inv*PidelV
     $     - ieheat*u(1,:,:)/(kappa_e*Te**1.5)
      s_u(8,9,:,:)=ieheat/(kappa_e*Te**2.5)
     $     *(1.5_r8*u(8,:,:)-0.5_r8*u(9,:,:))

      s_ux(8,1,:,:)=-two*eta_i0*kappa_i*Ti**2.5
     $     *((one+third)*vix(1,:,:)*vi_un(1,:,:)
     $     +third*viy(1,:,:)*vi_un(2,:,:)
     $     +vix(2,:,:)*vi_un(2,:,:) + vix(3,:,:)*vi_un(3,:,:))
      s_ux(8,2,:,:)=-two*eta_i0*kappa_i*Ti**2.5
     $     *(one+third)*vix(1,:,:)*n_inv
      s_ux(8,3,:,:)=-two*eta_i0*kappa_i*Ti**2.5*n_inv
     $     *(third*viy(1,:,:) + vix(2,:,:))
      s_ux(8,4,:,:)=-two*eta_i0*kappa_i*Ti**2.5*vix(3,:,:)*n_inv
      s_ux(8,8,:,:)=vi(1,:,:)

      s_uy(8,1,:,:)=-two*eta_i0*kappa_i*Ti**2.5
     $     *((one+third)*viy(2,:,:)*vi_un(2,:,:)
     $     + third*vix(2,:,:)*vi_un(1,:,:)
     $     + viy(1,:,:)*vi_un(1,:,:) + viy(3,:,:)*vi_un(3,:,:))
      s_uy(8,2,:,:)=-two*eta_i0*kappa_i*Ti**2.5*n_inv
     $     *(third*vix(2,:,:) + viy(1,:,:))
      s_uy(8,3,:,:)=-two*eta_i0*kappa_i*Ti**2.5
     $     *(one+third)*viy(2,:,:)*n_inv
      s_uy(8,4,:,:)=-two*eta_i0*kappa_i*Ti**2.5*viy(3,:,:)*n_inv
      s_uy(8,8,:,:)=vi(2,:,:)
c-----------------------------------------------------------------------
c     electron pressure.
c-----------------------------------------------------------------------
      fx_u(9,1,:,:) = gfac*u(9,:,:)*ve_un(1,:,:)
     $     - kface*BdotTe_un(1,:,:) - kperpe*Tex_un
     $     - (kpare_un - kperpe_un)*BdotTe(1,:,:) - kperpe_un*Tex
     $     + kperpe_bsq*Bsq_u(1,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*(Te_un*delVpar*bx + Te*delVpar*Ay_un(3,:,:)
     $     + Te*bx*delVpar_un)
      fx_u(9,2,:,:) = -fheat*Te*bx**2*invBsq
      fx_u(9,3,:,:) = -fheat*Te*bx*by*invBsq
      fx_u(9,4,:,:) = -fheat*Te*bx*bz*invBsq
      fx_u(9,5,:,:) = kface*BdotTe_bz(1,:,:)*ny_inv
     $     + kperpe_bsq*Bsq_u(2,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*bx*invBsq
     $     *(-(u(12,:,:)-u(4,:,:))*ny_inv - delVpar*Bsq_u(2,:,:))
      fx_u(9,6,:,:) = -kface*BdotTe_bz(1,:,:)*nx_inv    
     $     + kperpe_bsq*Bsq_u(3,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*bx*invBsq
     $     *((u(12,:,:)-u(4,:,:))*nx_inv - delVpar*Bsq_u(3,:,:))
      fx_u(9,7,:,:) = -kface
     $     *(BdotTe_bx(1,:,:)*ny_inv - BdotTe_by(1,:,:)*nx_inv)
     $     + kperpe_bsq*Bsq_u(4,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*(delVpar*ny_inv + bx*invBsq
     $     *((u(10,:,:)-u(2,:,:))*ny_inv - (u(11,:,:)-u(3,:,:))*nx_inv 
     $     - delVpar*Bsq_u(4,:,:)))
      fx_u(9,9,:,:) = gfac*ve(1,:,:) - kface
     $     *(BdotTe_Tex(1,:,:)*nx_inv + BdotTe_Tey(1,:,:)*ny_inv)
     $     - kperpe*nx_inv
     $     - (kpare_pe - kperpe_pe)*BdotTe(1,:,:) - kperpe_pe*Tex
     $     + fheat*n_inv*delVpar*bx
      fx_u(9,10,:,:) = Te*(gfac + fheat*bx**2*invBsq)
      fx_u(9,11,:,:) = fheat*Te*bx*by*invBsq
      fx_u(9,12,:,:) = fheat*Te*bx*bz*invBsq

      fx_ux(9,1,:,:) = -kface*BdotTe_uxn(1,:,:) - kperpe*Te_un
     $     + kperpe_bsq*Bsq_ux(1,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*delVpar_uxn*bx
      fx_ux(9,6,:,:) = -kface*BdotTe_bz(1,:,:)*n_inv
     $     + kperpe_bsq*Bsq_ux(3,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*bx*invBsq
     $     *((u(12,:,:)-u(4,:,:))*n_inv - delVpar*Bsq_ux(3,:,:))
      fx_ux(9,7,:,:) = kface*BdotTe_by(1,:,:)*n_inv
     $     + kperpe_bsq*Bsq_ux(4,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*bx*invBsq
     $     *(-(u(11,:,:)-u(3,:,:))*n_inv - delVpar*Bsq_ux(4,:,:))
      fx_ux(9,9,:,:) = -kface*BdotTe_Tex(1,:,:)*n_inv - kperpe*n_inv

      fx_uy(9,1,:,:) = -kface*BdotTe_uyn(1,:,:)
     $     + kperpe_bsq*Bsq_uy(1,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*(delVpar*A_un(3,:,:) + bx*delVpar_uyn)
      fx_uy(9,5,:,:) = kface*BdotTe_bz(1,:,:)*n_inv
     $     + kperpe_bsq*Bsq_uy(2,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*bx*invBsq
     $     *(-(u(12,:,:)-u(4,:,:))*n_inv - delVpar*Bsq_uy(2,:,:))
      fx_uy(9,7,:,:) = -kface*BdotTe_bx(1,:,:)*n_inv
     $     + kperpe_bsq*Bsq_uy(4,:,:)*(BdotTe(1,:,:) - Tex)
     $     + fheat*Te*(delVpar*n_inv + bx*invBsq
     $     *((u(10,:,:)-u(2,:,:))*n_inv - delVpar*Bsq_uy(4,:,:)))
      fx_uy(9,9,:,:) = -kface*BdotTe_Tey(1,:,:)*n_inv

      fy_u(9,1,:,:) = gfac*u(9,:,:)*ve_un(2,:,:)
     $     - kface*BdotTe_un(2,:,:) - kperpe*Tey_un
     $     - (kpare_un - kperpe_un)*BdotTe(2,:,:) - kperpe_un*Tey
     $     + kperpe_bsq*Bsq_u(1,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*(Te_un*delVpar*by - Te*delVpar*Ax_un(3,:,:)
     $     + Te*by*delVpar_un)
      fy_u(9,2,:,:) = -fheat*Te*bx*by*invBsq
      fy_u(9,3,:,:) = -fheat*Te*by**2*invBsq
      fy_u(9,4,:,:) = -fheat*Te*by*bz*invBsq
      fy_u(9,5,:,:) = kface*BdotTe_bz(2,:,:)*ny_inv
     $     + kperpe_bsq*Bsq_u(2,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*by*invBsq
     $     *(-(u(12,:,:)-u(4,:,:))*ny_inv - delVpar*Bsq_u(2,:,:))
      fy_u(9,6,:,:) = -kface*BdotTe_bz(2,:,:)*nx_inv    
     $     + kperpe_bsq*Bsq_u(3,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*by*invBsq
     $     *((u(12,:,:)-u(4,:,:))*nx_inv - delVpar*Bsq_u(3,:,:))
      fy_u(9,7,:,:) = -kface
     $     *(BdotTe_bx(2,:,:)*ny_inv - BdotTe_by(2,:,:)*nx_inv)
     $     + kperpe_bsq*Bsq_u(4,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*(-delVpar*nx_inv + by*invBsq
     $     *((u(10,:,:)-u(2,:,:))*ny_inv - (u(11,:,:)-u(3,:,:))*nx_inv 
     $     - delVpar*Bsq_u(4,:,:)))
      fy_u(9,9,:,:) = gfac*ve(2,:,:) - kface
     $     *(BdotTe_Tex(2,:,:)*nx_inv + BdotTe_Tey(2,:,:)*ny_inv)
     $     - kperpe*ny_inv
     $     - (kpare_pe - kperpe_pe)*BdotTe(2,:,:) - kperpe_pe*Tey
     $     + fheat*n_inv*delVpar*by
      fy_u(9,10,:,:) = fheat*Te*bx*by*invBsq
      fy_u(9,11,:,:) = Te*(gfac + fheat*by**2*invBsq)
      fy_u(9,12,:,:) = fheat*Te*by*bz*invBsq

      fy_ux(9,1,:,:) = -kface*BdotTe_uxn(2,:,:)
     $     + kperpe_bsq*Bsq_ux(1,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*(-delVpar*A_un(3,:,:) + by*delVpar_uxn)
      fy_ux(9,6,:,:) = -kface*BdotTe_bz(2,:,:)*n_inv
     $     + kperpe_bsq*Bsq_ux(3,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*by*invBsq
     $     *((u(12,:,:)-u(4,:,:))*n_inv - delVpar*Bsq_ux(3,:,:))
      fy_ux(9,7,:,:) = kface*BdotTe_by(2,:,:)*n_inv
     $     + kperpe_bsq*Bsq_ux(4,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*(-delVpar*n_inv + by*invBsq
     $     *(-(u(11,:,:)-u(3,:,:))*n_inv - delVpar*Bsq_ux(4,:,:)))
      fy_ux(9,9,:,:) = -kface*BdotTe_Tex(2,:,:)*n_inv

      fy_uy(9,1,:,:) = -kface*BdotTe_uyn(2,:,:) - kperpe*Te_un
     $     + kperpe_bsq*Bsq_uy(1,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*by*delVpar_uyn
      fy_uy(9,5,:,:) = kface*BdotTe_bz(2,:,:)*n_inv
     $     + kperpe_bsq*Bsq_uy(2,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*by*invBsq
     $     *(-(u(12,:,:)-u(4,:,:))*n_inv - delVpar*Bsq_uy(2,:,:))
      fy_uy(9,7,:,:) = -kface*BdotTe_bx(2,:,:)*n_inv
     $     + kperpe_bsq*Bsq_uy(4,:,:)*(BdotTe(2,:,:) - Tey)
     $     + fheat*Te*by*invBsq
     $     *((u(10,:,:)-u(2,:,:))*n_inv - delVpar*Bsq_uy(4,:,:))
      fy_uy(9,9,:,:) = -kface*BdotTe_Tey(2,:,:)*n_inv - kperpe*n_inv

      s_u(9,1,:,:)=ve_un(1,:,:)*ux(9,:,:) + ve_un(2,:,:)*uy(9,:,:)
     $     - eta_e0*kappa_e*Te**1.5*(2.5_r8*Te_un*PedelV
     $     + two*Te*((one+third)*(vex(1,:,:)*vex_un(1,:,:)
     $     + vey(2,:,:)*vey_un(2,:,:))
     $     + third*vex(2,:,:)*vey_un(1,:,:)
     $     + third*vey(1,:,:)*vex_un(2,:,:)
     $     + vey(1,:,:)*vey_un(1,:,:) + vex(2,:,:)*vex_un(2,:,:)
     $     + vex(3,:,:)*vex_un(3,:,:) + vey(3,:,:)*vey_un(3,:,:)))
     $     + (2.5_r8*ieheat*(u(8,:,:)-u(9,:,:))*Te 
     $     - 1.5_r8*Te_un*resist*delVsq)/(kappa_e*Te**2.5)
     $     + fheat*((u(10,:,:)-u(2,:,:))*BdotTe_un(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_un(2,:,:) 
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_un(3,:,:))
      s_u(9,2,:,:)=-two*resist/(kappa_e*Te**1.5)*(u(10,:,:)-u(2,:,:)) 
     $     - fheat*BdotTe(1,:,:)
      s_u(9,3,:,:)=-two*resist/(kappa_e*Te**1.5)*(u(11,:,:)-u(3,:,:))
     $     - fheat*BdotTe(2,:,:)
      s_u(9,4,:,:)=-two*resist/(kappa_e*Te**1.5)*(u(12,:,:)-u(4,:,:))
     $     - fheat*BdotTe(3,:,:)
      s_u(9,5,:,:)=-ny_inv*fheat*((u(10,:,:)-u(2,:,:))*BdotTe_bz(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_bz(2,:,:) 
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_bz(3,:,:))
      s_u(9,6,:,:)=nx_inv*fheat*((u(10,:,:)-u(2,:,:))*BdotTe_bz(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_bz(2,:,:) 
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_bz(3,:,:))
      s_u(9,7,:,:)=fheat*((u(10,:,:)-u(2,:,:))
     $     *(ny_inv*BdotTe_bx(1,:,:) - nx_inv*BdotTe_by(1,:,:))
     $     + (u(11,:,:)-u(3,:,:))
     $     *(ny_inv*BdotTe_bx(2,:,:) - nx_inv*BdotTe_by(2,:,:))
     $     + (u(12,:,:)-u(4,:,:))
     $     *(ny_inv*BdotTe_bx(3,:,:) - nx_inv*BdotTe_by(3,:,:)))
      s_u(9,8,:,:)=ieheat*u(1,:,:)/(kappa_e*Te**1.5)
      s_u(9,9,:,:)=-eta_e0*kappa_e*2.5_r8*Te**1.5*n_inv*PedelV
     $     + (ieheat*(0.5_r8*u(9,:,:) - 1.5_r8*u(8,:,:)) 
     $     - 1.5_r8*n_inv*resist*delVsq)/(kappa_e*Te**2.5)
     $     + fheat*((u(10,:,:)-u(2,:,:))
     $     *(nx_inv*BdotTe_Tex(1,:,:) + ny_inv*BdotTe_Tey(1,:,:))
     $     + (u(11,:,:)-u(3,:,:))
     $     *(nx_inv*BdotTe_Tex(2,:,:) + ny_inv*BdotTe_Tey(2,:,:))
     $     + (u(12,:,:)-u(4,:,:))
     $     *(nx_inv*BdotTe_Tex(3,:,:) + ny_inv*BdotTe_Tey(3,:,:)))
      s_u(9,10,:,:)=n_inv*ux(9,:,:)
     $     - two*eta_e0*kappa_e*Te**2.5
     $     *((one+third)*vex(1,:,:)*nx_inv
     $     + third*vex(2,:,:)*ny_inv + vey(1,:,:)*ny_inv)
     $     + two*resist/(kappa_e*Te**1.5)*(u(10,:,:)-u(2,:,:))
     $     + fheat*BdotTe(1,:,:)
      s_u(9,11,:,:)=n_inv*uy(9,:,:)
     $     - two*eta_e0*kappa_e*Te**2.5
     $     *((one+third)*vey(2,:,:)*ny_inv
     $     + third*vey(1,:,:)*nx_inv + vex(2,:,:)*nx_inv)
     $     + two*resist/(kappa_e*Te**1.5)*(u(11,:,:)-u(3,:,:))
     $     + fheat*BdotTe(2,:,:)
      s_u(9,12,:,:)=-two*eta_e0*kappa_e*Te**2.5
     $     *(vex(3,:,:)*nx_inv + vey(3,:,:)*ny_inv)
     $     + two*resist/(kappa_e*Te**1.5)*(u(12,:,:)-u(4,:,:))
     $     + fheat*BdotTe(3,:,:)

      s_ux(9,1,:,:)=-two*eta_e0*kappa_e*Te**2.5
     $     *((one+third)*vex(1,:,:)*ve_un(1,:,:)
     $     + third*vey(1,:,:)*ve_un(2,:,:)
     $     + vex(2,:,:)*ve_un(2,:,:) + vex(3,:,:)*ve_un(3,:,:))
     $     + fheat*((u(10,:,:)-u(2,:,:))*BdotTe_uxn(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_uxn(2,:,:) 
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_uxn(3,:,:))
      s_ux(9,6,:,:)=n_inv*fheat*((u(10,:,:)-u(2,:,:))*BdotTe_bz(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_bz(2,:,:) 
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_bz(3,:,:))
      s_ux(9,7,:,:)=-n_inv*fheat*((u(10,:,:)-u(2,:,:))*BdotTe_by(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_by(2,:,:)
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_by(3,:,:))
      s_ux(9,9,:,:)=ve(1,:,:)
     $     + n_inv*fheat*((u(10,:,:)-u(2,:,:))*BdotTe_Tex(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_Tex(2,:,:)
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_Tex(3,:,:))
      s_ux(9,10,:,:)=-two*eta_e0*kappa_e*Te**2.5
     $     *(one+third)*vex(1,:,:)*n_inv
      s_ux(9,11,:,:)=-two*eta_e0*kappa_e*Te**2.5*n_inv
     $     *(third*vey(1,:,:) + vex(2,:,:))
      s_ux(9,12,:,:)=-two*eta_e0*kappa_e*Te**2.5*vex(3,:,:)*n_inv

      s_uy(9,1,:,:)=-two*eta_e0*kappa_e*Te**2.5
     $     *((one+third)*vey(2,:,:)*ve_un(2,:,:)
     $     + third*vex(2,:,:)*ve_un(1,:,:)
     $     + vey(1,:,:)*ve_un(1,:,:) + vey(3,:,:)*ve_un(3,:,:))
     $     + fheat*((u(10,:,:)-u(2,:,:))*BdotTe_uyn(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_uyn(2,:,:) 
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_uyn(3,:,:))
      s_uy(9,5,:,:)=-n_inv*fheat*((u(10,:,:)-u(2,:,:))*BdotTe_bz(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_bz(2,:,:) 
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_bz(3,:,:))
      s_uy(9,7,:,:)=n_inv*fheat*((u(10,:,:)-u(2,:,:))*BdotTe_bx(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_bx(2,:,:)
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_bx(3,:,:))
      s_uy(9,9,:,:)=ve(2,:,:)
     $     + n_inv*fheat*((u(10,:,:)-u(2,:,:))*BdotTe_Tey(1,:,:)
     $     + (u(11,:,:)-u(3,:,:))*BdotTe_Tey(2,:,:)
     $     + (u(12,:,:)-u(4,:,:))*BdotTe_Tey(3,:,:))
      s_uy(9,10,:,:)=-two*eta_e0*kappa_e*Te**2.5*n_inv
     $     *(third*vex(2,:,:) + vey(1,:,:))
      s_uy(9,11,:,:)=-two*eta_e0*kappa_e*Te**2.5
     $     *(one+third)*vey(2,:,:)*n_inv
      s_uy(9,12,:,:)=-two*eta_e0*kappa_e*Te**2.5*vey(3,:,:)*n_inv
c-----------------------------------------------------------------------
c     x-Ampere's Law.
c-----------------------------------------------------------------------
      fx_u(10,1,:,:)=-Ay_un(2,:,:)
      fx_u(10,6,:,:)=-ny_inv
      fx_uy(10,1,:,:)=-A_un(2,:,:)
      fx_uy(10,6,:,:)=-n_inv

      fy_u(10,1,:,:)=Ay_un(1,:,:)
      fy_u(10,5,:,:)=ny_inv
      fy_uy(10,1,:,:)=A_un(1,:,:)
      fy_uy(10,5,:,:)=n_inv

      s_u(10,2,:,:)=-d_inv
      s_u(10,10,:,:)=d_inv
c-----------------------------------------------------------------------
c     y-Ampere's Law.
c-----------------------------------------------------------------------
      fx_u(11,1,:,:)=Ax_un(2,:,:)
      fx_u(11,6,:,:)=nx_inv
      fx_ux(11,1,:,:)=A_un(2,:,:)
      fx_ux(11,6,:,:)=n_inv

      fy_u(11,1,:,:)=-Ax_un(1,:,:)
      fy_u(11,5,:,:)=-nx_inv
      fy_ux(11,1,:,:)=-A_un(1,:,:)
      fy_ux(11,5,:,:)=-n_inv

      s_u(11,3,:,:)=-d_inv
      s_u(11,11,:,:)=d_inv
c-----------------------------------------------------------------------
c     z-Ampere's Law.
c-----------------------------------------------------------------------
      fx_u(12,1,:,:)=Ax_un(3,:,:)
      fx_u(12,7,:,:)=nx_inv
      fx_ux(12,1,:,:)=A_un(3,:,:)
      fx_ux(12,7,:,:)=n_inv

      fy_u(12,1,:,:)=Ay_un(3,:,:)
      fy_u(12,7,:,:)=ny_inv
      fy_uy(12,1,:,:)=A_un(3,:,:)
      fy_uy(12,7,:,:)=n_inv

      s_u(12,4,:,:)=-d_inv
      s_u(12,12,:,:)=d_inv
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
      USE XMHD_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      mass(2,10,:,:)=mass_r
      mass(3,11,:,:)=mass_r
      mass(4,12,:,:)=mass_r

      mass(5,5,:,:)=d_inv
      mass(5,10,:,:)=-mass_r
      mass(6,6,:,:)=d_inv
      mass(6,11,:,:)=-mass_r
      mass(7,7,:,:)=d_inv
      mass(7,12,:,:)=-mass_r
      mass(8,8,:,:)=one/(gamma-one)
      mass(9,9,:,:)=one/(gamma-one)
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
      USE XMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,eta
c-----------------------------------------------------------------------
c     special assignments.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("DispTest1a","DispTest1b","DispTest1c",
     $     "DispTest2a","DispTest2b","DispTest2c",
     $     "DispTest3a","DispTest3b","DispTest3c")
         ksi=twopi*x/kx
         eta=twopi*y/ky     
      CASE("DispTest4")
         ksi=lx*twopi*x/kx
         eta=lx*twopi*y/ky     
      CASE("GEM","uniform")
         ksi=lx*(x-half)
         eta=ly*(y-half)
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
      USE XMHD_mod
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
      USE XMHD_mod
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
