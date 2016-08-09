c-----------------------------------------------------------------------
c     file SolarLnMHD_2F.f
c     contains specifications for 2-pressure MHD model with gravity
c     and evolving ln(density).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. SolarLnMHD_mod.
c     1. SolarLnMHD_atmsph.
c     2. SolarLnMHD_psi.
c     3. SolarLnMHD_equil.
c     4. ChenShibata_eta.
c     5. ChenShibata_mu.
c     6. MAST_eta.
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
c     subprogram 0. SolarLnMHD_mod.
c     module declarations.
c-----------------------------------------------------------------------
c     u(1)  =log(n)
c     u(2)  =psi
c     u(3)  =bz 
c     u(4:6)=nv
c     u(7)  =jz
c     u(8)  =pe
c     u(9)  =pi
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE SolarLnMHD_mod
      USE transport_mod
      USE extra_mod
      IMPLICIT NONE

      LOGICAL :: source=.FALSE.,cylinder=.FALSE.
      CHARACTER(18) :: init_type=".",eta_case=".",kappa_case="."
      REAL, PARAMETER :: Zeff=1.
      INTEGER :: cyl_fac=0
      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,qe=1.602e-19,
     $     me=9.109e-31,mp=1.673e-27,ep0=8.854e-12,mu0=4.e-7*pi,
     $     k_B=1.381e-23,g0=2.74e2,chod_const=0.1
      REAL(r8) :: R0=6.955e8 ! Solar radius = 6.955e8(m)
      REAL(r8) :: eta=0.,eta_chrsp=0.,eta_vbl=0.,j_c=1.e8,etavac=1.,
     $     mu=0.,mu_vbl=0.,kappa_prp=0.,kappa_min=0.,kappa_max=1.e8,
     $     ieheat=3._r8,n0=1.e15,b0=1.e-3,bz0=0.,beta0=1.,lx=0.,
     $     ly=0.,x_curve=1.,y_curve=1.,rad0=1.,x0=0.,c_psi=0.,
     $     c_psi_e=0.,h_psi=0.,t_e=0.,ke_norm=1.,ki_norm=1.,xe_norm=1.,
     $     xi_norm=1.,gamma_fac=1.,gravity=1.,v_chod_norm=1.,
     $     etac_norm=1.,etas_norm=1.,hyper_eta=0.,htr=1.,wtr=1.,
     $     Tph=1.,Tco=1.,p0=1.,y_eta=1.,y0_eta=0.,mu_min=0.,epsilon=0.,
     $     rad_fac=1.,T0=1.0, c_rho=0., c_T = 0., 
     $     jc_norm = 1., eta_anm_norm = 1.,
     $     hlfw = 0., lambda = 0., ddiff = 0.
c      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: heat_blnc_radlos

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. SolarLnMHD_atmsph.
c     computes atmospherically stratified density and pressure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE SolarLnMHD_atmsph(x,y,rho,p,px,py,Tx,Ty)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: rho,p,px,py,Tx,Ty

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: w,Tmpr
c-----------------------------------------------------------------------
c     compute density and pressure.
c-----------------------------------------------------------------------
      w = (y-htr)/wtr       
c Elena htr,wtr?
      Tmpr = Tph + half*(Tco-Tph)*(one + TANH(w))
c Elena Tph - photospheric temperature, Tco - coronal temperature     
      p = half*p0*EXP(-half*gravity*wtr/Tco
     $     *(w - half*(Tco-Tph)/Tph*(LOG(EXP(-two*w)+Tco/Tph)-9.352)))
c Elena see Lee etal. 2014    
      px = 0
      py = -p*half*gravity/Tco
     $     *(one + (Tco-Tph)*EXP(-two*w)/(Tph*EXP(-two*w)+Tco))

      rho = p/Tmpr

      Tx = 0
      Ty = half*(Tco-Tph)/(wtr*COSH(w)**2)
      
!     htr = 0.5    ! htr*R0=3Mm
!     wtr = 0.1    ! wtr*R0=0.75Mm
!     T0 = b0*b0/(n0*mu0*k_B)
!     Tph = 5000/T0 ! Photospheric temperature (K)
!     Tco = 1.e6/T0 ! Coronal temperature (K)
!     T = Tph+0.5*(Tco-Tph)*(1.+tanh(w))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SolarLnMHD_atmsph
c-----------------------------------------------------------------------
c     subprogram 2. SolarLnMHD_psi.
c     sets psi.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION SolarLnMHD_psi(x,y)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r,SolarLnMHD_psi
c-----------------------------------------------------------------------
c     compute psi.
c-----------------------------------------------------------------------
      r = SQRT(x**2+(y-h_psi)**2)/rad0            
      SolarLnMHD_psi= -c_psi*LOG(((x+0.3)**2+(y+0.3)**2)
     $     * ((x-0.3)**2+(y+0.3)**2)
     $     /( ((x+1.5)**2+(y+0.3)**2) * ((x-1.5)**2+(y+0.3)**2) ))
     $     + half*rad0*LOG(x**2+(y+h_psi)**2)
      WHERE(r.LE.one)
         SolarLnMHD_psi=SolarLnMHD_psi
     $        + half*rad0*(-r**2 + half*(r**2-one)**2)
      ELSEWHERE
         SolarLnMHD_psi=SolarLnMHD_psi - rad0*(half + LOG(r))
      END WHERE
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION SolarLnMHD_psi
c-----------------------------------------------------------------------
c     subprogram 3. SolarLnMHD_equil.
c     computes equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE SolarLnMHD_equil(x,y,u, ux, uy, deriv)
c Elena: Equilibrium initial state      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u, ux, uy
      LOGICAL, INTENT(IN) :: deriv

      REAL(r8), DIMENSION(1,1) :: psitmp,xt,yt
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: rsq,r,p,px,py,Tx,Ty
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      CALL RANDOM_NUMBER(u)
      u = 1.e-20*u

      SELECT CASE(init_type)
c Elena begin: h_psi - width of the current sheet            
      CASE("CurrentSheet-hlf")
c        u(1,:,:)=LOG(one + c_rho*(one - (TANH(x/h_psi))**2))       
c Harris sheet
        u(1,:,:)=LOG(one + 1./beta0*(one - (TANH(x/h_psi))**2))
        u(2,:,:)= h_psi*LOG(COSH(x/h_psi))
c        u(8,:,:)=0.25_r8*beta0*(one + c_T*(one - (TANH(x/h_psi))**2))
c     $          *EXP(u(1,:,:))
c        u(8,:,:)=0.25_r8*beta0*(EXP(u(1,:,:)))**(gamma)
        
c        u(3,:,:)= SQRT(bz0**2 + one - (TANH(x/h_psi))**2
c     $    + 4.*0.25_r8*beta0 - 4.*u(8,:,:))-bz0
        u(3,:,:)=0.
c        u(3,:,:)= SQRT(bz0**2 + one - (TANH(x/h_psi))**2)-bz0
        u(7,:,:)=1._r8/h_psi/COSH(x/h_psi)/COSH(x/h_psi)
c        u(8,:,:)=0.25_r8*beta0
c        u(9,:,:)=0.25_r8*beta0
        u(8,:,:)=0.25_r8*beta0*(1.+1./beta0*(one - (TANH(x/h_psi))**2))
        u(9,:,:)=0.25_r8*beta0*(1.+1./beta0*(one - (TANH(x/h_psi))**2)) 
c        u(9,:,:)=0.25_r8*beta0*(one + c_T*(one - (TANH(x/h_psi))**2))
c     $          *EXP(u(1,:,:)) 
c        u(9,:,:)=0.25_r8*beta0*(EXP(u(1,:,:)))**(gamma)
        IF(.NOT. deriv)RETURN
c        ux(1,:,:)= -EXP(-(x/h_psi)**2)*2.*x/(one 
c     $  + EXP(-(x/h_psi)**2))/h_psi**2
c        ux(2,:,:)= TANH(x/h_psi)
c        ux(3,:,:)= 0.6*0.5_r8/SQRT(bz0**2 + one - TANH(x/h_psi)**2)
c     $  *(-2.*TANH(x/h_psi)/COSH(x/h_psi)**2/h_psi)
c        ux(7,:,:) = -2.*TANH(x/h_psi)/h_psi**2/COSH(x/h_psi)**2     
c Elena end
      CASE("TwoFR")
        u(1,:,:) = LOG(one)
        u(2,:,:) = - lambda*LOG(COSH(y/lambda)+hlfw*COS(x/lambda))
        u(3,:,:) = SQRT(bz0*bz0 + (1. - hlfw**2)/
     $             (COSH(y/lambda) + hlfw*COS(x/lambda))**2) - bz0 
        u(7,:,:) = hlfw/lambda*(COS(x/lambda)*(COSH(y/lambda)+
     $             hlfw*COS(x/lambda)) + hlfw*SIN(x/lambda)**2)*
     $             (COSH(y/lambda)+hlfw*COS(x/lambda))**(-2) +
     $             (-COSH(y/lambda)*(COSH(y/lambda)+hlfw*COS(x/lambda))+
     $             SINH(y/lambda)**2)/lambda*
     $             (COSH(y/lambda)+hlfw*COS(x/lambda))**(-2)
        u(8,:,:) = 0.25_r8*beta0
        u(9,:,:) = u(8,:,:)
      CASE("Chen-Shibata","Chen-Shibata-hlf")
        r=SQRT(x**2+(y-h_psi)**2)/rad0
        xt = zero
        yt = h_psi
        psitmp = SolarLnMHD_psi(xt,yt)

        u(1,:,:)=LOG(one)
        u(2,:,:)=SolarLnMHD_psi(x,y)

        WHERE(y.GT.0.5)
           rsq = MAX(0.,3.+4./rad0*(u(2,:,:)-psitmp(1,1)+0.25*rad0))
           rsq = MIN(1., 2.-SQRT(rsq))
           u(3,:,:)=SQRT(bz0**2 + MAX(0.,10./3.-8.*rsq
     $          + 6.*rsq**2 - (4./3.)*rsq**3)) - bz0
        END WHERE

        WHERE(r.LE.one)
           u(7,:,:)=u(7,:,:)-4.*(1-r**2)/rad0
        END WHERE
        u(8,:,:)=0.25_r8*beta0*EXP(u(1,:,:))
        u(9,:,:)=0.25_r8*beta0*EXP(u(1,:,:))

      CASE("CS-stratified")
        r=SQRT(x**2+(y-h_psi)**2)/rad0
        xt = zero
        yt = h_psi
        psitmp = SolarLnMHD_psi(xt,yt)

        CALL SolarLnMHD_atmsph(x,y,u(1,:,:),p,px,py,Tx,Ty)
        u(1,:,:)=LOG(u(1,:,:))
        u(2,:,:)=SolarLnMHD_psi(x,y)

        WHERE(y.GT.0.5)
           rsq = MAX(0.,3.+4./rad0*(u(2,:,:)-psitmp(1,1)+0.25*rad0))
           rsq = MIN(1., 2.-SQRT(rsq))
           u(3,:,:)=SQRT(bz0**2 + MAX(0.,10./3.-8.*rsq
     $          + 6.*rsq**2 - (4./3.)*rsq**3)) - bz0
        END WHERE

        WHERE(r.LE.one)
           u(7,:,:)=u(7,:,:)-4.*(1-r**2)/rad0
        END WHERE

      CASE("CurrentSheet")      
        CALL SolarLnMHD_atmsph(x,y,u(1,:,:),p,px,py,Tx,Ty)
        rsq = h_psi*(y+rad0)/rad0      
        u(1,:,:) = LOG(u(1,:,:))
        u(2,:,:) = rsq*LOG(COSH(x/rsq))
        u(3,:,:) = SQRT(bz0**2 + one 
     $       - (one + x**2/(y+rad0)**2)*TANH(x/rsq)**2 
     $       - (h_psi/rad0)**2*LOG(COSH(x/rsq))**2
     $       + two*h_psi*x/rad0/(y+rad0)*LOG(COSH(x/rsq))*TANH(x/rsq)) 
     $       - bz0
        u(7,:,:) = (one + x**2/(y+one)**2)/(rsq*COSH(x/rsq)**2)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SolarLnMHD_equil
c-----------------------------------------------------------------------
c     subprogram 4. ChenShibata_eta.
c     computes resistive layers
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ChenShibata_eta(x,y,jtot,eta_vbl,eta_chrsp,eta_local)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y,jtot
      REAL(r8), INTENT(IN) :: eta_vbl,eta_chrsp
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: eta_local

c-----------------------------------------------------------------------
c     compute resistivity.
c-----------------------------------------------------------------------

      eta_local = eta + eta_chrsp*EXP(-((y-y0_eta)/y_eta)**2)
     $     + eta_vbl*EXP(-((x-half*lx)/0.125)**2)
     $     + eta_vbl*EXP(-((y-ly)/0.5)**2)

      WHERE(jtot.GE.j_c .AND. jtot.LE.(two*j_c))
         eta_local = eta_local
     $        + etavac*half*(one-COS(pi*(jtot-j_c)/j_c))
      ELSEWHERE(jtot > (two*j_c))
         eta_local = eta_local + etavac
      END WHERE
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ChenShibata_eta
c-----------------------------------------------------------------------
c     subprogram 5. ChenShibata_mu.
c     computes viscous layers
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ChenShibata_mu(x,y,rho,mu_vbl,mu_local)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y,rho
      REAL(r8), INTENT(IN) :: mu_vbl
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: mu_local

c-----------------------------------------------------------------------
c     compute viscosity.
c-----------------------------------------------------------------------

c Elena: original expression
      mu_local = rho*(mu + mu_vbl*EXP(-((x-half*lx)/0.1)**2)
     $     + mu_vbl*EXP(-((y-ly)/0.4)**2)) + mu_min

c Elena: for the case with MAST initial conditions
c      mu_local = rho*(mu + mu_vbl*EXP(-((x-lx)/0.1)**2)
c     $     + mu_vbl*EXP(-((y-ly)/0.1)**2)) + mu_min
 
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ChenShibata_mu
c-----------------------------------------------------------------------
c     subprogram 6. MAST_eta.
c     computes resistive layers
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE MAST_eta(x,y,eta_vbl,eta_local)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), INTENT(IN) :: eta_vbl 
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: eta_local

c-----------------------------------------------------------------------
c     compute resistivity.
c-----------------------------------------------------------------------

      eta_local = eta  
     $     + eta_vbl*EXP(-((x-lx)/0.1)**2)
     $     + eta_vbl*EXP(-((y-ly)/0.1)**2)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE MAST_eta
c-----------------------------------------------------------------------
c     subprogram 7. MAST_ddiff.
c     computes diffusive layers
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE MAST_ddiff(x,y,ddiff,ddiff_local)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), INTENT(IN) :: ddiff
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ddiff_local

c-----------------------------------------------------------------------
c     compute mass diffusion coefficient
c-----------------------------------------------------------------------

      ddiff_local = ddiff*EXP(-((x-lx)/0.1)**2)
     $            + ddiff*EXP(-((y-ly)/0.1)**2)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE MAST_ddiff

      END MODULE SolarLnMHD_mod

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
      USE SolarLnMHD_mod
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
      NAMELIST/SolarLnMHD_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,init_type,source,cylinder,eta,eta_case,
     $     eta_chrsp,eta_vbl,j_c,etavac,mu,mu_vbl,mu_min,kappa_case,
     $     kappa_prp,kappa_min,kappa_max,ieheat,n0,b0,beta0,lx,ly,
     $     x_curve,y_curve,rad0,x0,c_psi,y0_eta,y_eta,c_psi_e,h_psi,t_e,
     $     hyper_eta,bz0,p0,R0,epsilon, c_rho, c_T, hlfw, lambda, ddiff
c Elena: source? cylinder?(geometry?) eta? - resistivity
c Elena: mu - viscosity? mu_vbl? kappa? 
c Elena: n0,b0,beta0 - normalizations?
c Elena: all other parameters what they mean     
c-----------------------------------------------------------------------
c     Sample namelist.
c-----------------------------------------------------------------------
c$$$&SolarLnMHD_list
c$$$
c$$$	cylinder=f
c$$$
c$$$	init_type="Chen-Shibata"
c$$$	source=f
c$$$
c$$$  R0=6.955e8
c$$$  lx=24.0
c$$$  ly=18.0
c$$$  x_curve=5.
c$$$  y_curve=0.
c$$$
c$$$  beta0=0.01
c$$$
c$$$  rad0=0.5
c$$$  x0=0.0 or 3.9
c$$$  c_psi=2.5628
c$$$  c_psi_e=11. or -25.
c$$$  h_psi=2
c$$$
c$$$	b0=1.
c$$$	n0=1.
c$$$  bz0=0.
c$$$  p0=1.
c$$$
c$$$	eta_case="uniform" 
c$$$	eta=1.e-4  
c$$$	eta_chrsp=1.e-2
c$$$	eta_vbl=1.e-2
c$$$    j_c=0.5
c$$$	etavac=0.
c$$$
c$$$    kappa_case="isotropic"
c$$$	kappa_prp=1.e-3
c$$$	kappa_min=1.e-6
c$$$	kappa_max=1.
c$$$    ieheat=3.
c$$$
c$$$	mu=1.e-5
c$$$    mu_vbl=1.
c$$$    mu_min=5.e-6
c$$$/
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=SolarLnMHD_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="SolarLnMHD"

      nqty=9
c Elena: nqty - number of variables      
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
      USE SolarLnMHD_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c      REAL(r8), DIMENSION(1,1) :: rho_init, p_init, alpha
c      REAL(r8), DIMENSION(SIZE(u,2), SIZE(u,3)) :: rho_init,p_init,alpha

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
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta_chrsp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta_vbl,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(y0_eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(y_eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(j_c,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(etavac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu_vbl,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu_min,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_prp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_min,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_max,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ieheat,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(R0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(p0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bz0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(x_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(y_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rad0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(x0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(c_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(c_psi_e,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(h_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(t_e,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(hyper_eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(c_rho,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(c_T,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(hlfw,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ddiff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(7)=.TRUE.
      adapt_qty(7)=.FALSE.
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
      tnorm=R0*SQRT(mu0*n0*mp*Zeff)/b0
      SELECT CASE(init_type)
      CASE("Chen-Shibata","Chen-Shibata-hlf","CurrentSheet-hlf","MAST")
         gravity=0.
      CASE DEFAULT
         gravity=g0*tnorm**2/R0     
      END SELECT

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
      rad_fac = n0*n0*R0/(b0*b0/mu0)/(b0/SQRT(mu0*mp*n0)) 
c      jc_norm = 2.14*iqe*R0*SQRT(mu0*n0/mp)
c Elena: reducing critical value of current
      jc_norm = 0.0011*2.14*qe*R0*SQRT(mu0*n0/mp)
      eta_anm_norm = 1.4e-3*(3.e8)*mp/(42.85*qe*b0*R0)

      SELECT CASE(init_type)
      CASE("Chen-Shibata","Chen-Shibata-hlf","CurrentSheet-hlf")
        T0 = b0*b0/(n0*mu0*k_B)
        Tco = 1.e6/T0  ! Coronal temperature (K)
      CASE("CS-stratified","CurrentSheet")
        htr = 2.5e6/R0 ! htr=2.5Mm/R0=0.5
        wtr = 5.e5/R0  ! wtr=0.5Mm/R0=0.1
        T0 = b0*b0/(n0*mu0*k_B)
        Tph = 5000/T0  ! Photospheric temp (K)
        Tco = 1.e6/T0  ! Coronal temperature (K)
      CASE("TwoFR")
        lambda = lambda/2./pi
      END SELECT

c      SELECT CASE(init_type)
c      CASE("CurrentSheet-hlf")
c         rho_init = one + c_rho*(one - (TANH(x/h_psi))**2) 
c         p_init = 0.25*beta0
c         CALL transport_radloss(rad_fac,T0,rho_init,p_init,alpha,
c     $     heat_blnc_radlos)
c      DEFAULT
c         rho_init = one
c         p_init = 0.25*beta0
c      CALL transport_radloss(rad_fac,T0,rho_init,p_init,alpha,
c     $     heat_blnc_radlos)
c      END SELECT
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
      USE SolarLnMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     SolarLnMHD initial conditions.
c-----------------------------------------------------------------------
      CALL SolarLnMHD_equil(x,y,u,ux,uy,.FALSE.)
      SELECT CASE(init_type)
      CASE("CurrentSheet-hlf")
      u(2,:,:) = u(2,:,:)+ h_psi*epsilon
     $        *EXP(-(y-half*ly)**2/(two*h_psi)**2)
     $        *EXP(-x**2/(half*h_psi)**2)

      CASE("CurrentSheet")
         u(2,:,:) = u(2,:,:) + c_psi*EXP(-x**2/(half*h_psi)**2)
     $        *EXP(-(y-3.*htr)**2/(two*h_psi)**2)
      CASE("TwoFR")
         u(2,:,:) = u(2,:,:) 
     $             + epsilon*SIN(x/lx/lambda)*COS(y*pi/ly/2.)
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
      USE SolarLnMHD_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------

c Elena: boundary conditions by default
c robin means a boundary conditions in equation form
      top%bc_type="robin"
      top%static=.FALSE.
      bottom%bc_type="robin"
      bottom%static=.FALSE.
      left%bc_type="robin"
      left%static=.FALSE.
      right%bc_type="robin"
      right%static=.FALSE.

c Elena: 1 - left, 4 - bottom, 2 - top, 3 - right
      SELECT CASE(init_type)
      CASE("CurrentSheet-hlf")
         edge_order=(/1,4,2,3/)
         top%bc_type="periodic"   
               
         bottom%bc_type="periodic"
         
         left%bc_type(1:3)="zeroflux" ! symmetric
         left%bc_type(5)="zeroflux"   ! symmetric
         left%bc_type(7:9)="zeroflux" ! only x,z-momentum asymmetric
         left%static(4)=.TRUE.
         left%static(6)=.TRUE.
         
         right%static(3:9)=.TRUE.
         right%bc_type(1)="natural"

      CASE("TwoFR")
         top%static(3:9)=.TRUE.
         top%bc_type(1)="natural"          

         left%bc_type(1:3)="zeroflux"
         left%bc_type(5)="zeroflux"
         left%bc_type(7:9)="zeroflux"
         left%static(4)=.TRUE.
         left%static(6)=.TRUE.
         
         right%bc_type(1:3)="zeroflux"
         right%bc_type(5)="zeroflux"
         right%bc_type(7:9)="zeroflux"
         right%static(4)=.TRUE.
         right%static(6)=.TRUE.

         bottom%bc_type(1:4)="zeroflux"
         bottom%bc_type(7:9)="zeroflux"
         bottom%static(5)=.TRUE.
         bottom%static(6)=.TRUE.
              
      CASE("MAST")
         edge_order=(/2,4,1,3/)

c         top%bc_type(1)="natural"
         top%static(1)=.TRUE.
         top%static(3:9)=.TRUE.

         left%bc_type(1:3)="zeroflux"
         left%bc_type(5)="zeroflux"
         left%bc_type(7:9)="zeroflux"
         left%static(4)=.TRUE.
         left%static(6)=.TRUE.

c         right%bc_type(1)="natural"
         right%static(1)=.TRUE.
         right%static(3:9)=.TRUE.

         bottom%bc_type(1:4)="zeroflux"
         bottom%bc_type(7:9)="zeroflux"
         bottom%static(5)=.TRUE.
         bottom%static(6)=.TRUE. 

      CASE("Chen-Shibata") !boundary condition type
         top%static(1)=.TRUE.
         right%static(1)=.TRUE.
         left%static(1)=.TRUE.
         top%static(3:9)=.TRUE.
         right%static(3:9)=.TRUE.
         left%static(3:9)=.TRUE.
         bottom%static(1:9)=.TRUE.
      CASE("Chen-Shibata-hlf") !boundary condition type
         edge_order=(/4,1,2,3/)

         bottom%static=.TRUE.
         bottom%bc_type(1)="natural"         

         left%bc_type(1:3)="zeroflux" ! symmetric
         left%bc_type(5)="zeroflux"   ! symmetric
         left%bc_type(7:9)="zeroflux" ! only x,z-momentum asymmetric
         left%static(4)=.TRUE.
         left%static(6)=.TRUE.

         top%static(1)=.TRUE.
         top%static(3:9)=.TRUE.

         right%static(1)=.TRUE.
         right%static(3:9)=.TRUE.
      CASE("CS-stratified") !boundary condition type
         edge_order=(/1,4,2,3/)

         bottom%static(1)=.FALSE.
         bottom%static(2:9)=.TRUE.

         left%bc_type(1:3)="zeroflux" ! symmetric
         left%bc_type(5)="zeroflux"   ! symmetric
         left%bc_type(7:9)="zeroflux" ! only x,z-momentum asymmetric
         left%static(4)=.TRUE.
         left%static(6)=.TRUE.

         top%static(1)=.TRUE.
         top%static(3:9)=.TRUE.

         right%static(1)=.TRUE.
         right%static(3:9)=.TRUE.
      CASE("CurrentSheet") !boundary condition type
         edge_order=(/1,4,2,3/)

         bottom%static(1:2)=.FALSE.
         bottom%static(3:9)=.TRUE.

         left%bc_type(1:3)="zeroflux" ! symmetric
         left%bc_type(5)="zeroflux"   ! symmetric
         left%bc_type(7:9)="zeroflux" ! only x,z-momentum asymmetric
         left%static(4)=.TRUE.
         left%static(6)=.TRUE.

         top%static(1)=.TRUE.
         top%static(3:9)=.TRUE.

         right%static(1:2)=.FALSE.
         right%static(3:9)=.TRUE.

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
c Elena: nhat - normal      
      USE SolarLnMHD_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      INTEGER :: i
      REAL(r8) :: tfrac
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci,psi_bc,
     $     rho0,pt0,px0,py0,Tx0,Ty0
c-----------------------------------------------------------------------
c     cylindrical to cartesian relationships:
c     1: z --> x
c     2: r --> y
c     3: phi --> z
c-----------------------------------------------------------------------
c     variables: n, psi,  bz, nvx, nvy, nvz,  jz,   p
c                1,   2,   3,   4,   5,   6,   7,   8
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
c Elena begin      
      CASE("CurrentSheet-hlf")
         SELECT CASE(lrtb)
         CASE("left")
            c(4,:,:)=u(4,:,:)                                    !0=n*vx
            c(6,:,:)=u(6,:,:)   
         CASE("right")                               
            c(3,:,:)=nhat(1,:,:)*ux(3,:,:)+nhat(2,:,:)*uy(3,:,:) !0=grad_n(u)
            c(4,:,:)=u(4,:,:)                                    !0=V_x
            c(7,:,:)=u(7,:,:)
            DO i=5,6                                  
               c(i,:,:)=nhat(1,:,:)*(ux(i,:,:)-u(i,:,:)*ux(1,:,:))
     $         +nhat(2,:,:)*(uy(i,:,:)-u(i,:,:)*uy(1,:,:))              
            ENDDO      
c thermal isolator gradT_e = 0. gradT_i =0.           
            DO i=8,9
               c(i,:,:)=nhat(1,:,:)*(ux(i,:,:)-u(i,:,:)*ux(1,:,:))
     $         +nhat(2,:,:)*(uy(i,:,:)-u(i,:,:)*uy(1,:,:)) 
            ENDDO 
         END SELECT     
c Elena end             
      CASE("TwoFR")
         SELECT CASE(lrtb)
         CASE("left", "right")
            c(4,:,:)=u(4,:,:)    !0=n*vx
            c(6,:,:)=u(6,:,:)    !0=n*vz
         CASE("top")
            c(3,:,:)=nhat(1,:,:)*ux(3,:,:)+nhat(2,:,:)*uy(3,:,:) !0=grad_n(bz)
c            c(4,:,:)=nhat(1,:,:)*(ux(4,:,:)-u(4,:,:)*ux(1,:,:))
c     $         +nhat(2,:,:)*(uy(4,:,:)-u(4,:,:)*uy(1,:,:))
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
c            c(6,:,:)=nhat(1,:,:)*(ux(6,:,:)-u(6,:,:)*ux(1,:,:))
c     $         +nhat(2,:,:)*(uy(6,:,:)-u(6,:,:)*uy(1,:,:))
            c(7,:,:)=u(7,:,:) 
c thermal isolator gradT_e = 0. gradT_i =0.           
            DO i=8,9
               c(i,:,:)=nhat(1,:,:)*(ux(i,:,:)-u(i,:,:)*ux(1,:,:))
     $         +nhat(2,:,:)*(uy(i,:,:)-u(i,:,:)*uy(1,:,:))
            ENDDO
         CASE("bottom")
            c(5,:,:)=u(5,:,:)    !0=n*vy
            c(6,:,:)=u(6,:,:)    !0=n*vz               
      END SELECT

        CASE("MAST")
         SELECT CASE(lrtb)
         CASE("top")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:) !0=grad_n(u)
            c(3,:,:)=nhat(1,:,:)*ux(3,:,:)+nhat(2,:,:)*uy(3,:,:) !0=grad_n(bz)
            c(4,:,:)=uy(4,:,:)-u(4,:,:)*uy(1,:,:)     !d(vx)/dy=0 
            c(5,:,:)=u(5,:,:)                         !vy*u1=0
            c(6,:,:)=uy(6,:,:)-u(6,:,:)*uy(1,:,:)     ! d(vz)/dy=0
            c(7,:,:)=u(7,:,:) !0=jz
c thermal isolator gradT_e = 0. gradT_i =0.           
            DO i=8,9
               c(i,:,:)=nhat(1,:,:)*(ux(i,:,:)-u(i,:,:)*ux(1,:,:))
     $         +nhat(2,:,:)*(uy(i,:,:)-u(i,:,:)*uy(1,:,:))
            ENDDO
         CASE("left")
            c(4,:,:)=u(4,:,:)    !0=n*vx
            c(6,:,:)=u(6,:,:)    !0=n*vz
         CASE("right")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:) !0=grad_n(u)
            c(3,:,:)=nhat(1,:,:)*ux(3,:,:)+nhat(2,:,:)*uy(3,:,:)    !grad_n(bz)=0
            c(4,:,:)=u(4,:,:)  !u1*vx=0
            c(5,:,:)=ux(5,:,:)-u(5,:,:)*ux(1,:,:) !0=d/dx(vy)
            c(6,:,:)=ux(6,:,:)-u(6,:,:)*ux(1,:,:) !0=d/dx(vz)
            c(7,:,:)=u(7,:,:)    !jz=0
c thermal isolator gradT_e = 0. gradT_i =0.            
            DO i=8,9
               c(i,:,:)=nhat(1,:,:)*(ux(i,:,:)-u(i,:,:)*ux(1,:,:))
     $         +nhat(2,:,:)*(uy(i,:,:)-u(i,:,:)*uy(1,:,:))
            ENDDO
         CASE("bottom")
            c(5,:,:)=u(5,:,:)    !0=n*vy
            c(6,:,:)=u(6,:,:)    !0=n*vz               
      END SELECT

      CASE("Chen-Shibata") !BC equations for rhs
         SELECT CASE(lrtb)
         CASE("left","right","top")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:) ! 0=grad_n(u)
            DO i=3,9
              c(i,:,:)=nhat(1,:,:)*ux(i,:,:)+nhat(2,:,:)*uy(i,:,:)
            ENDDO
         CASE("bottom")
            psi_bc=SolarLnMHD_psi(x,y)
            IF(t.LE.t_e)THEN
              tfrac=t/t_e
            ELSE
              tfrac=1.
            ENDIF
            WHERE(ABS(x-x0).LE.0.3)
                psi_bc=psi_bc-c_psi_e*half*(one+COS(pi*(x-x0)/0.3))
     $            *half*(tfrac+1./twopi*SIN(twopi*(tfrac-half)))
            END WHERE
            c(1,:,:)=uy(1,:,:)              !0=dy(ln(n))
            c(2,:,:)=u(2,:,:)-psi_bc          !0=Psi-Psi_0-Psi_e*t/t_e
            c(3,:,:)=uy(3,:,:)              !0=dy(B_z)
            c(4,:,:)=u(4,:,:)               !0=n*v_x
            c(5,:,:)=uy(5,:,:)-u(5,:,:)*uy(1,:,:) !0=n*dy(vy)=dy(n*vy)-n*vy*dy(ln(n))
            c(6,:,:)=u(6,:,:)               !0=n*v_z
            c(7,:,:)=uy(7,:,:)               !0=dy(j_z)
            c(8,:,:)=uy(8,:,:)-u(8,:,:)*uy(1,:,:) !0=dy(T)/n=dy(p)-pdy(ln(n))
            c(9,:,:)=uy(9,:,:)-u(9,:,:)*uy(1,:,:) !0=dy(T)/n=dy(p)-pdy(ln(n))
         END SELECT
      CASE("Chen-Shibata-hlf") !BC equations for rhs
         SELECT CASE(lrtb)
         CASE("left")
            c(4,:,:)=u(4,:,:)                                    ! 0=n*vx
            c(6,:,:)=u(6,:,:)                                    ! 0=n*vz
         CASE("right","top")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:) ! 0=grad_n(u)
            c(3,:,:)=nhat(1,:,:)*ux(3,:,:)+nhat(2,:,:)*uy(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=nhat(1,:,:)*ux(3,:,:)+nhat(2,:,:)*uy(3,:,:)
            c(6,:,:)=u(6,:,:)
            DO i=7,9
              c(i,:,:)=nhat(1,:,:)*ux(i,:,:)+nhat(2,:,:)*uy(i,:,:)
            ENDDO
         CASE("bottom")
            psi_bc=SolarLnMHD_psi(x,y)
            IF(t.LE.t_e)THEN
              tfrac=t/t_e
            ELSE
              tfrac=1.
            ENDIF
            WHERE(ABS(x-x0).LE.0.3)
                psi_bc=psi_bc-c_psi_e*half*(one+COS(pi*(x-x0)/0.3))
     $            *half*(tfrac+1./twopi*SIN(twopi*(tfrac-half)))
            END WHERE
c            c(1,:,:)=uy(1,:,:)              !0=dy(ln(n))
            c(2,:,:)=u(2,:,:)-psi_bc        !0=Psi-Psi_0-Psi_e*t/t_e
            c(3,:,:)=uy(3,:,:)              !0=dy(B_z)
            c(4,:,:)=u(4,:,:)               !0=n*v_x
            c(5,:,:)=u(5,:,:)               !0=n*v_y
            c(6,:,:)=u(6,:,:)               !0=n*v_z
            c(7,:,:)=uy(7,:,:)              !0=dy(j_z)
            c(8,:,:)=uy(8,:,:)-u(8,:,:)*uy(1,:,:) !0=dy(T)/n=dy(p)-pdy(ln(n))
            c(9,:,:)=uy(9,:,:)-u(9,:,:)*uy(1,:,:) !0=dy(T)/n=dy(p)-pdy(ln(n))
         END SELECT
      CASE("CS-stratified") !BC equations for rhs
         SELECT CASE(lrtb)
         CASE("left")
            c(4,:,:)=u(4,:,:)                                  ! 0=n*vx
            c(6,:,:)=u(6,:,:)                                  ! 0=n*vz
         CASE("right","top")
            c(1,:,:)=nhat(1,:,:)*ux(1,:,:)+nhat(2,:,:)*uy(1,:,:) ! 0=grad_n(u)
            c(3,:,:)=nhat(1,:,:)*ux(3,:,:)+nhat(2,:,:)*uy(3,:,:) ! 0=grad_n(u)
            c(4,:,:)=u(4,:,:)                                  !0=n*v_x 
            c(5,:,:)=u(5,:,:)                                  !0=n*v_y
            c(6,:,:)=u(6,:,:)
            DO i=7,9
              c(i,:,:)=nhat(1,:,:)*ux(i,:,:)+nhat(2,:,:)*uy(i,:,:)
            ENDDO
         CASE("bottom")
            psi_bc=SolarLnMHD_psi(x,y)
            IF(t.LE.t_e)THEN
              tfrac=t/t_e
            ELSE
              tfrac=1.
            ENDIF
            WHERE(ABS(x-x0).LE.0.3)
                psi_bc=psi_bc-c_psi_e*half*(one+COS(pi*(x-x0)/0.3))
     $            *half*(tfrac+1./twopi*sin(twopi*(tfrac-half)))
            END WHERE
            c(2,:,:)=u(2,:,:)-psi_bc        !0=Psi-(Psi_0-Psi_e*t/t_e)
            c(3,:,:)=uy(3,:,:)              !0=dy(B_z)
            c(4,:,:)=uy(4,:,:)-u(4,:,:)*uy(1,:,:) !0=dy(v_x)   
            c(5,:,:)=u(5,:,:)              !0=n*v_y
            c(6,:,:)=uy(6,:,:)-u(6,:,:)*uy(1,:,:) !0=dy(v_z)
            c(7,:,:)=uy(7,:,:)              !0=dy(j_z)

            CALL SolarLnMHD_atmsph(x,y,rho0,pt0,px0,py0,Tx0,Ty0)
            c(8,:,:)=uy(8,:,:)+py0 - (u(8,:,:)+pt0)*uy(1,:,:) 
     $           - Ty0*EXP(u(1,:,:)) !0=dy(delta T)
            c(9,:,:)=uy(9,:,:)+py0 - (u(9,:,:)+pt0)*uy(1,:,:) 
     $           - Ty0*EXP(u(1,:,:)) !0=dy(delta T)
         END SELECT
      CASE("CurrentSheet") !BC equations for rhs
         SELECT CASE(lrtb)
         CASE("left")
            c(4,:,:)=u(4,:,:)                                  ! 0=n*vx
            c(6,:,:)=u(6,:,:)                                  ! 0=n*vz
         CASE("right")
            c(3,:,:)=nhat(1,:,:)*ux(3,:,:)+nhat(2,:,:)*uy(3,:,:) ! 0=grad_n(u)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            DO i=7,9
              c(i,:,:)=nhat(1,:,:)*ux(i,:,:)+nhat(2,:,:)*uy(i,:,:)
            ENDDO
         CASE("top")
            c(1,:,:)=uy(1,:,:) ! 0=grad_n(u)
            c(3,:,:)=uy(3,:,:) ! 0=grad_n(u)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=uy(7,:,:)

            CALL SolarLnMHD_atmsph(x,y,rho0,pt0,px0,py0,Tx0,Ty0)
            c(8,:,:)=uy(8,:,:)+py0 - (u(8,:,:)+pt0)*uy(1,:,:) 
     $           - Ty0*EXP(u(1,:,:)) !0=dy(delta T)
            c(9,:,:)=uy(9,:,:)+py0 - (u(9,:,:)+pt0)*uy(1,:,:) 
     $           - Ty0*EXP(u(1,:,:)) !0=dy(delta T)
         CASE("bottom")
            c(3,:,:)=uy(3,:,:)              !0=dy(B_z)
            c(4,:,:)=uy(4,:,:)-u(4,:,:)*uy(1,:,:) !0=dy(v_x)
            c(5,:,:)=u(5,:,:)              !0=n*v_y
            c(6,:,:)=uy(6,:,:)-u(6,:,:)*uy(1,:,:) !0=dy(v_z)
            c(7,:,:)=uy(7,:,:)              !0=dy(j_z)

            CALL SolarLnMHD_atmsph(x,y,rho0,pt0,px0,py0,Tx0,Ty0)
            c(8,:,:)=uy(8,:,:)+py0 - (u(8,:,:)+pt0)*uy(1,:,:) 
     $           - Ty0*EXP(u(1,:,:)) !0=dy(delta T)
            c(9,:,:)=uy(9,:,:)+py0 - (u(9,:,:)+pt0)*uy(1,:,:) 
     $           - Ty0*EXP(u(1,:,:)) !0=dy(delta T)
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
      USE SolarLnMHD_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci,
     $     rho0,pt0,px0,py0,Tx0,Ty0
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
c Elena begin      
      CASE("CurrentSheet-hlf") 
         SELECT CASE(lrtb)
         CASE("left")
            c_u(4,4,:,:)=one
            c_u(6,6,:,:)=one   
         CASE("right")    
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)                           
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(5,5,:,:)=nhat(1,:,:)
            c_uy(5,5,:,:)=nhat(2,:,:)
            c_ux(5,1,:,:)=-nhat(1,:,:)*u(5,:,:)
            c_uy(5,1,:,:)=-nhat(2,:,:)*u(5,:,:)
            c_u(7,7,:,:)=one                                  
            c_u(6,6,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(6,6,:,:)=nhat(1,:,:)
            c_uy(6,6,:,:)=nhat(2,:,:)
            c_ux(6,1,:,:)=-nhat(1,:,:)*u(6,:,:)
            c_uy(6,1,:,:)=-nhat(2,:,:)*u(6,:,:)
            c_u(8,8,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(8,8,:,:)=nhat(1,:,:)
            c_uy(8,8,:,:)=nhat(2,:,:)
            c_ux(8,1,:,:)=-nhat(1,:,:)*u(8,:,:)
            c_uy(8,1,:,:)=-nhat(2,:,:)*u(8,:,:)   
            c_u(9,9,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(9,9,:,:)=nhat(1,:,:)
            c_uy(9,9,:,:)=nhat(2,:,:)
            c_ux(9,1,:,:)=-nhat(1,:,:)*u(9,:,:)
            c_uy(9,1,:,:)=-nhat(2,:,:)*u(9,:,:)   
         END SELECT     
c Elena end
      CASE("TwoFR")
         SELECT CASE(lrtb)
         CASE("left","right")
            c_u(4,4,:,:)=one
            c_u(6,6,:,:)=one
         CASE("top")
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)
c            c_u(4,4,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
c            c_ux(4,4,:,:)=nhat(1,:,:)
c            c_uy(4,4,:,:)=nhat(2,:,:)
c            c_ux(4,1,:,:)=-nhat(1,:,:)*u(4,:,:)
c            c_uy(4,1,:,:)=-nhat(2,:,:)*u(4,:,:) 
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
c            c_u(6,6,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
c            c_ux(6,6,:,:)=nhat(1,:,:)
c            c_uy(6,6,:,:)=nhat(2,:,:)
c            c_ux(6,6,:,:)=-nhat(1,:,:)*u(6,:,:)
c            c_uy(6,6,:,:)=-nhat(2,:,:)*u(6,:,:)
            c_u(7,7,:,:)=one
            c_u(8,8,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(8,8,:,:)=nhat(1,:,:)
            c_uy(8,8,:,:)=nhat(2,:,:)
            c_ux(8,1,:,:)=-nhat(1,:,:)*u(8,:,:)
            c_uy(8,1,:,:)=-nhat(2,:,:)*u(8,:,:)
            c_u(9,9,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(9,9,:,:)=nhat(1,:,:)
            c_uy(9,9,:,:)=nhat(2,:,:)
            c_ux(9,1,:,:)=-nhat(1,:,:)*u(9,:,:)
            c_uy(9,1,:,:)=-nhat(2,:,:)*u(9,:,:)
         CASE("bottom")
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
         END SELECT

        CASE("MAST")
         SELECT CASE(lrtb)
         CASE("left")
            c_u(4,4,:,:)=one
            c_u(6,6,:,:)=one
         CASE("right")
            c_ux(1,1,:,:)=nhat(1,:,:)
            c_uy(1,1,:,:)=nhat(2,:,:) 
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=-ux(1,:,:)
            c_ux(5,1,:,:)=-u(5,:,:)
            c_ux(5,5,:,:)=one
            c_u(6,6,:,:)=-ux(1,:,:)
            c_ux(6,1,:,:)=-u(6,:,:)
            c_ux(6,6,:,:)=one
            c_u(7,7,:,:)=one
            c_u(8,8,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(8,8,:,:)=nhat(1,:,:)
            c_uy(8,8,:,:)=nhat(2,:,:)
            c_ux(8,1,:,:)=-nhat(1,:,:)*u(8,:,:)
            c_uy(8,1,:,:)=-nhat(2,:,:)*u(8,:,:)
            c_u(9,9,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(9,9,:,:)=nhat(1,:,:)
            c_uy(9,9,:,:)=nhat(2,:,:)
            c_ux(9,1,:,:)=-nhat(1,:,:)*u(9,:,:)
            c_uy(9,1,:,:)=-nhat(2,:,:)*u(9,:,:)
         CASE("top")
            c_ux(1,1,:,:)=nhat(1,:,:)
            c_uy(1,1,:,:)=nhat(2,:,:)
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)
            c_u(4,4,:,:)=-uy(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)
            c_uy(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=-uy(1,:,:)
            c_uy(6,1,:,:)=-u(6,:,:)
            c_uy(6,6,:,:)=one
            c_u(7,7,:,:)=one
            c_u(8,8,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(8,8,:,:)=nhat(1,:,:)
            c_uy(8,8,:,:)=nhat(2,:,:)
            c_ux(8,1,:,:)=-nhat(1,:,:)*u(8,:,:)
            c_uy(8,1,:,:)=-nhat(2,:,:)*u(8,:,:)
            c_u(9,9,:,:)=-nhat(1,:,:)*ux(1,:,:) - nhat(2,:,:)*uy(1,:,:)
            c_ux(9,9,:,:)=nhat(1,:,:)
            c_uy(9,9,:,:)=nhat(2,:,:)
            c_ux(9,1,:,:)=-nhat(1,:,:)*u(9,:,:)
            c_uy(9,1,:,:)=-nhat(2,:,:)*u(9,:,:)
         CASE("bottom")
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
         END SELECT

      CASE("Chen-Shibata") !derivatives of BC equations
         SELECT CASE(lrtb)
         CASE("left","right","top")
            c_ux(1,1,:,:)=nhat(1,:,:)
            c_uy(1,1,:,:)=nhat(2,:,:)
            DO i=3,9
              c_ux(i,i,:,:)=nhat(1,:,:)
              c_uy(i,i,:,:)=nhat(2,:,:)
            ENDDO
         CASE("bottom")
            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_uy(5,5,:,:)=one
            c_u(5,5,:,:)=-uy(1,:,:)
            c_uy(5,1,:,:)=-u(5,:,:)
            c_u(6,6,:,:)=one
            c_uy(7,7,:,:)=one
            c_u(8,8,:,:)=-uy(1,:,:)
            c_uy(8,1,:,:)=-u(8,:,:)
            c_uy(8,8,:,:)=one
            c_u(9,9,:,:)=-uy(1,:,:)
            c_uy(9,1,:,:)=-u(9,:,:)
            c_uy(9,9,:,:)=one
         END SELECT
      CASE("Chen-Shibata-hlf") !derivatives of BC equations
         SELECT CASE(lrtb)
         CASE("left")
            c_u(4,4,:,:)=one
            c_u(6,6,:,:)=one
         CASE("right","top")
            c_ux(1,1,:,:)=nhat(1,:,:)
            c_uy(1,1,:,:)=nhat(2,:,:)
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)
            c_u(4,4,:,:)=one
            c_ux(5,5,:,:)=nhat(1,:,:)
            c_uy(5,5,:,:)=nhat(2,:,:)
            c_u(6,6,:,:)=one
            DO i=7,9
              c_ux(i,i,:,:)=nhat(1,:,:)
              c_uy(i,i,:,:)=nhat(2,:,:)
            ENDDO
         CASE("bottom")
c            c_uy(1,1,:,:)=one
            c_u(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
            c_uy(7,7,:,:)=one
            c_u(8,8,:,:)=-uy(1,:,:)
            c_uy(8,1,:,:)=-u(8,:,:)
            c_uy(8,8,:,:)=one
            c_u(9,9,:,:)=-uy(1,:,:)
            c_uy(9,1,:,:)=-u(9,:,:)
            c_uy(9,9,:,:)=one
         END SELECT
      CASE("CS-stratified") !derivatives of BC equations
         SELECT CASE(lrtb)
         CASE("left")
            c_u(4,4,:,:)=one
            c_u(6,6,:,:)=one
         CASE("right","top")
            c_ux(1,1,:,:)=nhat(1,:,:)
            c_uy(1,1,:,:)=nhat(2,:,:)
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
            DO i=7,9
              c_ux(i,i,:,:)=nhat(1,:,:)
              c_uy(i,i,:,:)=nhat(2,:,:)
            ENDDO
         CASE("bottom")
            c_u(2,2,:,:)=one
            c_uy(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_u(4,4,:,:)=-uy(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)            
            c_u(5,5,:,:)=one
            c_uy(6,6,:,:)=one
            c_u(6,6,:,:)=-uy(1,:,:)
            c_uy(6,1,:,:)=-u(6,:,:)
            c_uy(7,7,:,:)=one

            CALL SolarLnMHD_atmsph(x,y,rho0,pt0,px0,py0,Tx0,Ty0)
            c_u(8,1,:,:)=-Ty0*EXP(u(1,:,:))
            c_u(8,8,:,:)=-uy(1,:,:)
            c_uy(8,1,:,:)=-(u(8,:,:)+pt0)
            c_uy(8,8,:,:)=one
            c_u(9,1,:,:)=-Ty0*EXP(u(1,:,:))
            c_u(9,9,:,:)=-uy(1,:,:)
            c_uy(9,1,:,:)=-(u(9,:,:)+pt0)
            c_uy(9,9,:,:)=one
         END SELECT
      CASE("CurrentSheet") !derivatives of BC equations
         SELECT CASE(lrtb)
         CASE("left")
            c_u(4,4,:,:)=one
            c_u(6,6,:,:)=one
         CASE("right")
            c_ux(3,3,:,:)=nhat(1,:,:)
            c_uy(3,3,:,:)=nhat(2,:,:)
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
            DO i=7,9
              c_ux(i,i,:,:)=nhat(1,:,:)
              c_uy(i,i,:,:)=nhat(2,:,:)
            ENDDO
         CASE("top")
            c_uy(1,1,:,:)=one
            c_uy(3,3,:,:)=one
            c_u(4,4,:,:)=one
            c_u(5,5,:,:)=one
            c_u(6,6,:,:)=one
            c_uy(7,7,:,:)=one

            CALL SolarLnMHD_atmsph(x,y,rho0,pt0,px0,py0,Tx0,Ty0)
            c_u(8,1,:,:)=-Ty0*EXP(u(1,:,:))
            c_u(8,8,:,:)=-uy(1,:,:)
            c_uy(8,1,:,:)=-(u(8,:,:)+pt0)
            c_uy(8,8,:,:)=one
            c_u(9,1,:,:)=-Ty0*EXP(u(1,:,:))
            c_u(9,9,:,:)=-uy(1,:,:)
            c_uy(9,1,:,:)=-(u(9,:,:)+pt0)
            c_uy(9,9,:,:)=one
         CASE("bottom")
            c_uy(3,3,:,:)=one
            c_uy(4,4,:,:)=one
            c_u(4,4,:,:)=-uy(1,:,:)
            c_uy(4,1,:,:)=-u(4,:,:)            
            c_u(5,5,:,:)=one
            c_uy(6,6,:,:)=one
            c_u(6,6,:,:)=-uy(1,:,:)
            c_uy(6,1,:,:)=-u(6,:,:)
            c_uy(7,7,:,:)=one

            CALL SolarLnMHD_atmsph(x,y,rho0,pt0,px0,py0,Tx0,Ty0)
            c_u(8,1,:,:)=-Ty0*EXP(u(1,:,:))
            c_u(8,8,:,:)=-uy(1,:,:)
            c_uy(8,1,:,:)=-(u(8,:,:)+pt0)
            c_uy(8,8,:,:)=one
            c_u(9,1,:,:)=-Ty0*EXP(u(1,:,:))
            c_u(9,9,:,:)=-uy(1,:,:)
            c_uy(9,1,:,:)=-(u(9,:,:)+pt0)
            c_uy(9,9,:,:)=one
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
      USE SolarLnMHD_mod
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
c Elena begin
      CASE("CurrentSheet-hlf")
         SELECT CASE(lrtb)
         CASE("right")
            mass(2,2,:,:)=one
         END SELECT
c Elena end      
      CASE("TwoFR")
         SELECT CASE(lrtb)
         CASE("top")
            mass(2,2,:,:)=one
         END SELECT
      CASE("MAST")
         SELECT CASE(lrtb)
         CASE("top", "right")
            mass(2,2,:,:)=one
         END SELECT
      CASE("Chen-Shibata") !lhs
         SELECT CASE(lrtb)
         CASE("bottom")
         CASE("top")
            mass_x(2,2,:,:)=nhat(1,:,:)      !nhat.grad(Psi)=...=0
            mass_y(2,2,:,:)=nhat(2,:,:)      !nhat.grad(Psi)=...=0
         CASE("left")
            mass_x(2,2,:,:)=nhat(1,:,:)
            mass_y(2,2,:,:)=nhat(2,:,:)
         CASE("right")
            mass_x(2,2,:,:)=nhat(1,:,:)
            mass_y(2,2,:,:)=nhat(2,:,:)
         END SELECT
      CASE("Chen-Shibata-hlf") !lhs
         SELECT CASE(lrtb)
         CASE("bottom")
         CASE("top")
            mass_x(2,2,:,:)=nhat(1,:,:)      !nhat.grad(Psi)=...=0
            mass_y(2,2,:,:)=nhat(2,:,:)      !nhat.grad(Psi)=...=0
         CASE("right")
            mass_x(2,2,:,:)=nhat(1,:,:)
            mass_y(2,2,:,:)=nhat(2,:,:)
         END SELECT
      CASE("CS-stratified") !lhs
         SELECT CASE(lrtb)
         CASE("bottom")
            mass_y(1,1,:,:)=one            !dt(dy(ln n))=...=0
         CASE("top")
            mass_x(2,2,:,:)=nhat(1,:,:)      !nhat.grad(Psi)=...=0
            mass_y(2,2,:,:)=nhat(2,:,:)      !nhat.grad(Psi)=...=0
         CASE("right")
            mass_x(2,2,:,:)=nhat(1,:,:)
            mass_y(2,2,:,:)=nhat(2,:,:)
         END SELECT
      CASE("CurrentSheet") !lhs
         SELECT CASE(lrtb)
         CASE("bottom")
            mass_y(1,1,:,:)=one            !dt(dy(ln n))=...=0
            mass_y(2,2,:,:)=one
         CASE("top")
            mass_y(2,2,:,:)=one      !nhat.grad(Psi)=...=0
         CASE("right")
            mass_x(1,1,:,:)=nhat(1,:,:)
            mass_y(1,1,:,:)=nhat(2,:,:)
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
      USE SolarLnMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Tix,Tiy,Tex,Tey,
     $     kperpi,kfaci,kperpe,kface,heat_exch,
     $     r_fac,r_faci,j1,j2,j3,jtot,eta_local,b1,b2,Bsq,rho0,pt0,px0,
     $     py0,Tx0,Ty0,rho,rho_inv,rhox_inv,rhoy_inv,gx_vec,gy_vec,
     $     r_sphr,mu_local, rad_loss, alpha, ddiff_local, 
     $     heat_blnc_radlos
      REAL(r8), DIMENSION(SIZE(u,2), SIZE(u,3)) :: rho_init, p_init
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: vi,vix,viy,BdotTi,
     $     BdotTe
      REAL(r8), DIMENSION(4,SIZE(u,2),SIZE(u,3)) :: coeff
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
      SELECT CASE(init_type)
      CASE DEFAULT
         gx_vec=0.
         gy_vec=gravity
      END SELECT
c-----------------------------------------------------------------------
c     IF(init_type==CS-stratified), get stratified atmosphere
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("CS-stratified","CurrentSheet")
         CALL SolarLnMHD_atmsph(x,y,rho0,pt0,px0,py0,Tx0,Ty0)      
      CASE DEFAULT
         rho0=zero
         pt0=zero
         px0=zero
         py0=zero
         Tx0=zero
         Ty0=zero
      END SELECT
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
c Elena: (p_e+p_0) = n*(T_e+T_0)
c Elena: (p_i+p_0) = n*(T_i+T_0)
      Tex = rho_inv*(ux(8,:,:)+px0 - (u(8,:,:)+pt0)*ux(1,:,:)) - Tx0
      Tey = rho_inv*(uy(8,:,:)+py0 - (u(8,:,:)+pt0)*uy(1,:,:)) - Ty0
      Tix = rho_inv*(ux(9,:,:)+px0 - (u(9,:,:)+pt0)*ux(1,:,:)) - Tx0
      Tiy = rho_inv*(uy(9,:,:)+py0 - (u(9,:,:)+pt0)*uy(1,:,:)) - Ty0
c-----------------------------------------------------------------------
c     magnetic fields and currents.
c-----------------------------------------------------------------------
c Elena b1 = -dpsi/dy = B_x b2 = dpsi/dx = B_y
c Elena bz = bz0 + bz
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2 + (u(3,:,:) + bz0)**2
c Components of current j1 = dB_z/dy j2 = - dB_z/dx j3 = u(7)      
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
     $        u(8,:,:)+pt0,u(9,:,:)+pt0,jtot,jc_norm,eta_anm_norm,
     $        0._r8,etas_norm,0._r8,0._r8,0._r8,eta,
     $        etavac,coeff,eta_local)
      CASE("IAT-anomalous")
         CALL transport_seteta("IAT-anomalous",r_sphr,rho,
     $        u(8,:,:)+pt0,u(9,:,:)+pt0,jtot,jc_norm,eta_anm_norm,
     $        0._r8,etas_norm,0._r8,0._r8,0._r8,eta,
     $        etavac, coeff,eta_local)
      CASE("spitzer-chodura")
         CALL transport_seteta(eta_case,r_sphr,rho,u(8,:,:)+pt0,
     $        u(9,:,:)+pt0, jtot,jc_norm,eta_anm_norm,
     $        chod_const,etas_norm,etac_norm,v_chod_norm,0._r8,eta,
     $        etavac,coeff,eta_local)
      CASE("Chen-Shibata")
         CALL ChenShibata_eta(x,y,jtot,eta_vbl,eta_chrsp,eta_local)
      CASE("MAST")
         CALL MAST_eta(x,y,eta_vbl,eta_local)
      CASE DEFAULT
         eta_local=eta
      END SELECT
      
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT(b1,b2,u(3,:,:)+bz0,Tex,Tey,BdotTe)
         CALL transport_BdotT(b1,b2,u(3,:,:)+bz0,Tix,Tiy,BdotTi)
         CALL transport_kbrag(rho,u(8,:,:)+pt0,one,Bsq,ke_norm,
     $        zero,xe_norm,zero,kappa_min,kappa_max,kperpe,kface)
         CALL transport_kbrag(rho,u(9,:,:)+pt0,zero,Bsq,zero,
     $        ki_norm,zero,xi_norm,kappa_min,kappa_max,kperpi,kfaci)
      CASE DEFAULT
         BdotTe=0
         BdotTi=0
         kface=0
         kfaci=0
c In our case we use scalar for heat conduction        
         IF(kappa_case == "scalar")THEN
            kperpe=kappa_prp
            kperpi=kappa_prp
         ELSE
            kperpe=kappa_prp*rho
            kperpi=kappa_prp*rho
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     heat exchange.
c-----------------------------------------------------------------------
      CALL transport_hexch(ieheat,ke_norm,rho,u(9,:,:)+pt0,u(8,:,:)+pt0,
     $     heat_exch)
c-----------------------------------------------------------------------
c     radiative losses
c-----------------------------------------------------------------------

      SELECT CASE(init_type)
      CASE("CurrentSheet-hlf","MAST")
         CALL transport_radloss(rad_fac,T0,rho,u(8,:,:)+pt0,alpha,
     $        rad_loss)
      CASE("Chen-Shibata-hlf")
         rad_loss=0.0
      END SELECT
c-----------------------------------------------------------------------
c     constant heating function based on electron initial density and 
c     temperature profiles
c-----------------------------------------------------------------------

      SELECT CASE(init_type)
      CASE("CurrentSheet-hlf")
c        rho_init = one + c_rho*(one - (TANH(x/h_psi))**2) 
c        p_init = 0.25*beta0*(EXP(u(1,:,:)))**(gamma)
        rho_init = one
        p_init = 0.25*beta0
        CALL transport_radloss(rad_fac,T0,rho_init,p_init,alpha,
     $     heat_blnc_radlos)
      CASE ("Chen-Shibata-hlf")
         heat_blnc_radlos=0.0    
      CASE DEFAULT
        rho_init = one
        p_init = 0.25*beta0
        CALL transport_radloss(rad_fac,T0,rho_init,p_init,alpha,
     $     heat_blnc_radlos)
      END SELECT
c-----------------------------------------------------------------------
c     viscous boundary layer.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Chen-Shibata-hlf","CS-stratified","CurrentSheet", "MAST")
         CALL ChenShibata_mu(x,y,rho,mu_vbl,mu_local)
      CASE DEFAULT
         mu_local = mu*rho + mu_min
      END SELECT

      SELECT CASE(init_type)
      CASE("MAST")
         CALL MAST_ddiff(x,y,ddiff,ddiff_local)
      CASE DEFAULT
         ddiff_local=ddiff
      END SELECT
c-----------------------------------------------------------------------
c     velocities and their gradients.
c-----------------------------------------------------------------------
c Elena: three components of ion velocity vi(1)=u(4)/n and so on
      DO i=1,3
         vi(i,:,:)=u(i+3,:,:)*rho_inv
         vix(i,:,:)=ux(i+3,:,:)*rho_inv - vi(i,:,:)*ux(1,:,:)
         viy(i,:,:)=uy(i+3,:,:)*rho_inv - vi(i,:,:)*uy(1,:,:)
      ENDDO
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
c Elena r_fac=1 for Cartesian
      fx(1,:,:) = r_fac*(vi(1,:,:) - ddiff_local*ux(1,:,:))
      fy(1,:,:) = r_fac*(vi(2,:,:) - ddiff_local*uy(1,:,:))
      s(1,:,:) = -r_fac*(vi(1,:,:)*ux(1,:,:) + vi(2,:,:)*uy(1,:,:)  
     $           - ddiff_local*ux(1,:,:)*ux(1,:,:) 
     $           - ddiff_local*uy(1,:,:)*uy(1,:,:))
c-----------------------------------------------------------------------
c     poloidal magnetic flux equation.
c-----------------------------------------------------------------------
      fx(2,:,:) = r_fac*hyper_eta*ux(7,:,:)
      fy(2,:,:) = r_fac*hyper_eta*uy(7,:,:)
      
      s(2,:,:) = r_fac*(vi(2,:,:)*b1 - vi(1,:,:)*b2
     $     + eta_local*j3) + cyl_fac*r_faci*hyper_eta*u(7,:,:)
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
c B_z = bz0+u(3)
      fx(3,:,:) = (u(3,:,:)+bz0)*vi(1,:,:) - vi(3,:,:)*b1
     $     + eta_local*j2

      fy(3,:,:) = (u(3,:,:)+bz0)*vi(2,:,:) - vi(3,:,:)*b2
     $     - eta_local*j1
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(4,:,:) = r_fac*(u(4,:,:)*vi(1,:,:)
     $     - two*mu_local*vix(1,:,:) 
     $     + half*u(3,:,:)**2 + u(3,:,:)*bz0 + u(8,:,:) + u(9,:,:))
         
      fy(4,:,:) = r_fac*(u(4,:,:)*vi(2,:,:) 
     $     - mu_local*(viy(1,:,:) + vix(2,:,:)))

      s(4,:,:) = -r_fac*(j3*b2 + gx_vec*(rho-rho0))
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(5,:,:) = r_fac*(u(5,:,:)*vi(1,:,:) 
     $     - mu_local*(vix(2,:,:) + viy(1,:,:)))
      
      fy(5,:,:) = r_fac*(u(5,:,:)*vi(2,:,:) + half*u(3,:,:)**2
     $     + u(3,:,:)*bz0 - two*mu_local*viy(2,:,:))

c Elena: why uy(8) and uy(9) goes to source     
      
      s(5,:,:) = r_fac*(j3*b1 - uy(8,:,:) - uy(9,:,:) 
     $     - gy_vec*(rho-rho0)) + cyl_fac*(u(6,:,:)*vi(3,:,:) 
     $     - half*u(3,:,:)**2 - u(3,:,:)*bz0 
     $     - two*mu_local*vi(2,:,:)*r_faci)
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx(6,:,:) = r_fac*(u(6,:,:)*vi(1,:,:) - mu_local*vix(3,:,:))
      
      fy(6,:,:) = r_fac*(u(6,:,:)*vi(2,:,:) - mu_local*viy(3,:,:))

      s(6,:,:) = r_fac*(j1*b2 - j2*b1) - cyl_fac*(u(6,:,:)*vi(2,:,:) 
     $     + r_faci*(mu_local*vi(3,:,:)))
c-----------------------------------------------------------------------
c     out-of-plane electron momentum (current) equation.
c-----------------------------------------------------------------------
      fx(7,:,:)=b2
      fy(7,:,:)=-b1
      s(7,:,:)=j3
c-----------------------------------------------------------------------
c     electron pressure equation.
c-----------------------------------------------------------------------
      IF(eta_case == "Chen-Shibata" .AND. (init_type == "CurrentSheet"
     $   .OR. init_type == "Chen-Shibata-hlf"))
     $     CALL ChenShibata_eta(x,y,jtot,zero,eta_chrsp,eta_local)
  
      IF(eta_case == "MAST") CALL MAST_eta(x,y,zero,eta_local) 

      fx(8,:,:)=r_fac*(gamma_fac*(u(8,:,:)+pt0)*vi(1,:,:)
     $     - kface*BdotTe(1,:,:) - kperpe*Tex)
      fy(8,:,:)=r_fac*(gamma_fac*(u(8,:,:)+pt0)*vi(2,:,:)
     $     - kface*BdotTe(2,:,:) - kperpe*Tey)
      s(8,:,:)=r_fac*(vi(1,:,:)*(ux(8,:,:)+px0) 
     $     + vi(2,:,:)*(uy(8,:,:)+py0) + eta_local*jtot**2 
     $     + hyper_eta*(ux(7,:,:)**2 + uy(7,:,:)**2) - heat_exch
     $     - rad_loss + heat_blnc_radlos)
c-----------------------------------------------------------------------
c     ion pressure equation.
c-----------------------------------------------------------------------
      IF(init_type == "CurrentSheet" .OR. init_type=="Chen-Shibata-hlf"
     $   .OR. init_type == "MAST")
     $     CALL ChenShibata_mu(x,y,rho,zero,mu_local)

      fx(9,:,:)=r_fac*(gamma_fac*(u(9,:,:)+pt0)*vi(1,:,:)
     $     - kfaci*BdotTi(1,:,:) - kperpi*Tix)
      fy(9,:,:)=r_fac*(gamma_fac*(u(9,:,:)+pt0)*vi(2,:,:)
     $     - kfaci*BdotTi(2,:,:) - kperpi*Tiy)
      s(9,:,:)=r_fac*(vi(1,:,:)*(ux(9,:,:)+px0) 
     $     + vi(2,:,:)*(uy(9,:,:)+py0)
     $     + mu_local*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:)
     $     + vix(3,:,:)**2 + viy(3,:,:)**2) + heat_exch)
     $     + cyl_fac*two*mu_local*vi(2,:,:)**2*r_faci
c-----------------------------------------------------------------------
c     initial equilibrium source term.
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("CurrentSheet-hlf")
            CALL SolarLnMHD_equil(x,y,u0,u0x,u0y,.TRUE.)
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
      USE SolarLnMHD_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci,Tix,Tiy,
     $     Ti_un,Tex,Tey,Te_un,kperpi,kfaci,kperpe,kface,kpari_un,
     $     kpari_p,kperpi_un,kperpi_p,kperpi_bsq,kpare_un,
     $     kpare_p,kperpe_un,kperpe_p,kperpe_bsq,heat_exch,hexch_un,
     $     hexch_pi,hexch_pe,j1,j2,j3,jtot,b1,b2,Bsq,rho0,pt0,px0,py0,
     $     Tx0,Ty0,eta_local,eta_rho,eta_p,eta_j,j_u3,j_ux3,j_uy3,j_u7,
     $     rho,rho_inv,rhox_inv,rhoy_inv,gx_vec,gy_vec,r_sphr,mu_local,
     $     rad_loss, rloss_un, rloss_pe, alpha, eta_pi, ddiff_local
      REAL(r8), DIMENSION(3,SIZE(x,1),SIZE(x,2)) :: vi,vix,viy,
     $     BdotTe,BdotTe_b1,BdotTe_b2,BdotTe_b3,BdotTe_Tx,BdotTe_Ty,
     $     BdotTi,BdotTi_b1,BdotTi_b2,BdotTi_b3,BdotTi_Tx,BdotTi_Ty
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
      SELECT CASE(init_type)
      CASE DEFAULT
         gx_vec=0.
         gy_vec=gravity
      END SELECT
c-----------------------------------------------------------------------
c     IF(init_type==CS-stratified), get stratified atmosphere
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("CS-stratified","CurrentSheet")
         CALL SolarLnMHD_atmsph(x,y,rho0,pt0,px0,py0,Tx0,Ty0)
      CASE DEFAULT
         rho0=zero
         pt0=zero
         px0=zero
         py0=zero
         Tx0=zero
         Ty0=zero
      END SELECT
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
      Tex = rho_inv*(ux(8,:,:)+px0 - (u(8,:,:)+pt0)*ux(1,:,:)) - Tx0
      Tey = rho_inv*(uy(8,:,:)+py0 - (u(8,:,:)+pt0)*uy(1,:,:)) - Ty0
      Te_un=-(u(8,:,:)+pt0)*rho_inv
      Tix = rho_inv*(ux(9,:,:)+px0 - (u(9,:,:)+pt0)*ux(1,:,:)) - Tx0
      Tiy = rho_inv*(uy(9,:,:)+py0 - (u(9,:,:)+pt0)*uy(1,:,:)) - Ty0
      Ti_un=-(u(9,:,:)+pt0)*rho_inv
c-----------------------------------------------------------------------
c     magnetic fields and currents.
c-----------------------------------------------------------------------
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2 + (u(3,:,:) + bz0)**2
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
         CALL transport_seteta_u("spitzer-chodura",r_sphr,rho,
     $        u(8,:,:)+pt0,u(9,:,:)+pt0,jtot,jc_norm,eta_anm_norm,
     $        0._r8,etas_norm,0._r8,0._r8,0._r8,eta,
     $        etavac,eta_local,eta_rho,eta_p,eta_pi,eta_j)
         eta_rho = eta_rho*rho
      CASE("spitzer-chodura")
         CALL transport_seteta_u(eta_case,r_sphr,rho,u(8,:,:)+pt0,
     $        u(9,:,:)+pt0,jtot,jc_norm,eta_anm_norm,
     $        chod_const,etas_norm,etac_norm,v_chod_norm,0._r8,
     $        eta,etavac,eta_local,eta_rho,eta_p,eta_pi,eta_j)
         eta_rho = eta_rho*rho
         WHERE(eta_local == etavac)
            eta_rho=0.
            eta_p=0.
            eta_j=0.
         END WHERE
      CASE("Chen-Shibata")
         CALL ChenShibata_eta(x,y,jtot,eta_vbl,eta_chrsp,eta_local)
         eta_rho=0.
         eta_p=0.
         eta_j=0.
         WHERE(jtot.GE.j_c .AND. jtot.LE.(two*j_c))
            eta_j = etavac*half*pi/j_c*SIN(pi*(jtot-j_c)/j_c)
         END WHERE
      CASE("MAST")
         CALL MAST_eta(x,y,eta_vbl,eta_local)
         eta_rho=0.
         eta_p=0.
         eta_j=0.
      CASE("IAT-anomalous")
         CALL transport_seteta_u("IAT-anomalous",r_sphr,rho,
     $        u(8,:,:)+pt0,u(9,:,:)+pt0,jtot,jc_norm,eta_anm_norm,
     $        0._r8,etas_norm,0._r8,0._r8,0._r8,eta,
     $        etavac,eta_local,eta_rho,eta_p,eta_pi,eta_j)
         eta_rho = eta_rho*rho
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
         CALL transport_BdotT_u(b1,b2,u(3,:,:)+bz0,Tex,Tey,BdotTe,
     $        BdotTe_b1,BdotTe_b2,BdotTe_b3,BdotTe_Tx,BdotTe_Ty)
         CALL transport_kbrag_u(rho,u(8,:,:)+pt0,one,Bsq,ke_norm,
     $        zero,xe_norm,zero,kappa_min,kappa_max,kperpe,kface,
     $        kpare_un,kpare_p,kperpe_un,kperpe_p,kperpe_bsq)
         CALL transport_BdotT_u(b1,b2,u(3,:,:)+bz0,Tix,Tiy,BdotTi,
     $        BdotTi_b1,BdotTi_b2,BdotTi_b3,BdotTi_Tx,BdotTi_Ty)
         CALL transport_kbrag_u(rho,u(9,:,:)+pt0,zero,Bsq,zero,
     $        ki_norm,zero,xi_norm,kappa_min,kappa_max,kperpi,kfaci,
     $        kpari_un,kpari_p,kperpi_un,kperpi_p,kperpi_bsq)
         kpare_un = kpare_un*rho
         kperpe_un = kperpe_un*rho
         kpari_un = kpari_un*rho
         kperpi_un = kperpi_un*rho
      CASE DEFAULT
         BdotTe=0
         BdotTe_b1=0
         BdotTe_b2=0
         BdotTe_b3=0
         BdotTe_Tx=0
         BdotTe_Ty=0
         kface=0
         kpare_un=0
         kpare_p=0
         kperpe_p=0
         kperpe_bsq=0
         BdotTi=0
         BdotTi_b1=0
         BdotTi_b2=0
         BdotTi_b3=0
         BdotTi_Tx=0
         BdotTi_Ty=0
         kfaci=0
         kpari_un=0
         kpari_p=0
         kperpi_p=0
         kperpi_bsq=0
         IF(kappa_case == "scalar")THEN
            kperpe=kappa_prp
            kperpe_un=0
            kperpi=kappa_prp
            kperpi_un=0
         ELSE
            kperpe=kappa_prp*rho
            kperpe_un=kappa_prp*rho
            kperpi=kappa_prp*rho
            kperpi_un=kappa_prp*rho
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     heat exchange.
c-----------------------------------------------------------------------
      CALL transport_hexch_u(ieheat,ke_norm,rho,u(9,:,:)+pt0,
     $     u(8,:,:)+pt0,heat_exch,hexch_un,hexch_pi,hexch_pe)
      hexch_un = hexch_un*rho
c-----------------------------------------------------------------------
c     radiative losses.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("CurrentSheet-hlf","MAST")
         CALL transport_radloss_u(rad_fac,T0,rho,u(8,:,:)+pt0,
     &     rad_loss, rloss_un, rloss_pe)
      CASE("Chen-Shibata-hlf")
         rad_loss = 0.0
         rloss_un = 0.0
         rloss_pe = 0.0
      END SELECT
      rloss_un = rloss_un*rho
c-----------------------------------------------------------------------
c     viscous boundary layer.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Chen-Shibata-hlf","CS-stratified","CurrentSheet","MAST")
         CALL ChenShibata_mu(x,y,rho,mu_vbl,mu_local)
      CASE DEFAULT
         mu_local = mu*rho + mu_min
      END SELECT

      SELECT CASE(init_type)
      CASE("MAST")
         CALL MAST_ddiff(x,y,ddiff,ddiff_local)
      CASE DEFAULT
         ddiff_local=ddiff
      END SELECT
c-----------------------------------------------------------------------
c     velocities and their derivatives.
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:) = u(i+3,:,:)*rho_inv
         vix(i,:,:) = ux(i+3,:,:)*rho_inv - vi(i,:,:)*ux(1,:,:)
         viy(i,:,:) = uy(i+3,:,:)*rho_inv - vi(i,:,:)*uy(1,:,:)
      ENDDO
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
      fx_u(1,1,:,:) = -r_fac*vi(1,:,:)
      fx_u(1,4,:,:) = r_fac*rho_inv
      fx_ux(1,1,:,:)= -r_fac*ddiff_local
      fy_u(1,1,:,:) = -r_fac*vi(2,:,:)
      fy_u(1,5,:,:) = r_fac*rho_inv
      fy_uy(1,1,:,:)= -r_fac*ddiff_local
      s_u(1,1,:,:) = r_fac*(vi(1,:,:)*ux(1,:,:) + vi(2,:,:)*uy(1,:,:))
      s_u(1,4,:,:) = r_fac*rhox_inv
      s_u(1,5,:,:) = r_fac*rhoy_inv
      s_ux(1,1,:,:) = -r_fac*(vi(1,:,:) - 2.*ddiff_local)
      s_uy(1,1,:,:) = -r_fac*(vi(2,:,:) - 2.*ddiff_local)
c-----------------------------------------------------------------------
c     poloidal magnetic flux equation.
c-----------------------------------------------------------------------
      fx_ux(2,7,:,:) = r_fac*hyper_eta
      
      fy_uy(2,7,:,:) = r_fac*hyper_eta
      
      s_u(2,1,:,:)=r_fac*(-vi(2,:,:)*b1 + vi(1,:,:)*b2 + eta_rho*j3)
      s_u(2,2,:,:)=-cyl_fac*vi(2,:,:)
      s_u(2,3,:,:)=r_fac*eta_j*j_u3*j3
      s_u(2,4,:,:)=-r_fac*b2*rho_inv
      s_u(2,5,:,:)=r_fac*b1*rho_inv
      s_u(2,7,:,:)=r_fac*(eta_j*j_u7*j3 + eta_local)
     $     + cyl_fac*r_faci*hyper_eta
      s_u(2,8,:,:)=r_fac*eta_p*j3
      s_u(2,9,:,:)=r_fac*eta_pi*j3
      s_ux(2,2,:,:)=-r_fac*vi(1,:,:)
      s_ux(2,3,:,:)=r_fac*eta_j*j_ux3*j3
      s_uy(2,2,:,:)=-r_fac*vi(2,:,:)
      s_uy(2,3,:,:)=r_fac*eta_j*j_uy3*j3
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=-(u(3,:,:)+bz0)*vi(1,:,:) + vi(3,:,:)*b1 
     $     + eta_rho*j2
      fx_u(3,2,:,:)=cyl_fac*r_faci*vi(3,:,:)
      fx_u(3,3,:,:)=vi(1,:,:) + eta_j*j_u3*j2
      fx_u(3,4,:,:)=(u(3,:,:)+bz0)*rho_inv
      fx_u(3,6,:,:)=-rho_inv*b1
      fx_u(3,7,:,:)=eta_j*j_u7*j2
      fx_u(3,8,:,:)= eta_p*j2
      fx_u(3,9,:,:)= eta_pi*j2
      fx_ux(3,3,:,:)=-eta_local + eta_j*j_ux3*j2
      fx_uy(3,2,:,:)=vi(3,:,:)
      fx_uy(3,3,:,:)=eta_j*j_uy3*j2
         
      fy_u(3,1,:,:)=-(u(3,:,:)+bz0)*vi(2,:,:) + vi(3,:,:)*b2 
     $     - eta_rho*j1
      fy_u(3,3,:,:)=vi(2,:,:) - eta_local*cyl_fac*r_faci - eta_j*j_u3*j1
      fy_u(3,5,:,:)=(u(3,:,:)+bz0)*rho_inv
      fy_u(3,6,:,:)=-rho_inv*b2
      fy_u(3,7,:,:)=-eta_j*j_u7*j1
      fy_u(3,8,:,:)=-eta_p*j1
      fy_u(3,9,:,:)=-eta_pi*j1
      fy_ux(3,2,:,:)=-vi(3,:,:)
      fy_ux(3,3,:,:)=-eta_j*j_ux3*j1
      fy_uy(3,3,:,:)=-eta_local - eta_j*j_uy3*j1
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:)=-r_fac*(u(4,:,:)*vi(1,:,:) - two*mu_min*vix(1,:,:))
      fx_u(4,3,:,:)=r_fac*(u(3,:,:)+bz0)
      fx_u(4,4,:,:)=r_fac*two*(vi(1,:,:) - mu_local*rhox_inv)
      fx_u(4,8,:,:)=r_fac
      fx_u(4,9,:,:)=r_fac
      fx_ux(4,1,:,:)= r_fac*two*mu_local*vi(1,:,:)
      fx_ux(4,4,:,:)=-r_fac*two*mu_local*rho_inv

      fy_u(4,1,:,:)=-r_fac*(u(4,:,:)*vi(2,:,:)
     $     - mu_min*(viy(1,:,:) + vix(2,:,:)))
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
     $     - mu_min*(vix(2,:,:) + viy(1,:,:))) 
      fx_u(5,4,:,:)=r_fac*(vi(2,:,:) - mu_local*rhoy_inv)
      fx_u(5,5,:,:)=r_fac*(vi(1,:,:) - mu_local*rhox_inv)
      fx_ux(5,1,:,:)= r_fac*mu_local*vi(2,:,:)
      fx_ux(5,5,:,:)=-r_fac*mu_local*rho_inv
      fx_uy(5,1,:,:)= r_fac*mu_local*vi(1,:,:)
      fx_uy(5,4,:,:)=-r_fac*mu_local*rho_inv
      
      fy_u(5,1,:,:)=-r_fac*(u(5,:,:)*vi(2,:,:)
     $     - two*mu_min*viy(2,:,:)) 
      fy_u(5,3,:,:)= r_fac*(u(3,:,:)+bz0)
      fy_u(5,5,:,:)= r_fac*two*(vi(2,:,:) - mu_local*rhoy_inv)
      fy_uy(5,1,:,:)= r_fac*two*mu_local*vi(2,:,:)
      fy_uy(5,5,:,:)=-r_fac*two*mu_local*rho_inv

      s_u(5,1,:,:)=-r_fac*gy_vec*rho - cyl_fac*(u(6,:,:)*vi(3,:,:) 
     $     - two*r_faci*mu_min*vi(2,:,:))
      s_u(5,2,:,:)=-cyl_fac*j3
      s_u(5,3,:,:)=-cyl_fac*(u(3,:,:)+bz0)
      s_u(5,5,:,:)=-cyl_fac*r_faci*two*mu_local*rho_inv
      s_u(5,6,:,:)=cyl_fac*two*vi(3,:,:)
      s_u(5,7,:,:)=r_fac*b1
      s_uy(5,2,:,:)=-r_fac*j3
      s_uy(5,8,:,:)=-r_fac
      s_uy(5,9,:,:)=-r_fac
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------      
      fx_u(6,1,:,:)=-r_fac*(u(6,:,:)*vi(1,:,:) - mu_min*vix(3,:,:)) 
      fx_u(6,4,:,:)=r_fac*u(6,:,:)*rho_inv
      fx_u(6,6,:,:)=r_fac*(vi(1,:,:) - mu_local*rhox_inv)
      fx_ux(6,1,:,:)= r_fac*mu_local*vi(3,:,:)
      fx_ux(6,6,:,:)=-r_fac*mu_local*rho_inv
      
      fy_u(6,1,:,:)=-r_fac*(u(6,:,:)*vi(2,:,:) - mu_min*viy(3,:,:)) 
      fy_u(6,5,:,:)=r_fac*u(6,:,:)*rho_inv
      fy_u(6,6,:,:)=r_fac*(vi(2,:,:) - mu_local*rhoy_inv)
      fy_uy(6,1,:,:)= r_fac*mu_local*vi(3,:,:)
      fy_uy(6,6,:,:)=-r_fac*mu_local*rho_inv
      
      s_u(6,1,:,:)=cyl_fac*(u(6,:,:)*vi(2,:,:) 
     $     + r_faci*mu_min*vi(3,:,:))
      s_u(6,2,:,:)=cyl_fac*j2
      s_u(6,3,:,:)=cyl_fac*b2
      s_u(6,5,:,:)=-cyl_fac*vi(3,:,:)
      s_u(6,6,:,:)=-cyl_fac*(vi(2,:,:) + r_faci*mu_local*rho_inv)
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
c     electron pressure equation.
c-----------------------------------------------------------------------
      IF(eta_case == "Chen-Shibata" .AND. (init_type == "CurrentSheet"
     $   .OR. init_type == "Chen-Shibata-hlf"))
     $     CALL ChenShibata_eta(x,y,jtot,zero,eta_chrsp,eta_local)

      IF(eta_case == "MAST") CALL MAST_eta(x,y,zero,eta_local) 

      fx_u(8,1,:,:)=r_fac*(-gamma_fac*(u(8,:,:)+pt0)*vi(1,:,:)
     $     + (kface*BdotTe_Tx(1,:,:) + kperpe)*(Tex+Tx0)
     $     + kface*BdotTe_Ty(1,:,:)*(Tey+Ty0)
     $     - (kpare_un - kperpe_un)*BdotTe(1,:,:) - kperpe_un*Tex)
      fx_u(8,2,:,:)=-cyl_fac*(-kface*BdotTe_b1(1,:,:)
     $     + two*b1*kperpe_bsq*(BdotTe(1,:,:) - Tex))
      fx_u(8,3,:,:)=r_fac*(-kface*BdotTe_b3(1,:,:)
     $     + two*(u(3,:,:)+bz0)*kperpe_bsq*(BdotTe(1,:,:) - Tex))
      fx_u(8,4,:,:)=r_fac*gamma_fac*(u(8,:,:)+pt0)*rho_inv
      fx_u(8,8,:,:)=r_fac*(gamma_fac*vi(1,:,:)
     $     - (kface*BdotTe_Tx(1,:,:) + kperpe)*rhox_inv
     $     - kface*BdotTe_Ty(1,:,:)*rhoy_inv
     $     - (kpare_p - kperpe_p)*BdotTe(1,:,:) - kperpe_p*Tex)

      fx_ux(8,1,:,:)=-r_fac*(kface*BdotTe_Tx(1,:,:) + kperpe)*Te_un
      fx_ux(8,2,:,:)=r_fac*(-kface*BdotTe_b2(1,:,:)
     $     + two*b2*kperpe_bsq*(BdotTe(1,:,:) - Tex))
      fx_ux(8,8,:,:)=-r_fac*(kface*BdotTe_Tx(1,:,:) + kperpe)*rho_inv

      fx_uy(8,1,:,:)=-r_fac*kface*BdotTe_Ty(1,:,:)*Te_un
      fx_uy(8,2,:,:)=-r_fac*(-kface*BdotTe_b1(1,:,:)
     $     + two*b1*kperpe_bsq*(BdotTe(1,:,:) - Tex))
      fx_uy(8,8,:,:)=-r_fac*kface*BdotTe_Ty(1,:,:)*rho_inv

      fy_u(8,1,:,:)=r_fac*(-gamma_fac*(u(8,:,:)+pt0)*vi(2,:,:)
     $     + (kface*BdotTe_Ty(2,:,:) + kperpe)*(Tey+Ty0)
     $     + kface*BdotTe_Tx(2,:,:)*(Tex+Tx0)
     $     - (kpare_un - kperpe_un)*BdotTe(2,:,:) - kperpe_un*Tey)
      fy_u(8,2,:,:)=-cyl_fac*(-kface*BdotTe_b1(2,:,:)
     $     + two*b1*kperpe_bsq*(BdotTe(2,:,:) - Tey))
      fy_u(8,3,:,:)=r_fac*(-kface*BdotTe_b3(2,:,:)
     $     + two*(u(3,:,:)+bz0)*kperpe_bsq*(BdotTe(2,:,:) - Tey))
      fy_u(8,5,:,:)=r_fac*gamma_fac*(u(8,:,:)+pt0)*rho_inv
      fy_u(8,8,:,:)=r_fac*(gamma_fac*vi(2,:,:)
     $     - (kface*BdotTe_Ty(2,:,:) + kperpe)*rhoy_inv
     $     - kface*BdotTe_Tx(2,:,:)*rhox_inv
     $     - (kpare_p - kperpe_p)*BdotTe(2,:,:) - kperpe_p*Tey)

      fy_ux(8,1,:,:)=-r_fac*kface*BdotTe_Tx(2,:,:)*Te_un
      fy_ux(8,2,:,:)=r_fac*(-kface*BdotTe_b2(2,:,:)
     $     + two*b2*kperpe_bsq*(BdotTe(2,:,:) - Tey))
      fy_ux(8,8,:,:)=-r_fac*kface*BdotTe_Tx(2,:,:)*rho_inv

      fy_uy(8,1,:,:)=-r_fac*(kface*BdotTe_Ty(2,:,:) + kperpe)*Te_un
      fy_uy(8,2,:,:)=-r_fac*(-kface*BdotTe_b1(2,:,:)
     $     + two*b1*kperpe_bsq*(BdotTe(2,:,:) - Tey))
      fy_uy(8,8,:,:)=-r_fac*(kface*BdotTe_Ty(2,:,:) + kperpe)*rho_inv

      s_u(8,1,:,:)=r_fac*(-vi(1,:,:)*(ux(8,:,:)+px0)
     $     - vi(2,:,:)*(uy(8,:,:)+py0) - hexch_un - rloss_un
     $     + eta_rho*jtot**2)
      s_u(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_u3
      s_u(8,4,:,:)=r_fac*(ux(8,:,:)+px0)*rho_inv
      s_u(8,5,:,:)=r_fac*(uy(8,:,:)+py0)*rho_inv 
      s_u(8,7,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_u7
      s_u(8,8,:,:)=r_fac*(eta_p*jtot**2 - hexch_pe - rloss_pe)
      s_u(8,9,:,:)=r_fac*(-hexch_pi+eta_pi*jtot**2)

      s_ux(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_ux3
      s_ux(8,7,:,:)=r_fac*hyper_eta*two*ux(7,:,:)
      s_ux(8,8,:,:)=r_fac*vi(1,:,:)

      s_uy(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot)*jtot*j_uy3
      s_uy(8,7,:,:)=r_fac*hyper_eta*two*uy(7,:,:)
      s_uy(8,8,:,:)=r_fac*vi(2,:,:)
c-----------------------------------------------------------------------
c     ion pressure equation.
c-----------------------------------------------------------------------
      IF(init_type == "CurrentSheet" .OR. init_type=="Chen-Shibata-hlf"
     $   .OR. init_type == "MAST")
     $     CALL ChenShibata_mu(x,y,rho,zero,mu_local)
      
      fx_u(9,1,:,:)=r_fac*(-gamma_fac*(u(9,:,:)+pt0)*vi(1,:,:)
     $     + (kfaci*BdotTi_Tx(1,:,:) + kperpi)*(Tix+Tx0)
     $     + kfaci*BdotTi_Ty(1,:,:)*(Tiy+Ty0)
     $     - (kpari_un - kperpi_un)*BdotTi(1,:,:) - kperpi_un*Tix)
      fx_u(9,2,:,:)=-cyl_fac*(-kfaci*BdotTi_b1(1,:,:)
     $     + two*b1*kperpi_bsq*(BdotTi(1,:,:) - Tix))
      fx_u(9,3,:,:)=r_fac*(-kfaci*BdotTi_b3(1,:,:)
     $     + two*(u(3,:,:)+bz0)*kperpi_bsq*(BdotTi(1,:,:) - Tix))
      fx_u(9,4,:,:)=r_fac*gamma_fac*(u(9,:,:)+pt0)*rho_inv
      fx_u(9,9,:,:)=r_fac*(gamma_fac*vi(1,:,:)
     $     - (kfaci*BdotTi_Tx(1,:,:) + kperpi)*rhox_inv
     $     - kfaci*BdotTi_Ty(1,:,:)*rhoy_inv
     $     - (kpari_p - kperpi_p)*BdotTi(1,:,:) - kperpi_p*Tix)

      fx_ux(9,1,:,:)=-r_fac*(kfaci*BdotTi_Tx(1,:,:) + kperpi)*Ti_un
      fx_ux(9,2,:,:)=r_fac*(-kfaci*BdotTi_b2(1,:,:)
     $     + two*b2*kperpi_bsq*(BdotTi(1,:,:) - Tix))
      fx_ux(9,9,:,:)=-r_fac*(kfaci*BdotTi_Tx(1,:,:) + kperpi)*rho_inv

      fx_uy(9,1,:,:)=-r_fac*kfaci*BdotTi_Ty(1,:,:)*Ti_un
      fx_uy(9,2,:,:)=-r_fac*(-kfaci*BdotTi_b1(1,:,:)
     $     + two*b1*kperpi_bsq*(BdotTi(1,:,:) - Tix))
      fx_uy(9,9,:,:)=-r_fac*kfaci*BdotTi_Ty(1,:,:)*rho_inv

      fy_u(9,1,:,:)=r_fac*(-gamma_fac*(u(9,:,:)+pt0)*vi(2,:,:)
     $     + (kfaci*BdotTi_Ty(2,:,:) + kperpi)*(Tiy+Ty0)
     $     + kfaci*BdotTi_Tx(2,:,:)*(Tix+Tx0)
     $     - (kpari_un - kperpi_un)*BdotTi(2,:,:) - kperpi_un*Tiy)
      fy_u(9,2,:,:)=-cyl_fac*(-kfaci*BdotTi_b1(2,:,:)
     $     + two*b1*kperpi_bsq*(BdotTi(2,:,:) - Tiy))
      fy_u(9,3,:,:)=r_fac*(-kfaci*BdotTi_b3(2,:,:)
     $     + two*(u(3,:,:)+bz0)*kperpi_bsq*(BdotTi(2,:,:) - Tiy))
      fy_u(9,5,:,:)=r_fac*gamma_fac*(u(9,:,:)+pt0)*rho_inv
      fy_u(9,9,:,:)=r_fac*(gamma_fac*vi(2,:,:)
     $     - (kfaci*BdotTi_Ty(2,:,:) + kperpi)*rhoy_inv
     $     - kfaci*BdotTi_Tx(2,:,:)*rhox_inv
     $     - (kpari_p - kperpi_p)*BdotTi(2,:,:) - kperpi_p*Tiy)

      fy_ux(9,1,:,:)=-r_fac*kfaci*BdotTi_Tx(2,:,:)*Ti_un
      fy_ux(9,2,:,:)=r_fac*(-kfaci*BdotTi_b2(2,:,:)
     $     + two*b2*kperpi_bsq*(BdotTi(2,:,:) - Tiy))
      fy_ux(9,9,:,:)=-r_fac*kfaci*BdotTi_Tx(2,:,:)*rho_inv

      fy_uy(9,1,:,:)=-r_fac*(kfaci*BdotTi_Ty(2,:,:) + kperpi)*Ti_un
      fy_uy(9,2,:,:)=-r_fac*(-kfaci*BdotTi_b1(2,:,:)
     $     + two*b1*kperpi_bsq*(BdotTi(2,:,:) - Tiy))
      fy_uy(9,9,:,:)=-r_fac*(kfaci*BdotTi_Ty(2,:,:) + kperpi)*rho_inv

      s_u(9,1,:,:)=r_fac*(-vi(1,:,:)*(ux(9,:,:)+px0)
     $     - vi(2,:,:)*(uy(9,:,:)+py0)
     $     - (mu_local+mu_min)*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:)
     $     + vix(3,:,:)**2 + viy(3,:,:)**2) + hexch_un)
     $     - cyl_fac*r_faci*(mu_local+mu_min)*two*vi(2,:,:)**2
      s_u(9,4,:,:)=r_fac*((ux(9,:,:)+px0)*rho_inv 
     $     + mu_local*(4._r8*vix(1,:,:)*rhox_inv
     $     + two*(vix(2,:,:) + viy(1,:,:))*rhoy_inv))
      s_u(9,5,:,:)=r_fac*((uy(9,:,:)+py0)*rho_inv 
     $     + mu_local*(4._r8*viy(2,:,:)*rhoy_inv
     $     + two*(vix(2,:,:) + viy(1,:,:))*rhox_inv))
     $     + cyl_fac*r_faci*4._r8*mu_local*vi(2,:,:)*rho_inv
      s_u(9,6,:,:)=r_fac*two
     $     *mu_local*(vix(3,:,:)*rhox_inv + viy(3,:,:)*rhoy_inv)
      s_u(9,8,:,:)=r_fac*hexch_pe
      s_u(9,9,:,:)=r_fac*hexch_pi

      s_ux(9,1,:,:)= -r_fac*two*mu_local*(two*vix(1,:,:)*vi(1,:,:)  
     $     + vi(2,:,:)*(vix(2,:,:) + viy(1,:,:))
     $     + vix(3,:,:)*vi(3,:,:))
      s_ux(9,4,:,:)=r_fac*mu_local*4._r8*vix(1,:,:)*rho_inv
      s_ux(9,5,:,:)=r_fac*mu_local*two*(vix(2,:,:) + viy(1,:,:))*rho_inv
      s_ux(9,6,:,:)=r_fac*two*mu_local*vix(3,:,:)*rho_inv
      s_ux(9,9,:,:)=r_fac*vi(1,:,:)

      s_uy(9,1,:,:)= -r_fac*two*mu_local*(two*viy(2,:,:)*vi(2,:,:)  
     $     + vi(1,:,:)*(viy(1,:,:) + vix(2,:,:))
     $     + viy(3,:,:)*vi(3,:,:))
      s_uy(9,4,:,:)=r_fac*mu_local*two*(viy(1,:,:) + vix(2,:,:))*rho_inv
      s_uy(9,5,:,:)=r_fac*mu_local*4._r8*viy(2,:,:)*rho_inv
      s_uy(9,6,:,:)=r_fac*two*mu_local*viy(3,:,:)*rho_inv
      s_uy(9,9,:,:)=r_fac*vi(2,:,:)
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
      USE SolarLnMHD_mod
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
      mass(9,9,:,:)=r_fac/(gamma-one)
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
      USE SolarLnMHD_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: ksi,phi
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: x,y
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("TwoFR","MAST")
        x=lx*(x_curve*ksi**3 + ksi)/(x_curve + one)
        y=ly*(y_curve*phi**2 + phi)/(y_curve + one)
      CASE("Chen-Shibata")
        x=lx*0.5*(x_curve*(ksi-0.5)**3 + (ksi-0.5))
     $        /(.125*x_curve + .5)
        y=ly*(y_curve*phi**2 + phi)/(y_curve+one)
      CASE("Chen-Shibata-hlf","CS-stratified", "CurrentSheet-hlf")
c        x=lx*0.5*(x_curve*ksi**3 + ksi)/(x_curve + 1.)
        x=lx*0.5*(x_curve*ksi**5 + ksi)/(x_curve + 1.)
c        y=ly*(y_curve*phi**2 + phi)/(y_curve+one)
c Elena: y-grid modification to resolve CS near y=0.74
        y=ly*((phi-0.272)**3+y_curve*phi+0.272**3)/
     $       ((1.-0.272)**3+y_curve+0.272**3)
      CASE("CurrentSheet")
        y=ly*(y_curve*phi**2 + phi)/(y_curve+one)
        x=lx*0.5*(x_curve*ksi**3 + ksi)/(x_curve + 1.)
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
      USE SolarLnMHD_mod
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
      USE SolarLnMHD_mod
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
