c-----------------------------------------------------------------------
c     file XMHDhlx.f.
c     contains specifications for MHD equations with helical symmetry.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. XMHDhlx_mod.
c     1. XMHDhlx_input.
c     2. XMHDhlx_init.
c     3. XMHDhlx_boundary.
c     4. XMHDhlx_bound_rhs.
c     5. XMHDhlx_bound_drdu.
c     6. XMHDhlx_init_special.
c     7. XMHDhlx_rhs.
c     8. XMHDhlx_drdu.
c     9. XMHDhlx_mass.
c     10. XMHDhlx_equil.
c     11. XMHDhlx_grid.
c-----------------------------------------------------------------------
c     subprogram 0. XMHDhlx_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE XMHDhlx_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: source=.FALSE.
      CHARACTER(16), PRIVATE :: init_type=" "
      REAL(r8), PARAMETER, PRIVATE :: mass_r=5.44617e-4
      REAL(r8), PRIVATE :: di=0,eta=0,mu=0,nu=0,eta_0=0,mu_0=0,nu_0=0,
     $     t0_init=0,t0_decay=0,epsilon=0,Tconst=0,beta_e=0,beta0=0,
     $     r_s=0,Rinv=0,heatfdg=1.,n0=1.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. XMHDhlx_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_input(nx,ny,np,nq,nbx,xperiodic,
     $     yperiodic,nqty,dt,dtmax,tmax,nstep,gr_curve)

      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic
      INTEGER, INTENT(OUT) :: nqty
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax,gr_curve

      INTEGER :: myios
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/kink_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,di,eta,mu,nu,eta_0,mu_0,nu_0,
     $     t0_init,t0_decay,beta0,beta_e,Tconst,r_s,Rinv,
     $     epsilon,source,init_type
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=kink_list,IOSTAT=myios)

      nqty=9

      xperiodic=.FALSE.
      yperiodic=.TRUE.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_input
c-----------------------------------------------------------------------
c     subprogram 2. XMHDhlx_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_init(xpi,ypi,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: power
      REAL(r8), DIMENSION(0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) :: gfac,pfac
      REAL(r8), DIMENSION(SIZE(ui,1),0:SIZE(ui,2)-1,0:SIZE(ui,3)-1) :: 
     $     ux,uy
c-----------------------------------------------------------------------
c     add perturbation to an equilibrium
c-----------------------------------------------------------------------
      gfac=1/SQRT(1+(xpi*Rinv)**2)
      CALL XMHDhlx_equil(xpi,ypi,ui,ux,uy,.FALSE.)

      SELECT CASE(init_type)
      CASE("m1_axial","sawtooth1","sawtooth2")
         power = 2
         pfac = (2.5*xpi)**(2*power)

         ui(2,:,:)=ui(2,:,:) - epsilon*COS(ypi)
     $        *.5*xpi*EXP(-pfac)*ui(1,:,:)
         ui(3,:,:)=ui(3,:,:) + epsilon*SIN(ypi)
     $        *gfac*xpi*(1. - power*pfac)*EXP(-pfac)*ui(1,:,:)
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_init
c-----------------------------------------------------------------------
c     subprogram 3. XMHDhlx_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_init_special(static)

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(7)=.TRUE.
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
      IF(source)
     $     CALL program_stop("Do not run XMHDhlx with source term!")
      IF(di == 0 .OR. beta_e == 0)
     $     CALL program_stop("Neither di nor beta_e can be zero!")
      SELECT CASE(init_type)
      CASE("m1_axial")
         IF(r_s == .5)THEN
            Tconst = (300. - 3180.*Rinv**2 + 1800.*Rinv**4)
     $           /(3840.*beta0 - 2035. - 1195.*Rinv**2 + 1034.*Rinv**4)
         ELSE
            CALL program_stop("Non-zero density gradient at r=1 !")
         ENDIF
      CASE("sawtooth1")
         n0=0.9
      CASE("sawtooth2")
         n0=0.8
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_init_special
c-----------------------------------------------------------------------
c     subprogram 4. XMHDhlx_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_boundary(left,right,top,bottom,nqty)

      INTEGER, INTENT(IN) :: nqty
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      left%bc_type(1)="polar"
      left%bc_type(2:3)="robin"
      left%static(2:3)=.TRUE.
      left%bc_type(4:9)="polar"
      SELECT CASE(init_type)
      CASE("m1_axial","sawtooth1","sawtooth2")
         right%bc_type(1)="natural"
         right%static(2)=.TRUE.
      END SELECT
      right%static(3:4)=.TRUE.
      right%bc_type(5)="robin"
      right%a(5,5)=1.
      right%static(6:9)=.TRUE.
      top%bc_type = "periodic"
      bottom%bc_type = "periodic"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_boundary
c-----------------------------------------------------------------------
c     subprogram 5. XMHDhlx_bound_rhs.
c     computes rhs for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_bound_rhs(t,edge_type,x,y,u,ux,uy,c)

      REAL(r8), INTENT(IN) :: t
      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: c

      REAL(r8) :: gfac,hfac,Ezb,Teb,Jzb
      REAL(r8), DIMENSION(SIZE(x)) :: eta_local,nu_e,Btot,btau,Jpar,
     $     n_inv
c-----------------------------------------------------------------------
c     define auxiliary variables.
c-----------------------------------------------------------------------
      gfac=1/SQRT(1. + Rinv**2)
      hfac = Rinv**2*gfac**2

      n_inv = 1./u(1,:)
      Btot = SQRT(uy(5,:)**2 + (gfac*ux(5,:))**2 + (gfac*u(6,:))**2)
      btau = gfac*ux(5,:)/Btot
      Jpar = - (uy(5,:)*uy(6,:) + gfac**2*ux(5,:)*ux(6,:)
     $     - gfac/di*u(6,:)*(u(4,:) - u(1,:)*u(7,:)))/Btot

      eta_local = di*u(1,:)**1.5/(1.96*nu*u(9,:)**1.5)
      nu_e = 2.*mass_r*di/nu*u(1,:)**2.5/SQRT(u(9,:))
c-----------------------------------------------------------------------
c     set boundary conditions.
c-----------------------------------------------------------------------
      c=0
      SELECT CASE(init_type)
      CASE("m1_axial")
         Teb = beta0*beta_e + Tconst*(r_s**2 - 1.)
         Jzb = - (r_s**2 - 3.*Rinv**2 + 2.*(Rinv*r_s)**2)
      CASE("sawtooth1","sawtooth2")
         Teb = beta_e*(beta0
     $        - (2.*(Rinv*gfac)**2 + Rinv**4 + 32.*Rinv**2
     $        - 16.*Rinv**2*(gfac + 1./gfac) 
     $        + 4.*(1. - 8.*Rinv**4)*LOG(gfac))/(32.*Rinv**6))
     $        /(n0 + (one - n0)*EXP(-1./r_s**4))
         Jzb = - 0.5*gfac**5*(2. + Rinv**2)*(2. + 2.*Rinv**2 - 1/gfac)
      END SELECT
      Ezb = di/(1.96*nu*Teb**1.5)*Jzb
      IF(edge_type == "left")THEN
         c(2,:) = u(2,:)
         c(3,:) = u(3,:)
      ENDIF
      IF(edge_type == "right")THEN
         c(2,:) = u(2,:)
         c(3,:) = ux(3,:)*u(1,:) - ux(1,:)*u(3,:)
         c(4,:) = ux(4,:)*u(1,:) - ux(1,:)*u(4,:)
         c(5,:) = Ezb
         c(6,:) = Rinv*Ezb + di*n_inv*(uy(6,:)*u(6,:) + uy(9,:)/gfac**2
     $        + nu_e*2.*Rinv*gfac*(ux(7,:) - hfac*u(7,:)))
     $        - n_inv*u(2,:)*u(6,:)
     $        + eta_local*(1.96*ux(6,:) + .96/gfac*Jpar*btau)
         c(7,:) = ux(7,:)
         c(8,:) = u(8,:) - Teb*u(1,:)*(1. - beta_e)/beta_e
         c(9,:) = u(9,:) - Teb*u(1,:)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_bound_rhs
c-----------------------------------------------------------------------
c     subprogram 6. XMHDhlx_bound_drdu.
c     computes drdu for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_bound_drdu(t,edge_type,x,y,u,ux,uy,
     $     c_u,c_ux,c_uy)

      REAL(r8), INTENT(IN) :: t
      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy

      REAL(r8) :: gfac,hfac,Teb,Ezb,Jzb
      REAL(r8), DIMENSION(SIZE(x)) :: 
     $     eta_local,eta_local_u1,eta_local_u9,nu_e,nu_e_u1,nu_e_u9,
     $     Btot,btau,btau_u6,btau_ux5,btau_uy5,Jpar,Jpar_u1,Jpar_u4,
     $     Jpar_u6,Jpar_u7,Jpar_ux5,Jpar_ux6,Jpar_uy5,Jpar_uy6,n_inv
c-----------------------------------------------------------------------
c     define auxiliary variables.
c-----------------------------------------------------------------------
      gfac=1/SQRT(1+Rinv**2)
      hfac = Rinv**2*gfac**2

      SELECT CASE(init_type)
      CASE("m1_axial")
         Teb = beta0*beta_e + Tconst*(r_s**2 - 1.)
         Jzb = - (r_s**2 - 3.*Rinv**2 + 2.*(Rinv*r_s)**2)
      CASE("sawtooth1","sawtooth2")
         Teb = beta_e*(beta0
     $        - (2.*(Rinv*gfac)**2 + Rinv**4 + 32.*Rinv**2
     $        - 16.*Rinv**2*(gfac + 1./gfac) 
     $        + 4.*(1. - 8.*Rinv**4)*LOG(gfac))/(32.*Rinv**6))
     $        /(n0 + (one - n0)*EXP(-1./r_s**4))
         Jzb = - 0.5*gfac**5*(2. + Rinv**2)*(2. + 2.*Rinv**2 - 1/gfac)
      END SELECT

      Ezb = di/(1.96*nu*Teb**1.5)*Jzb

      n_inv = 1./u(1,:)
      Btot = SQRT(uy(5,:)**2 + (gfac*ux(5,:))**2 + (gfac*u(6,:))**2)
      btau = gfac*ux(5,:)/Btot
      Jpar = - (uy(5,:)*uy(6,:) + gfac**2*ux(5,:)*ux(6,:)
     $     - gfac/di*u(6,:)*(u(4,:) - u(1,:)*u(7,:)))/Btot

      btau_u6 = - gfac**3*ux(5,:)*u(6,:)/Btot**3
      btau_ux5 = gfac/Btot*(1 - (gfac*ux(5,:)/Btot)**2)
      btau_uy5 = - gfac*ux(5,:)*uy(5,:)/Btot**3

      Jpar_u1 = - gfac/di*u(6,:)*u(7,:)/Btot
      Jpar_u4 = gfac/di*u(6,:)/Btot
      Jpar_u6 = (gfac/di*(u(4,:) - u(1,:)*u(7,:)) 
     $     - gfac**2*u(6,:)*Jpar/Btot)/Btot
      Jpar_u7 = - gfac/di*u(6,:)*u(1,:)/Btot
      Jpar_ux5 = - gfac**2*(ux(6,:) + ux(5,:)*Jpar/Btot)/Btot
      Jpar_ux6 = - gfac**2*ux(5,:)/Btot
      Jpar_uy5 = - (uy(6,:) + uy(5,:)*Jpar/Btot)/Btot
      Jpar_uy6 = - uy(5,:)/Btot

      eta_local = di*u(1,:)**1.5/(1.96*nu*u(9,:)**1.5)
      eta_local_u1 = 1.5*eta_local/u(1,:)
      eta_local_u9 = - 1.5*eta_local/u(9,:)

      nu_e = 2.*mass_r*di/nu*u(1,:)**2.5/SQRT(u(9,:))
      nu_e_u1 = 2.5*nu_e/u(1,:)
      nu_e_u9 = - 0.5*nu_e/u(9,:)
c-----------------------------------------------------------------------
c     set drdu for boundary conditions.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      IF(edge_type == "left")THEN
         c_u(2,2,:) = 1.
         c_u(3,3,:) = 1.
      ENDIF
      IF(edge_type == "right")THEN
         
         c_u(2,2,:) = 1.

         c_u(3,1,:) = ux(3,:)
         c_u(3,3,:) = -ux(1,:)
         c_ux(3,1,:) = -u(3,:)
         c_ux(3,3,:) = u(1,:)

         c_u(4,1,:) = ux(4,:)
         c_u(4,4,:) = -ux(1,:)
         c_ux(4,1,:) = -u(4,:)
         c_ux(4,4,:) = u(1,:)
         
         c_u(6,1,:) = - di*2.*Rinv*gfac*n_inv
     $        *(nu_e*n_inv - nu_e_u1)*(ux(7,:) - hfac*u(7,:)) 
     $        - di*n_inv**2*uy(9,:)/gfac**2
     $        + (u(2,:) - di*uy(6,:))*u(6,:)*n_inv**2
     $        - eta_local*.96/gfac*Jpar_u1*btau
     $        + eta_local_u1*(1.96*ux(6,:) + .96/gfac*Jpar*btau)
         c_u(6,2,:) = - n_inv*u(6,:)
         c_u(6,4,:) = eta_local*.96/gfac*Jpar_u4*btau
         c_u(6,6,:) = - (u(2,:) - di*uy(6,:))*n_inv
     $        + eta_local*.96/gfac*(Jpar_u6*btau + Jpar*btau_u6)
         c_u(6,7,:) = - nu_e*di*2.*Rinv*gfac*hfac*n_inv
     $        + eta_local*.96/gfac*Jpar_u7*btau
         c_u(6,9,:) =  nu_e_u9*di*2.*Rinv*gfac*n_inv
     $        *(ux(7,:) - hfac*u(7,:))
     $        + eta_local_u9*(1.96*ux(6,:) + .96/gfac*Jpar*btau)
         c_ux(6,5,:) = eta_local*.96/gfac
     $        *(Jpar_ux5*btau + Jpar*btau_ux5)
         c_ux(6,6,:) = eta_local*(1.96 + .96/gfac*Jpar_ux6*btau)
         c_ux(6,7,:) = nu_e*di*2.*Rinv*gfac*n_inv
         c_uy(6,5,:) = eta_local*.96/gfac
     $        *(Jpar_uy5*btau + Jpar*btau_uy5)
         c_uy(6,6,:) = di*u(6,:)*n_inv 
     $        + eta_local*.96/gfac*Jpar_uy6*btau
         c_uy(6,9,:) = di*n_inv/gfac**2
         
         c_ux(7,7,:) = 1.

         c_u(8,1,:) = - Teb*(1.- beta_e)/beta_e
         c_u(8,8,:) = 1.
         c_u(9,1,:) = - Teb
         c_u(9,9,:) = 1.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_bound_drdu
c-----------------------------------------------------------------------
c     subprogram 7. XMHDhlx_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE XMHDhlx_rhs(t,x,y,u,ux,uy,fx,fy,s,first)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: 
     $     gfac,hfac,kfac,mfac,lfac,
     $     Vr,Vrx,Vry,Vtau,Vtaux,Vtauy,Ve,Vex,Vey,Btot,divVi,br,btau,be,
     $     Ti,Tix,Tiy,Te,Tex,Tey,n_inv,nx_inv,ny_inv,Jpar,
     $     eta_local,mu_c,mu_i,kap_pari,kap_hati,kap_prpi,heatfl,nu_e,
     $     kap_pare,kap_hate,kap_prpe,Te_source
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     calculate auxiliary variables.
c-----------------------------------------------------------------------
      gfac = 1./SQRT(1. + (x*Rinv)**2)
      hfac = x*Rinv**2*gfac**2
      kfac = x**2*gfac**3
      mfac = 1./(x*gfac**3)
      lfac = 1./(x*gfac)
      n_inv = 1./u(1,:,:)
      nx_inv = - ux(1,:,:)/u(1,:,:)**2
      ny_inv = - uy(1,:,:)/u(1,:,:)**2
c-----------------------------------------------------------------------
c     calculate auxiliary temperature variables.
c-----------------------------------------------------------------------
      Ti = u(8,:,:)*n_inv
      Tix = ux(8,:,:)*n_inv + u(8,:,:)*nx_inv
      Tiy = uy(8,:,:)*n_inv + u(8,:,:)*ny_inv
      Te = u(9,:,:)*n_inv
      Tex = ux(9,:,:)*n_inv + u(9,:,:)*nx_inv
      Tey = uy(9,:,:)*n_inv + u(9,:,:)*ny_inv
c-----------------------------------------------------------------------
c     calculate auxiliary B-field variables.
c-----------------------------------------------------------------------
      Btot = SQRT((uy(5,:,:)/x)**2 + (gfac*ux(5,:,:))**2 
     $     + (gfac*u(6,:,:))**2)
      br = -uy(5,:,:)/(x*Btot)
      btau = gfac*ux(5,:,:)/Btot
      be = gfac*u(6,:,:)/Btot
      Jpar = br/x*uy(6,:,:) - gfac*btau*ux(6,:,:) 
     $     + be/di*(u(4,:,:) - u(1,:,:)*u(7,:,:))
c-----------------------------------------------------------------------
c     calculate auxiliary velocity variables.
c-----------------------------------------------------------------------
      Vr = u(2,:,:)*n_inv
      Vrx = ux(2,:,:)*n_inv + u(2,:,:)*nx_inv
      Vry = uy(2,:,:)*n_inv + u(2,:,:)*ny_inv
      Vtau = u(3,:,:)*n_inv
      Vtaux = ux(3,:,:)*n_inv + u(3,:,:)*nx_inv
      Vtauy = uy(3,:,:)*n_inv + u(3,:,:)*ny_inv
      Ve = u(4,:,:)*n_inv
      Vex = ux(4,:,:)*n_inv + u(4,:,:)*nx_inv
      Vey = uy(4,:,:)*n_inv + u(4,:,:)*ny_inv
      divVi = Vr/x + Vrx + Vtauy*lfac
c-----------------------------------------------------------------------
c     define diffusion coefficients.
c-----------------------------------------------------------------------
      mu_c = 0.96*di*mu*Ti**2.5
      mu_i = 1.2*di/mu*u(1,:,:)**2.5/SQRT(u(8,:,:))
      kap_pari = 3.9*di*mu*Ti**2.5
      kap_hati = 2.5*di*u(8,:,:)
      kap_prpi = 2.*di/mu*u(1,:,:)**2.5/SQRT(u(8,:,:))

      heatfl = heatfdg*3.*u(1,:,:)**2.5/(di*nu*u(9,:,:)**1.5)
      eta_local = di/(1.96*nu*Te**1.5)
      nu_e = 2.*mass_r*di/nu*u(1,:,:)**2.5/SQRT(u(9,:,:))
      kap_pare = 3.2*di*nu*Te**2.5
      kap_hate = - 2.5*di*u(9,:,:)
      kap_prpe = 4.7*di/nu*u(1,:,:)**2.5/SQRT(u(9,:,:))
c-----------------------------------------------------------------------
c     define source terms.
c-----------------------------------------------------------------------
      Te_source = 0
      SELECT CASE(init_type)
      CASE("sawtooth1","sawtooth2")
         Te_source = Tconst*EXP(-x**4/r_s**4)
      END SELECT
c-----------------------------------------------------------------------
c     Helical MHD model.
c     Density Eq.
c-----------------------------------------------------------------------
      fx(1,:,:) = x*u(2,:,:)
      fy(1,:,:) = u(3,:,:)/gfac
c-----------------------------------------------------------------------
c     r-Momentum Eq.
c-----------------------------------------------------------------------
      fx(2,:,:) = x**2*gfac*Vr**2*u(1,:,:) - mu_c*x**2*gfac*divVi
     $     - 0.25*mu_i*x*(x*gfac*Vrx - Vtauy - gfac*Vr)

      fy(2,:,:) = x*Vr*Vtau*u(1,:,:)
     $     - 0.25*mu_i*(Vry/gfac + x*Vtaux - gfac**2*Vtau)

      s(2,:,:) = - kfac*u(6,:,:)*ux(6,:,:) 
     $     - (x*gfac)**2*ux(5,:,:)*(u(4,:,:)-u(1,:,:)*u(7,:,:))/di
     $     + x*gfac**3*Vr**2*u(1,:,:)
     $     + x*gfac**3*(u(3,:,:) - Rinv*x*u(4,:,:))*(Vtau - Rinv*x*Ve)
     $     - x**2*gfac*(ux(8,:,:) + ux(9,:,:)) 
     $     + mu_i*x*gfac**2*(2.*Rinv*Vey - 0.5*Rinv**2*x*gfac*Vr)
     $     - mu_c*x*gfac**3*(2. - (Rinv*x)**2)*divVi
     $     + nu_e*2.*Rinv*x*gfac**2*uy(7,:,:)
c-----------------------------------------------------------------------
c     tau-Momentum Eq.
c-----------------------------------------------------------------------
      fx(3,:,:) = x**2*gfac*Vr*Vtau*u(1,:,:)
     $     - 0.25*mu_i*x*(Vry + x*gfac*Vtaux - gfac**3*Vtau)

      fy(3,:,:) = x*Vtau**2*u(1,:,:) + 0.5*x*gfac**2*u(6,:,:)**2
     $     + x*(u(8,:,:) + u(9,:,:)) - mu_c*x*divVi
     $     - 0.25*mu_i*(Vtauy/gfac - x*Vrx 
     $     + gfac**2*(1. - (Rinv*x)**2)*Vr)
      
      s(3,:,:) = - x*gfac*uy(5,:,:)*(u(4,:,:) - u(1,:,:)*u(7,:,:))/di
     $     + 2*Rinv*kfac*u(2,:,:)*Ve - mu_i*(2.*Rinv*kfac
     $     *(Vex - hfac*Ve) + (Rinv*gfac)**2*kfac*Vtau)
     $     - nu_e*2.*Rinv*kfac*(ux(7,:,:) - hfac*u(7,:,:))
c-----------------------------------------------------------------------
c     e-Momentum Eq.
c-----------------------------------------------------------------------
      fx(4,:,:) = x/gfac*u(2,:,:)*Ve + u(6,:,:)*uy(5,:,:) 
     $     - mu_i*x*(0.5*Rinv*gfac*Vtau + (Vex - hfac*Ve)/gfac)
     $     - nu_e*x/gfac*(ux(7,:,:) - hfac*u(7,:,:))

      fy(4,:,:) = u(3,:,:)*Ve/gfac**2 - u(6,:,:)*ux(5,:,:) 
     $     + mu_i*(0.5*Rinv*Vr - Vey*mfac) - nu_e*mfac*uy(7,:,:)
c-----------------------------------------------------------------------
c     e-component of Ohm's Law
c-----------------------------------------------------------------------
      fx(5,:,:) = - nu_e*di*x/gfac*n_inv*(ux(7,:,:) - hfac*u(7,:,:))
      
      fy(5,:,:) = - nu_e*di*n_inv*mfac*uy(7,:,:)

      s(5,:,:) = - x*Vr*ux(5,:,:) - Vtau*uy(5,:,:)/gfac 
     $     + di*(ux(5,:,:)*uy(6,:,:) - uy(5,:,:)*ux(6,:,:))/u(1,:,:)
     $     + eta_local*x/gfac*(1.96*(u(4,:,:) - u(1,:,:)*u(7,:,:))/di
     $     - .96*Jpar*be)
     $     - nu_e*di*x/gfac*(nx_inv*(ux(7,:,:) - hfac*u(7,:,:))
     $     + ny_inv*lfac**2*uy(7,:,:))
c-----------------------------------------------------------------------
c     e-component of curl of Ohm's Law+
c-----------------------------------------------------------------------
      fx(6,:,:) = x*gfac**2*u(6,:,:)*Vr + gfac*u(7,:,:)*uy(5,:,:)
     $     - di*gfac**2*u(6,:,:)*uy(6,:,:)/u(1,:,:)
     $     - eta_local*x*gfac*(1.96*gfac*ux(6,:,:) + .96*Jpar*btau)
     $     - nu_e*di*n_inv*2.*Rinv*x*gfac**3*(ux(7,:,:) - hfac*u(7,:,:))

      fy(6,:,:) = gfac*u(6,:,:)*Vtau - gfac*u(7,:,:)*ux(5,:,:)
     $     + di*gfac**2*u(6,:,:)*ux(6,:,:)/u(1,:,:)
     $     - eta_local*(1.96*uy(6,:,:)/x - .96*Jpar*br)
     $     - nu_e*di*n_inv*2.*Rinv*gfac/x*uy(7,:,:)

      s(6,:,:) = di*(uy(9,:,:)*nx_inv - ux(9,:,:)*ny_inv)
c-----------------------------------------------------------------------
c     e-component of Ampere's Law
c-----------------------------------------------------------------------
      fx(7,:,:) = di*x*gfac**2*ux(5,:,:)

      fy(7,:,:) = di/x*uy(5,:,:)

      s(7,:,:) = 2.*di*Rinv*x*gfac**4*u(6,:,:)
     $     + x*gfac*(u(4,:,:) - u(1,:,:)*u(7,:,:))
c-----------------------------------------------------------------------
c     ion pressure equation
c-----------------------------------------------------------------------
      fx(8,:,:) = 1.5*x*u(8,:,:)*Vr 
     $     - (kap_pari - kap_prpi)*br*(x*br*Tix + btau/gfac*Tiy) 
     $     - kap_hati*be*Tiy/gfac - kap_prpi*x*Tix

      fy(8,:,:) = (1.5*u(8,:,:)*Vtau
     $     - (kap_pari - kap_prpi)*btau*(br*Tix + btau*lfac*Tiy)
     $     + kap_hati*be*Tix - kap_prpi*lfac*Tiy)/gfac

      s(8,:,:) = - u(8,:,:)*x*divVi + mu_c*x*divVi**2
     $     + mu_i*(0.25*x*(Vtaux + Vry*lfac - gfac**2/x*Vtau)**2
     $     + 0.25*x*(Vrx - Vtauy*lfac - gfac**2/x*Vr)**2
     $     + x*(Vex + Rinv*gfac**2*Vtau - hfac*Ve)**2
     $     + x*(Vey*lfac - Rinv*gfac**2*Vr)**2
     $     + 0.5*Rinv*gfac**2
     $     *(x*Vtau*Vex - x*hfac*Vtau*Ve - Vr*Vey/gfac)
     $     - 0.25*x*hfac*Vr*(divVi + hfac*Vr))
     $     + x*heatfl*(u(9,:,:) - u(8,:,:))
c-----------------------------------------------------------------------
c     electron pressure equation
c-----------------------------------------------------------------------
      fx(9,:,:) = 1.5*x*u(9,:,:)*Vr - 1.5*di*Te*uy(6,:,:)
     $     - (kap_pare - kap_prpe)*br*(x*br*Tex + btau/gfac*Tey) 
     $     - kap_hate*be*Tey/gfac - kap_prpe*x*Tex

      fy(9,:,:) = 1.5*di*Te*ux(6,:,:) + (1.5*u(9,:,:)*Vtau
     $     - (kap_pare - kap_prpe)*btau*(br*Tex + btau*lfac*Tey)
     $     + kap_hate*be*Tex - kap_prpe*lfac*Tey)/gfac

      s(9,:,:) = - u(9,:,:)*(Vr + x*Vrx + Vtauy/gfac 
     $     + di*(ux(6,:,:)*ny_inv - uy(6,:,:)*nx_inv))
     $     + x*heatfl*(u(8,:,:) - u(9,:,:))
     $     + nu_e*(x*(ux(7,:,:) - hfac*u(7,:,:))
     $     *(ux(7,:,:) - hfac*u(7,:,:) + 2.*Rinv*gfac**2
     $     *n_inv*(u(3,:,:) + di*gfac*ux(6,:,:)))
     $     + uy(7,:,:)*gfac*mfac*(uy(7,:,:) - 2.*Rinv*gfac**3
     $     *n_inv*(x*u(2,:,:) - di*uy(6,:,:))))
     $     + eta_local*x*(1.96*((uy(6,:,:)/x)**2 + (gfac*ux(6,:,:))**2
     $     + (u(4,:,:) - u(1,:,:)*u(7,:,:))**2/di**2) - .96*Jpar**2)
ccc     $     + 1.5*x*Te_source*u(1,:,:)
c-----------------------------------------------------------------------
c     if necessary, add equilibrium source term
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE DEFAULT
            CALL program_stop("No sources for init_type = "//init_type)
         END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_rhs
c-----------------------------------------------------------------------
c     subprogram 8. XMHDhlx_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_drdu(t,x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: 
     $     gfac,hfac,kfac,mfac,lfac,
     $     Vr,Vrx,Vry,Vtau,Vtaux,Vtauy,Ve,Vex,Vey,Btot,divVi,
     $     Vr_u1,Vrx_u1,Vry_u1,Vtau_u1,Vtaux_u1,Vtauy_u1,
     $     Ve_u1,Vex_u1,Vey_u1,
     $     Btot_u6,Btot_ux5,Btot_uy5,divVi_u1,divVi_u2,divVi_u3,
     $     divVi_ux1,divVi_ux2,divVi_uy1,divVi_uy3,
     $     Ti,Tix,Tiy,Ti_u1,Tix_u1,Tiy_u1,
     $     Te,Tex,Tey,Te_u1,Tex_u1,Tey_u1,
     $     br,btau,be,br_u6,br_ux5,br_uy5,btau_u6,btau_ux5,btau_uy5,
     $     be_u6,be_ux5,be_uy5,Jpar,Jpar_u6,Jpar_ux5,Jpar_uy5,
     $     eta_local,eta_local_u1,eta_local_u9,mu_c,mu_c_u1,mu_c_u8,
     $     mu_i,mu_i_u1,mu_i_u8,kap_pari,kap_pari_u1,kap_pari_u8,
     $     kap_hati,kap_prpi,kap_prpi_u1,kap_prpi_u8,
     $     kap_pare,kap_pare_u1,kap_pare_u9,
     $     kap_hate,kap_prpe,kap_prpe_u1,kap_prpe_u9,
     $     heatfl,heatfl_u1,heatfl_u9,nu_e,nu_e_u1,nu_e_u9,
     $     n_inv,nx_inv,ny_inv,Te_source
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
      gfac = 1/SQRT(1+(x*Rinv)**2)
      hfac = x*Rinv**2*gfac**2
      kfac = x**2*gfac**3
      mfac = 1./(x*gfac**3)
      lfac = 1./(x*gfac)
      n_inv = 1./u(1,:,:)
      nx_inv = - ux(1,:,:)/u(1,:,:)**2
      ny_inv = - uy(1,:,:)/u(1,:,:)**2
c-----------------------------------------------------------------------
c     calculate auxiliary temperature variables.
c-----------------------------------------------------------------------
      Ti = u(8,:,:)*n_inv
      Tix = ux(8,:,:)*n_inv + u(8,:,:)*nx_inv
      Tiy = uy(8,:,:)*n_inv + u(8,:,:)*ny_inv
      Te = u(9,:,:)*n_inv
      Tex = ux(9,:,:)*n_inv + u(9,:,:)*nx_inv
      Tey = uy(9,:,:)*n_inv + u(9,:,:)*ny_inv
c-----------------------------------------------------------------------
c     calculate auxiliary temperature derivatives.
c-----------------------------------------------------------------------
      Ti_u1 = - Ti*n_inv
      Tix_u1 = - (ux(8,:,:)*n_inv + 2.*u(8,:,:)*nx_inv)*n_inv
      Tiy_u1 = - (uy(8,:,:)*n_inv + 2.*u(8,:,:)*ny_inv)*n_inv

      Te_u1 = - Te*n_inv
      Tex_u1 = - (ux(9,:,:)*n_inv + 2.*u(9,:,:)*nx_inv)*n_inv
      Tey_u1 = - (uy(9,:,:)*n_inv + 2.*u(9,:,:)*ny_inv)*n_inv
c-----------------------------------------------------------------------
c     calculate auxiliary B-field variables.
c-----------------------------------------------------------------------
      Btot = SQRT((uy(5,:,:)/x)**2 + (gfac*ux(5,:,:))**2 
     $     + (gfac*u(6,:,:))**2)
      br = - uy(5,:,:)/(x*Btot)
      btau = gfac*ux(5,:,:)/Btot
      be = gfac*u(6,:,:)/Btot
      Jpar = br/x*uy(6,:,:) - gfac*btau*ux(6,:,:) 
     $     + be/di*(u(4,:,:) - u(1,:,:)*u(7,:,:))
c-----------------------------------------------------------------------
c     calculate auxiliary B-field derivatives.
c-----------------------------------------------------------------------
      Btot_u6 = gfac**2*u(6,:,:)/Btot
      Btot_ux5 = gfac**2*ux(5,:,:)/Btot
      Btot_uy5 = uy(5,:,:)/(x**2*Btot)
      br_u6 = uy(5,:,:)*Btot_u6/(x*Btot**2)
      br_ux5 = uy(5,:,:)*Btot_ux5/(x*Btot**2)
      br_uy5 = (uy(5,:,:)*Btot_uy5/Btot - 1.)/(x*Btot)
      btau_u6 = - gfac*ux(5,:,:)*Btot_u6/Btot**2
      btau_ux5 = gfac*(1. - ux(5,:,:)*Btot_ux5/Btot)/Btot
      btau_uy5 = - gfac*ux(5,:,:)*Btot_uy5/Btot**2
      be_u6 = gfac*(1. - u(6,:,:)*Btot_u6/Btot)/Btot
      be_ux5 = - gfac*u(6,:,:)*Btot_ux5/Btot**2
      be_uy5 = - gfac*u(6,:,:)*Btot_uy5/Btot**2

      Jpar_u6 = br_u6/x*uy(6,:,:) - gfac*btau_u6*ux(6,:,:) 
     $     + be_u6/di*(u(4,:,:) - u(1,:,:)*u(7,:,:))
      Jpar_ux5 = br_ux5/x*uy(6,:,:) - gfac*btau_ux5*ux(6,:,:) 
     $     + be_ux5/di*(u(4,:,:) - u(1,:,:)*u(7,:,:))
      Jpar_uy5 = br_uy5/x*uy(6,:,:) - gfac*btau_uy5*ux(6,:,:) 
     $     + be_uy5/di*(u(4,:,:) - u(1,:,:)*u(7,:,:))
c-----------------------------------------------------------------------
c     calculate auxiliary velocity variables.
c-----------------------------------------------------------------------
      Vr = u(2,:,:)*n_inv
      Vrx = ux(2,:,:)*n_inv + u(2,:,:)*nx_inv
      Vry = uy(2,:,:)*n_inv + u(2,:,:)*ny_inv
      Vtau = u(3,:,:)*n_inv
      Vtaux = ux(3,:,:)*n_inv + u(3,:,:)*nx_inv
      Vtauy = uy(3,:,:)*n_inv + u(3,:,:)*ny_inv
      Ve = u(4,:,:)*n_inv
      Vex = ux(4,:,:)*n_inv + u(4,:,:)*nx_inv
      Vey = uy(4,:,:)*n_inv + u(4,:,:)*ny_inv
      divVi = Vr/x + Vrx + Vtauy*lfac
c-----------------------------------------------------------------------
c     calculate auxiliary velocity derivatives.
c-----------------------------------------------------------------------
      Vr_u1 = - Vr*n_inv
      Vrx_u1 = - (ux(2,:,:)*n_inv + 2.*u(2,:,:)*nx_inv)*n_inv
      Vry_u1 = - (uy(2,:,:)*n_inv + 2.*u(2,:,:)*ny_inv)*n_inv
      Vtau_u1 = - Vtau*n_inv
      Vtaux_u1 = - (ux(3,:,:)*n_inv + 2.*u(3,:,:)*nx_inv)*n_inv
      Vtauy_u1 = - (uy(3,:,:)*n_inv + 2.*u(3,:,:)*ny_inv)*n_inv
      Ve_u1 = - Ve*n_inv
      Vex_u1 = - (ux(4,:,:)*n_inv + 2.*u(4,:,:)*nx_inv)*n_inv
      Vey_u1 = - (uy(4,:,:)*n_inv + 2.*u(4,:,:)*ny_inv)*n_inv

      divVi_u1 = Vr_u1/x + Vrx_u1 + Vtauy_u1*lfac
      divVi_u2 = n_inv/x + nx_inv
      divVi_u3 = ny_inv*lfac
      divVi_ux1 = Vr_u1
      divVi_ux2 = n_inv
      divVi_uy1 = Vtau_u1*lfac
      divVi_uy3 = n_inv*lfac      
c-----------------------------------------------------------------------
c     define diffusion coefficients
c-----------------------------------------------------------------------
      mu_c = 0.96*di*mu*Ti**2.5
      mu_c_u1 = - 2.5*mu_c/u(1,:,:)
      mu_c_u8 = 2.5*mu_c/u(8,:,:)

      mu_i = 1.2*di/mu*u(1,:,:)**2.5/SQRT(u(8,:,:))
      mu_i_u1 = 2.5*mu_i/u(1,:,:)
      mu_i_u8 = - 0.5*mu_i/u(8,:,:)

      kap_pari = 3.9*di*mu*Ti**2.5
      kap_pari_u1 = - 2.5*kap_pari/u(1,:,:)
      kap_pari_u8 = 2.5*kap_pari/u(8,:,:)

      kap_hati = 2.5*di*u(8,:,:)

      kap_prpi = 2.*di/mu*u(1,:,:)**2.5/SQRT(u(8,:,:))
      kap_prpi_u1 = 2.5*kap_prpi/u(1,:,:)
      kap_prpi_u8 = - 0.5*kap_prpi/u(8,:,:)

      heatfl = heatfdg*3.*u(1,:,:)**2.5/(di*nu*u(9,:,:)**1.5)
      heatfl_u1 = 2.5*heatfl/u(1,:,:)
      heatfl_u9 = - 1.5*heatfl/u(9,:,:)

      eta_local = di/(1.96*nu*Te**1.5)
      eta_local_u1 = 1.5*eta_local/u(1,:,:)
      eta_local_u9 = - 1.5*eta_local/u(9,:,:)

      nu_e = 2.*mass_r*di/nu*u(1,:,:)**2.5/SQRT(u(9,:,:))
      nu_e_u1 = 2.5*nu_e/u(1,:,:)
      nu_e_u9 = - 0.5*nu_e/u(9,:,:)

      kap_pare = 3.2*di*nu*Te**2.5
      kap_pare_u1 = - 2.5*kap_pare/u(1,:,:)
      kap_pare_u9 = 2.5*kap_pare/u(9,:,:)

      kap_hate = - 2.5*di*u(9,:,:)

      kap_prpe = 4.7*di/nu*u(1,:,:)**2.5/SQRT(u(9,:,:))
      kap_prpe_u1 = 2.5*kap_prpe/u(1,:,:)
      kap_prpe_u9 = - 0.5*kap_prpe/u(9,:,:)
c-----------------------------------------------------------------------
c     define source terms.
c-----------------------------------------------------------------------
      Te_source = 0
      SELECT CASE(init_type)
      CASE("sawtooth1","sawtooth2")
         Te_source = Tconst*EXP(-x**4/r_s**4)
      END SELECT
c-----------------------------------------------------------------------
c     Helical MHD model.
c     Density Eq.
c-----------------------------------------------------------------------
      fx_u(1,2,:,:) = x
      fy_u(1,3,:,:) = 1./gfac
c-----------------------------------------------------------------------
c     r-Momentum Eq.
c-----------------------------------------------------------------------
      fx_u(2,1,:,:) = - x**2*gfac*Vr**2 
     $     - mu_c*x**2*gfac*divVi_u1 - mu_c_u1*x**2*gfac*divVi
     $     - 0.25*mu_i*x*(x*gfac*Vrx_u1 - Vtauy_u1 - gfac*Vr_u1)
     $     - 0.25*mu_i_u1*x*(x*gfac*Vrx - Vtauy - gfac*Vr)
      fx_u(2,2,:,:) = gfac*(2*x**2*Vr - mu_c*x**2*divVi_u2
     $     - 0.25*mu_i*x*(x*nx_inv - n_inv))
      fx_u(2,3,:,:) = - mu_c*x**2*gfac*divVi_u3
     $     + 0.25*mu_i*x*ny_inv
      fx_u(2,8,:,:) =  - mu_c_u8*x**2*gfac*divVi
     $     - 0.25*mu_i_u8*x*(x*gfac*Vrx - Vtauy - gfac*Vr)

      fx_ux(2,1,:,:) = - gfac*x**2*(mu_c*divVi_ux1 + 0.25*mu_i*Vr_u1)
      fx_ux(2,2,:,:) = - gfac*x**2*(mu_c*divVi_ux2 + 0.25*mu_i*n_inv)
      fx_uy(2,1,:,:) = - mu_c*x**2*gfac*divVi_uy1 + 0.25*mu_i*x*Vtau_u1
      fx_uy(2,3,:,:) = - mu_c*x**2*gfac*divVi_uy3 + 0.25*mu_i*x*n_inv

      fy_u(2,1,:,:) = - x*Vr*Vtau
     $     - 0.25*mu_i*(Vry_u1/gfac + x*Vtaux_u1 - gfac**2*Vtau_u1)
     $     - 0.25*mu_i_u1*(Vry/gfac + x*Vtaux - gfac**2*Vtau)
      fy_u(2,2,:,:) = x*Vtau - 0.25*mu_i*ny_inv/gfac
      fy_u(2,3,:,:) = x*Vr - 0.25*mu_i*(x*nx_inv - gfac**2*n_inv)
      fy_u(2,8,:,:) = - 0.25*mu_i_u8*(Vry/gfac + x*Vtaux - gfac**2*Vtau)
      fy_ux(2,1,:,:) = - 0.25*mu_i*x*Vtau_u1
      fy_ux(2,3,:,:) = - 0.25*mu_i*x*n_inv
      fy_uy(2,1,:,:) = - 0.25*mu_i*Vr_u1/gfac
      fy_uy(2,2,:,:) = - 0.25*mu_i*n_inv/gfac
      
      s_u(2,1,:,:) = (x*gfac)**2*ux(5,:,:)*u(7,:,:)/di
     $     - x*gfac**3*Vr**2
     $     + x*gfac**3*(u(3,:,:) - Rinv*x*u(4,:,:))
     $     *(Vtau_u1 - Rinv*x*Ve_u1)
     $     + mu_i*x*gfac**2*(2*Rinv*Vey_u1 - 0.5*Rinv**2*x*gfac*Vr_u1)
     $     + mu_i_u1*x*gfac**2*(2.*Rinv*Vey - 0.5*Rinv**2*x*gfac*Vr)
     $     - mu_c*x*gfac**3*(2.-(Rinv*x)**2)*divVi_u1 
     $     - mu_c_u1*x*gfac**3*(2. - (Rinv*x)**2)*divVi
     $     + nu_e_u1*2.*Rinv*x*gfac**2*uy(7,:,:)
      s_u(2,2,:,:) = x*gfac**3*(2.*Vr - 0.5*mu_i*Rinv**2*x*n_inv
     $     - mu_c*(2. - (Rinv*x)**2)*divVi_u2)
      s_u(2,3,:,:) = x*gfac**3*(2.*(u(3,:,:) - Rinv*x*u(4,:,:))*n_inv
     $     - mu_c*(2.-(Rinv*x)**2)*divVi_u3)
      s_u(2,4,:,:) = - (x*gfac)**2*ux(5,:,:)/di
     $     - 2.*Rinv*kfac*(u(3,:,:) - Rinv*x*u(4,:,:))*n_inv
     $     + mu_i*2*Rinv*x*gfac**2*ny_inv
      s_u(2,6,:,:) = - kfac*ux(6,:,:)
      s_u(2,7,:,:) = (x*gfac)**2*ux(5,:,:)*u(1,:,:)/di
      s_u(2,8,:,:) = mu_i_u8*(2.*Rinv*x*gfac**2*Vey 
     $     - 0.5*(Rinv*x)**2*gfac**3*Vr)
     $     - mu_c_u8*x*gfac**3*(2.-(Rinv*x)**2)*divVi
      s_u(2,9,:,:) = nu_e_u9*2.*Rinv*x*gfac**2*uy(7,:,:)
      s_ux(2,1,:,:) = - mu_c*x*gfac**3*(2.-(Rinv*x)**2)*divVi_ux1
      s_ux(2,2,:,:) = - mu_c*x*gfac**3*(2.-(Rinv*x)**2)*divVi_ux2
      s_ux(2,5,:,:) = - (x*gfac)**2*(u(4,:,:)-u(1,:,:)*u(7,:,:))/di 
      s_ux(2,6,:,:) = - kfac*u(6,:,:)
      s_ux(2,8,:,:) = - x**2*gfac
      s_ux(2,9,:,:) = - x**2*gfac
      s_uy(2,1,:,:) = mu_i*2*Rinv*x*gfac**2*Ve_u1
     $     - mu_c*x*gfac**3*(2.-(Rinv*x)**2)*divVi_uy1
      s_uy(2,3,:,:) = - mu_c*x*gfac**3*(2.-(Rinv*x)**2)*divVi_uy3
      s_uy(2,4,:,:) = mu_i*2.*Rinv*x*gfac**2*n_inv
      s_uy(2,7,:,:) = nu_e*2.*Rinv*x*gfac**2
c-----------------------------------------------------------------------
c     tau-Momentum Eq.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:) = - x**2*gfac*Vr*Vtau
     $     - 0.25*mu_i*x*(Vry_u1 + x*gfac*Vtaux_u1 - gfac**3*Vtau_u1)
     $     - 0.25*mu_i_u1*x*(Vry + x*gfac*Vtaux - gfac**3*Vtau)
      fx_u(3,2,:,:) = x**2*gfac*Vtau
     $     - 0.25*mu_i*x*ny_inv
      fx_u(3,3,:,:) = x**2*gfac*Vr
     $     - 0.25*mu_i*x*(x*gfac*nx_inv - gfac**3*n_inv)
      fx_u(3,8,:,:) = - 0.25*mu_i_u8*x
     $     *(Vry + x*gfac*Vtaux - gfac**3*Vtau)
      fx_ux(3,1,:,:) = - 0.25*mu_i*x**2*gfac*Vtau_u1
      fx_ux(3,3,:,:) = - 0.25*mu_i*x**2*gfac*n_inv
      fx_uy(3,1,:,:) = - 0.25*mu_i*x*Vr_u1
      fx_uy(3,2,:,:) = - 0.25*mu_i*x*n_inv

      fy_u(3,1,:,:) = - x*Vtau**2 - mu_c*x*divVi_u1 - mu_c_u1*x*divVi
     $     - 0.25*mu_i*(Vtauy_u1/gfac - x*Vrx_u1 
     $     + gfac**2*(1. - (Rinv*x)**2)*Vr_u1)
     $     - 0.25*mu_i_u1*(Vtauy/gfac - x*Vrx 
     $     + gfac**2*(1. - (Rinv*x)**2)*Vr)
      fy_u(3,2,:,:) = - mu_c*x*divVi_u2
     $     + 0.25*mu_i*(x*nx_inv - gfac**2*(1. - (Rinv*x)**2)*n_inv)
      fy_u(3,3,:,:) = 2*x*Vtau - mu_c*x*divVi_u3
     $     - 0.25*mu_i*ny_inv/gfac
      fy_u(3,6,:,:) = x*gfac**2*u(6,:,:)
      fy_u(3,8,:,:) = x  - mu_c_u8*x*divVi
     $     - 0.25*mu_i_u8*(Vtauy/gfac - x*Vrx 
     $     + gfac**2*(1. - (Rinv*x)**2)*Vr)
      fy_u(3,9,:,:) = x
      fy_ux(3,1,:,:) = - mu_c*x*divVi_ux1 + 0.25*mu_i*x*Vr_u1
      fy_ux(3,2,:,:) = - mu_c*x*divVi_ux2 + 0.25*mu_i*x*n_inv
      fy_uy(3,1,:,:) = - mu_c*x*divVi_uy1 - 0.25*mu_i*Vtau_u1/gfac
      fy_uy(3,3,:,:) = - mu_c*x*divVi_uy3 - 0.25*mu_i*n_inv/gfac
      
      s_u(3,1,:,:) = x*gfac*uy(5,:,:)*u(7,:,:)/di
     $     + 2.*Rinv*kfac*u(2,:,:)*Ve_u1
     $     - mu_i*(2*Rinv*kfac*(Vex_u1 - hfac*Ve_u1)
     $     + (Rinv*gfac)**2*kfac*Vtau_u1)
     $     - mu_i_u1*(2.*Rinv*kfac
     $     *(Vex - hfac*Ve) + (Rinv*gfac)**2*kfac*Vtau)
     $     - nu_e_u1*2.*Rinv*kfac*(ux(7,:,:) - hfac*u(7,:,:))
      s_u(3,2,:,:) = 2.*Rinv*kfac*Ve
      s_u(3,3,:,:) = - mu_i*(Rinv*gfac)**2*kfac*n_inv
      s_u(3,4,:,:) = - x*gfac*uy(5,:,:)/di 
     $     + 2.*Rinv*kfac*u(2,:,:)*n_inv
     $     - mu_i*2*Rinv*kfac*(nx_inv - hfac*n_inv)
      s_u(3,7,:,:) = x*gfac*u(1,:,:)*uy(5,:,:)/di
     $     + nu_e*2.*Rinv*hfac*kfac
      s_u(3,8,:,:) = - mu_i_u8*(2.*Rinv*kfac
     $     *(Vex - hfac*Ve) + (Rinv*gfac)**2*kfac*Vtau)
      s_u(3,9,:,:) = - nu_e_u9*2.*Rinv*kfac*(ux(7,:,:) - hfac*u(7,:,:))
      s_ux(3,1,:,:) = - mu_i*2.*Rinv*kfac*Ve_u1
      s_ux(3,4,:,:) = - mu_i*2.*Rinv*kfac*n_inv
      s_ux(3,7,:,:) = - nu_e*2.*Rinv*kfac
      s_uy(3,5,:,:) = - x*gfac*(u(4,:,:) - u(1,:,:)*u(7,:,:))/di
c-----------------------------------------------------------------------
c     e-Momentum Eq.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:) = - x/gfac*Vr*Ve - mu_i*x*(0.5*Rinv*gfac*Vtau_u1
     $     + (Vex_u1 - hfac*Ve_u1)/gfac)
     $     - mu_i_u1*x*(0.5*Rinv*gfac*Vtau + (Vex - hfac*Ve)/gfac)
     $     - nu_e_u1*x/gfac*(ux(7,:,:) - x*(Rinv*gfac)**2*u(7,:,:))
      fx_u(4,2,:,:) = x/gfac*Ve
      fx_u(4,3,:,:) = - mu_i*0.5*Rinv*x*gfac*n_inv
      fx_u(4,4,:,:) = x/gfac*Vr
     $     - mu_i*x/gfac*(nx_inv - hfac*n_inv)
      fx_u(4,6,:,:) = uy(5,:,:)
      fx_u(4,7,:,:) = nu_e*Rinv**2*x**2*gfac
      fx_u(4,8,:,:) = - mu_i_u8*x*(0.5*Rinv*gfac*Vtau
     $     + (Vex - hfac*Ve)/gfac)
      fx_u(4,9,:,:) = - nu_e_u9*x/gfac*(ux(7,:,:) - hfac*u(7,:,:))
      fx_ux(4,1,:,:) = - mu_i*x/gfac*Ve_u1
      fx_ux(4,4,:,:) = - mu_i*x/gfac*n_inv
      fx_ux(4,7,:,:) = - nu_e*x/gfac
      fx_uy(4,5,:,:) = u(6,:,:)

      fy_u(4,1,:,:) = - Vtau*Ve/gfac**2 - nu_e_u1*mfac*uy(7,:,:)
     $     + mu_i*(0.5*Rinv*Vr_u1 - Vey_u1*mfac)
     $     + mu_i_u1*(0.5*Rinv*Vr - Vey*mfac)
      fy_u(4,2,:,:) = mu_i*0.5*Rinv*n_inv
      fy_u(4,3,:,:) = Ve/gfac**2
      fy_u(4,4,:,:) = Vtau/gfac**2 - mu_i*mfac*ny_inv
      fy_u(4,6,:,:) = - ux(5,:,:)
      fy_u(4,8,:,:) = mu_i_u8*(0.5*Rinv*Vr - Vey*mfac)
      fy_u(4,9,:,:) = - nu_e_u9*mfac*uy(7,:,:)
      fy_ux(4,5,:,:) = - u(6,:,:)
      fy_uy(4,1,:,:) = - mu_i*mfac*Ve_u1
      fy_uy(4,4,:,:) = - mu_i*mfac*n_inv
      fy_uy(4,7,:,:) = - nu_e*mfac
c-----------------------------------------------------------------------
c     e-component of Ohm's Law
c-----------------------------------------------------------------------
      fx_u(5,1,:,:) = di*x/gfac*n_inv*(ux(7,:,:) - hfac*u(7,:,:))
     $     *(nu_e*n_inv - nu_e_u1)
      fx_u(5,7,:,:) = nu_e*di*Rinv**2*x**2*gfac*n_inv
      fx_u(5,9,:,:) = - nu_e_u9*di*x/gfac*n_inv
     $     *(ux(7,:,:) - hfac*u(7,:,:))
      fx_ux(5,7,:,:) = - nu_e*di*x/gfac*n_inv

      fy_u(5,1,:,:) = di*n_inv*mfac*uy(7,:,:)*(nu_e*n_inv - nu_e_u1)
      fy_u(5,9,:,:) = - nu_e_u9*di*n_inv*mfac*uy(7,:,:)
      fy_uy(5,7,:,:) = - nu_e*di*n_inv*mfac

      s_u(5,1,:,:) = - x*Vr_u1*ux(5,:,:) - Vtau_u1*uy(5,:,:)/gfac 
     $     - di*(ux(5,:,:)*uy(6,:,:) - uy(5,:,:)*ux(6,:,:))/u(1,:,:)**2
     $     - eta_local*x/gfac*u(7,:,:)*(1.96 - .96*be**2)/di
     $     + eta_local_u1*x/gfac
     $     *(1.96*(u(4,:,:) - u(1,:,:)*u(7,:,:))/di - .96*Jpar*be)
     $     + di*x/gfac*(nx_inv*(ux(7,:,:) - hfac*u(7,:,:))
     $     + ny_inv*lfac**2*uy(7,:,:))*(2.*nu_e*n_inv - nu_e_u1)
      s_u(5,2,:,:) = - x*n_inv*ux(5,:,:)
      s_u(5,3,:,:) = - n_inv*uy(5,:,:)/gfac
      s_u(5,4,:,:) = eta_local*x/gfac*(1.96 - .96*be**2)/di
      s_u(5,6,:,:) = - .96*eta_local*x/gfac*(Jpar_u6*be + Jpar*be_u6)
      s_u(5,7,:,:) = - eta_local*x/gfac*u(1,:,:)*(1.96 - .96*be**2)/di
     $     + nu_e*di*x**2*Rinv**2*gfac*nx_inv
      s_u(5,9,:,:) = eta_local_u9*x/gfac
     $     *(1.96*(u(4,:,:) - u(1,:,:)*u(7,:,:))/di - .96*Jpar*be)
     $     - nu_e_u9*di*x/gfac*(nx_inv*(ux(7,:,:) - hfac*u(7,:,:)) 
     $     + ny_inv*lfac**2*uy(7,:,:))
      s_ux(5,1,:,:) = nu_e*di*x/gfac*n_inv**2
     $     *(ux(7,:,:) - hfac*u(7,:,:))
      s_ux(5,5,:,:) = - x*Vr + di*uy(6,:,:)/u(1,:,:)
     $     - .96*eta_local*x/gfac*(Jpar_ux5*be + Jpar*be_ux5)
      s_ux(5,6,:,:) = - di*uy(5,:,:)/u(1,:,:)
     $     + .96*eta_local*x*btau*be
      s_ux(5,7,:,:) = - nu_e*di*x/gfac*nx_inv
      s_uy(5,1,:,:) = nu_e*di*mfac*n_inv**2*uy(7,:,:)
      s_uy(5,5,:,:) = - Vtau/gfac - di*ux(6,:,:)/u(1,:,:)
     $     - .96*eta_local*x/gfac*(Jpar_uy5*be + Jpar*be_uy5)
      s_uy(5,6,:,:) = di*ux(5,:,:)/u(1,:,:)
     $     - .96*eta_local/gfac*br*be
      s_uy(5,7,:,:) = - nu_e*di*mfac*ny_inv
c-----------------------------------------------------------------------
c     e-component of curl of Ohm's Law+
c-----------------------------------------------------------------------
      fx_u(6,1,:,:) = x*gfac**2*u(6,:,:)*Vr_u1
     $     + di*gfac**2*u(6,:,:)*uy(6,:,:)/u(1,:,:)**2
     $     + .96*eta_local*x*gfac*be*btau*u(7,:,:)/di
     $     - eta_local_u1*x*gfac*(1.96*gfac*ux(6,:,:) + .96*Jpar*btau)
     $     + di*n_inv*2.*Rinv*x*gfac**3*(ux(7,:,:) - hfac*u(7,:,:))
     $     *(nu_e*n_inv - nu_e_u1)
      fx_u(6,2,:,:) = x*gfac**2*u(6,:,:)*n_inv
      fx_u(6,4,:,:) = - .96*eta_local*x*gfac*be*btau/di
      fx_u(6,6,:,:) = x*gfac**2*Vr - di*gfac**2*uy(6,:,:)/u(1,:,:)
     $     - .96*eta_local*x*gfac*(Jpar_u6*btau + Jpar*btau_u6)
      fx_u(6,7,:,:) = gfac*uy(5,:,:)
     $     + .96*eta_local*x*gfac*be*btau*u(1,:,:)/di
     $     + nu_e*di*n_inv*2.*Rinv**3*gfac**2*kfac
      fx_u(6,9,:,:) = - nu_e_u9*di*n_inv*2.*Rinv*x*gfac**3
     $     *(ux(7,:,:) - hfac*u(7,:,:))
     $     - eta_local_u9*x*gfac*(1.96*gfac*ux(6,:,:) + .96*Jpar*btau)
      fx_ux(6,5,:,:) =  - .96*eta_local*x*gfac
     $     *(Jpar_ux5*btau + Jpar*btau_ux5)
      fx_ux(6,6,:,:) = - eta_local*x*gfac**2*(1.96 - .96*btau**2)
      fx_ux(6,7,:,:) = - nu_e*di*n_inv*2.*Rinv*x*gfac**3
      fx_uy(6,5,:,:) = gfac*u(7,:,:) - .96*eta_local*x*gfac
     $     *(Jpar_uy5*btau + Jpar*btau_uy5)
      fx_uy(6,6,:,:) = - di*gfac**2*u(6,:,:)/u(1,:,:)
     $      - .96*eta_local*gfac*br*btau

      fy_u(6,1,:,:) = gfac*u(6,:,:)*Vtau_u1
     $     - di*gfac**2*u(6,:,:)*ux(6,:,:)/u(1,:,:)**2
     $     - .96*eta_local*br*be*u(7,:,:)/di
     $     - eta_local_u1*(1.96*uy(6,:,:)/x - .96*Jpar*br)
     $     + di*n_inv*2.*Rinv*gfac/x*uy(7,:,:)*(nu_e*n_inv - nu_e_u1)
      fy_u(6,3,:,:) = gfac*u(6,:,:)*n_inv
      fy_u(6,4,:,:) = .96*eta_local*br*be/di
      fy_u(6,6,:,:) = gfac*Vtau + di*gfac**2*ux(6,:,:)/u(1,:,:)
     $     + .96*eta_local*(Jpar_u6*br + Jpar*br_u6)
      fy_u(6,7,:,:) = - gfac*ux(5,:,:) - .96*eta_local*br*be*u(1,:,:)/di
      fy_u(6,9,:,:) = - nu_e_u9*di*n_inv*2.*Rinv*gfac/x*uy(7,:,:)
     $     - eta_local_u9*(1.96*uy(6,:,:)/x - .96*Jpar*br)
      fy_ux(6,5,:,:) = - gfac*u(7,:,:)
     $     + .96*eta_local*(Jpar_ux5*br + Jpar*br_ux5)
      fy_ux(6,6,:,:) = di*gfac**2*u(6,:,:)/u(1,:,:)
     $     - .96*eta_local*gfac*br*btau
      fy_uy(6,5,:,:) = .96*eta_local*(Jpar_uy5*br + Jpar*br_uy5)
      fy_uy(6,6,:,:) = - eta_local/x*(1.96 - .96*br**2)
      fy_uy(6,7,:,:) = - nu_e*di*n_inv*2.*Rinv*gfac/x

      s_u(6,1,:,:) = - 2.*di*n_inv*(uy(9,:,:)*nx_inv - ux(9,:,:)*ny_inv)
      s_ux(6,1,:,:) = - di*uy(9,:,:)*n_inv**2
      s_ux(6,9,:,:) = - di*ny_inv
      s_uy(6,1,:,:) = di*ux(9,:,:)*n_inv**2
      s_uy(6,9,:,:) = di*nx_inv
c-----------------------------------------------------------------------
c     e-component of Ampere's Law
c-----------------------------------------------------------------------
      fx_ux(7,5,:,:) = di*x*gfac**2

      fy_uy(7,5,:,:) = di/x

      s_u(7,1,:,:) = - x*gfac*u(7,:,:)
      s_u(7,4,:,:) = x*gfac
      s_u(7,6,:,:) = 2*di*Rinv*x*gfac**4
      s_u(7,7,:,:) = - x*gfac*u(1,:,:)
c-----------------------------------------------------------------------
c     ion pressure equation
c-----------------------------------------------------------------------
      fx_u(8,1,:,:) = 1.5*x*u(8,:,:)*Vr_u1
     $     - (kap_pari - kap_prpi)*br*(x*br*Tix_u1 + btau/gfac*Tiy_u1)
     $     - (kap_pari_u1 - kap_prpi_u1)*br*(x*br*Tix + btau/gfac*Tiy)
     $     - kap_hati*be*Tiy_u1/gfac
     $     - x*(kap_prpi_u1*Tix + kap_prpi*Tix_u1)
      fx_u(8,2,:,:) = 1.5*x*u(8,:,:)*n_inv
      fx_u(8,6,:,:) = - (kap_pari - kap_prpi)*(2.*x*br*br_u6*Tix 
     $     + btau*br_u6*Tiy/gfac + br*btau_u6*Tiy/gfac) 
     $     - kap_hati*be_u6*Tiy/gfac
      fx_u(8,8,:,:) = 1.5*x*Vr
     $     - (kap_pari - kap_prpi)*br*(x*br*nx_inv + btau/gfac*ny_inv)
     $     - (kap_pari_u8 - kap_prpi_u8)*br*(x*br*Tix + btau/gfac*Tiy)
     $     - 2.5*di*be/gfac*(u(8,:,:)*ny_inv + Tiy)
     $     - x*(kap_prpi_u8*Tix + kap_prpi*nx_inv)
      fx_ux(8,1,:,:) = - ((kap_pari - kap_prpi)*br**2 + kap_prpi)
     $     *x*Ti_u1
      fx_ux(8,5,:,:) = - (kap_pari - kap_prpi)
     $     *(br_ux5*(x*br*Tix + btau/gfac*Tiy)
     $     + br*(x*br_ux5*Tix + btau_ux5/gfac*Tiy))
     $     - kap_hati*be_ux5*Tiy/gfac
      fx_ux(8,8,:,:) = - ((kap_pari - kap_prpi)*br**2 + kap_prpi)
     $     *x*n_inv
      fx_uy(8,1,:,:) = - ((kap_pari - kap_prpi)*br*btau + kap_hati*be)
     $     *Ti_u1/gfac
      fx_uy(8,5,:,:) = - (kap_pari - kap_prpi)
     $     *(br_uy5*(x*br*Tix + btau/gfac*Tiy)
     $     + br*(x*br_uy5*Tix + btau_uy5/gfac*Tiy))
     $     - kap_hati*be_uy5*Tiy/gfac
      fx_uy(8,8,:,:) = - ((kap_pari - kap_prpi)*br*btau + kap_hati*be)
     $     *n_inv/gfac

      fy_u(8,1,:,:) = (1.5*u(8,:,:)*Vtau_u1
     $     - (kap_pari - kap_prpi)*btau*(br*Tix_u1 + btau*lfac*Tiy_u1) 
     $     - (kap_pari_u1 - kap_prpi_u1)*btau
     $     *(br*Tix + btau*lfac*Tiy) + kap_hati*be*Tix_u1
     $     - (kap_prpi_u1*Tiy + kap_prpi*Tiy_u1)*lfac)/gfac
      fy_u(8,3,:,:) = 1.5/gfac*u(8,:,:)*n_inv
      fy_u(8,6,:,:) = - ((kap_pari - kap_prpi)*(btau_u6*br*Tix 
     $     + 2.*btau*btau_u6*lfac*Tiy + btau*br_u6*Tix) 
     $     - kap_hati*be_u6*Tix)/gfac
      fy_u(8,8,:,:) = (1.5*Vtau
     $     - (kap_pari - kap_prpi)*btau*(br*nx_inv + btau*lfac*ny_inv)
     $     - (kap_pari_u8 - kap_prpi_u8)*btau*(br*Tix + btau*lfac*Tiy)
     $     + 2.5*di*be*(u(8,:,:)*nx_inv + Tix)
     $     - (kap_prpi_u8*Tiy + kap_prpi*ny_inv)*lfac)/gfac
      fy_ux(8,1,:,:) = - ((kap_pari - kap_prpi)*btau*br - kap_hati*be)
     $     *Ti_u1/gfac
      fy_ux(8,5,:,:) = - ((kap_pari - kap_prpi)
     $     *(btau_ux5*(br*Tix + btau*lfac*Tiy)
     $     + btau*(br_ux5*Tix + btau_ux5*lfac*Tiy))
     $     - kap_hati*be_ux5*Tix)/gfac
      fy_ux(8,8,:,:) = - ((kap_pari - kap_prpi)*btau*br - kap_hati*be)
     $     *n_inv/gfac
      fy_uy(8,1,:,:) = - ((kap_pari - kap_prpi)*btau**2 + kap_prpi)
     $     *mfac*gfac*Ti_u1
      fy_uy(8,5,:,:) = - ((kap_pari - kap_prpi)
     $     *(btau_uy5*(br*Tix + btau*lfac*Tiy)
     $     + btau*(br_uy5*Tix + btau_uy5*lfac*Tiy))
     $     - kap_hati*be_uy5*Tix)/gfac
      fy_uy(8,8,:,:) = - ((kap_pari - kap_prpi)*btau**2 + kap_prpi)
     $     *mfac*gfac*n_inv

      s_u(8,1,:,:) = - u(8,:,:)*x*divVi_u1
     $     + 2.*mu_c*x*divVi*divVi_u1 + mu_c_u1*x*divVi**2
     $     - mu_i*(0.25*x*hfac*
     $     (Vr_u1*(divVi + 2.*hfac*Vr) + Vr*divVi_u1)
     $     - 0.5*(x*Vtaux + Vry/gfac - gfac**2*Vtau)
     $     *(Vtaux_u1 + Vry_u1*lfac - gfac**2/x*Vtau_u1)
     $     - 0.5*(x*Vrx - Vtauy/gfac - gfac**2*Vr)
     $     *(Vrx_u1 - Vtauy_u1*lfac - gfac**2/x*Vr_u1)
     $     - 2.*x*(Vex + Rinv*gfac**2*Vtau - hfac*Ve)
     $     *(Vex_u1 + Rinv*gfac**2*Vtau_u1 - hfac*Ve_u1)
     $     - 2.*x*(Vey*lfac - Rinv*gfac**2*Vr)
     $     *(Vey_u1*lfac - Rinv*gfac**2*Vr_u1)
     $     - 0.5*Rinv*gfac**2*(x*(Vtau_u1*Vex + Vtau*Vex_u1) 
     $     - x*hfac*(Vtau_u1*Ve + Vtau*Ve_u1) 
     $     - (Vr_u1*Vey + Vr*Vey_u1)/gfac))
     $     + mu_i_u1*(0.25*x*(Vtaux + Vry*lfac - gfac**2/x*Vtau)**2
     $     + 0.25*x*(Vrx - Vtauy*lfac - gfac**2/x*Vr)**2
     $     + x*(Vex + Rinv*gfac**2*Vtau - hfac*Ve)**2
     $     + x*(Vey*lfac - Rinv*gfac**2*Vr)**2
     $     + 0.5*Rinv*gfac**2
     $     *(x*Vtau*Vex - x*hfac*Vtau*Ve - Vr*Vey/gfac)
     $     - 0.25*x*hfac*Vr*(divVi + hfac*Vr))
     $     + x*heatfl_u1*(u(9,:,:) - u(8,:,:))
      s_u(8,2,:,:) = - u(8,:,:)*x*divVi_u2 + 2.*mu_c*x*divVi*divVi_u2
     $     - mu_i*(0.25*x*hfac
     $     *(n_inv*(divVi + 2.*hfac*Vr) + Vr*divVi_u2)
     $     - 0.5*(x*Vtaux + Vry/gfac - gfac**2*Vtau)*ny_inv*lfac
     $     - 0.5*(x*Vrx - Vtauy/gfac - gfac**2*Vr)
     $     *(nx_inv - gfac**2/x*n_inv)
     $     + 2.*(Vey - Rinv*x*gfac**3*Vr)*Rinv*gfac*n_inv
     $     + 0.5*Rinv*gfac*n_inv*Vey)
      s_u(8,3,:,:) = - u(8,:,:)*x*divVi_u3 + 2.*mu_c*x*divVi*divVi_u3
     $     - mu_i*(0.25*x*hfac*Vr*divVi_u3
     $     - 0.5*(x*Vtaux + Vry/gfac - gfac**2*Vtau)
     $     *(nx_inv - gfac**2/x*n_inv)
     $     + 0.5*(x*Vrx - Vtauy/gfac - gfac**2*Vr)*ny_inv*lfac 
     $     - 2.*x*(Vex + Rinv*gfac**2*Vtau - hfac*Ve)
     $     *Rinv*gfac**2*n_inv 
     $     - 0.5*Rinv*x*gfac**2*n_inv*(Vex - hfac*Ve))
      s_u(8,4,:,:) = mu_i*(2.*x*(Vex + Rinv*gfac**2*Vtau - hfac*Ve)
     $     *(nx_inv - hfac*n_inv)
     $     + 2.*(Vey/gfac - Rinv*x*gfac**2*Vr)*ny_inv*lfac
     $     + 0.5*Rinv*gfac**2
     $     *(x*Vtau*nx_inv - x*hfac*Vtau*n_inv - Vr*ny_inv/gfac))
      s_u(8,8,:,:) = - x*divVi + mu_c_u8*x*divVi**2
     $     + mu_i_u8*(0.25*x*(Vtaux + Vry*lfac - gfac**2/x*Vtau)**2
     $     + 0.25*x*(Vrx - Vtauy*lfac - gfac**2/x*Vr)**2
     $     + x*(Vex + Rinv*gfac**2*Vtau - hfac*Ve)**2
     $     + x*(Vey*lfac - Rinv*gfac**2*Vr)**2
     $     + 0.5*Rinv*gfac**2
     $     *(x*Vtau*Vex - x*hfac*Vtau*Ve - Vr*Vey/gfac)
     $     - 0.25*x*hfac*Vr*(divVi + hfac*Vr))
     $     - x*heatfl
      s_u(8,9,:,:) = x*(heatfl - heatfl_u9*(u(8,:,:) - u(9,:,:)))

      s_ux(8,1,:,:) = - u(8,:,:)*x*divVi_ux1 + 2.*mu_c*x*divVi*divVi_ux1
     $     - mu_i*(0.25*x*hfac*Vr*divVi_ux1
     $     - 0.5*(x*Vtaux + Vry/gfac - gfac**2*Vtau)*Vtau_u1
     $     - 0.5*(x*Vrx - Vtauy/gfac - gfac**2*Vr)*Vr_u1
     $     - 2.*x*(Vex + Rinv*gfac**2*Vtau - hfac*Ve)*Ve_u1
     $     - 0.5*Rinv*gfac**2*x*Vtau*Ve_u1)
      s_ux(8,2,:,:) = - u(8,:,:)*x*divVi_ux2 + 2.*mu_c*x*divVi*divVi_ux2
     $     - mu_i*(0.25*x*hfac*Vr*divVi_ux2
     $     - 0.5*(x*Vrx - Vtauy/gfac - gfac**2*Vr)*n_inv)
      s_ux(8,3,:,:) = mu_i*0.5*(x*Vtaux + Vry/gfac - gfac**2*Vtau)*n_inv
      s_ux(8,4,:,:) = mu_i*(2.*x*(Vex + Rinv*gfac**2*Vtau - hfac*Ve)
     $     *n_inv + 0.5*Rinv*gfac**2*x*Vtau*n_inv)

      s_uy(8,1,:,:) = - u(8,:,:)*x*divVi_uy1 + 2.*mu_c*x*divVi*divVi_uy1
     $     - mu_i*(0.25*x*hfac*Vr*divVi_uy1
     $     - 0.5*(x*Vtaux + Vry/gfac - gfac**2*Vtau)*Vr_u1*lfac
     $     + 0.5*(x*Vrx - Vtauy/gfac - gfac**2*Vr)*Vtau_u1*lfac
     $     - 2.*(Vey/gfac - Rinv*x*gfac**2*Vr)*Ve_u1*lfac
     $     + 0.5*Rinv*gfac*Vr*Ve_u1)
      s_uy(8,2,:,:) = mu_i*0.5*(x*Vtaux + Vry/gfac - gfac**2*Vtau)
     $     *n_inv*lfac
      s_uy(8,3,:,:) = - u(8,:,:)*x*divVi_uy3 + 2.*mu_c*x*divVi*divVi_uy3
     $     - mu_i*(0.25*x*hfac*Vr*divVi_uy3
     $     + 0.5*(x*Vrx - Vtauy/gfac - gfac**2*Vr)*n_inv*lfac)
      s_uy(8,4,:,:) = mu_i*(2.*(Vey/gfac - Rinv*x*gfac**2*Vr)
     $     *n_inv*lfac - 0.5*Rinv*gfac*Vr*n_inv)
c-----------------------------------------------------------------------
c     electron pressure equation
c-----------------------------------------------------------------------
      fx_u(9,1,:,:) = 1.5*x*u(9,:,:)*Vr_u1 - 1.5*di*Te_u1*uy(6,:,:)
     $     - (kap_pare - kap_prpe)*br*(x*br*Tex_u1 + btau/gfac*Tey_u1)
     $     - (kap_pare_u1 - kap_prpe_u1)*br*(x*br*Tex + btau/gfac*Tey)
     $     - kap_hate*be*Tey_u1/gfac
     $     - x*(kap_prpe_u1*Tex + kap_prpe*Tex_u1)
      fx_u(9,2,:,:) = 1.5*x*u(9,:,:)*n_inv
      fx_u(9,6,:,:) = - (kap_pare - kap_prpe)*(2.*x*br*br_u6*Tex 
     $     + btau*br_u6*Tey/gfac + br*btau_u6*Tey/gfac) 
     $     - kap_hate*be_u6*Tey/gfac
      fx_u(9,9,:,:) = 1.5*x*Vr - 1.5*di*n_inv*uy(6,:,:)
     $     - (kap_pare - kap_prpe)*br*(x*br*nx_inv + btau/gfac*ny_inv)
     $     - (kap_pare_u9 - kap_prpe_u9)*br*(x*br*Tex + btau/gfac*Tey)
     $     + 2.5*di*be/gfac*(u(9,:,:)*ny_inv + Tey)
     $     - x*(kap_prpe_u9*Tex + kap_prpe*nx_inv)
      fx_ux(9,1,:,:) = - ((kap_pare - kap_prpe)*br**2 + kap_prpe)
     $     *x*Te_u1
      fx_ux(9,5,:,:) = - (kap_pare - kap_prpe)
     $     *(br_ux5*(x*br*Tex + btau/gfac*Tey)
     $     + br*(x*br_ux5*Tex + btau_ux5/gfac*Tey))
     $     - kap_hate*be_ux5*Tey/gfac
      fx_ux(9,9,:,:) = - ((kap_pare - kap_prpe)*br**2 + kap_prpe)
     $     *x*n_inv
      fx_uy(9,1,:,:) = - ((kap_pare - kap_prpe)*br*btau + kap_hate*be)
     $     *Te_u1/gfac
      fx_uy(9,5,:,:) = - (kap_pare - kap_prpe)
     $     *(br_uy5*(x*br*Tex + btau/gfac*Tey)
     $     + br*(x*br_uy5*Tex + btau_uy5/gfac*Tey))
     $     - kap_hate*be_uy5*Tey/gfac
      fx_uy(9,6,:,:) = - 1.5*di*Te
      fx_uy(9,9,:,:) = - ((kap_pare - kap_prpe)*br*btau + kap_hate*be)
     $     *n_inv/gfac

      fy_u(9,1,:,:) = 1.5*di*Te_u1*ux(6,:,:) + (1.5*u(9,:,:)*Vtau_u1 
     $     - (kap_pare - kap_prpe)*btau*(br*Tex_u1 + btau*lfac*Tey_u1) 
     $     - (kap_pare_u1 - kap_prpe_u1)*btau
     $     *(br*Tex + btau*lfac*Tey) + kap_hate*be*Tex_u1
     $     - (kap_prpe_u1*Tey + kap_prpe*Tey_u1)*lfac)/gfac
      fy_u(9,3,:,:) = 1.5/gfac*u(9,:,:)*n_inv
      fy_u(9,6,:,:) = - ((kap_pare - kap_prpe)*(btau_u6*br*Tex 
     $     + 2.*btau*btau_u6*lfac*Tey + btau*br_u6*Tex) 
     $     - kap_hate*be_u6*Tex)/gfac
      fy_u(9,9,:,:) = 1.5*di*n_inv*ux(6,:,:) + (1.5*Vtau 
     $     - (kap_pare - kap_prpe)*btau*(br*nx_inv + btau*lfac*ny_inv)
     $     - (kap_pare_u9 - kap_prpe_u9)*btau*(br*Tex + btau*lfac*Tey)
     $     - 2.5*di*be*(u(9,:,:)*nx_inv + Tex)
     $     - (kap_prpe_u9*Tey + kap_prpe*ny_inv)*lfac)/gfac
      fy_ux(9,1,:,:) = - ((kap_pare - kap_prpe)*btau*br - kap_hate*be)
     $     *Te_u1/gfac
      fy_ux(9,5,:,:) = - ((kap_pare - kap_prpe)
     $     *(btau_ux5*(br*Tex + btau*lfac*Tey)
     $     + btau*(br_ux5*Tex + btau_ux5*lfac*Tey))
     $     - kap_hate*be_ux5*Tex)/gfac
      fy_ux(9,9,:,:) = - ((kap_pare - kap_prpe)*btau*br - kap_hate*be)
     $     *n_inv/gfac
      fy_uy(9,1,:,:) = - ((kap_pare - kap_prpe)*btau**2 + kap_prpe)
     $     *mfac*gfac*Te_u1
      fy_uy(9,5,:,:) = - ((kap_pare - kap_prpe)
     $     *(btau_uy5*(br*Tex + btau*lfac*Tey)
     $     + btau*(br_uy5*Tex + btau_uy5*lfac*Tey))
     $     - kap_hate*be_uy5*Tex)/gfac
      fy_ux(9,6,:,:) = 1.5*di*Te
      fy_uy(9,9,:,:) = - ((kap_pare - kap_prpe)*btau**2 + kap_prpe)
     $     *mfac*gfac*n_inv

      s_u(9,1,:,:) = - u(9,:,:)*(Vr_u1 + x*Vrx_u1 + Vtauy_u1/gfac 
     $     + 2.*di/u(1,:,:)*(uy(6,:,:)*nx_inv - ux(6,:,:)*ny_inv))
     $     + x*heatfl_u1*(u(8,:,:) - u(9,:,:))
     $     - nu_e*2.*Rinv*n_inv**2
     $     *(x*gfac**2*(ux(7,:,:) - hfac*u(7,:,:))
     $     *(u(3,:,:) + di*gfac*ux(6,:,:))
     $     - uy(7,:,:)*(u(2,:,:) - di/x*uy(6,:,:)))
     $     + nu_e_u1*(x*(ux(7,:,:) - hfac*u(7,:,:))
     $     *(ux(7,:,:) - hfac*u(7,:,:) + 2.*Rinv*gfac**2
     $     *n_inv*(u(3,:,:) + di*gfac*ux(6,:,:)))
     $     + uy(7,:,:)*gfac*mfac*(uy(7,:,:) - 2.*Rinv*gfac**3
     $     *n_inv*(x*u(2,:,:) - di*uy(6,:,:))))
     $     - eta_local*x*u(7,:,:)/di
     $     *(3.92*(u(4,:,:) - u(1,:,:)*u(7,:,:))/di - 1.92*Jpar*be)
     $     + eta_local_u1*x
     $     *(1.96*((uy(6,:,:)/x)**2 + (gfac*ux(6,:,:))**2
     $     + (u(4,:,:) - u(1,:,:)*u(7,:,:))**2/di**2) - .96*Jpar**2)
ccc     $     + 1.5*x*Te_source
      s_u(9,2,:,:) = - u(9,:,:)*(n_inv + x*nx_inv)
     $     - nu_e*2.*Rinv*gfac*n_inv*uy(7,:,:)
      s_u(9,3,:,:) = - u(9,:,:)*ny_inv/gfac
     $     + nu_e*2.*Rinv*x*gfac**2*n_inv*(ux(7,:,:) - hfac*u(7,:,:))
      s_u(9,4,:,:) = eta_local*x/di
     $     *(3.92*(u(4,:,:) - u(1,:,:)*u(7,:,:))/di - 1.92*Jpar*be)
      s_u(9,6,:,:) = - 1.92*eta_local*x*Jpar*Jpar_u6
      s_u(9,7,:,:) = - nu_e*2.*x*hfac*(ux(7,:,:) - hfac*u(7,:,:) 
     $     + Rinv*gfac**2*n_inv*(u(3,:,:) + di*gfac*ux(6,:,:)))
     $     - eta_local*x*u(1,:,:)/di
     $     *(3.92*(u(4,:,:) - u(1,:,:)*u(7,:,:))/di - 1.92*Jpar*be)
      s_u(9,8,:,:) = x*heatfl
      s_u(9,9,:,:) = - (Vr + x*Vrx + Vtauy/gfac 
     $     - di*(uy(6,:,:)*nx_inv - ux(6,:,:)*ny_inv)) 
     $     - x*(heatfl - heatfl_u9*(u(8,:,:) - u(9,:,:)))
     $     + nu_e_u9*(x*(ux(7,:,:) - hfac*u(7,:,:))
     $     *(ux(7,:,:) - hfac*u(7,:,:) + 2.*Rinv*gfac**2
     $     *n_inv*(u(3,:,:) + di*gfac*ux(6,:,:)))
     $     + uy(7,:,:)*gfac*mfac*(uy(7,:,:) - 2.*Rinv*gfac**3
     $     *n_inv*(x*u(2,:,:) - di*uy(6,:,:))))
     $     + eta_local_u9*x
     $     *(1.96*((uy(6,:,:)/x)**2 + (gfac*ux(6,:,:))**2
     $     + (u(4,:,:) - u(1,:,:)*u(7,:,:))**2/di**2) - .96*Jpar**2)
      s_ux(9,1,:,:) = - u(9,:,:)*(x*Vr_u1 + di*uy(6,:,:)*n_inv**2)
      s_ux(9,2,:,:) = - u(9,:,:)*x*n_inv
      s_ux(9,5,:,:) = - 1.92*eta_local*x*Jpar*Jpar_ux5
      s_ux(9,6,:,:) = - di*u(9,:,:)*ny_inv
     $     + nu_e*di*2.*Rinv*x*gfac**3*n_inv*(ux(7,:,:) - hfac*u(7,:,:))
     $     + eta_local*x*gfac*(3.92*gfac*ux(6,:,:) + 1.92*Jpar*btau)
      s_ux(9,7,:,:) = nu_e*2.*x*(ux(7,:,:) - hfac*u(7,:,:) 
     $     + Rinv*gfac**2*n_inv*(u(3,:,:) + di*gfac*ux(6,:,:)))
      s_uy(9,1,:,:) = - u(9,:,:)*(Vtau_u1/gfac - di*ux(6,:,:)*n_inv**2)
      s_uy(9,3,:,:) = - u(9,:,:)*n_inv/gfac
      s_uy(9,5,:,:) = - 1.92*eta_local*x*Jpar*Jpar_uy5
      s_uy(9,6,:,:) = di*u(9,:,:)*nx_inv
     $     + nu_e*di*2.*Rinv*gfac/x*n_inv*uy(7,:,:)
     $     + eta_local*(3.92*uy(6,:,:)/x - 1.92*Jpar*br)
      s_uy(9,7,:,:) = nu_e*2.*gfac*mfac*(uy(7,:,:) 
     $     - Rinv*gfac**3*n_inv*(x*u(2,:,:) - di*uy(6,:,:)))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_drdu
c-----------------------------------------------------------------------
c     subprogram 9. XMHDhlx_mass.
c     computes mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_mass(x,y,mass)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      gfac=1/SQRT(1+(x*Rinv)**2)

      mass(1,1,:,:)=x
      mass(2,2,:,:)=x**2*gfac
      mass(3,3,:,:)=x**2*gfac
      mass(4,4,:,:)=x/gfac
      mass(5,5,:,:)=x
      mass(6,6,:,:)=x*gfac**2
      mass(6,5,:,:)=-2*x*Rinv*gfac**4
      mass(8,8,:,:)=1.5*x
      mass(9,9,:,:)=1.5*x
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_mass
c-----------------------------------------------------------------------
c     subprogram 10. XMHDhlx_equil.
c     computes equilibrium for Helical RMHD model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_equil(x,y,u,ux,uy,derivs)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs

      REAL(r8) :: gfac_s,psi0,fac1,fac2
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      SELECT CASE(init_type)
      CASE("m1_axial")
         IF(r_s <= 0)CALL program_stop("r_s has to be positive!")
c-----------------------------------------------------------------------
c     additional constant definitions.
c-----------------------------------------------------------------------
         gfac=1./(1. + (x*Rinv)**2)
         gfac_s=(1. + (Rinv*r_s)**2/3.)
         fac1 = 1. - (Rinv*r_s)**2
         fac2 = 2. + r_s**2
         psi0=half/gfac_s
c-----------------------------------------------------------------------
c     ion pressure.
c-----------------------------------------------------------------------
         u(8,:,:) = (1. - beta_e)*(beta0 + .5*x**2*(-.3*Rinv**4*x**8 
     $        - 5./8.*Rinv**2*fac1*x**6
     $        - 1./3.*(fac1**2 - 2.*Rinv**2*fac2)*x**4
     $        + .75*fac1*fac2*x**2 - .5*fac2**2))
         IF(MINVAL(u(8,:,:)) < min_eps)
     $        CALL program_stop("Negative initial pressure!")
c-----------------------------------------------------------------------
c     electron pressure.
c-----------------------------------------------------------------------
         u(9,:,:) = u(8,:,:)/(1. - beta_e)*beta_e
c-----------------------------------------------------------------------
c     density.
c-----------------------------------------------------------------------
         u(1,:,:) = u(9,:,:)/(beta0*beta_e + Tconst*(r_s**2 - x**2))
c-----------------------------------------------------------------------
c     magnetic flux.
c-----------------------------------------------------------------------
         u(5,:,:)=(x**4*Rinv**2/3. + x**2*(1. - (r_s*Rinv)**2)/2. 
     $        - r_s**2)*x**2/4.
c-----------------------------------------------------------------------
c     Be.
c-----------------------------------------------------------------------
         u(6,:,:)=1/(gfac*Rinv) - Rinv*x**2/(2*gfac)*(x**2-r_s**2)
c-----------------------------------------------------------------------
c     electron e-velocity.
c-----------------------------------------------------------------------
         u(7,:,:)=-di*SQRT(gfac)/u(1,:,:)
     $        *((1. + 2.*(Rinv*x)**2)*(x**2 - r_s**2) + x**2/gfac - 2)
      CASE("sawtooth1","sawtooth2")
         IF(r_s <= 0)CALL program_stop("r_s has to be positive!")
c-----------------------------------------------------------------------
c     ion pressure.
c-----------------------------------------------------------------------
         gfac=1./SQRT(1. + (x*Rinv)**2)
         u(8,:,:) = (1. - beta_e)*(beta0
     $        - (2.*(Rinv*x*gfac)**2 + (Rinv*x)**4 + 32.*Rinv**2
     $        - 16.*Rinv**2*(gfac + 1./gfac) 
     $        + 4.*(1. - 8.*Rinv**4)*LOG(gfac))/(32.*Rinv**6))
         IF(MINVAL(u(8,:,:)) < min_eps)
     $        CALL program_stop("Negative initial pressure!")
c-----------------------------------------------------------------------
c     electron pressure.
c-----------------------------------------------------------------------
         u(9,:,:) = u(8,:,:)/(1. - beta_e)*beta_e
c-----------------------------------------------------------------------
c     density.
c-----------------------------------------------------------------------
         u(1,:,:) = n0 + (one - n0)*EXP(-x**4/r_s**4)
c-----------------------------------------------------------------------
c     ion r-momentum.
c-----------------------------------------------------------------------
ccc         u(2,:,:) = di*u(1,:,:)**2.5/(1.96*nu*u(9,:,:)**1.5)
ccc     $        *x*Rinv**2*gfac**6/(1. + .5*Rinv**2*x**4*gfac**3)
ccc     $        *(2.*(1. + x**2*(2. + (x*Rinv)**2)*(Rinv**2 + 1./gfac))
ccc     $        + gfac**2*(2./gfac**3 - x**2)
ccc     $        *(4./gfac**3 - 4.*x**2 + Rinv**2*x**4)
ccc     $        /(4. + Rinv**2*x**6*gfac**6))
c-----------------------------------------------------------------------
c     ion tau-momentum.
c-----------------------------------------------------------------------
         u(3,:,:) = -di*(1. - beta_e)*Rinv*x
     $        *0.5*gfac**3*(8.-x**2*(4./gfac-2.*x**2-Rinv**2*(8.+x**4)))
     $        /(4./gfac - x**2)
c-----------------------------------------------------------------------
c     ion e-momentum.
c-----------------------------------------------------------------------
         u(4,:,:) = -di*(1. - beta_e)
     $        *0.5*gfac**3*(8.-x**2*(4./gfac-2.*x**2-Rinv**2*(8.+x**4)))
     $        /(4./gfac - x**2)
c-----------------------------------------------------------------------
c     magnetic flux.
c-----------------------------------------------------------------------
         u(5,:,:) = x**4/16._r8
c-----------------------------------------------------------------------
c     Be.
c-----------------------------------------------------------------------
         u(6,:,:) = 1/(Rinv*gfac)
c-----------------------------------------------------------------------
c     electron e-velocity.
c-----------------------------------------------------------------------
         u(7,:,:) = (di*gfac**3*(2./gfac - x**2*(1. + 1/gfac**2)/2._r8)
     $        + u(4,:,:))/u(1,:,:)
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_equil
c-----------------------------------------------------------------------
c     subprogram 11. XMHDhlx_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE XMHDhlx_grid(x,y,ksi,eta)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta

      REAL(r8) :: rs
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      rs=1/(1+((1-r_s)/r_s)**third)
      rs=rs+1.5*gr_curve*(r_s-rs)
      ksi=((x-rs)**3+5*gr_curve*x+rs**3)/(1+5*gr_curve+3*rs*(rs-1))
      eta=twopi*y
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE XMHDhlx_grid
      END MODULE XMHDhlx_mod
