c-----------------------------------------------------------------------
c     file frc.f.
c     contains specifications for compressible MHD model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. frc_mod.
c     1. frc_equil.
c     2. frc_wall.
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
c     subprogram 0. frc_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE frc_mod
      USE extra_mod
      USE transport_mod
      IMPLICIT NONE

      LOGICAL :: source=.FALSE.,cylinder=.FALSE.,
     $     if_kinvisc=.FALSE.,flux_inflow=.FALSE.

      CHARACTER(16) :: equilfile=".",
     $     interp="bicube",eta_case=".",kappa_case="."
      INTEGER :: cyl_fac=0
      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,qe=1.602e-19,
     $     me=9.109e-31,mi=3.345e-27,mp=1.673e-27,ep0=8.854e-12,
     $     mu0=4.e-7*pi,chod_const=0.1,kappa_min=0.,kappa_max=1.e8

      REAL(r8) :: eta=0.,r_eta=1.e10,etavac=1.,v_chod_norm=1.,
     $     etac_norm=1.,etas_norm=1.,mu=0.,kappa_par=0.,kappa_perp=0.,
     $     ddiff=0.,rhomin=0.,pmin=0.,pmax=1.,Tmin=1.,n0=1.e20,
     $     T0=1.,b0=1.,L0=1.,bx=0.,lx=0.,zmin=-1.,gamma_fac=1.,
     $     ke_norm=1.,xe_norm=1.,ki_norm=1.,xi_norm=1.,gr_curve=1.,
     $     awall,bwall,omwall,rphase=0

      TYPE(bicube_type) :: equil_bc
      REAL(r8), DIMENSION(:,:,:), POINTER :: equil,equilxy

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. frc_equil.
c     computes equilibrium for magnetic reconnection.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE frc_equil(x,y,u,ux,uy,deriv)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: deriv

      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: eq
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
c-----------------------------------------------------------------------
c     compute equilibrium,eq(1) = flux, eq(2) = pressure, and eq(3) =
c     jphi
c-----------------------------------------------------------------------
      CALL extra_interp(interp,x,y,equil_bc,equil,equilxy,eq)
c-----------------------------------------------------------------------
c     apply cosine smoothing to floor pressure (pmin).
c-----------------------------------------------------------------------
      WHERE(eq(2,:,:) < pmin)
         eq(2,:,:)=pmin
      ELSEWHERE(eq(2,:,:) >= pmin .AND. eq(2,:,:) < two*pmin)
         eq(2,:,:)=pmin 
     $        + half*pmin*(one - COS(pi*(eq(2,:,:)-pmin)/pmin))
      END WHERE
c-----------------------------------------------------------------------
c     apply cosine smoothing to floor density (rhomin).
c-----------------------------------------------------------------------
      u(1,:,:) = SQRT(eq(2,:,:))*b0**2/(two*mu0*n0*qe*T0)
      WHERE(u(1,:,:) < rhomin)
         u(1,:,:) = rhomin
      ELSEWHERE(u(1,:,:) >= rhomin .AND. u(1,:,:) < two*rhomin)
         u(1,:,:) = rhomin
     $        + half*rhomin*(one - COS(pi*(u(1,:,:)-rhomin)/rhomin))
      END WHERE

      u(2,:,:) = -eq(1,:,:)/y
      WHERE(y == 0)u(2,:,:) = zero
      u(3,:,:) = eq(2,:,:)
      u(6,:,:) = eq(3,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE frc_equil
c-----------------------------------------------------------------------
c     subprogram 2. frc_wall.
c     computes moving wall quantities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE frc_wall(t,rwall,vwall)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), INTENT(OUT) :: rwall,vwall
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      rwall=awall*COS(omwall*t+rphase)+bwall
      vwall=-omwall*awall*SIN(omwall*t+rphase)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE frc_wall
      END MODULE frc_mod
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
      USE frc_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: physics_type
      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty,nqty_schur
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep,exit_flag
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: myios
      REAL(r8) :: rinit=1,rstag=1,tstag=1
      REAL(r8) :: rwall,vwall
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/frc_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,source,cylinder,eta,eta_case,
     $     r_eta,etavac,mu,if_kinvisc,kappa_par,kappa_perp,
     $     ddiff,rhomin,pmin,pmax,b0,n0,T0,L0,bx,lx,
     $     zmin,equilfile,interp,flux_inflow,
     $     kappa_case,gr_curve,rinit,rstag,tstag,rphase
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1p,2(a,e10.3))
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=frc_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      nqty=6
      nqty_schur=0
      physics_type="frc"
      cylinder=.TRUE.
c-----------------------------------------------------------------------
c     compute moving wall parameters.
c-----------------------------------------------------------------------
      awall=(rinit-rstag)/2
      bwall=(rinit+rstag)/2
      omwall=pi/tstag
      rphase=rphase*pi
c-----------------------------------------------------------------------
c     diagnose wall parameters.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL frc_wall(zero,rwall,vwall)
         WRITE(*,10)"awall =  ",awall,", bwall =  ",bwall
         WRITE(*,10)"omwall = ",omwall,", rphase = ",rphase
         WRITE(*,10)"rwall =  ",rwall,", vwall =  ",vwall
      ENDIF
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
      USE frc_mod
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
      CALL MPI_Bcast(equilfile,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(interp,16,MPI_CHARACTER,0,comm,ierr)
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
      CALL MPI_Bcast(kappa_par,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_perp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ddiff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rhomin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(T0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(L0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(zmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(awall,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bwall,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(omwall,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rphase,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      IF(cylinder)cyl_fac=1
c-----------------------------------------------------------------------
c     compute normalizations in MKS units.
c     v_chod_norm normalizes the ve/vs term in Chodura resistivity.
c     etac_norm absorbs the constants in front of Chodura resistivity.
c     etas_norm absorbs the constants in front of Spitzer resistivity.
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

      CALL extra_read_marklin(equilfile,interp,pmax,
     $     equil_bc,equil,equilxy)

      Tmin=pmin/rhomin
      gamma_fac=gamma/(gamma-one)
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
      USE frc_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ux,uy
c-----------------------------------------------------------------------
c     frc initial conditions.
c-----------------------------------------------------------------------
      u=zero
      CALL frc_equil(x,y,u,ux,uy,.FALSE.)
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
      USE frc_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     set boundary condition type.
c-----------------------------------------------------------------------
      top%bc_type="robin"
      top%bc_type(6)="natural"
      bottom%bc_type="robin"
      left%bc_type="robin"
      right%bc_type="robin"
c-----------------------------------------------------------------------
c     set static.
c-----------------------------------------------------------------------
      top%static=.TRUE.
      top%static(2)=.FALSE.
      bottom%static=.TRUE.
      left%static=.TRUE.
      right%static=.TRUE.
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
      USE frc_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      REAL(r8) :: rwall,vwall
c-----------------------------------------------------------------------
c     zero c.
c-----------------------------------------------------------------------
      c=0
c-----------------------------------------------------------------------
c     wall boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(lrtb)
      CASE("top")
         CALL frc_wall(t,rwall,vwall)
         c(1,:,:)=uy(1,:,:)
         c(3,:,:)=u(1,:,:)*uy(3,:,:)-uy(1,:,:)*u(3,:,:)
         c(4,:,:)=u(4,:,:)
         c(5,:,:)=u(5,:,:)-u(1,:,:)*vwall
c-----------------------------------------------------------------------
c     axis boundary conditions.
c-----------------------------------------------------------------------
      CASE("bottom")
         c(1,:,:)=uy(1,:,:)
         c(2,:,:)=u(2,:,:)
         c(3,:,:)=uy(3,:,:)
         c(4,:,:)=uy(4,:,:)
         c(5,:,:)=u(5,:,:)
         c(6,:,:)=u(6,:,:)
c-----------------------------------------------------------------------
c     end boundary conditions.
c-----------------------------------------------------------------------
      CASE("left","right")
         c(1,:,:)=u(1,:,:)-rhomin
         c(3,:,:)=u(3,:,:)-Tmin*u(1,:,:)
         c(4,:,:)=u(4,:,:)
         c(5,:,:)=u(5,:,:)
         c(6,:,:)=u(6,:,:)
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
      USE frc_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy

      REAL(r8) :: rwall,vwall
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
      c_uxy=0
c-----------------------------------------------------------------------
c     wall boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(lrtb)
      CASE("top")
         CALL frc_wall(t,rwall,vwall)
         c_uy(1,1,:,:)=one
         c_u(3,1,:,:)=uy(3,:,:)
         c_u(3,3,:,:)=-uy(1,:,:)
         c_uy(3,1,:,:)=-u(3,:,:)
         c_uy(3,3,:,:)=u(1,:,:)
         c_u(4,4,:,:)=one
         c_u(5,1,:,:)=-vwall
         c_u(5,5,:,:)=one
c-----------------------------------------------------------------------
c     axis boundary conditions.
c-----------------------------------------------------------------------
      CASE("bottom")
         c_uy(1,1,:,:)=one
         c_u(2,2,:,:)=one
         c_uy(3,3,:,:)=one
         c_uy(4,4,:,:)=one
         c_u(5,5,:,:)=one
         c_u(6,6,:,:)=one
c-----------------------------------------------------------------------
c     end boundary conditions.
c-----------------------------------------------------------------------
      CASE("left","right")
         c_u(1,1,:,:)=one
         c_u(3,1,:,:)=-Tmin
         c_u(3,3,:,:)=one
         c_u(4,4,:,:)=one
         c_u(5,5,:,:)=one
         c_u(6,6,:,:)=one
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
      USE frc_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y

      INTEGER :: iqty
c-----------------------------------------------------------------------
c     zero mass matrices.
c-----------------------------------------------------------------------
      mass=0
      mass_x=0
      mass_y=0
c-----------------------------------------------------------------------
c     set mass matrix for top -A_phi variable.
c-----------------------------------------------------------------------
      IF(lrtb == "top")mass(2,2,:,:)=1
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
      USE frc_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i,iqty
      REAL(r8) :: rwall,vwall
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Tix,Tiy,kperp,kfac,
     $     eta_local,r_fac,r_faci,visc,nil,b1,b2,Bsq,n_inv,nx_inv,
     $     ny_inv,massfac,massfacy
      REAL(r8), DIMENSION(2,SIZE(u,2),SIZE(u,3)) :: vi,vix,viy
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: BdotT
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     compute moving wall parameters.
c-----------------------------------------------------------------------
      CALL frc_wall(t,rwall,vwall)
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
      ny_inv=-uy(1,:,:)*n_inv**2/rwall
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
      CALL transport_seteta(eta_case,one-y,u(1,:,:),u(3,:,:),u(6,:,:),
     $     chod_const,etas_norm,etac_norm,v_chod_norm,one-r_eta,eta,
     $     etavac,eta_local)
c-----------------------------------------------------------------------
c     temperature gradients.
c-----------------------------------------------------------------------
      Tix = ux(3,:,:)*n_inv + u(3,:,:)*nx_inv
      Tiy = uy(3,:,:)*n_inv/rwall + u(3,:,:)*ny_inv
c-----------------------------------------------------------------------
c     magnetic fields.
c-----------------------------------------------------------------------
      nil=0
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))/rwall
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
         viy(i,:,:)=uy(i+3,:,:)*n_inv/rwall + u(i+3,:,:)*ny_inv
      ENDDO
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
      fx(1,:,:)=r_fac*(u(4,:,:) - ddiff*ux(1,:,:))
      fy(1,:,:)=r_fac*(u(5,:,:) - ddiff*uy(1,:,:)/rwall)
c-----------------------------------------------------------------------
c     equation for -Aphi.
c-----------------------------------------------------------------------
      s(2,:,:) = vi(2,:,:)*(b1+bx/rwall**2)
     $     - vi(1,:,:)*b2 + eta_local*u(6,:,:)
c-----------------------------------------------------------------------
c     pressure equation.
c-----------------------------------------------------------------------
      fx(3,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi(1,:,:)
     $     - kfac*BdotT(1,:,:) - kperp*Tix)
      fy(3,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi(2,:,:)
     $     - kfac*BdotT(2,:,:) - kperp*Tiy)
      s(3,:,:)=r_fac*(vi(1,:,:)*ux(3,:,:) + vi(2,:,:)*uy(3,:,:)/rwall
     $     + eta_local*u(6,:,:)**2
     $     + visc*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 + viy(1,:,:)**2
     $     + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:))
     $     + cyl_fac*two*visc*(vi(2,:,:)*r_faci/rwall)**2)
c-----------------------------------------------------------------------
c     axial momentum equation.
c-----------------------------------------------------------------------
      fx(4,:,:)=r_fac*(u(4,:,:)*vi(1,:,:)
     $     - two*visc*vix(1,:,:) + u(3,:,:))

      fy(4,:,:)=r_fac*(u(4,:,:)*vi(2,:,:) 
     $     - visc*(viy(1,:,:) + vix(2,:,:)))
      
      s(4,:,:)=-r_fac*u(6,:,:)*b2
c-----------------------------------------------------------------------
c     radial momentum equation.
c-----------------------------------------------------------------------
      fx(5,:,:)=r_fac*(u(5,:,:)*vi(1,:,:)
     $     - visc*(vix(2,:,:) + viy(1,:,:)))
      
      fy(5,:,:)=r_fac*(u(5,:,:)*vi(2,:,:)
     $     - two*visc*viy(2,:,:))

      s(5,:,:)=r_fac*(u(6,:,:)*(b1+bx/rwall**2) - uy(3,:,:)/rwall
     $     - cyl_fac*two*visc*vi(2,:,:)*(r_faci/rwall)**2)
c-----------------------------------------------------------------------
c     current equation.
c-----------------------------------------------------------------------
      fx(6,:,:)=b2
      fy(6,:,:)=-b1
      s(6,:,:)=u(6,:,:)
c-----------------------------------------------------------------------
c     moving wall terms.
c-----------------------------------------------------------------------
      DO iqty=1,5
         SELECT CASE(iqty)
         CASE(1,4,5)
            massfac=y
            massfacy=1
         CASE(2)
            massfac=1
            massfacy=0
         CASE(3)
            massfac=y/(gamma-1)
            massfacy=1/(gamma-1)
         END SELECT
         fy(iqty,:,:)=fy(iqty,:,:)-massfac*u(iqty,:,:)*y*vwall
         s(iqty,:,:)=s(iqty,:,:)
     $        -(2*massfac+y*massfacy)*u(iqty,:,:)*vwall/rwall
      ENDDO
      fy=fy/rwall
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
      USE frc_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i,iqty
      REAL(r8) :: rwall,vwall
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: r_fac,r_faci,Tix,Tiy,
     $     Tix_un,Tiy_un,Ti_un,kperp,kfac,kpar_u1,kperp_bsq,kperp_u1,
     $     kperp_u3,kpar_u3,eta_local,visc,visc_rho,eta_u1,eta_u3,
     $     eta_u6,nil,b1,b2,Bsq,n_inv,nx_inv,ny_inv,massfac,massfacy
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
c     compute moving wall parameters.
c-----------------------------------------------------------------------
      CALL frc_wall(t,rwall,vwall)
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
      ny_inv=-uy(1,:,:)*n_inv**2/rwall
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
      CALL transport_seteta_u(eta_case,one-y,u(1,:,:),u(3,:,:),u(6,:,:),
     $     chod_const,etas_norm,etac_norm,v_chod_norm,one-r_eta,eta,
     $     etavac,eta_local,eta_u1,eta_u3,eta_u6) 
c-----------------------------------------------------------------------
c     temperature gradients and derivatives.
c-----------------------------------------------------------------------
      Tix = ux(3,:,:)*n_inv + u(3,:,:)*nx_inv
      Tiy = uy(3,:,:)*n_inv/rwall + u(3,:,:)*ny_inv
      Ti_un=-u(3,:,:)*n_inv**2
      Tix_un=-(two*u(3,:,:)*nx_inv + ux(3,:,:)*n_inv)*n_inv
      Tiy_un=-(two*u(3,:,:)*ny_inv + uy(3,:,:)*n_inv/rwall)*n_inv
c-----------------------------------------------------------------------
c     magnetic fields.
c-----------------------------------------------------------------------
      nil=0
      nil3=0
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))/rwall
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
         viy(i,:,:)=uy(i+3,:,:)*n_inv/rwall + u(i+3,:,:)*ny_inv

         vi_un(i,:,:)=-vi(i,:,:)*n_inv
         vix_un(i,:,:) = -(vi(i,:,:)*nx_inv + vix(i,:,:)*n_inv)
         viy_un(i,:,:) = -(vi(i,:,:)*ny_inv + viy(i,:,:)*n_inv)
      ENDDO
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
      fx_u(1,4,:,:)=r_fac
      fy_u(1,5,:,:)=r_fac
      fx_ux(1,1,:,:)=-r_fac*ddiff
      fy_uy(1,1,:,:)=-r_fac*ddiff/rwall
c-----------------------------------------------------------------------
c     equation for -Aphi.
c-----------------------------------------------------------------------
      s_u(2,1,:,:)=vi_un(2,:,:)*(b1+bx/rwall**2) - vi_un(1,:,:)*b2
     $     + eta_u1*u(6,:,:)
      s_u(2,2,:,:)=-vi(2,:,:)*cyl_fac*r_faci/rwall
      s_u(2,3,:,:)=eta_u3*u(6,:,:)
      s_u(2,4,:,:)=-n_inv*b2
      s_u(2,5,:,:)=n_inv*(b1+bx/rwall**2)
      s_u(2,6,:,:)=eta_local + u(6,:,:)*eta_u6
      s_ux(2,2,:,:)=-vi(1,:,:)
      s_uy(2,2,:,:)=-vi(2,:,:)/rwall
c-----------------------------------------------------------------------
c     pressure equation.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi_un(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*Tix_un 
     $     - kfac*BdotT_Ty(1,:,:)*Tiy_un
     $     - (kpar_u1 - kperp_u1)*BdotT(1,:,:) - kperp_u1*Tix)
      fx_u(3,2,:,:)=cyl_fac*(kfac*BdotT_b1(1,:,:)
     $     - two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))/rwall
      fx_u(3,3,:,:)=r_fac*(gamma_fac*vi(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*nx_inv
     $     - kfac*BdotT_Ty(1,:,:)*ny_inv
     $     - (kpar_u3 - kperp_u3)*BdotT(1,:,:) - kperp_u3*Tix)
      fx_u(3,4,:,:)=r_fac*gamma_fac*u(3,:,:)*n_inv

      fx_ux(3,1,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*Ti_un
      fx_ux(3,2,:,:)=r_fac*(-kfac*BdotT_b2(1,:,:)
     $     + two*b2*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_ux(3,3,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*n_inv

      fx_uy(3,1,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*Ti_un/rwall
      fx_uy(3,2,:,:)=r_fac*(kfac*BdotT_b1(1,:,:)
     $     - two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))/rwall
      fx_uy(3,3,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*n_inv/rwall

      fy_u(3,1,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi_un(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*Tiy_un
     $     - kfac*BdotT_Tx(2,:,:)*Tix_un
     $     - (kpar_u1 - kperp_u1)*BdotT(2,:,:) - kperp_u1*Tiy)
      fy_u(3,2,:,:)=cyl_fac*(kfac*BdotT_b1(2,:,:)
     $     - two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))/rwall
      fy_u(3,3,:,:)=r_fac*(gamma_fac*vi(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*ny_inv
     $     - kfac*BdotT_Tx(2,:,:)*nx_inv
     $     - (kpar_u3 - kperp_u3)*BdotT(2,:,:) - kperp_u3*Tiy)
      fy_u(3,5,:,:)=r_fac*gamma_fac*u(3,:,:)*n_inv

      fy_ux(3,1,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*Ti_un
      fy_ux(3,2,:,:)=r_fac*(-kfac*BdotT_b2(2,:,:)
     $     + two*b2*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_ux(3,3,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*n_inv

      fy_uy(3,1,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*Ti_un/rwall
      fy_uy(3,2,:,:)=r_fac*(kfac*BdotT_b1(2,:,:)
     $     - two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))/rwall
      fy_uy(3,3,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*n_inv/rwall

      s_u(3,1,:,:)=r_fac*(vi_un(1,:,:)*ux(3,:,:) 
     $     + vi_un(2,:,:)*uy(3,:,:)/rwall
     $     + eta_u1*u(6,:,:)**2
     $     + two*visc*(two*vix(1,:,:)*vix_un(1,:,:) 
     $     + two*viy(2,:,:)*viy_un(2,:,:) 
     $     + viy(1,:,:)*viy_un(1,:,:) + vix(2,:,:)*vix_un(2,:,:)
     $     + viy(1,:,:)*vix_un(2,:,:) + vix(2,:,:)*viy_un(1,:,:))
     $     + visc_rho*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:))
     $     + cyl_fac*two*(r_faci/rwall)**2
     $     *(visc*two*vi(2,:,:)*vi_un(2,:,:)
     $     + visc_rho*vi(2,:,:)**2))

      s_u(3,3,:,:)=r_fac*eta_u3*u(6,:,:)**2
      s_u(3,4,:,:)=r_fac*(ux(3,:,:)*n_inv + two*visc*
     $     (two*vix(1,:,:)*nx_inv + (viy(1,:,:) + vix(2,:,:))*ny_inv))
      s_u(3,5,:,:)=r_fac*(uy(3,:,:)*n_inv/rwall + two*visc*
     $     (two*viy(2,:,:)*ny_inv + (vix(2,:,:) + viy(1,:,:))*nx_inv)
     $     + cyl_fac*4._r8*visc*(r_faci/rwall)**2*vi(2,:,:)*n_inv)
      s_u(3,6,:,:)=r_fac*(two*eta_local*u(6,:,:) + eta_u6*u(6,:,:)**2)

      s_ux(3,1,:,:)=r_fac*two*visc*(two*vix(1,:,:)*vi_un(1,:,:)  
     $     + vi_un(2,:,:)*(vix(2,:,:) + viy(1,:,:)))
      s_ux(3,3,:,:)=r_fac*vi(1,:,:)
      s_ux(3,4,:,:)=r_fac*visc*4._r8*vix(1,:,:)*n_inv
      s_ux(3,5,:,:)=r_fac*visc*two*n_inv*(vix(2,:,:) + viy(1,:,:))

      s_uy(3,1,:,:)=r_fac*two*visc*(two*viy(2,:,:)*vi_un(2,:,:)  
     $     + vi_un(1,:,:)*(viy(1,:,:) + vix(2,:,:)))/rwall
      s_uy(3,3,:,:)=r_fac*vi(2,:,:)/rwall
      s_uy(3,4,:,:)=r_fac*visc*two*n_inv*(viy(1,:,:) + vix(2,:,:))/rwall
      s_uy(3,5,:,:)=r_fac*visc*4._r8*viy(2,:,:)*n_inv/rwall
c-----------------------------------------------------------------------
c     axial momentum equation.
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
      fy_uy(4,1,:,:)=-r_fac*visc*vi_un(1,:,:)/rwall
      fy_uy(4,4,:,:)=-r_fac*visc*n_inv/rwall

      s_u(4,6,:,:)=-r_fac*b2
      s_ux(4,2,:,:)=-r_fac*u(6,:,:)
c-----------------------------------------------------------------------
c     radial momentum equation.
c-----------------------------------------------------------------------
      fx_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(1,:,:) 
     $     - visc*(vix_un(2,:,:) + viy_un(1,:,:))
     $     - visc_rho*(vix(2,:,:) + viy(1,:,:)))
      fx_u(5,4,:,:)=r_fac*(vi(2,:,:) - visc*ny_inv)
      fx_u(5,5,:,:)=r_fac*(vi(1,:,:) - visc*nx_inv)
      fx_ux(5,1,:,:)=-r_fac*visc*vi_un(2,:,:)
      fx_ux(5,5,:,:)=-r_fac*visc*n_inv
      fx_uy(5,1,:,:)=-r_fac*visc*vi_un(1,:,:)/rwall
      fx_uy(5,4,:,:)=-r_fac*visc*n_inv/rwall
      
      fy_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(2,:,:) 
     $     - two*visc*viy_un(2,:,:) - two*visc_rho*viy(2,:,:))
      fy_u(5,5,:,:)=r_fac*two*(vi(2,:,:) - visc*ny_inv)
      fy_uy(5,1,:,:)=-r_fac*visc*two*vi_un(2,:,:)/rwall
      fy_uy(5,5,:,:)=-r_fac*visc*two*n_inv/rwall

      s_u(5,1,:,:)=-cyl_fac*two*r_faci/rwall**2*(visc*vi_un(2,:,:)
     $     + visc_rho*vi(2,:,:))
      s_u(5,2,:,:)=-cyl_fac*u(6,:,:)/rwall
      s_u(5,5,:,:)=-cyl_fac*two*visc*r_faci/rwall**2*n_inv
      s_u(5,6,:,:)=r_fac*(b1+bx/rwall**2)
      s_uy(5,2,:,:)=-r_fac*u(6,:,:)/rwall
      s_uy(5,3,:,:)=-r_fac/rwall
c-----------------------------------------------------------------------
c     current equation.
c-----------------------------------------------------------------------
      fx_ux(6,2,:,:)=one
      fy_u(6,2,:,:)=cyl_fac*r_faci/rwall
      fy_uy(6,2,:,:)=one/rwall
      s_u(6,6,:,:)=one
c-----------------------------------------------------------------------
c     moving wall terms.
c-----------------------------------------------------------------------
      DO iqty=1,5
         SELECT CASE(iqty)
         CASE(1,4,5)
            massfac=y
            massfacy=1
         CASE(2)
            massfac=1
            massfacy=0
         CASE(3)
            massfac=y/(gamma-1)
            massfacy=1/(gamma-1)
         END SELECT
         fy_u(iqty,iqty,:,:)=fy_u(iqty,iqty,:,:)-massfac*y*vwall
         s_u(iqty,iqty,:,:)=s_u(iqty,iqty,:,:)
     $        -(2*massfac+y*massfacy)*vwall/rwall
      ENDDO
      fy_u=fy_u/rwall
      fy_ux=fy_ux/rwall
      fy_uy=fy_uy/rwall
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
      USE frc_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac
c-----------------------------------------------------------------------
c        cylindrical to cartesian relationships:
c        1: z --> x
c        2: r --> y
c        3: phi --> z
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     modify mass matrices.
c-----------------------------------------------------------------------
      r_fac=y**cyl_fac
      mass(1,1,:,:)=r_fac
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
      USE frc_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,etag

      REAL(r8) :: fac,a
c-----------------------------------------------------------------------
c     set grid according to init_type.
c-----------------------------------------------------------------------
      a=zmin/L0
      IF(gr_curve == 0.)THEN
         ksi=lx*x+a
         etag=y
      ELSE
         fac=(2.1_r8 - a)/lx
         ksi=((x - fac)**3 + gr_curve*x + fac**3)*lx/
     $        ((one - fac)**3 + gr_curve + fac**3) + a
         etag=(y**2 + gr_curve*y)/(one + gr_curve)
      ENDIF
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
      USE frc_mod
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
      USE frc_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     deallocate appropriate arrays.
c-----------------------------------------------------------------------
      CALL extra_equil_dealloc(interp,equil_bc,equil,equilxy)
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
