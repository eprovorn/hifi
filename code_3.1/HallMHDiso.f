c-----------------------------------------------------------------------
c     file HallMHDiso.f.
c     contains specifications for two-fluid isothermal model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. HallMHDiso_mod.
c     1. HallMHDiso_input.
c     2. HallMHDiso_init.
c     3. HallMHDiso_init_special.
c     4. HallMHDiso_boundary.
c     5. HallMHDiso_rhs.
c     6. HallMHDiso_drdu.
c     7. HallMHDiso_equil.
c     8. HallMHDiso_grid.
c-----------------------------------------------------------------------
c     subprogram 0. HallMHDiso_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE HallMHDiso_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: eta_spitzer=.FALSE.,source=.FALSE.
      CHARACTER(16), PRIVATE :: init_type=" "
      REAL(r8), PRIVATE :: di=0,eta=0,mu=0,nu=0,kappa_par=0,
     $     kappa_perp=0,Dn=0,lx=0,ly=0,mach=0,alfven=0,lambda_psi=0,
     $     lambda_phi=0,epsilon=0,Bguide=0,beta_e=0,beta0=0,
     $     kx=0,ky=0,ksq=0,grad_rho=0

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. HallMHDiso_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDiso_input(nx,ny,np,nq,nbx,xperiodic,
     $     yperiodic,nqty,dt,dtmax,tmax,nstep,gr_curve)

      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic
      INTEGER, INTENT(OUT) :: nqty
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax,gr_curve

      INTEGER :: myios
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/HallMHD_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,di,eta,mu,nu,kappa_par,kappa_perp,
     $     eta_spitzer,Dn,lx,ly,mach,alfven,lambda_psi,lambda_phi,
     $     epsilon,Bguide,beta_e,beta0,source,init_type
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=HallMHD_list,IOSTAT=myios)

      nqty=7
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDiso_input
c-----------------------------------------------------------------------
c     subprogram 2. HallMHDiso_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDiso_init(xpi,ypi,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      REAL(r8), DIMENSION(SIZE(ui,1),0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1)
     $     :: uy,ux
      REAL(r8), DIMENSION(0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) :: cosx,cosy,
     $     sinx,siny
c-----------------------------------------------------------------------
c     magnetic reconnection, initial conditions.
c-----------------------------------------------------------------------
      CALL HallMHDiso_equil(xpi,ypi,ui,ux,uy,.FALSE.)
      sinx=SIN(kx*xpi)
      siny=SIN(ky*ypi)
      cosx=COS(kx*xpi)
      cosy=COS(ky*ypi)
      SELECT CASE(init_type)
      CASE("Craig")
         ui(1,:,:)=ui(1,:,:)+epsilon*sinx*siny/kx
         ui(2,:,:)=ui(2,:,:)+epsilon*cosx/kx
         ui(7,:,:)=ui(7,:,:)+di*epsilon*cosx*kx
      CASE("Alan")
         ui(2,:,:)=ui(2,:,:)-epsilon*cosx*cosy*alfven/ksq
         ui(4,:,:)=ui(4,:,:)-epsilon*cosx*siny*mach*ky/ksq
         ui(5,:,:)=ui(5,:,:)+epsilon*sinx*cosy*mach*kx/ksq
         ui(7,:,:)=ui(7,:,:)-di*epsilon*cosx*cosy*alfven
c-----------------------------------------------------------------------
c     abort if not recognized.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT

      RETURN
      END SUBROUTINE HallMHDiso_init
c-----------------------------------------------------------------------
c     subprogram 3. HallMHDiso_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDiso_init_special(static)

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
      CALL MPI_Bcast(eta_spitzer,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_par,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_perp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Dn,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta_e,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_phi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Bguide,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alfven,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      grad_rho=0
      SELECT CASE(init_type)
      CASE("Craig")
         kx=pi
         ky=pi
      CASE("Alan")
         kx=twopi/lx
         ky=pi
      END SELECT
      ksq=kx**2+ky**2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDiso_init_special
c-----------------------------------------------------------------------
c     subprogram 4. HallMHDiso_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDiso_boundary(left,right,top,bottom,nqty)

      INTEGER, INTENT(IN) :: nqty
      TYPE(edge_type) :: left,right,top,bottom

      INTEGER :: k,j
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Craig")
         left%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         right%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         top%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         bottom%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         left%b=0
         right%b=0
         top%b=0
         bottom%b=0
      CASE("Alan")
         left%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         right%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         top%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         bottom%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         left%b=0
         right%b=0
         top%b=0
         bottom%b=0
         top%bc_type(1)="natural"
         bottom%bc_type(1)="natural"
         top%bc_type(7)="natural"
         bottom%bc_type(7)="natural"
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDiso_boundary
c-----------------------------------------------------------------------
c     subprogram 5. HallMHDiso_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE HallMHDiso_rhs(x,y,u,ux,uy,fx,fy,s,first)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8), DIMENSION(SIZE(u,1),SIZE(x,1),SIZE(x,2)) :: u0,u0y,u0x,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     two-fluid isothermal formulation.
c-----------------------------------------------------------------------
      fx(1,:,:)=u(4,:,:)
      
      fy(1,:,:)=u(5,:,:)
      
      fx(2,:,:)=-eta*ux(2,:,:)
         
      fy(2,:,:)=-eta*uy(2,:,:)
         
      s(2,:,:)=-(u(4,:,:)*ux(2,:,:)+u(5,:,:)*uy(2,:,:))/u(1,:,:)
     $     -di*(ux(3,:,:)*uy(2,:,:)-uy(3,:,:)*ux(2,:,:))/u(1,:,:)
     $     +di*beta0*grad_rho/u(1,:,:)

      fx(3,:,:)=u(3,:,:)*u(4,:,:)/u(1,:,:)+u(7,:,:)*uy(2,:,:)/u(1,:,:)
     $     -eta*ux(3,:,:)-two*di*beta0*grad_rho
     $     *u(3,:,:)/(u(1,:,:)*ux(2,:,:))
         
      fy(3,:,:)=u(3,:,:)*u(5,:,:)/u(1,:,:)-u(7,:,:)*ux(2,:,:)/u(1,:,:)
     $     -eta*uy(3,:,:)

      s(3,:,:)=di*u(3,:,:)*(ux(3,:,:)*uy(1,:,:)-uy(3,:,:)*ux(1,:,:))
     $     /u(1,:,:)**2
     $     -grad_rho*u(3,:,:)*u(7,:,:)/u(1,:,:)**2

      fx(4,:,:)=u(4,:,:)**2/u(1,:,:)+two*beta0*u(1,:,:)
     $     +half*(ux(2,:,:)**2-uy(2,:,:)**2+u(3,:,:)**2)
     $     -mu*(ux(4,:,:)/u(1,:,:)-u(4,:,:)*ux(1,:,:)/u(1,:,:)**2)
         
      fy(4,:,:)=u(4,:,:)*u(5,:,:)/u(1,:,:)+uy(2,:,:)*ux(2,:,:)
     $     -mu*(uy(4,:,:)/u(1,:,:)-u(4,:,:)*uy(1,:,:)/u(1,:,:)**2)
      
      s(4,:,:)=grad_rho*u(4,:,:)*u(6,:,:)/u(1,:,:)**2

      fx(5,:,:)=u(4,:,:)*u(5,:,:)/u(1,:,:)+uy(2,:,:)*ux(2,:,:)
     $     -mu*(ux(5,:,:)/u(1,:,:)-u(5,:,:)*ux(1,:,:)/u(1,:,:)**2)
      
      fy(5,:,:)=u(5,:,:)**2/u(1,:,:)+two*beta0*u(1,:,:)
     $     +half*(-ux(2,:,:)**2+uy(2,:,:)**2+u(3,:,:)**2)
     $     -mu*(uy(5,:,:)/u(1,:,:)-u(5,:,:)*uy(1,:,:)/u(1,:,:)**2)
      
      s(5,:,:)=grad_rho*(u(5,:,:)*u(6,:,:)/u(1,:,:)**2
     $     -two*beta0*u(3,:,:)/ux(2,:,:))

      fx(6,:,:)=u(4,:,:)*u(6,:,:)/u(1,:,:)+u(3,:,:)*uy(2,:,:)
     $     -mu*(ux(6,:,:)/u(1,:,:)-u(6,:,:)*ux(1,:,:)/u(1,:,:)**2)
      
      fy(6,:,:)=u(5,:,:)*u(6,:,:)/u(1,:,:)-u(3,:,:)*ux(2,:,:)
     $     -mu*(uy(6,:,:)/u(1,:,:)-u(6,:,:)*uy(1,:,:)/u(1,:,:)**2)
      
      s(6,:,:)=grad_rho*(u(6,:,:)/u(1,:,:))**2

      fx(7,:,:)=di*ux(2,:,:)
      fy(7,:,:)=di*uy(2,:,:)
      s(7,:,:)=u(6,:,:)-u(7,:,:)

      IF(source .AND. first)THEN
         CALL HallMHDiso_equil(x,y,u0,u0x,u0y,.TRUE.)
         first=.FALSE.
         CALL HallMHDiso_rhs(x,y,u0,u0x,u0y,fx0,fy0,s0,first)
         fx=fx-fx0
         fy=fy-fy0
         s=s-s0
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDiso_rhs
c-----------------------------------------------------------------------
c     subprogram 6. HallMHDiso_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDiso_drdu(x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

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
c     two-fluid isothermal reconnection.
c-----------------------------------------------------------------------
      fx_u(1,4,:,:)=1
      fy_u(1,5,:,:)=1

      fx_ux(2,2,:,:)=-eta
      fy_uy(2,2,:,:)=-eta

      s_u(2,1,:,:)=(u(4,:,:)*ux(2,:,:)+u(5,:,:)*uy(2,:,:)
     $     -di*(ux(3,:,:)*uy(2,:,:)-uy(3,:,:)*ux(2,:,:)))
     $     /u(1,:,:)**2
     $     -di*beta0*grad_rho/u(1,:,:)**2
      s_u(2,4,:,:)=-ux(2,:,:)/u(1,:,:)
      s_u(2,5,:,:)=-uy(2,:,:)/u(1,:,:)
      s_ux(2,2,:,:)=(di*uy(3,:,:)-u(4,:,:))/u(1,:,:)
      s_ux(2,3,:,:)=-di*uy(2,:,:)/u(1,:,:)
      s_uy(2,2,:,:)=-(u(5,:,:)+di*ux(3,:,:))/u(1,:,:)
      s_uy(2,3,:,:)=di*ux(2,:,:)/u(1,:,:)

      fx_u(3,1,:,:)=-(u(3,:,:)*u(4,:,:)+u(7,:,:)*uy(2,:,:))
     $     /u(1,:,:)**2
     $     +two*di*beta0*grad_rho
     $     *u(3,:,:)/(ux(2,:,:)*u(1,:,:)**2)
      fx_u(3,3,:,:)=u(4,:,:)/u(1,:,:)
     $     -two*di*beta0*grad_rho/(u(1,:,:)*ux(2,:,:))

      fx_u(3,4,:,:)=u(3,:,:)/u(1,:,:)
      fx_u(3,7,:,:)=uy(2,:,:)/u(1,:,:)
      fx_ux(3,2,:,:)=two*di*beta0*grad_rho
     $     *u(3,:,:)/(u(1,:,:)*ux(2,:,:)**2)
      fx_ux(3,3,:,:)=-eta
      fx_uy(3,2,:,:)=u(7,:,:)/u(1,:,:)
         
      fy_u(3,1,:,:)=-(u(3,:,:)*u(5,:,:)-u(7,:,:)*ux(2,:,:))/u(1,:,:)**2
      fy_u(3,3,:,:)=u(5,:,:)/u(1,:,:)
      fy_u(3,5,:,:)=u(3,:,:)/u(1,:,:)
      fy_u(3,7,:,:)=-ux(2,:,:)/u(1,:,:)
      fy_ux(3,2,:,:)=-u(7,:,:)/u(1,:,:)
      fy_uy(3,3,:,:)=-eta
         
      s_u(3,1,:,:)=-two*di*u(3,:,:)*
     $     (ux(3,:,:)*uy(1,:,:)-uy(3,:,:)*ux(1,:,:))/u(1,:,:)**3
     $     +two*grad_rho*u(3,:,:)*u(7,:,:)/u(1,:,:)**3
      s_u(3,3,:,:)=di*(ux(3,:,:)*uy(1,:,:)-uy(3,:,:)*ux(1,:,:))
     $     /u(1,:,:)**2
     $     -grad_rho*u(7,:,:)/u(1,:,:)**2
      s_u(3,7,:,:)=-grad_rho*u(3,:,:)/u(1,:,:)**2
      s_ux(3,1,:,:)=-di*u(3,:,:)*uy(3,:,:)/u(1,:,:)**2
      s_ux(3,3,:,:)=di*u(3,:,:)*uy(1,:,:)/u(1,:,:)**2
      s_uy(3,1,:,:)=di*u(3,:,:)*ux(3,:,:)/u(1,:,:)**2
      s_uy(3,3,:,:)=-di*u(3,:,:)*ux(1,:,:)/u(1,:,:)**2

      fx_u(4,1,:,:)=two*beta0 + (-u(4,:,:)**2+mu*(ux(4,:,:) 
     $     -two*u(4,:,:)*ux(1,:,:)/u(1,:,:)))/u(1,:,:)**2
      fx_u(4,3,:,:)=u(3,:,:)
      fx_u(4,4,:,:)=two*u(4,:,:)/u(1,:,:)+mu*ux(1,:,:)/u(1,:,:)**2
      fx_ux(4,1,:,:)=mu*u(4,:,:)/u(1,:,:)**2
      fx_ux(4,2,:,:)=ux(2,:,:)
      fx_ux(4,4,:,:)=-mu/u(1,:,:)
      fx_uy(4,2,:,:)=-uy(2,:,:)

      fy_u(4,1,:,:)=(-u(4,:,:)*u(5,:,:)+mu*(uy(4,:,:) 
     $     -two*u(4,:,:)*uy(1,:,:)/u(1,:,:)))/u(1,:,:)**2
      fy_u(4,4,:,:)=u(5,:,:)/u(1,:,:)+mu*uy(1,:,:)/u(1,:,:)**2
      fy_u(4,5,:,:)=u(4,:,:)/u(1,:,:)
      fy_ux(4,2,:,:)=uy(2,:,:)
      fy_uy(4,1,:,:)=mu*u(4,:,:)/u(1,:,:)**2
      fy_uy(4,2,:,:)=ux(2,:,:)
      fy_uy(4,4,:,:)=-mu/u(1,:,:)

      s_u(4,1,:,:)=-two*grad_rho*u(4,:,:)*u(6,:,:)
     $     /u(1,:,:)**3
      s_u(4,4,:,:)=grad_rho*u(6,:,:)/u(1,:,:)**2
      s_u(4,6,:,:)=grad_rho*u(4,:,:)/u(1,:,:)**2

      fx_u(5,1,:,:)=(-u(4,:,:)*u(5,:,:)+mu*(ux(5,:,:) 
     $     -two*u(5,:,:)*ux(1,:,:)/u(1,:,:)))/u(1,:,:)**2
      fx_u(5,4,:,:)=u(5,:,:)/u(1,:,:)
      fx_u(5,5,:,:)=u(4,:,:)/u(1,:,:)+mu*ux(1,:,:)/u(1,:,:)**2
      fx_ux(5,1,:,:)=mu*u(5,:,:)/u(1,:,:)**2
      fx_ux(5,2,:,:)=uy(2,:,:)
      fx_ux(5,5,:,:)=-mu/u(1,:,:)
      fx_uy(5,2,:,:)=ux(2,:,:)
      
      fy_u(5,1,:,:)=two*beta0 + (-u(5,:,:)**2+mu*(uy(5,:,:) 
     $     -two*u(5,:,:)*uy(1,:,:)/u(1,:,:)))/u(1,:,:)**2
      fy_u(5,3,:,:)=u(3,:,:)
      fy_u(5,5,:,:)=two*u(5,:,:)/u(1,:,:)+mu*uy(1,:,:)/u(1,:,:)**2
      fy_ux(5,2,:,:)=-ux(2,:,:)
      fy_uy(5,1,:,:)=mu*u(5,:,:)/u(1,:,:)**2
      fy_uy(5,2,:,:)=uy(2,:,:)
      fy_uy(5,5,:,:)=-mu/u(1,:,:)
      
      s_u(5,1,:,:)=-two*grad_rho*u(5,:,:)*u(6,:,:)
     $     /u(1,:,:)**3
      s_u(5,5,:,:)=grad_rho*u(6,:,:)/u(1,:,:)**2
      s_u(5,3,:,:)=-two*beta0*grad_rho/ux(2,:,:)
      s_ux(5,2,:,:)=two*beta0*grad_rho*u(3,:,:)
     $     /ux(2,:,:)**2

      fx_u(6,1,:,:)=(-u(4,:,:)*u(6,:,:)+mu*(ux(6,:,:) 
     $     -two*u(6,:,:)*ux(1,:,:)/u(1,:,:)))/u(1,:,:)**2
      fx_u(6,3,:,:)=uy(2,:,:)
      fx_u(6,4,:,:)=u(6,:,:)/u(1,:,:)
      fx_u(6,6,:,:)=u(4,:,:)/u(1,:,:)+mu*ux(1,:,:)/u(1,:,:)**2
      fx_ux(6,1,:,:)=mu*u(6,:,:)/u(1,:,:)**2
      fx_ux(6,6,:,:)=-mu/u(1,:,:)
      fx_uy(6,2,:,:)=u(3,:,:)
      
      fy_u(6,1,:,:)=(-u(5,:,:)*u(6,:,:)+mu*(uy(6,:,:) 
     $     -two*u(6,:,:)*uy(1,:,:)/u(1,:,:)))/u(1,:,:)**2
      fy_u(6,3,:,:)=-ux(2,:,:)
      fy_u(6,5,:,:)=u(6,:,:)/u(1,:,:)
      fy_u(6,6,:,:)=u(5,:,:)/u(1,:,:)+mu*uy(1,:,:)/u(1,:,:)**2
      fy_ux(6,2,:,:)=-u(3,:,:)
      fy_uy(6,1,:,:)=mu*u(6,:,:)/u(1,:,:)**2
      fy_uy(6,6,:,:)=-mu/u(1,:,:)
      
      s_u(6,1,:,:)=-two*grad_rho*u(6,:,:)**2/u(1,:,:)**3
      s_u(6,6,:,:)=two*grad_rho*u(6,:,:)/u(1,:,:)**2

      fx_ux(7,2,:,:)=di
      fy_uy(7,2,:,:)=di
      s_u(7,6,:,:)=one
      s_u(7,7,:,:)=-one
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDiso_drdu
c-----------------------------------------------------------------------
c     subprogram 7. HallMHDiso_equil.
c     computes equilibrium for magnetic reconnection.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDiso_equil(x,y,u,ux,uy,derivs)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: coshy_psi,coshy_phi,
     $     sinhy_psi,sinhy_phi,coshx_psi,coshx_phi,sinhx_psi,sinhx_phi
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
      IF(init_type=="Craig")THEN
         u(1,:,:)=1
         u(2,:,:)=alfven*SIN(kx*x)*SIN(ky*y)/kx
         u(7,:,:)=di*u(2,:,:)*ksq
         u(4,:,:)=-mach*ky/kx*SIN(kx*x)*COS(ky*y)
         u(5,:,:)=mach*COS(kx*x)*SIN(ky*y)
         IF(.NOT. derivs)RETURN
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         ux(2,:,:)=alfven*COS(kx*x)*SIN(ky*y)
         ux(7,:,:)=di*ux(2,:,:)*ksq
         ux(4,:,:)=-mach*ky*COS(kx*x)*COS(ky*y)
         ux(5,:,:)=-mach*kx*SIN(kx*x)*SIN(ky*y)
         uy(2,:,:)=alfven*SIN(kx*x)*COS(ky*y)*ky/kx
         uy(7,:,:)=di*uy(2,:,:)*ksq
         uy(4,:,:)=mach*ky**2/kx*SIN(kx*x)*SIN(ky*y)
         uy(5,:,:)=mach*ky*COS(kx*x)*COS(ky*y)
      ELSEIF(init_type=="Alan")THEN
         u(1,:,:)=1
         u(2,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(4,:,:)=-mach*sinhy_phi/coshy_phi
         u(7,:,:)=-di*alfven/(lambda_psi*coshy_psi**2)
         IF(.NOT. derivs)RETURN
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         uy(2,:,:)=alfven*sinhy_psi/coshy_psi
         uy(4,:,:)=-mach/lambda_phi*(1-(sinhy_phi/coshy_phi)**2)
         uy(7,:,:)=2*di*alfven/lambda_psi**2*sinhy_psi/coshy_psi**3
      ELSE
         CALL program_stop("No initial condition specified")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDiso_equil
c-----------------------------------------------------------------------
c     subprogram 8. HallMHDiso_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDiso_grid(x,y,ksi,eta)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Craig")
         ksi=lx*(x-half)
         eta=lx*(y-half)
      CASE("Alan")
         ksi=lx*(x-half)
         eta=y-half
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDiso_grid
      END MODULE HallMHDiso_mod
