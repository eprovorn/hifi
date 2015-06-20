c-----------------------------------------------------------------------
c     file fivefield.f.
c     two-fluid incompressible model with density and electron
c     viscosity/inertia.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. fivefield_mod.
c     1. fivefield_input.
c     2. fivefield_init.
c     3. fivefield_init_special.
c     4. fivefield_boundary.
c     5. fivefield_bound_rhs.
c     6. fivefield_bound_drdu.
c     7. fivefield_rhs.
c     8. fivefield_drdu.
c     9. fivefield_mass.
c     10. fivefield_equil.
c     11. fivefield_grid.
c-----------------------------------------------------------------------
c     subprogram 0. fivefield_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE fivefield_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: source=.FALSE.
      CHARACTER(16), PRIVATE :: init_type=" "
      REAL(r8), PRIVATE :: di=0,eta=0,mu=0,nu=0,lx=0,ly=0,alfven=1,
     $     mach=0,lambda_psi=0,lambda_phi=0,epsilon=0,bound_eps=0,tau=0,
     $     Bguide=0,kx=0,ky=0,ksq=0,mass_r=5.44617e-4,de_sq=0

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fivefield_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_input(nx,ny,np,nq,nbx,xperiodic,
     $     yperiodic,nqty,dt,dtmax,tmax,nstep,gr_curve)

      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic
      INTEGER, INTENT(OUT) :: nqty
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax,gr_curve

      INTEGER :: myios
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/fivefield_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,di,eta,mu,nu,lx,ly,alfven,mach,
     $     lambda_psi,lambda_phi,epsilon,bound_eps,tau,Bguide,source,
     $     init_type
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=fivefield_list,IOSTAT=myios)

      nqty=8
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_input
c-----------------------------------------------------------------------
c     subprogram 2. fivefield_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_init(xpi,ypi,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      REAL(r8), DIMENSION(0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) :: cosx,cosy,
     $     sinx,siny
      REAL(r8), DIMENSION(SIZE(ui,1),SIZE(ui,2),SIZE(ui,3)) :: uix,uiy
c-----------------------------------------------------------------------
c     magnetic reconnection, initial conditions.
c-----------------------------------------------------------------------
      CALL fivefield_equil(xpi,ypi,ui,uix,uiy,.FALSE.)
      sinx=SIN(kx*xpi)
      siny=SIN(ky*ypi)
      cosx=COS(kx*xpi)
      cosy=COS(ky*ypi)
      SELECT CASE(init_type)
      CASE("GEM")
         ui(2,:,:)=ui(2,:,:)+epsilon*cosx*cosy
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_init
c-----------------------------------------------------------------------
c     subprogram 3. fivefield_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_init_special(static,ground)

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      ground(7)=.TRUE.
      static(6:8)=.TRUE.
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
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_phi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bound_eps,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(tau,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Bguide,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alfven,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      de_sq=mass_r*di**2

      SELECT CASE(init_type)
      CASE("GEM")
         epsilon=-.1_r8
         lambda_psi=.5_r8
         lx=25.6_r8
         ly=12.8_r8
         di=one
         kx=twopi/lx
         ky=pi/ly
      END SELECT
      ksq=kx**2+ky**2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_init_special
c-----------------------------------------------------------------------
c     subprogram 4. fivefield_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_boundary(left,right,top,bottom,nqty)

      INTEGER, INTENT(IN) :: nqty
      TYPE(edge_type) :: left,right,top,bottom

      INTEGER :: k,j
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("GEM")
         top%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         bottom%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         top%b=0
         bottom%b=0
         top%bc_type(1)="natural"
         bottom%bc_type(1)="natural"
         top%static(3:8)=.TRUE.
         bottom%static(3:8)=.TRUE.
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_boundary
c-----------------------------------------------------------------------
c     subprogram 5. fivefield_bound_rhs.
c     computes rhs for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_bound_rhs(edge_type,t,x,y,u,ux,uy,uxx,uyy,c)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy,uxx,uyy
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
      SELECT CASE(init_type)
      CASE("GEM")
         SELECT CASE(edge_type)
         CASE("top","bottom")
            c(3,:)=uy(3,:)
            c(4,:)=u(4,:)
            c(5,:)=uy(5,:)
            c(6:7,:)=u(6:7,:)
            c(8,:)=uy(8,:)
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_bound_rhs
c-----------------------------------------------------------------------
c     subprogram 6. fivefield_bound_drdu.
c     computes drdu for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_bound_drdu(edge_type,x,y,u,ux,uy,uxx,uyy,
     $     c_u,c_ux,c_uy,c_uxx,c_uyy)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy,uxx,uyy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy
c-----------------------------------------------------------------------
c     zero out output.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("GEM")
         SELECT CASE(edge_type)
         CASE("top","bottom")
            c_uy(3,3,:)=one
            c_u(4,4,:)=one
            c_uy(5,5,:)=one
            c_u(6,6,:)=one
            c_u(7,7,:)=one
            c_uy(8,8,:)=one
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_bound_drdu
c-----------------------------------------------------------------------
c     subprogram 7. fivefield_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE fivefield_rhs(x,y,u,ux,uy,fx,fy,s,
     $     first)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: fi,fix,fiy,F,ui,Wi
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     two-fluid incompressible formulation with density.
c-----------------------------------------------------------------------
      fi=u(2,:,:)+mass_r*di*u(6,:,:)
      fix=ux(2,:,:)+mass_r*di*ux(6,:,:)
      fiy=uy(2,:,:)+mass_r*di*uy(6,:,:)
      F=u(3,:,:)-mass_r*di*u(8,:,:)
      ui=u(4,:,:)+mass_r*u(8,:,:)
      Wi=u(5,:,:)+mass_r*u(6,:,:)

      fx(1,:,:)=u(7,:,:)*uy(1,:,:)
      
      fy(1,:,:)=-u(7,:,:)*ux(1,:,:)
      
      fx(2,:,:)=-fi*uy(7,:,:)-eta*ux(2,:,:)-di*nu*ux(6,:,:)/u(1,:,:)

      fy(2,:,:)=fi*ux(7,:,:)-eta*uy(2,:,:)-di*nu*uy(6,:,:)/u(1,:,:)

      s(2,:,:)=di*(fix*uy(3,:,:)-fiy*ux(3,:,:))/u(1,:,:)
     $     +di*nu*(ux(6,:,:)*ux(1,:,:)+uy(6,:,:)*uy(1,:,:))
     $     /u(1,:,:)**2
      
      fx(3,:,:)=-F*(uy(7,:,:)+di*uy(3,:,:)/u(1,:,:))
     $     +u(6,:,:)*uy(2,:,:)
     $     -eta*ux(3,:,:)+di*nu*ux(8,:,:)/u(1,:,:)
 
      fy(3,:,:)=F*(ux(7,:,:)+di*ux(3,:,:)/u(1,:,:))
     $     -u(6,:,:)*ux(2,:,:)
     $     -eta*uy(3,:,:)+di*nu*uy(8,:,:)/u(1,:,:)

      fx(4,:,:)=-ui*uy(7,:,:) + F*uy(3,:,:)/u(1,:,:)
     $     +(u(5,:,:)-u(6,:,:))*uy(2,:,:)/MAX(di,min_eps)
     $     -(mu*ux(4,:,:)+nu*ux(8,:,:))/u(1,:,:)

      fy(4,:,:)=ui*ux(7,:,:) - F*ux(3,:,:)/u(1,:,:)
     $     -(u(5,:,:)-u(6,:,:))*ux(2,:,:)/MAX(di,min_eps)
     $     -(mu*uy(4,:,:)+nu*uy(8,:,:))/u(1,:,:)
      
      fx(5,:,:)=-Wi*uy(7,:,:)-(mu*ux(5,:,:)+nu*ux(6,:,:))/u(1,:,:)
      
      fy(5,:,:)=Wi*ux(7,:,:)-(mu*uy(5,:,:)+nu*uy(6,:,:))/u(1,:,:)

      s(5,:,:)=(fix*uy(3,:,:)-fiy*ux(3,:,:))/u(1,:,:)
     $     +mu*(ux(5,:,:)*ux(1,:,:)+uy(5,:,:)*uy(1,:,:))
     $     /u(1,:,:)**2
     $     +nu*(ux(6,:,:)*ux(1,:,:)+uy(6,:,:)*uy(1,:,:))
     $     /u(1,:,:)**2
      
      fx(6,:,:)=di*ux(2,:,:)
      fy(6,:,:)=di*uy(2,:,:)
      s(6,:,:)=u(1,:,:)*(u(5,:,:)-u(6,:,:))
      
      fx(7,:,:)=ux(7,:,:)
      fy(7,:,:)=uy(7,:,:)
      s(7,:,:)=u(4,:,:)

      fx(8,:,:)=di*ux(3,:,:)/u(1,:,:)+ux(7,:,:)
      fy(8,:,:)=di*uy(3,:,:)/u(1,:,:)+uy(7,:,:)
      s(8,:,:)=u(8,:,:)

      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("GEM")
            CALL fivefield_equil(x,y,u0,u0x,u0y,.TRUE.)
            first=.FALSE.
            CALL fivefield_rhs(x,y,u0,u0x,u0y,fx0,fy0,s0,first)
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
      END SUBROUTINE fivefield_rhs
c-----------------------------------------------------------------------
c     subprogram 8. fivefield_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_drdu(x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: fi,fix,fiy,F,ui,Wi
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
c     prepare intermediate variables
c-----------------------------------------------------------------------
      fi=u(2,:,:)+mass_r*di*u(6,:,:)
      fix=ux(2,:,:)+mass_r*di*ux(6,:,:)
      fiy=uy(2,:,:)+mass_r*di*uy(6,:,:)
      F=u(3,:,:)-mass_r*di*u(8,:,:)
      ui=u(4,:,:)+mass_r*u(8,:,:)
      Wi=u(5,:,:)+mass_r*u(6,:,:)
c-----------------------------------------------------------------------
c     two-fluid incompressible reconnection.
c-----------------------------------------------------------------------
      fx_u(1,7,:,:)=uy(1,:,:)
      fx_uy(1,1,:,:)=u(7,:,:)
      
      fy_u(1,7,:,:)=-ux(1,:,:)
      fy_ux(1,1,:,:)=-u(7,:,:)
      
      fx_u(2,1,:,:)=di*nu*ux(6,:,:)/u(1,:,:)**2
      fx_u(2,2,:,:)=-uy(7,:,:)
      fx_u(2,6,:,:)=-mass_r*di*uy(7,:,:)
      fx_ux(2,2,:,:)=-eta
      fx_ux(2,6,:,:)=-di*nu/u(1,:,:)
      fx_uy(2,7,:,:)=-fi
      
      fy_u(2,1,:,:)=di*nu*uy(6,:,:)/u(1,:,:)**2
      fy_u(2,2,:,:)=ux(7,:,:)
      fy_u(2,6,:,:)=mass_r*di*ux(7,:,:)
      fy_ux(2,7,:,:)=fi
      fy_uy(2,2,:,:)=-eta
      fy_uy(2,6,:,:)=-di*nu/u(1,:,:)

      s_u(2,1,:,:)=-di*(fix*uy(3,:,:)-fiy*ux(3,:,:))/u(1,:,:)**2
     $     -2*di*nu
     $     *(ux(6,:,:)*ux(1,:,:)+uy(6,:,:)*uy(1,:,:))/u(1,:,:)**3
      s_ux(2,1,:,:)=di*nu*ux(6,:,:)/u(1,:,:)**2
      s_ux(2,2,:,:)=di*uy(3,:,:)/u(1,:,:)
      s_ux(2,3,:,:)=-di*fiy/u(1,:,:)
      s_ux(2,6,:,:)=de_sq*uy(3,:,:)/u(1,:,:)
     $     +di*nu*ux(1,:,:)/u(1,:,:)**2
      s_uy(2,1,:,:)=di*nu*uy(6,:,:)/u(1,:,:)**2
      s_uy(2,2,:,:)=-di*ux(3,:,:)/u(1,:,:)
      s_uy(2,3,:,:)=di*fix/u(1,:,:)
      s_uy(2,6,:,:)=-de_sq*ux(3,:,:)/u(1,:,:)
     $     +di*nu*uy(1,:,:)/u(1,:,:)**2
      
      fx_u(3,1,:,:)=di*F*uy(3,:,:)/u(1,:,:)**2
     $     -di*nu*ux(8,:,:)/u(1,:,:)**2
      fx_u(3,3,:,:)=-(uy(7,:,:)+di*uy(3,:,:)/u(1,:,:))
      fx_u(3,6,:,:)=uy(2,:,:)
      fx_u(3,8,:,:)=mass_r*di*uy(7,:,:)+de_sq*uy(3,:,:)/u(1,:,:) 
      fx_ux(3,3,:,:)=-eta
      fx_ux(3,8,:,:)=di*nu/u(1,:,:)
      fx_uy(3,2,:,:)=u(6,:,:)
      fx_uy(3,3,:,:)=-di*F/u(1,:,:)
      fx_uy(3,7,:,:)=-F
      
      fy_u(3,1,:,:)=-di*F*ux(3,:,:)/u(1,:,:)**2
     $     -di*nu*uy(8,:,:)/u(1,:,:)**2
      fy_u(3,3,:,:)=(ux(7,:,:)+di*ux(3,:,:)/u(1,:,:))
      fy_u(3,6,:,:)=-ux(2,:,:)
      fy_u(3,8,:,:)=-mass_r*di*ux(7,:,:)-de_sq*ux(3,:,:)/u(1,:,:) 
      fy_ux(3,2,:,:)=-u(6,:,:)
      fy_ux(3,3,:,:)=di*F/u(1,:,:)
      fy_ux(3,7,:,:)=F
      fy_uy(3,3,:,:)=-eta
      fy_uy(3,8,:,:)=di*nu/u(1,:,:)

      fx_u(4,1,:,:)=-F*uy(3,:,:)/u(1,:,:)**2
     $     +(mu*ux(4,:,:)+nu*ux(8,:,:))/u(1,:,:)**2
      fx_u(4,3,:,:)=uy(3,:,:)/u(1,:,:)
      fx_u(4,4,:,:)=-uy(7,:,:)
      fx_u(4,5,:,:)=uy(2,:,:)/MAX(di,min_eps)
      fx_u(4,6,:,:)=-uy(2,:,:)/MAX(di,min_eps)
      fx_u(4,8,:,:)=-mass_r*uy(7,:,:)-mass_r*di*uy(3,:,:)/u(1,:,:)
      fx_ux(4,4,:,:)=-mu/u(1,:,:)
      fx_ux(4,8,:,:)=-nu/u(1,:,:)
      fx_uy(4,2,:,:)=(u(5,:,:)-u(6,:,:))/MAX(di,min_eps)
      fx_uy(4,3,:,:)=F/u(1,:,:)
      fx_uy(4,7,:,:)=-ui
      
      fy_u(4,1,:,:)=F*ux(3,:,:)/u(1,:,:)**2
     $     +(mu*uy(4,:,:)+nu*uy(8,:,:))/u(1,:,:)**2
      fy_u(4,3,:,:)=-ux(3,:,:)/u(1,:,:)
      fy_u(4,4,:,:)=ux(7,:,:)
      fy_u(4,5,:,:)=-ux(2,:,:)/MAX(di,min_eps)
      fy_u(4,6,:,:)=ux(2,:,:)/MAX(di,min_eps)
      fy_u(4,8,:,:)=mass_r*ux(7,:,:)+mass_r*di*ux(3,:,:)/u(1,:,:)
      fy_ux(4,2,:,:)=-(u(5,:,:)-u(6,:,:))/MAX(di,min_eps)
      fy_ux(4,3,:,:)=-F/u(1,:,:)
      fy_ux(4,7,:,:)=ui
      fy_uy(4,4,:,:)=-mu/u(1,:,:)
      fy_uy(4,8,:,:)=-nu/u(1,:,:)
      
      fx_u(5,1,:,:)=(mu*ux(5,:,:)+nu*ux(6,:,:))/u(1,:,:)**2
      fx_u(5,5,:,:)=-uy(7,:,:)
      fx_u(5,6,:,:)=-mass_r*uy(7,:,:)
      fx_ux(5,5,:,:)=-mu/u(1,:,:)
      fx_ux(5,6,:,:)=-nu/u(1,:,:)
      fx_uy(5,7,:,:)=-Wi

      fy_u(5,1,:,:)=(mu*uy(5,:,:)+nu*uy(6,:,:))/u(1,:,:)**2
      fy_u(5,5,:,:)=ux(7,:,:)
      fy_u(5,6,:,:)=mass_r*ux(7,:,:)
      fy_ux(5,7,:,:)=Wi
      fy_uy(5,5,:,:)=-mu/u(1,:,:)
      fy_uy(5,6,:,:)=-nu/u(1,:,:)

      s_u(5,1,:,:)=-(fix*uy(3,:,:)-fiy*ux(3,:,:))/u(1,:,:)**2
     $     -2*mu*(ux(5,:,:)*ux(1,:,:)+uy(5,:,:)*uy(1,:,:))/u(1,:,:)**3
     $     -2*nu*(ux(6,:,:)*ux(1,:,:)+uy(6,:,:)*uy(1,:,:))/u(1,:,:)**3
      s_ux(5,1,:,:)=(mu*ux(5,:,:)+nu*ux(6,:,:))/u(1,:,:)**2
      s_ux(5,2,:,:)=uy(3,:,:)/u(1,:,:)
      s_ux(5,3,:,:)=-fiy/u(1,:,:)
      s_ux(5,5,:,:)=mu*ux(1,:,:)/u(1,:,:)**2
      s_ux(5,6,:,:)=mass_r*di*uy(3,:,:)/u(1,:,:)
     $     +nu*ux(1,:,:)/u(1,:,:)**2
      s_uy(5,1,:,:)=(mu*uy(5,:,:)+nu*uy(6,:,:))/u(1,:,:)**2
      s_uy(5,2,:,:)=-ux(3,:,:)/u(1,:,:)
      s_uy(5,3,:,:)=fix/u(1,:,:)
      s_uy(5,5,:,:)=mu*uy(1,:,:)/u(1,:,:)**2
      s_uy(5,6,:,:)=-mass_r*di*ux(3,:,:)/u(1,:,:)
     $     +nu*uy(1,:,:)/u(1,:,:)**2
      
      fx_ux(6,2,:,:)=di
      fy_uy(6,2,:,:)=di
      s_u(6,1,:,:)=u(5,:,:)-u(6,:,:)
      s_u(6,5,:,:)=u(1,:,:)
      s_u(6,6,:,:)=-u(1,:,:)
 
      fx_ux(7,7,:,:)=one
      fy_uy(7,7,:,:)=one
      s_u(7,4,:,:)=one
      
      fx_u(8,1,:,:)=-di*ux(3,:,:)/u(1,:,:)**2
      fx_ux(8,3,:,:)=di/u(1,:,:)
      fx_ux(8,7,:,:)=one
      fy_u(8,1,:,:)=-di*uy(3,:,:)/u(1,:,:)**2
      fy_uy(8,3,:,:)=di/u(1,:,:)
      fy_uy(8,7,:,:)=one
      s_u(8,8,:,:)=one
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_drdu
c-----------------------------------------------------------------------
c     subprogram 9. fivefield_mass.
c     computes mass matrix couplings.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_mass(mass)

      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      mass(2,6,:,:)=mass_r*di
      mass(3,8,:,:)=-mass_r*di
      mass(4,8,:,:)=mass_r
      mass(5,6,:,:)=mass_r
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_mass
c-----------------------------------------------------------------------
c     subprogram 10. fivefield_equil.
c     computes equilibrium for magnetic reconnection.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_equil(x,y,u,ux,uy,deriv)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: deriv

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
      SELECT CASE(init_type)
      CASE("GEM")
         u(1,:,:)=one/coshy_psi**2 + .2_r8
         u(2,:,:)=-lambda_psi*LOG(coshy_psi)
         u(3,:,:)=Bguide
         u(6,:,:)=di/(lambda_psi*u(1,:,:)*coshy_psi**2)
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=-two*TANH(y/lambda_psi)/(lambda_psi*coshy_psi**2)
         uy(2,:,:)=-TANH(y/lambda_psi)
         uy(6,:,:)=-u(6,:,:)*(two*TANH(y/lambda_psi)/lambda_psi
     $        +uy(1,:,:)/u(1,:,:))
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_equil
c-----------------------------------------------------------------------
c     subprogram 11. fivefield_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fivefield_grid(x,y,ksi,eta)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("GEM")
         ksi=lx*(x-half)
         eta=ly*(y-half)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fivefield_grid
      END MODULE fivefield_mod
