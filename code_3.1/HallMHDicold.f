c-----------------------------------------------------------------------
c     file HallMHDicold.f.
c     contains specifications for Hall MHD model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. HallMHDicold_mod.
c     1. HallMHDicold_input.
c     2. HallMHDicold_init.
c     3. HallMHDicold_init_special.
c     4. HallMHDicold_boundary.
c     5. HallMHDicold_bound_rhs.
c     6. HallMHDicold_bound_drdu.
c     7. HallMHDicold_rhs.
c     8. HallMHDicold_drdu.
c     9. HallMHDicold_mass.
c     10. HallMHDicold_equil.
c     11. HallMHDicold_grid.
c     12. HallMHDicold_BdotT.
c-----------------------------------------------------------------------
c     subprogram 0. HallMHDicold_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE HallMHDicold_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: xper=.FALSE.,eta_spitzer=.FALSE.,
     $     source=.FALSE.
      CHARACTER(16), PRIVATE :: init_type=" "
      REAL(r8), PARAMETER, PRIVATE :: gamma=5._r8/3._r8,Bmin=1e-16
      REAL(r8), PRIVATE :: di=0,eta=0,mu=0,nu=0,kappa_par=0,
     $     kappa_perp=0,Dn=0,lx=0,ly=0,mach=0,alfven=0,lambda_psi=0,
     $     lambda_phi=0,epsilon=0,Bguide=0,beta_e=0,beta0=0,
     $     kx=0,ky=0,ksq=0

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. HallMHDicold_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_input(nx,ny,np,nq,nbx,xperiodic,
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
      xper=xperiodic
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_input
c-----------------------------------------------------------------------
c     subprogram 2. HallMHDicold_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_init(xpi,ypi,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: ix,iy
      REAL(r8), DIMENSION(0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) :: cosx,cosy,
     $     sinx,siny,coshy_psi,tanhx_phi

      REAL(r8), DIMENSION(SIZE(ui,1),SIZE(ui,2),SIZE(ui,3)) :: uix,uiy
c-----------------------------------------------------------------------
c     magnetic reconnection, initial conditions.
c-----------------------------------------------------------------------
      CALL HallMHDicold_equil(xpi,ypi,ui,uix,uiy,.FALSE.)
      sinx=SIN(kx*xpi)
      siny=SIN(ky*ypi)
      cosx=COS(kx*xpi)
      cosy=COS(ky*ypi)
      coshy_psi=COSH(ypi/lambda_psi)
      tanhx_phi=TANH(xpi/lambda_phi)
      SELECT CASE(init_type)
      CASE("Harris")
         ui(1,:,:)=ui(1,:,:)-epsilon*cosx*cosy*alfven/ksq
         ui(3,:,:)=ui(3,:,:)-epsilon*cosx*siny*mach*ky/ksq
         ui(4,:,:)=ui(4,:,:)+epsilon*sinx*cosy*mach*kx/ksq
         ui(6,:,:)=ui(6,:,:)-di*epsilon*cosx*cosy*alfven
      CASE("GEM","history")
         ui(1,:,:)=ui(1,:,:) + epsilon*cosx*cosy
      CASE("history_exp")
         ui(1,:,:)=ui(1,:,:) + epsilon
     $        *EXP(-xpi**2/(two*lambda_psi)**2)
     $        *EXP(-ypi**2/(half*lambda_psi)**2)
      CASE("GEM_open")
         DO ix=0,SIZE(xpi,1)-1
            DO iy=0,SIZE(xpi,2)-1
               IF(xpi(ix,iy) <= lx/4 .AND. ypi(ix,iy) <= ly/2)THEN
                  ui(1,ix,iy)=ui(1,ix,iy)+epsilon*COS(kx*xpi(ix,iy))
     $                 *COS(ky*ypi(ix,iy))
               ENDIF
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     abort if not recognized.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT

      RETURN
      END SUBROUTINE HallMHDicold_init
c-----------------------------------------------------------------------
c     subprogram 3. HallMHDicold_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_init_special(static)

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(6)=.TRUE.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(eta_spitzer,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(xper,1,MPI_LOGICAL,0,comm,ierr)
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
         kx=twopi/lx
         ky=pi/ly
      CASE("history","history_exp")
         lambda_psi=.5_r8
         kx=twopi/lx
         ky=pi/ly
      END SELECT
      ksq=kx**2+ky**2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_init_special
c-----------------------------------------------------------------------
c     subprogram 4. HallMHDicold_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_boundary(left,right,top,bottom,nqty)

      INTEGER, INTENT(IN) :: nqty
      TYPE(edge_type) :: left,right,top,bottom

      INTEGER :: k,j
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Harris")
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
         top%bc_type(6:7)="natural"
         bottom%bc_type(6:7)="natural"
         top%static(3:5)=.TRUE.
         bottom%static(3:5)=.TRUE.
      CASE("GEM")
         top%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         top%b=0
         top%static(2:7)=.TRUE.

         bottom%static=.TRUE.
         left%static=.TRUE.
         right%static=.TRUE.
      CASE("GEM_open")
         top%static=.TRUE.
         right%static=.TRUE.
         left%static=.TRUE.
         bottom%static=.TRUE.
      CASE("history","history_exp")
         top%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         top%b=0
         top%static(2:7)=.TRUE.
         bottom%a=RESHAPE((/(1.,(0.,k=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         bottom%b=0
         bottom%static(2:7)=.TRUE.
         left%static=.TRUE.
         right%static=.TRUE.
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_boundary
c-----------------------------------------------------------------------
c     subprogram 5. HallMHDicold_bound_rhs.
c     computes rhs for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_bound_rhs(edge_type,x,y,u,ux,uy,uxx,uyy,c)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy,uxx,uyy
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Harris")
         SELECT CASE(edge_type)
         CASE("top","bottom")
            c(3,:)=uy(3,:)
            c(4,:)=u(4,:)
            c(5,:)=uy(5,:)
         END SELECT
      CASE("GEM")
         SELECT CASE(edge_type)
         CASE("top")
            c(2,:)=-u(2,:)*ux(2,:)+eta*uy(2,:)
            c(3,:)=uy(3,:)
            c(4,:)=u(4,:)
            c(5,:)=uy(5,:)
            c(6,:)=u(6,:)-u(5,:)
            c(7,:)=uy(7,:)
         CASE("bottom")
c-----------------------------------------------------------------------
c     even boundary conditions
c-----------------------------------------------------------------------
            c(1,:)=uy(1,:)
            c(3,:)=uy(3,:)
            c(5,:)=uy(5,:)
            c(6,:)=uy(6,:)
            c(7,:)=uy(7,:)
c-----------------------------------------------------------------------
c     odd boundary conditions
c-----------------------------------------------------------------------
            c(2,:)=u(2,:)
            c(4,:)=u(4,:)
         CASE("left","right")
c-----------------------------------------------------------------------
c     even boundary conditions
c-----------------------------------------------------------------------
            c(1,:)=ux(1,:)
            c(4,:)=ux(4,:)
            c(5,:)=ux(5,:)
            c(6,:)=ux(6,:)
            c(7,:)=ux(7,:)
c-----------------------------------------------------------------------
c     odd boundary conditions
c-----------------------------------------------------------------------
            c(2,:)=u(2,:)
            c(3,:)=u(3,:)
         END SELECT
      CASE("GEM_open")
         SELECT CASE(edge_type)
         CASE("top")
            c(1,:)=uyy(1,:)
            c(2,:)=uy(2,:)
            c(3,:)=u(3,:)
            c(4,:)=uy(4,:)
            c(5,:)=uy(5,:)
            c(6,:)=uy(6,:)
            c(7,:)=uy(7,:)
         CASE("bottom")
c-----------------------------------------------------------------------
c     even boundary conditions
c-----------------------------------------------------------------------
            c(1,:)=uy(1,:)
            c(3,:)=uy(3,:)
            c(5,:)=uy(5,:)
            c(6,:)=uy(6,:)
            c(7,:)=uy(7,:)
c-----------------------------------------------------------------------
c     odd boundary conditions
c-----------------------------------------------------------------------
            c(2,:)=u(2,:)
            c(4,:)=u(4,:)
         CASE("left")
c-----------------------------------------------------------------------
c     even boundary conditions
c-----------------------------------------------------------------------
            c(1,:)=ux(1,:)
            c(4,:)=ux(4,:)
            c(5,:)=ux(5,:)
            c(6,:)=ux(6,:)
            c(7,:)=ux(7,:)
c-----------------------------------------------------------------------
c     odd boundary conditions
c-----------------------------------------------------------------------
            c(2,:)=u(2,:)
            c(3,:)=u(3,:)
         CASE("right")
            c(1,:)=uxx(1,:)
            c(2,:)=ux(2,:)
            c(3,:)=ux(3,:)
            c(4,:)=u(4,:)
            c(5,:)=ux(5,:)
            c(6,:)=ux(6,:)
            c(7,:)=ux(7,:)
         END SELECT
      CASE("history","history_exp")
         SELECT CASE(edge_type)
         CASE("top","bottom")
            c(2,:)=uy(2,:)
            c(3,:)=uy(3,:)
            c(4,:)=u(4,:)
            c(5,:)=uy(5,:)
            c(6,:)=uy(6,:)
            c(7,:)=uy(7,:)
         CASE("left","right")
            c(1,:)=ux(1,:)
            c(2,:)=ux(2,:)
            c(3,:)=ux(3,:)
            c(4,:)=ux(4,:)
            c(5,:)=ux(5,:)
            c(6,:)=ux(6,:)
            c(7,:)=ux(7,:)
         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_bound_rhs
c-----------------------------------------------------------------------
c     subprogram 6. HallMHDicold_bound_drdu.
c     computes drdu for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_bound_drdu(edge_type,x,y,u,ux,uy,uxx,uyy,
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
      CASE("Harris")
         SELECT CASE(edge_type)
         CASE("top","bottom")
            c_uy(3,3,:)=1
            c_u(4,4,:)=1
            c_uy(5,5,:)=1
         END SELECT
      CASE("GEM")
         SELECT CASE(edge_type)
         CASE("top")
            c_u(2,2,:)=-ux(2,:)
            c_ux(2,2,:)=-u(2,:)
            c_uy(2,2,:)=eta
            c_uy(3,3,:)=one
            c_u(4,4,:)=one
            c_uy(5,5,:)=one
            c_u(6,5,:)=-one
            c_u(6,6,:)=one
            c_uy(7,7,:)=one
         CASE("bottom")
            c_uy(1,1,:)=one
            c_u(2,2,:)=one
            c_uy(3,3,:)=one
            c_u(4,4,:)=one
            c_uy(5,5,:)=one
            c_uy(6,6,:)=one
            c_uy(7,7,:)=one
         CASE("left","right")
            c_ux(1,1,:)=one
            c_u(2,2,:)=one
            c_u(3,3,:)=one
            c_ux(4,4,:)=one
            c_ux(5,5,:)=one
            c_ux(6,6,:)=one
            c_ux(7,7,:)=one
        END SELECT
      CASE("GEM_open")
         SELECT CASE(edge_type)
         CASE("top")
            c_uyy(1,1,:)=one
            c_uy(2,2,:)=one
            c_u(3,3,:)=one
            c_uy(4,4,:)=one
            c_uy(5,5,:)=one
            c_uy(6,6,:)=one
            c_uy(7,7,:)=one
         CASE("bottom")
            c_uy(1,1,:)=one
            c_u(2,2,:)=one
            c_uy(3,3,:)=one
            c_u(4,4,:)=one
            c_uy(5,5,:)=one
            c_uy(6,6,:)=one
            c_uy(7,7,:)=one
         CASE("left")
            c_ux(1,1,:)=one
            c_u(2,2,:)=one
            c_u(3,3,:)=one
            c_ux(4,4,:)=one
            c_ux(5,5,:)=one
            c_ux(6,6,:)=one
            c_ux(7,7,:)=one
         CASE("right")
            c_uxx(1,1,:)=one
            c_ux(2,2,:)=one
            c_ux(3,3,:)=one
            c_u(4,4,:)=one
            c_ux(5,5,:)=one
            c_ux(6,6,:)=one
            c_ux(7,7,:)=one
        END SELECT
      CASE("history","history_exp")
         SELECT CASE(edge_type)
         CASE("top","bottom")
            c_uy(2,2,:)=one
            c_uy(3,3,:)=one
            c_u(4,4,:)=one
            c_uy(5,5,:)=one
            c_uy(6,6,:)=one
            c_uy(7,7,:)=one
         CASE("left","right")
            c_ux(1,1,:)=one
            c_ux(2,2,:)=one
            c_ux(3,3,:)=one
            c_ux(4,4,:)=one
            c_ux(5,5,:)=one
            c_ux(6,6,:)=one
            c_ux(7,7,:)=one
        END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_bound_drdu
c-----------------------------------------------------------------------
c     subprogram 7. HallMHDicold_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE HallMHDicold_rhs(x,y,u,ux,uy,fx,fy,s,first)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: ix,iy
      REAL(r8) :: gamma_fac,Bsq
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: invBsq,
     $     kappa_prp,kappa_fac,mu_local,eta_local
      REAL(r8), DIMENSION(2,SIZE(u,2),SIZE(u,3)) :: BdotT
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0      
      u0=0
      u0x=0
      u0y=0
c-----------------------------------------------------------------------
c     prepare additional variables.
c-----------------------------------------------------------------------
      gamma_fac=gamma/(gamma-one)
      mu_local=mu
      IF((init_type=="history" .OR. init_type=="history_exp")
     $     .AND. .NOT. xper)
     $     mu_local=mu + 1.e-1
     $     *(EXP(-(20*(x-lx/2)/lx)**4) + EXP(-(20*(x+lx/2)/lx)**4))
c-----------------------------------------------------------------------
c     initial equilibrium source term.
c-----------------------------------------------------------------------
      IF((init_type == "history" .OR. init_type=="history_exp")
     $     .AND. source)CALL HallMHDicold_equil(x,y,u0,u0x,u0y,.TRUE.)
c-----------------------------------------------------------------------
c     if needed, make resistivity temperature dependent
c-----------------------------------------------------------------------
      eta_local=eta
      IF(eta_spitzer)eta_local=eta*(beta0/u(7,:,:))**1.5
c-----------------------------------------------------------------------
c     Cold ion Hall MHD formulation.
c-----------------------------------------------------------------------
      fx(1,:,:)=-eta_local*(ux(1,:,:)-u0x(1,:,:))
     $     -di*nu*(ux(6,:,:)-u0x(6,:,:))
      fy(1,:,:)=-eta_local*(uy(1,:,:)-u0y(1,:,:))
     $     -di*nu*(uy(6,:,:)-u0y(6,:,:))

      s(1,:,:)=-(u(3,:,:)*ux(1,:,:)+u(4,:,:)*uy(1,:,:))
     $     -di*(ux(2,:,:)*uy(1,:,:)-uy(2,:,:)*ux(1,:,:))

      fx(2,:,:)=(u(2,:,:)*u(3,:,:)+u(6,:,:)*uy(1,:,:))
     $     -eta_local*(ux(2,:,:)-u0x(2,:,:))
      
      fy(2,:,:)=(u(2,:,:)*u(4,:,:)-u(6,:,:)*ux(1,:,:))
     $     -eta_local*(uy(2,:,:)-u0y(2,:,:))

      fx(3,:,:)=u(3,:,:)**2 + u(7,:,:)
     $     +half*(ux(1,:,:)**2-uy(1,:,:)**2+u(2,:,:)**2)
     $     -mu_local*2*ux(3,:,:)
         
      fy(3,:,:)=u(3,:,:)*u(4,:,:) + uy(1,:,:)*ux(1,:,:)
     $     -mu_local*(ux(4,:,:)+uy(3,:,:))
      
      fx(4,:,:)=u(3,:,:)*u(4,:,:) + uy(1,:,:)*ux(1,:,:)
     $     -mu_local*(ux(4,:,:)+uy(3,:,:))
      
      fy(4,:,:)=u(4,:,:)**2 + u(7,:,:)
     $     +half*(-ux(1,:,:)**2+uy(1,:,:)**2+u(2,:,:)**2)
     $     -mu_local*2*uy(4,:,:)
      
      fx(5,:,:)=u(3,:,:)*u(5,:,:) + u(2,:,:)*uy(1,:,:)
     $     -mu_local*ux(5,:,:)-nu*(ux(6,:,:)-u0x(6,:,:))
      
      fy(5,:,:)=u(4,:,:)*u(5,:,:) - u(2,:,:)*ux(1,:,:)
     $     -mu_local*uy(5,:,:)-nu*(uy(6,:,:)-u0y(6,:,:))
      
      fx(6,:,:)=di*ux(1,:,:)
      fy(6,:,:)=di*uy(1,:,:)
      s(6,:,:)=u(5,:,:)-u(6,:,:)
c-----------------------------------------------------------------------
c     ion pressure equation.
c-----------------------------------------------------------------------
      IF(kappa_par > 0 .OR. kappa_perp > 0)THEN
         invBsq=one/(ux(1,:,:)**2+uy(1,:,:)**2+u(2,:,:)**2)
         BdotT(1,:,:)=uy(1,:,:)*invBsq*
     $        (uy(1,:,:)*ux(7,:,:) - ux(1,:,:)*uy(7,:,:))
         BdotT(2,:,:)=-ux(1,:,:)*invBsq*
     $        (uy(1,:,:)*ux(7,:,:) - ux(1,:,:)*uy(7,:,:))
         DO ix=1,SIZE(u,2)
            DO iy=1,SIZE(u,3)
               Bsq=ux(1,ix,iy)**2+uy(1,ix,iy)**2+u(2,ix,iy)**2
               kappa_prp(ix,iy)=kappa_par*kappa_perp/
     $              (Bsq*kappa_par + kappa_perp)
               IF(Bsq <= Bmin)THEN
                  BdotT(:,ix,iy)=zero
               ENDIF
            ENDDO
         ENDDO
         kappa_fac=kappa_par-kappa_prp
      ELSE
         kappa_prp=0
         kappa_fac=0
         BdotT=0
      ENDIF

      fx(7,:,:)=gamma_fac*u(7,:,:)*(u(3,:,:)-di*uy(2,:,:))
     $     -kappa_fac*BdotT(1,:,:) - kappa_prp*ux(7,:,:)
      fy(7,:,:)=gamma_fac*u(7,:,:)*(u(4,:,:)+di*ux(2,:,:))
     $     -kappa_fac*BdotT(2,:,:) - kappa_prp*uy(7,:,:)
      s(7,:,:)=(u(3,:,:) - di*uy(2,:,:))*ux(7,:,:)+
     $     (u(4,:,:) + di*ux(2,:,:))*uy(7,:,:)
     $     + eta_local*((u(5,:,:)-u(6,:,:)+u0(6,:,:))**2/di**2 
     $     + (ux(2,:,:)-u0x(2,:,:))**2 + (uy(2,:,:)-u0y(2,:,:))**2) 
     $     + nu*((ux(6,:,:)-u0x(6,:,:))**2 
     $     + (uy(6,:,:)-u0y(6,:,:))**2)
c-----------------------------------------------------------------------
c     initial equilibrium source term.
c-----------------------------------------------------------------------
      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("GEM","GEM_open")
            CALL HallMHDicold_equil(x,y,u0,u0x,u0y,.TRUE.)
            first=.FALSE.
            CALL HallMHDicold_rhs(x,y,u0,u0x,u0y,fx0,fy0,s0,first)
            fx=fx-fx0
            fy=fy-fy0
            s=s-s0
         END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_rhs
c-----------------------------------------------------------------------
c     subprogram 8. HallMHDicold_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_drdu(x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8) :: gamma_fac
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: kappa_prp,kappa_fac,
     $     kappa_u2,kappa_ux1,kappa_uy1,mu_local,eta_local,eta_u7
      REAL(r8), DIMENSION(2,SIZE(u,2),SIZE(u,3)) :: BdotT,BdotT_u2,
     $     BdotT_ux1,BdotT_ux7,BdotT_uy1,BdotT_uy7
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y
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
      u0=0
      u0x=0
      u0y=0
c-----------------------------------------------------------------------
c     prepare additional variables.
c-----------------------------------------------------------------------
      gamma_fac=gamma/(gamma-one)
      mu_local=mu
      IF((init_type=="history" .OR. init_type=="history_exp")
     $     .AND. .NOT. xper)THEN
         mu_local=mu + 1.e-1
     $        *(EXP(-(20*(x-lx/2)/lx)**4) + EXP(-(20*(x+lx/2)/lx)**4))
      ENDIF
c-----------------------------------------------------------------------
c     initial equilibrium source term.
c-----------------------------------------------------------------------
      IF((init_type=="history" .OR. init_type=="history_exp")
     $     .AND. source)THEN
         CALL HallMHDicold_equil(x,y,u0,u0x,u0y,.TRUE.)
      ENDIF
c-----------------------------------------------------------------------
c     if needed, make resistivity temperature dependent
c-----------------------------------------------------------------------
      eta_local=eta
      eta_u7=0
      IF(eta_spitzer)THEN
         eta_local=eta*(beta0/u(7,:,:))**1.5
         eta_u7=-1.5*eta*beta0**1.5/u(7,:,:)**2.5
      ENDIF
c-----------------------------------------------------------------------
c     Cold ion Hall MHD.
c-----------------------------------------------------------------------
      fx_u(1,7,:,:)=-eta_u7*(ux(1,:,:)-u0x(1,:,:))
      fx_ux(1,1,:,:)=-eta_local
      fx_ux(1,6,:,:)=-di*nu

      fy_u(1,7,:,:)=-eta_u7*(uy(1,:,:)-u0y(1,:,:))
      fy_uy(1,1,:,:)=-eta_local
      fy_uy(1,6,:,:)=-di*nu

      s_u(1,3,:,:)=-ux(1,:,:)
      s_u(1,4,:,:)=-uy(1,:,:)
      s_ux(1,1,:,:)=(di*uy(2,:,:)-u(3,:,:))
      s_ux(1,2,:,:)=-di*uy(1,:,:)
      s_uy(1,1,:,:)=-(u(4,:,:)+di*ux(2,:,:))
      s_uy(1,2,:,:)=di*ux(1,:,:)

      fx_u(2,2,:,:)=u(3,:,:)
      fx_u(2,3,:,:)=u(2,:,:)
      fx_u(2,6,:,:)=uy(1,:,:)
      fx_u(2,7,:,:)=-eta_u7*(ux(2,:,:)-u0x(2,:,:))
      fx_ux(2,2,:,:)=-eta_local
      fx_uy(2,1,:,:)=u(6,:,:)
 
      fy_u(2,2,:,:)=u(4,:,:)
      fy_u(2,4,:,:)=u(2,:,:)
      fy_u(2,6,:,:)=-ux(1,:,:)
      fy_u(2,7,:,:)=-eta_u7*(uy(2,:,:)-u0y(2,:,:))
      fy_ux(2,1,:,:)=-u(6,:,:)
      fy_uy(2,2,:,:)=-eta_local
         
      fx_u(3,2,:,:)=u(2,:,:)
      fx_u(3,3,:,:)=two*u(3,:,:)
      fx_u(3,7,:,:)=one
      fx_ux(3,1,:,:)=ux(1,:,:)
      fx_ux(3,3,:,:)=-2*mu_local
      fx_uy(3,1,:,:)=-uy(1,:,:)

      fy_u(3,3,:,:)=u(4,:,:)
      fy_u(3,4,:,:)=u(3,:,:)
      fy_ux(3,1,:,:)=uy(1,:,:)
      fy_ux(3,4,:,:)=-mu_local
      fy_uy(3,1,:,:)=ux(1,:,:)
      fy_uy(3,3,:,:)=-mu_local

      fx_u(4,3,:,:)=u(4,:,:)
      fx_u(4,4,:,:)=u(3,:,:)
      fx_ux(4,1,:,:)=uy(1,:,:)
      fx_ux(4,4,:,:)=-mu_local
      fx_uy(4,1,:,:)=ux(1,:,:)
      fx_uy(4,3,:,:)=-mu_local
      
      fy_u(4,2,:,:)=u(2,:,:)
      fy_u(4,4,:,:)=two*u(4,:,:)
      fy_u(4,7,:,:)=one
      fy_ux(4,1,:,:)=-ux(1,:,:)
      fy_uy(4,1,:,:)=uy(1,:,:)
      fy_uy(4,4,:,:)=-2*mu_local
      
      fx_u(5,2,:,:)=uy(1,:,:)
      fx_u(5,3,:,:)=u(5,:,:)
      fx_u(5,5,:,:)=u(3,:,:)
      fx_ux(5,5,:,:)=-mu_local
      fx_ux(5,6,:,:)=-nu
      fx_uy(5,1,:,:)=u(2,:,:)
      
      fy_u(5,2,:,:)=-ux(1,:,:)
      fy_u(5,4,:,:)=u(5,:,:)
      fy_u(5,5,:,:)=u(4,:,:)
      fy_ux(5,1,:,:)=-u(2,:,:)
      fy_uy(5,5,:,:)=-mu_local
      fy_uy(5,6,:,:)=-nu
      
      fx_ux(6,1,:,:)=di
      fy_uy(6,1,:,:)=di
      s_u(6,5,:,:)=one
      s_u(6,6,:,:)=-one
c-----------------------------------------------------------------------
c     ion pressure equation.
c-----------------------------------------------------------------------
      IF(kappa_par > 0 .OR. kappa_perp > 0)THEN
         CALL HallMHDicold_BdotT(u,ux,uy,kappa_prp,kappa_u2,kappa_ux1,
     $        kappa_uy1,BdotT,BdotT_u2,BdotT_ux1,BdotT_ux7,BdotT_uy1,
     $        BdotT_uy7)
         kappa_fac=kappa_par-kappa_prp
      ELSE
         kappa_prp=0
         kappa_fac=0
         kappa_u2=0
         kappa_ux1=0
         kappa_uy1=0
         BdotT=0
      ENDIF
         
      fx_u(7,2,:,:)=-kappa_fac*BdotT_u2(1,:,:) 
     $     + kappa_u2*(BdotT(1,:,:) - ux(7,:,:))
      fx_u(7,3,:,:)=gamma_fac*u(7,:,:)
      fx_u(7,7,:,:)=gamma_fac*(u(3,:,:)-di*uy(2,:,:))

      fx_ux(7,1,:,:)=-kappa_fac*BdotT_ux1(1,:,:) 
     $     + kappa_ux1*(BdotT(1,:,:) - ux(7,:,:))
      fx_ux(7,7,:,:)=-kappa_fac*BdotT_ux7(1,:,:)
     $     - kappa_prp

      fx_uy(7,1,:,:)=-kappa_fac*BdotT_uy1(1,:,:)
     $     + kappa_uy1*(BdotT(1,:,:) - ux(7,:,:))
      fx_uy(7,2,:,:)=-di*gamma_fac*u(7,:,:)
      fx_uy(7,7,:,:)=-kappa_fac*BdotT_uy7(1,:,:)

      fy_u(7,2,:,:)=-kappa_fac*BdotT_u2(2,:,:)
     $     + kappa_u2*(BdotT(2,:,:) - uy(7,:,:))
      fy_u(7,4,:,:)=gamma_fac*u(7,:,:)
      fy_u(7,7,:,:)=gamma_fac*(u(4,:,:)+di*ux(2,:,:))

      fy_ux(7,1,:,:)=-kappa_fac*BdotT_ux1(2,:,:)
     $     + kappa_ux1*(BdotT(2,:,:) - uy(7,:,:))
      fy_ux(7,2,:,:)=di*gamma_fac*u(7,:,:)
      fy_ux(7,7,:,:)=-kappa_fac*BdotT_ux7(2,:,:)

      fy_uy(7,1,:,:)=-kappa_fac*BdotT_uy1(2,:,:)
     $     + kappa_uy1*(BdotT(2,:,:) - uy(7,:,:))
      fy_uy(7,7,:,:)=-kappa_fac*BdotT_uy7(2,:,:)
     $     - kappa_prp

      s_u(7,3,:,:)=ux(7,:,:)
      s_u(7,4,:,:)=uy(7,:,:)
      s_u(7,5,:,:)=2*eta_local*(u(5,:,:)-u(6,:,:)+u0(6,:,:))
     $     /di**2
      s_u(7,6,:,:)=-2*eta_local*(u(5,:,:)-u(6,:,:)+u0(6,:,:))
     $     /di**2
      s_u(7,7,:,:)=eta_u7
     $     *((u(5,:,:)-u(6,:,:)+u0(6,:,:))**2/di**2 
     $     + (ux(2,:,:)-u0x(2,:,:))**2 + (uy(2,:,:)-u0y(2,:,:))**2)

      s_ux(7,2,:,:)=di*uy(7,:,:) + 2*eta_local
     $     *(ux(2,:,:)-u0x(2,:,:))
      s_ux(7,6,:,:)=2*nu*(ux(6,:,:)-u0x(6,:,:))
      s_ux(7,7,:,:)=u(3,:,:) - di*uy(2,:,:)

      s_uy(7,2,:,:)=-di*ux(7,:,:) + 2*eta_local
     $     *(uy(2,:,:)-u0y(2,:,:))
      s_uy(7,6,:,:)=2*nu*(uy(6,:,:)-u0y(6,:,:))
      s_uy(7,7,:,:)=u(4,:,:) + di*ux(2,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_drdu
c-----------------------------------------------------------------------
c     subprogram 9. HallMHDicold_mass.
c     computes mass matrix couplings.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_mass(mass)

      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      mass(7,7,:,:)=one/(gamma-one)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_mass
c-----------------------------------------------------------------------
c     subprogram 10. HallMHDicold_equil.
c     computes equilibrium for magnetic reconnection.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_equil(x,y,u,ux,uy,deriv)
      
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
      CASE("Harris")
         u(1,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(3,:,:)=-mach*sinhy_phi/coshy_phi
         u(6,:,:)=-di*alfven/(lambda_psi*coshy_psi**2)
         u(7,:,:)=beta0
      CASE("GEM","GEM_open")
         u(1,:,:)=-lambda_psi*LOG(coshy_psi)
         u(2,:,:)=Bguide
         u(5,:,:)=-di/(lambda_psi*coshy_psi**2)
         u(7,:,:)=half
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=-TANH(y/lambda_psi)
         uy(5,:,:)=-u(5,:,:)*two*TANH(y/lambda_psi)/lambda_psi
      CASE("history","history_exp")
         u(1,:,:)=-lambda_psi*LOG(coshy_psi)
         u(2,:,:)=SQRT(one/coshy_psi**2 + Bguide**2)
         u(6,:,:)=di/(lambda_psi*coshy_psi**2)
         u(7,:,:)=beta0
         IF(.NOT. deriv)RETURN
         uy(1,:,:)=-TANH(y/lambda_psi)
         uy(2,:,:)=-TANH(y/lambda_psi)
     $        /(u(2,:,:)*lambda_psi*coshy_psi**2)
         uy(6,:,:)=-u(6,:,:)*two*TANH(y/lambda_psi)/lambda_psi
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_equil
c-----------------------------------------------------------------------
c     subprogram 11. HallMHDicold_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_grid(x,y,ksi,eta)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("Harris")
         ksi=lx*(x-half)
         eta=y-half
      CASE("GEM")
         ksi=lx*half*(x**4+4*gr_curve*x)/(one+4*gr_curve)
         eta=ly*half*(y**4+gr_curve*y)/(one+gr_curve)
      CASE("GEM_open")
         ksi=two*lx*x
         eta=two*ly*(y**2+gr_curve*y)/(one+gr_curve)
      CASE("history")
         IF(xper)THEN
            ksi=lx*((x-.5)**3 + 2*gr_curve*(x-.5))/(.25 + 2*gr_curve)
         ELSE
            ksi=lx*half*TANH(4.5*x-2.25)/TANH(2.25)
         ENDIF
         eta=ly*((y-.5)**3+gr_curve*(y-.5))/(.25+gr_curve)
      CASE("history_exp")
         ksi=lx*((x-.5)**3+2*gr_curve*(x-.5))/(.25+2*gr_curve)
         eta=ly*((y-.5)**3+gr_curve*(y-.5))/(.25+gr_curve)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_grid
c-----------------------------------------------------------------------
c     subprogram 12. HallMHDicold_BdotT.
c     computes derivatives of b(b*gradT) with respect to variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE HallMHDicold_BdotT(u,ux,uy,kappa_prp,kappa_u2,
     $     kappa_ux1,kappa_uy1,BdotT,BdotT_u2,BdotT_ux1,BdotT_ux7,
     $     BdotT_uy1,BdotT_uy7)
      
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: kappa_prp,kappa_u2,
     $     kappa_ux1,kappa_uy1
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: BdotT,BdotT_u2,
     $    BdotT_ux1,BdotT_ux7,BdotT_uy1,BdotT_uy7
      
      INTEGER :: ix,iy
      REAL(r8) :: Bsq
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: invBsq
c-----------------------------------------------------------------------
c     ion pressure equation.
c-----------------------------------------------------------------------
      invBsq=one/(ux(1,:,:)**2 + uy(1,:,:)**2 + u(2,:,:)**2)

      BdotT(1,:,:)=uy(1,:,:)*invBsq*
     $     (uy(1,:,:)*ux(7,:,:) - ux(1,:,:)*uy(7,:,:))
      BdotT(2,:,:)=-ux(1,:,:)*invBsq*
     $     (uy(1,:,:)*ux(7,:,:) - ux(1,:,:)*uy(7,:,:))
      BdotT_u2(1,:,:)=-2*u(2,:,:)*uy(1,:,:)*invBsq**2
     $     *(uy(1,:,:)*ux(7,:,:) - ux(1,:,:)*uy(7,:,:))
      BdotT_ux1(1,:,:)=uy(1,:,:)*invBsq**2
     $     *(ux(1,:,:)**2*uy(7,:,:) - 2*ux(1,:,:)*uy(1,:,:)*ux(7,:,:)
     $     - uy(1,:,:)**2*uy(7,:,:) - u(2,:,:)**2*uy(7,:,:))
      BdotT_ux7(1,:,:)=uy(1,:,:)**2*invBsq
      BdotT_uy1(1,:,:)=invBsq**2*(ux(1,:,:)
     $     *(uy(1,:,:)**2*uy(7,:,:) + 2*ux(1,:,:)*uy(1,:,:)*ux(7,:,:)
     $     - ux(1,:,:)**2*uy(7,:,:)) 
     $     + u(2,:,:)**2*(2*uy(1,:,:)*ux(7,:,:) - ux(1,:,:)*uy(7,:,:)))
      BdotT_uy7(1,:,:)=-ux(1,:,:)*uy(1,:,:)*invBsq
      
      BdotT_u2(2,:,:)=-2*u(2,:,:)*ux(1,:,:)*invBsq**2
     $     *(ux(1,:,:)*uy(7,:,:) - uy(1,:,:)*ux(7,:,:))
      BdotT_ux1(2,:,:)=invBsq**2*(uy(1,:,:)
     $     *(ux(1,:,:)**2*ux(7,:,:) + 2*ux(1,:,:)*uy(1,:,:)*uy(7,:,:)
     $     - uy(1,:,:)**2*ux(7,:,:)) 
     $     + u(2,:,:)**2*(2*ux(1,:,:)*uy(7,:,:) - uy(1,:,:)*ux(7,:,:)))
      BdotT_ux7(2,:,:)=-ux(1,:,:)*uy(1,:,:)*invBsq
      BdotT_uy1(2,:,:)=ux(1,:,:)*invBsq**2
     $     *(uy(1,:,:)**2*ux(7,:,:) - 2*ux(1,:,:)*uy(1,:,:)*uy(7,:,:)
     $     - ux(1,:,:)**2*ux(7,:,:) - u(2,:,:)**2*ux(7,:,:))
      BdotT_uy7(2,:,:)=ux(1,:,:)**2*invBsq

      DO ix=1,SIZE(u,2)
         DO iy=1,SIZE(u,3)
            Bsq=ux(1,ix,iy)**2+uy(1,ix,iy)**2+u(2,ix,iy)**2
            kappa_prp(ix,iy)=kappa_par*kappa_perp/
     $           (Bsq*kappa_par + kappa_perp)
            kappa_u2(ix,iy)=-2*kappa_par**2*kappa_perp*u(2,ix,iy)
     $           /(Bsq*kappa_par + kappa_perp)**2
            kappa_ux1(ix,iy)=-2*kappa_par**2*kappa_perp*ux(1,ix,iy)
     $           /(Bsq*kappa_par + kappa_perp)**2
            kappa_uy1(ix,iy)=-2*kappa_par**2*kappa_perp*uy(1,ix,iy)
     $           /(Bsq*kappa_par + kappa_perp)**2
            IF(Bsq <= Bmin)THEN
               BdotT(:,ix,iy)=zero
               BdotT_u2(:,ix,iy)=zero
               BdotT_ux1(:,ix,iy)=zero
               BdotT_ux7(:,ix,iy)=zero
               BdotT_uy1(:,ix,iy)=zero
               BdotT_uy7(:,ix,iy)=zero
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE HallMHDicold_BdotT
      END MODULE HallMHDicold_mod
