c-----------------------------------------------------------------------
c     file epmhd.f.
c     contains specifications for electron-positron MHD.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. epmhd_mod.
c     1. epmhd_input.
c     2. epmhd_init.
c     3. epmhd_init_special.
c     4. epmhd_boundary.
c     5. epmhd_bound_rhs.
c     6. epmhd_bound_drdu.
c     7. epmhd_rhs.
c     8. epmhd_drdu.
c     9. epmhd_mass.
c     10. epmhd_equil.
c     11. epmhd_grid.
c-----------------------------------------------------------------------
c     subprogram 0. epmhd_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE epmhd_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: source=.FALSE.
      CHARACTER(16), PRIVATE :: init_type=" "
      REAL(r8), PRIVATE :: eta=0,mu=0,de_sq=0,
     $     lambda_psi=0,lambda_phi=0,epsilon=0,tau=1,bound_eps=0,
     $     mach=0,alfven=1,lx=0,ly=0,kx=0,ky=0,ksq=0

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. epmhd_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $     dt,dtmax,tmax,nstep,gr_curve)

      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic
      INTEGER, INTENT(OUT) :: nqty
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax,gr_curve

      INTEGER :: myios
c-----------------------------------------------------------------------
c     epmhd namelist.
c-----------------------------------------------------------------------
      NAMELIST/epmhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,gr_curve,eta,mu,de_sq,tau,source,lx,ly,
     $     lambda_psi,lambda_phi,epsilon,bound_eps,mach,alfven,init_type
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=epmhd_list,IOSTAT=myios)

      IF(myios /= 0)CALL program_stop("job_input incorrect")

      nqty=4
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_input
c-----------------------------------------------------------------------
c     subprogram 2. epmhd_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_init(xpi,ypi,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: ix,nmodes=1
      REAL(r8), DIMENSION(SIZE(ui,1),SIZE(ui,2),SIZE(ui,3)) :: ux,uy
      REAL(r8), DIMENSION(0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) :: cosx,sinx,
     $     cosy,siny
c-----------------------------------------------------------------------
c     initial equilibrium.
c-----------------------------------------------------------------------
      CALL epmhd_equil(xpi,ypi,ui,ux,uy,.FALSE.)
      cosx=COS(kx*xpi)
      sinx=SIN(kx*xpi)
      cosy=COS(ky*ypi)
      siny=SIN(ky*ypi)
c-----------------------------------------------------------------------
c     add perturbation
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("islands")
         ui(1,:,:)=ui(1,:,:)+epsilon*cosx*cosy
      CASE("Harris")
         ui(1,:,:)=ui(1,:,:)-epsilon*cosx*cosy*alfven/ksq
         ui(2,:,:)=ui(2,:,:)+epsilon*cosx*cosy*mach
         ui(3,:,:)=ui(3,:,:)+epsilon*cosx*cosy*alfven
         ui(4,:,:)=ui(4,:,:)-epsilon*cosx*cosy*mach/ksq
      CASE("Harris_local")
         ui(1,:,:)=ui(1,:,:) + epsilon
     $        *EXP(-xpi**2/(2*lambda_psi)**2)
     $        *EXP(-ypi**2/(lambda_psi/2)**2)
ccc         DO ix=nmodes,nmodes
ccc            ui(1,:,:)=ui(1,:,:) - epsilon
ccc     $           *EXP(-ypi**2/(lambda_psi/2)**2)*COS(ix*kx*xpi)
ccc         ENDDO
      END SELECT
      RETURN
      END SUBROUTINE epmhd_init
c-----------------------------------------------------------------------
c     subprogram 3. epmhd_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_init_special(static,ground)

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      ground(4)=.TRUE.
      static(3:4)=.TRUE.
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
      CALL MPI_Bcast(de_sq,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_phi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bound_eps,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(tau,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alfven,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("islands")
         lambda_psi=one/twopi
         kx=pi
         ky=half*pi
         alfven=.2
      CASE("Harris")
         kx=twopi/lx
         ky=pi/ly
      CASE("Harris_local")
         kx=pi/lx
         ky=half*pi/ly
      CASE("Bhimsen","Bhimsen_mod")
         ky=twopi/lx
      END SELECT
      ksq=kx**2+ky**2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_init_special
c-----------------------------------------------------------------------
c     subprogram 4. epmhd_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_boundary(left,right,top,bottom,nqty)

      INTEGER, INTENT(IN) :: nqty
      TYPE(edge_type) :: left,right,top,bottom

      INTEGER :: i,j
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("islands")
         top%a=RESHAPE((/(1.,(0.,i=1,nqty),j=1,nqty-1),1./),
     $        (/nqty,nqty/))
         top%b=0
         left%static=.TRUE.
         right%static=.TRUE.
         bottom%static=.TRUE.
      CASE("Harris_open")
         left%static=.TRUE.
         right%a=0
         right%b(1,1)=1
         right%static(2:4)=.TRUE.
         top%a=0
         top%b(1,1)=1
         top%static(2:4)=.TRUE.
         bottom%static=.TRUE.
      CASE("Harris_local")
         left%static=.TRUE.
         right%static=.TRUE.
         top%b(1,1)=1
         top%a=0
         top%static(2:4)=.TRUE.
         bottom%static=.TRUE.
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_boundary
c-----------------------------------------------------------------------
c     subprogram 5. epmhd_bound_rhs.
c     computes rhs for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_bound_rhs(edge_type,t,x,y,u,ux,uy,uxx,uyy,c)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy,uxx,uyy
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: c

      REAL(r8) :: dksi_dt,ky
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
      SELECT CASE(init_type)
      CASE("islands")
         SELECT CASE(edge_type)
         CASE("bottom")
            c(1,:)=uy(1,:)
            c(2,:)=u(2,:)
            c(3,:)=uy(3,:)
            c(4,:)=u(4,:)
         CASE("left","right")
            c(1,:)=ux(1,:)
            c(2,:)=u(2,:)
            c(3,:)=ux(3,:)
            c(4,:)=u(4,:)
         END SELECT
      CASE("Harris")
         SELECT CASE(edge_type)
         CASE("top","bottom")
            c(2,:)=u(2,:)
            c(3,:)=u(3,:)
         END SELECT
      CASE("Harris_open")
         SELECT CASE(edge_type)
         CASE("top")
            c(2,:)=u(2,:)-uxx(4,:)
            c(3,:)=u(3,:)-uxx(1,:)-uyy(1,:)
            c(4,:)=uy(4,:)
         CASE("bottom")
            c(1,:)=uy(1,:)
            c(2,:)=u(2,:)
            c(3,:)=uy(3,:)
            c(4,:)=u(4,:)
         CASE("right")
            c(2,:)=u(2,:)-uyy(4,:)
            c(3,:)=u(3,:)-uxx(1,:)-uyy(1,:)
            c(4,:)=ux(4,:)
         CASE("left")
            c(1,:)=ux(1,:)
            c(2,:)=u(2,:)
            c(3,:)=ux(3,:)
            c(4,:)=u(4,:)
         END SELECT
      CASE("Harris_local")
         SELECT CASE(edge_type)
         CASE("top")
            c(2,:)=u(2,:)
            c(3,:)=uyy(3,:)
            c(4,:)=uy(4,:)
         CASE("bottom")
            c(1,:)=uy(1,:)
            c(2,:)=u(2,:)
            c(3,:)=uy(3,:)
            c(4,:)=u(4,:)
         CASE("left","right")
            c(1,:)=ux(1,:)
            c(2,:)=u(2,:)
            c(3,:)=ux(3,:)
            c(4,:)=u(4,:)
         END SELECT
      CASE("Bhimsen")
         ky=twopi/lx
         dksi_dt=bound_eps/tau**2*EXP(-t/tau)
         SELECT CASE(edge_type)
         CASE("left")
            c(1,:)=-dksi_dt*t*SIN(ky*y)
            c(4,:)=-one/ky*dksi_dt*(one-t/tau)*COS(ky*y)
         CASE("right")
            c(1,:)=-dksi_dt*t*SIN(ky*y)
            c(4,:)=one/ky*dksi_dt*(one-t/tau)*COS(ky*y)
         END SELECT
         c(2,:)=u(2,:)
         c(3,:)=u(3,:)
      CASE("Bhimsen_mod")
         ky=twopi/lx
         dksi_dt=bound_eps*EXP(-t/tau)*t/(two*tau**3)
         SELECT CASE(edge_type)
         CASE("left")
            c(1,:)=-dksi_dt*t*SIN(ky*y)
            c(4,:)=-one/ky*dksi_dt*(two-t/tau)*COS(ky*y)
         CASE("right")
            c(1,:)=-dksi_dt*t*SIN(ky*y)
            c(4,:)=one/ky*dksi_dt*(two-t/tau)*COS(ky*y)
         END SELECT
         c(2,:)=u(2,:)
         c(3,:)=u(3,:)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_bound_rhs
c-----------------------------------------------------------------------
c     subprogram 6. epmhd_bound_drdu.
c     computes drdu for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_bound_drdu(edge_type,t,x,y,u,ux,uy,uxx,uyy,
     $     c_u,c_ux,c_uy,c_uxx,c_uyy)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), INTENT(IN) :: t
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
      SELECT CASE(init_type)
      CASE("islands")
         SELECT CASE(edge_type)
         CASE("bottom")
            c_uy(1,1,:)=one
            c_u(2,2,:)=one
            c_uy(3,3,:)=one
            c_u(4,4,:)=one
         CASE("left","right")
            c_ux(1,1,:)=one
            c_u(2,2,:)=one
            c_ux(3,3,:)=one
            c_u(4,4,:)=one
         END SELECT
      CASE("Harris")
         SELECT CASE(edge_type)
         CASE("top","bottom")
            c_u(2,2,:)=one
            c_u(3,3,:)=one
         END SELECT
      CASE("Harris_open")
         SELECT CASE(edge_type)
         CASE("top")
            c_u(2,2,:)=one
            c_uxx(2,4,:)=-one
            c_u(3,3,:)=one
            c_uxx(3,1,:)=-one
            c_uyy(3,1,:)=-one
            c_uy(4,4,:)=one
         CASE("bottom")
            c_uy(1,1,:)=one
            c_u(2,2,:)=one
            c_uy(3,3,:)=one
            c_u(4,4,:)=one
         CASE("right")
            c_u(2,2,:)=one
            c_uyy(2,4,:)=-one
            c_u(3,3,:)=one
            c_uxx(3,1,:)=-one
            c_uyy(3,1,:)=-one
            c_ux(4,4,:)=one
         CASE("left")
            c_ux(1,1,:)=one
            c_u(2,2,:)=one
            c_ux(3,3,:)=one
            c_u(4,4,:)=one
         END SELECT
      CASE("Harris_local")
         SELECT CASE(edge_type)
         CASE("top")
            c_u(2,2,:)=one
            c_uyy(3,3,:)=one
            c_uy(4,4,:)=one
         CASE("bottom")
            c_uy(1,1,:)=one
            c_u(2,2,:)=one
            c_uy(3,3,:)=one
            c_u(4,4,:)=one
         CASE("left","right")
            c_ux(1,1,:)=one
            c_u(2,2,:)=one
            c_ux(3,3,:)=one
            c_u(4,4,:)=one
         END SELECT
      CASE("Bhimsen","Bhimsen_mod")
         c_u(2,2,:)=one
         c_u(3,3,:)=one
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_bound_drdu
c-----------------------------------------------------------------------
c     subprogram 7. epmhd_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE epmhd_rhs(t,x,y,u,ux,uy,fx,fy,s,first)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8) :: d_eta,d2
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: eta_local 
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0y,u0x,
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
      eta_local=eta
      IF(init_type=="Harris_open")THEN
         d_eta=10*eta
         d2=10
         eta_local=eta + d_eta
     $        *(EXP(-(d2*(x-lx)/lx)**4) + EXP(-(d2*(x+lx)/lx)**4))
      ENDIF
c-----------------------------------------------------------------------
c     electron-positron MHD
c-----------------------------------------------------------------------
      fx(1,:,:)=-uy(4,:,:)*(u(1,:,:) - half*de_sq*u(3,:,:)) 
     $     + half*de_sq*mu*ux(3,:,:)
      
      fy(1,:,:)=ux(4,:,:)*(u(1,:,:) - half*de_sq*u(3,:,:))
     $     + half*de_sq*mu*uy(3,:,:)

      s(1,:,:)=eta_local*u(3,:,:)

      fx(2,:,:)=-uy(4,:,:)*u(2,:,:) + half*uy(1,:,:)*u(3,:,:)
     $     -mu*ux(2,:,:)
      
      fy(2,:,:)=ux(4,:,:)*u(2,:,:) - half*ux(1,:,:)*u(3,:,:)
     $     -mu*uy(2,:,:)
      
      fx(3,:,:)=-ux(1,:,:)
      fy(3,:,:)=-uy(1,:,:)
      s(3,:,:)=-u(3,:,:)
      
      fx(4,:,:)=-ux(4,:,:)
      fy(4,:,:)=-uy(4,:,:)
      s(4,:,:)=-u(2,:,:)

      IF(source .AND. first)THEN
         CALL epmhd_equil(x,y,u0,u0x,u0y,.TRUE.)
         first=.FALSE.
         CALL epmhd_rhs(t,x,y,u0,u0x,u0y,fx0,fy0,s0,first)
         fx=fx-fx0
         fy=fy-fy0
         s=s-s0
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_rhs
c-----------------------------------------------------------------------
c     subprogram 8. epmhd_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_drdu(t,x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8) :: d_eta,d2
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: eta_local 
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
      eta_local=eta
      IF(init_type=="Harris_open")THEN
         d_eta=10*eta
         d2=10
         eta_local=eta + d_eta
     $        *(EXP(-(d2*(x-lx)/lx)**4) + EXP(-(d2*(x+lx)/lx)**4))
      ENDIF
c-----------------------------------------------------------------------
c     reduced MHD + electron inertia + sound radius model.
c-----------------------------------------------------------------------
      fx_u(1,1,:,:)=-uy(4,:,:)
      fx_u(1,3,:,:)=half*de_sq*uy(4,:,:)
      fx_ux(1,3,:,:)=half*de_sq*mu
      fx_uy(1,4,:,:)=-(u(1,:,:) - half*de_sq*u(3,:,:))

      fy_u(1,1,:,:)=ux(4,:,:)
      fy_u(1,3,:,:)=-half*de_sq*ux(4,:,:)
      fy_ux(1,4,:,:)=u(1,:,:) - half*de_sq*u(3,:,:)
      fy_uy(1,3,:,:)=half*de_sq*mu

      s_u(1,3,:,:)=eta_local
      
      fx_u(2,2,:,:)=-uy(4,:,:)
      fx_u(2,3,:,:)=half*uy(1,:,:)
      fx_ux(2,2,:,:)=-mu
      fx_uy(2,1,:,:)=half*u(3,:,:)
      fx_uy(2,4,:,:)=-u(2,:,:)

      fy_u(2,2,:,:)=ux(4,:,:)
      fy_u(2,3,:,:)=-half*ux(1,:,:)
      fy_ux(2,1,:,:)=-half*u(3,:,:)
      fy_ux(2,4,:,:)=u(2,:,:)
      fy_uy(2,2,:,:)=-mu

      fx_ux(3,1,:,:)=-one
      fy_uy(3,1,:,:)=-one
      s_u(3,3,:,:)=-one
      
      fx_ux(4,4,:,:)=-one
      fy_uy(4,4,:,:)=-one
      s_u(4,2,:,:)=-one
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_drdu
c-----------------------------------------------------------------------
c     subprogram 9. epmhd_mass.
c     computes mass matrix couplings.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_mass(mass)

      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      mass(1,3,:,:)=-half*de_sq
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_mass
c-----------------------------------------------------------------------
c     subprogram 10. epmhd_equil.
c     computes equilibrium for magnetic reconnection.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_equil(x,y,u,ux,uy,derivs)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: coshy_psi,coshy_phi,
     $     sinhy_psi,sinhy_phi,coshx_psi,coshx_phi,sinhx_psi,cosx_psi,
     $     sinx_psi
c-----------------------------------------------------------------------
c     compute cosh.
c-----------------------------------------------------------------------
      coshy_psi=COSH(y/lambda_psi)
      coshy_phi=COSH(y/lambda_phi)
      coshx_psi=COSH(x/lambda_psi)
      coshx_phi=COSH(x/lambda_phi)
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      SELECT CASE(init_type)
      CASE("Harris")
         u(1,:,:)=alfven*lambda_psi*log(coshy_psi)
         u(2,:,:)=mach/(lambda_phi*coshy_phi**2)
         u(3,:,:)=alfven/(lambda_psi*coshy_psi**2)
         u(4,:,:)=mach*lambda_phi*log(coshy_phi)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhy_psi=SINH(y/lambda_psi)
         sinhy_phi=SINH(y/lambda_phi)
         uy(1,:,:)=alfven*sinhy_psi/coshy_psi
         uy(2,:,:)=-2*mach/lambda_phi**2*sinhy_phi/coshy_phi**3
         uy(3,:,:)=-2*alfven/lambda_psi**2*sinhy_psi/coshy_psi**3
         uy(4,:,:)=mach*sinhy_phi/coshy_phi
      CASE("Harris_open")
         coshx_psi=COSH(x/(lambda_psi*lx/ly))
         u(1,:,:)=lambda_psi*(log(coshy_psi)-alfven*log(coshx_psi))
         u(3,:,:)=one/lambda_psi*
     $        (one/coshy_psi**2 - alfven/(lx/ly*coshx_psi)**2)
c-----------------------------------------------------------------------
c           compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhx_psi=SINH(x/(lambda_psi*lx/ly))
         sinhy_psi=SINH(y/lambda_psi)
         ux(1,:,:)=-alfven*sinhx_psi/(lx/ly*coshx_psi)
         ux(3,:,:)=2*alfven/lambda_psi**2*sinhx_psi/(lx/ly*coshx_psi)**3
         uy(1,:,:)=sinhy_psi/coshy_psi
         uy(3,:,:)=-2/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE("Harris_local")
         u(1,:,:)=lambda_psi*log(coshy_psi)
         u(3,:,:)=one/(lambda_psi*coshy_psi**2)
c-----------------------------------------------------------------------
c           compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinhy_psi=SINH(y/lambda_psi)
         uy(1,:,:)=sinhy_psi/coshy_psi
         uy(3,:,:)=-2/lambda_psi**2*sinhy_psi/coshy_psi**3
      CASE("islands")
         cosx_psi=COS(x/lambda_psi)
         u(1,:,:)=-lambda_psi*LOG(coshy_psi+alfven*cosx_psi)
         u(3,:,:)=(alfven**2-one)/lambda_psi
     $        /(coshy_psi+alfven*cosx_psi)**2
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         sinx_psi=SIN(x/lambda_psi)
         sinhy_psi=SINH(y/lambda_psi)
         ux(1,:,:)=alfven*sinx_psi/(coshy_psi+alfven*cosx_psi)
         uy(1,:,:)=-sinhy_psi/(coshy_psi+alfven*cosx_psi)
         ux(3,:,:)=u(3,:,:)*two*alfven/lambda_psi*sinx_psi
     $        /(coshy_psi+alfven*cosx_psi)
         uy(3,:,:)=-u(3,:,:)*two/lambda_psi*sinhy_psi
     $        /(coshy_psi+alfven*cosx_psi)
      CASE("Bhimsen","Bhimsen_mod")
         u(1,:,:)=x**2*half
         u(3,:,:)=one
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         IF(.NOT. derivs)RETURN
         ux(1,:,:)=x
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_equil
c-----------------------------------------------------------------------
c     subprogram 11. epmhd_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE epmhd_grid(x,y,ksi,eta)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("islands")
         ksi=(x**2+.2*x)/1.2
         eta=(y**2+.5*y)/1.5
      CASE("Harris")
         ksi=lx*(x-.5)
         eta=ly*(y-.5)*((y-.5)**2+.125)/.375
      CASE("Harris_open")
         ksi=lx*x
ccc         ksi=lx*(TANH(3*x-1.5)+TANH(1.5))/(2*TANH(1.5))
         eta=ly*(TANH(4*y-2.)+TANH(2.))/(2*TANH(2.))
ccc         eta=ly*((y-.5)**3+gr_curve*(y-.5))/(.25+gr_curve)
      CASE("Harris_local")
         ksi=lx*x*(gr_curve*x**3 + 1/gr_curve)/(gr_curve + 1/gr_curve)
         eta=ly*(TANH(gr_curve*(y-1))+TANH(gr_curve))/TANH(gr_curve)
      CASE("Bhimsen","Bhimsen_mod")
         ksi=two*(x-half)
         eta=lx*y
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE epmhd_grid
      END MODULE epmhd_mod
