c-----------------------------------------------------------------------
c     file beltrami.f.
c     Beltrami's equation for adaptive grid generation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. beltrami_mod.
c     1. beltrami_init.
c     2. beltrami_boundary.
c     3. beltrami_edge_rhs.
c     4. beltrami_edge_drdu.
c     5. beltrami_rhs.
c     6. beltrami_drdu.
c     7. beltrami_phi.
c     8. beltrami_griderr.
c     9. beltrami_gmat1.
c     10. beltrami_gmat2.
c     11. beltrami_gmat_fit.
c     12. beltrami_getdata.
c     13. beltrami_dealloc.
c     14. beltrami_diagnose.
c-----------------------------------------------------------------------
c     subprogram 0. beltrami_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE beltrami_mod
      USE jacobi_mod
      USE bicube_mod
      IMPLICIT NONE

      TYPE :: beltr_type
      LOGICAL :: xperiodic,yperiodic
      INTEGER :: nx,ny,np
      REAL(r8), DIMENSION(:,:), POINTER :: mass
      REAL(r8), DIMENSION(:,:,:), POINTER :: gfit
      TYPE(jacobi_type) :: basis
      TYPE(bicube_type) :: griderr
      END TYPE beltr_type

      LOGICAL :: bel_diagnose=.FALSE.
      LOGICAL, DIMENSION(:), POINTER :: adapt_qty 
      REAL(r8) :: bel_phifac=1.
      TYPE(beltr_type) :: beltr

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. beltrami_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_init(x,y,u)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: u
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      u(1,:,:)=x
      u(2,:,:)=y
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_init
c-----------------------------------------------------------------------
c     subprogram 2. beltrami_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_boundary(left,right,top,bottom)

      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     set periodicity and static flags.
c-----------------------------------------------------------------------
      IF(beltr%xperiodic)THEN
         left%bc_type(1)="robin"
         right%bc_type(1)="robin"
         left%static(1)=.TRUE.
         right%static(1)=.TRUE.
         left%bc_type(2)="periodic"
         right%bc_type(2)="periodic"
         top%bc_type="robin"
         bottom%bc_type="robin"
         top%static=.TRUE.
         bottom%static=.TRUE.
      ELSEIF(beltr%yperiodic)THEN
         left%bc_type="robin"
         right%bc_type="robin"
         left%static=.TRUE.
         right%static=.TRUE.
         top%bc_type(1)="periodic"
         bottom%bc_type(1)="periodic"
         top%bc_type(2)="robin"
         bottom%bc_type(2)="robin"
         top%static(2)=.TRUE.
         bottom%static(2)=.TRUE.
      ELSE
         left%bc_type="robin"
         right%bc_type="robin"
         top%bc_type="robin"
         bottom%bc_type="robin"
         left%static=.TRUE.
         right%static=.TRUE.
         top%static=.TRUE.
         bottom%static=.TRUE.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_boundary
c-----------------------------------------------------------------------
c     subprogram 3. beltrami_edge_rhs.
c     computes rhs for boundary terms.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_edge_rhs(lrtb,x,y,u,ux,uy,c)

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     specify rhs.
c-----------------------------------------------------------------------
      c=0
      SELECT CASE(lrtb)
      CASE("left")
         c(1,:,:)=u(1,:,:)
         c(2,:,:)=ux(2,:,:)
      CASE("right")
         c(1,:,:)=u(1,:,:)-one
         c(2,:,:)=ux(2,:,:)
      CASE("top")
         c(1,:,:)=uy(1,:,:)
         c(2,:,:)=u(2,:,:)-one
      CASE("bottom")
         c(1,:,:)=uy(1,:,:)
         c(2,:,:)=u(2,:,:)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_edge_rhs
c-----------------------------------------------------------------------
c     subprogram 4. beltrami_edge_drdu.
c     computes drdu for boundary terms.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_edge_drdu(lrtb,x,y,u,ux,uy,
     $     c_u,c_ux,c_uy)

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy
c-----------------------------------------------------------------------
c     specify drdu.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      SELECT CASE(lrtb)
      CASE("left","right")
         c_u(1,1,:,:)=one
         c_ux(2,2,:,:)=one
      CASE("top","bottom")
         c_uy(1,1,:,:)=one
         c_u(2,2,:,:)=one
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_edge_drdu
c-----------------------------------------------------------------------
c     subprogram 5. beltrami_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_rhs(x,y,u,ux,uy,fx,fy,s)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s

      REAL(r8), DIMENSION(4,SIZE(x,1),SIZE(x,2)) :: gmat,gxmat,gymat
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: jac,fac,gdet
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      fac = 0
      jac = ux(1,:,:)*uy(2,:,:)-ux(2,:,:)*uy(1,:,:)
      CALL beltrami_gmat2(u(1,:,:),u(2,:,:),gmat,gxmat,gymat)
      gdet = ABS(gmat(4,:,:))
      WHERE(gdet == 0 .OR. jac == 0)
         fac = HUGE(fac)
      ELSEWHERE
         fac = one/(SQRT(gdet)*jac)
      END WHERE
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     fluxes.
c-----------------------------------------------------------------------
      fx(1,:,:) = fac*(gmat(1,:,:)*uy(2,:,:)**2
     $     + gmat(2,:,:)*uy(1,:,:)**2
     $     - two*gmat(3,:,:)*uy(1,:,:)*uy(2,:,:))
      fy(1,:,:) = -fac*(gmat(1,:,:)*ux(2,:,:)*uy(2,:,:)
     $     + gmat(2,:,:)*ux(1,:,:)*uy(1,:,:)
     $     - gmat(3,:,:)*(ux(1,:,:)*uy(2,:,:)+uy(1,:,:)*ux(2,:,:)))
      fx(2,:,:) = fy(1,:,:)
      fy(2,:,:) = fac*(gmat(1,:,:)*ux(2,:,:)**2
     $     + gmat(2,:,:)*ux(1,:,:)**2
     $     - two*gmat(3,:,:)*ux(1,:,:)*ux(2,:,:))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_rhs
c-----------------------------------------------------------------------
c     subprogram 6. beltrami_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_drdu(x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8), DIMENSION(4,SIZE(x,1),SIZE(x,2)) :: gmat,gxmat,gymat
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gdet,
     $     jac,fac,facx,facy,invjac
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: fx,fy
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      jac = ux(1,:,:)*uy(2,:,:)-ux(2,:,:)*uy(1,:,:)
      CALL beltrami_gmat2(u(1,:,:),u(2,:,:),gmat,gxmat,gymat)

      WHERE(jac == 0)
         invjac = HUGE(jac)
      ELSEWHERE
         invjac = one/jac
      ENDWHERE

      gdet = ABS(gmat(4,:,:))
      WHERE(gdet == 0)
         fac=HUGE(gdet)
         facx=HUGE(gdet)
         facy=HUGE(gdet)
      ELSEWHERE
         fac=invjac/SQRT(gdet)
         facx=-gxmat(4,:,:)/(two*gdet)
         facy=-gymat(4,:,:)/(two*gdet)
      END WHERE

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
c     fluxes.
c-----------------------------------------------------------------------
      fx(1,:,:) = fac*(gmat(1,:,:)*uy(2,:,:)**2
     $     + gmat(2,:,:)*uy(1,:,:)**2
     $     - two*gmat(3,:,:)*uy(1,:,:)*uy(2,:,:))
      fy(1,:,:) = -fac*(gmat(1,:,:)*ux(2,:,:)*uy(2,:,:)
     $     + gmat(2,:,:)*ux(1,:,:)*uy(1,:,:)
     $     - gmat(3,:,:)*(ux(1,:,:)*uy(2,:,:)+uy(1,:,:)*ux(2,:,:)))
      fx(2,:,:) = fy(1,:,:)
      fy(2,:,:) = fac*(gmat(1,:,:)*ux(2,:,:)**2
     $     + gmat(2,:,:)*ux(1,:,:)**2
     $     - two*gmat(3,:,:)*ux(1,:,:)*ux(2,:,:))
c-----------------------------------------------------------------------
c     derivatives of fx1.
c-----------------------------------------------------------------------
      fx_u(1,1,:,:) = fx(1,:,:)*facx 
     $     + fac*(gxmat(1,:,:)*uy(2,:,:)**2 + gxmat(2,:,:)*uy(1,:,:)**2
     $     - two*gxmat(3,:,:)*uy(1,:,:)*uy(2,:,:))
      fx_u(1,2,:,:) = fx(1,:,:)*facy
     $     + fac*(gymat(1,:,:)*uy(2,:,:)**2 + gymat(2,:,:)*uy(1,:,:)**2
     $     - two*gymat(3,:,:)*uy(1,:,:)*uy(2,:,:))
      fx_ux(1,1,:,:) = -fx(1,:,:)*uy(2,:,:)*invjac
      fx_ux(1,2,:,:) = fx(1,:,:)*uy(1,:,:)*invjac
      fx_uy(1,1,:,:) = fx(1,:,:)*ux(2,:,:)*invjac
     $     + two*fac*(gmat(2,:,:)*uy(1,:,:) - gmat(3,:,:)*uy(2,:,:))
      fx_uy(1,2,:,:) = -fx(1,:,:)*ux(1,:,:)*invjac
     $     + two*fac*(gmat(1,:,:)*uy(2,:,:) - gmat(3,:,:)*uy(1,:,:))
c-----------------------------------------------------------------------
c     derivatives of fy2.
c-----------------------------------------------------------------------
      fy_u(2,1,:,:) = fy(2,:,:)*facx
     $     + fac*(gxmat(1,:,:)*ux(2,:,:)**2 + gxmat(2,:,:)*ux(1,:,:)**2
     $     - two*gxmat(3,:,:)*ux(1,:,:)*ux(2,:,:))
      fy_u(2,2,:,:) = fy(2,:,:)*facy
     $     + fac*(gymat(1,:,:)*ux(2,:,:)**2 + gymat(2,:,:)*ux(1,:,:)**2
     $     - two*gymat(3,:,:)*ux(1,:,:)*ux(2,:,:))
      fy_ux(2,1,:,:) = -fy(2,:,:)*uy(2,:,:)*invjac
     $     + two*fac*(gmat(2,:,:)*ux(1,:,:) - gmat(3,:,:)*ux(2,:,:))
      fy_ux(2,2,:,:) = fy(2,:,:)*uy(1,:,:)*invjac
     $     + two*fac*(gmat(1,:,:)*ux(2,:,:) - gmat(3,:,:)*ux(1,:,:))
      fy_uy(2,1,:,:) = fy(2,:,:)*ux(2,:,:)*invjac
      fy_uy(2,2,:,:) = -fy(2,:,:)*ux(1,:,:)*invjac
c-----------------------------------------------------------------------
c     derivatives of fx2.
c-----------------------------------------------------------------------
      fx_u(2,1,:,:) = fx(2,:,:)*facx
     $     - fac*(gxmat(1,:,:)*ux(2,:,:)*uy(2,:,:)
     $     + gxmat(2,:,:)*ux(1,:,:)*uy(1,:,:)
     $     - gxmat(3,:,:)*(ux(1,:,:)*uy(2,:,:)+uy(1,:,:)*ux(2,:,:)))
      fx_u(2,2,:,:) = fx(2,:,:)*facy
     $     - fac*(gymat(1,:,:)*ux(2,:,:)*uy(2,:,:)
     $     + gymat(2,:,:)*ux(1,:,:)*uy(1,:,:)
     $     - gymat(3,:,:)*(ux(1,:,:)*uy(2,:,:)+uy(1,:,:)*ux(2,:,:)))
      fx_ux(2,1,:,:) = -fx(2,:,:)*uy(2,:,:)*invjac
     $     - fac*(gmat(2,:,:)*uy(1,:,:) - gmat(3,:,:)*uy(2,:,:))
      fx_uy(2,1,:,:) = fx(2,:,:)*ux(2,:,:)*invjac
     $     - fac*(gmat(2,:,:)*ux(1,:,:) - gmat(3,:,:)*ux(2,:,:))
      fx_ux(2,2,:,:) = fx(2,:,:)*uy(1,:,:)*invjac
     $     - fac*(gmat(1,:,:)*uy(2,:,:) - gmat(3,:,:)*uy(1,:,:))
      fx_uy(2,2,:,:) = -fx(2,:,:)*ux(1,:,:)*invjac
     $     - fac*(gmat(1,:,:)*ux(2,:,:) - gmat(3,:,:)*ux(1,:,:))
c-----------------------------------------------------------------------
c     derivatives of fy1.
c-----------------------------------------------------------------------
      fy_u(1,1,:,:) = fx_u(2,1,:,:)
      fy_u(1,2,:,:) = fx_u(2,2,:,:)
      fy_ux(1,1,:,:) = fx_ux(2,1,:,:)
      fy_uy(1,1,:,:) = fx_uy(2,1,:,:)
      fy_ux(1,2,:,:) = fx_ux(2,2,:,:)
      fy_uy(1,2,:,:) = fx_uy(2,2,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_drdu
c-----------------------------------------------------------------------
c     subprogram 7. beltrami_phi.
c     computes monitor function phi and its derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_phi(x,y,phi,phix,phiy)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: phi,phix,phiy

      INTEGER :: ix,iy
c-----------------------------------------------------------------------
c     compute phi based on griderr.
c-----------------------------------------------------------------------
      DO ix=0,SIZE(x,1)-1
         DO iy=0,SIZE(x,2)-1
            CALL bicube_eval(beltr%griderr,x(ix,iy),y(ix,iy))
            phi(ix,iy)=beltr%griderr%f(3)
            phix(ix,iy)=beltr%griderr%f(1)
            phiy(ix,iy)=beltr%griderr%f(2)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_phi
c-----------------------------------------------------------------------
c     subprogram 8. beltrami_griderr.
c     computes spatial truncation errors.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_griderr(alloc,uu,tinv)

      LOGICAL, INTENT(IN) :: alloc
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(r8), DIMENSION(0:,0:), INTENT(IN), OPTIONAL :: tinv

      INTEGER :: ix,jx,iy,jy,nx,ny,np,iqty,nqty
      REAL(r8) :: maxerr,maxnorm
      REAL(r8), DIMENSION(0:beltr%np,0:beltr%np) :: u,delu,mass
      REAL(r8), DIMENSION(beltr%nx,beltr%ny) :: norm,errx_max,erry_max,
     $     errxy_max
      REAL(r8), DIMENSION(SIZE(uu,1),beltr%nx,beltr%ny) :: errx,erry
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      nx=beltr%nx
      ny=beltr%ny
      np=beltr%np
      nqty=SIZE(uu,1)
      mass=beltr%mass
c-----------------------------------------------------------------------
c     compute norms and errors.
c-----------------------------------------------------------------------
      errx=0
      erry=0
      DO iqty=1,nqty
         IF(.NOT. adapt_qty(iqty))CYCLE
         jy=0
         DO iy=1,ny
            jx=0
            DO ix=1,nx
               u=uu(iqty,jx:jx+np,jy:jy+np)
               IF(PRESENT(tinv))u=MATMUL(MATMUL(tinv,u),TRANSPOSE(tinv))
               norm(ix,iy)=SUM(u*MATMUL(MATMUL(mass,u),mass))
               delu=0
               delu(np-1,:)=u(np-1,:)
               errx(iqty,ix,iy)
     $              =SUM(delu*MATMUL(MATMUL(mass,delu),mass))
               delu=0
               delu(:,np-1)=u(:,np-1)
               erry(iqty,ix,iy)
     $              =SUM(delu*MATMUL(MATMUL(mass,delu),mass))
               jx=jx+np+1
            ENDDO
            jy=jy+np+1
         ENDDO
c-----------------------------------------------------------------------
c     eliminate zero norms.
c-----------------------------------------------------------------------
         maxnorm=MAXVAL(norm)
         WHERE(norm <= maxnorm*min_eps**.25)
            norm = one/min_eps
         END WHERE
         WHERE(norm <= min_eps)
            norm = one/min_eps
         END WHERE
c-----------------------------------------------------------------------
c     normalize.
c-----------------------------------------------------------------------
         errx(iqty,:,:)=SQRT(ABS(errx(iqty,:,:)/norm))
         erry(iqty,:,:)=SQRT(ABS(erry(iqty,:,:)/norm))
      ENDDO
c-----------------------------------------------------------------------
c     compute maxima.
c-----------------------------------------------------------------------
      errx_max=(MAXVAL(errx,1))**0.25
      erry_max=(MAXVAL(erry,1))**0.25
      errxy_max=MAX(errx_max,erry_max)
c-----------------------------------------------------------------------
c     rescale and fit to bicubic splines.
c-----------------------------------------------------------------------
      IF(alloc)CALL bicube_alloc(beltr%griderr,nx,ny,3)
      beltr%griderr%xs=(/(ix,ix=0,nx)/)/REAL(nx,r8)
      beltr%griderr%ys=(/(iy,iy=0,ny)/)/REAL(ny,r8)
      beltr%griderr%f=0
      beltr%griderr%fs=0
      beltr%griderr%fs(1:nx,1:ny,1)=errx_max
      beltr%griderr%fs(1:nx,1:ny,2)=erry_max
      beltr%griderr%fs(1:nx,1:ny,3)=errxy_max

      maxerr=MAXVAL(beltr%griderr%fs(1:nx,1:ny,3))
      beltr%griderr%fs(1:nx,1:ny,:)=one + bel_phifac/maxerr
     $     *beltr%griderr%fs(1:nx,1:ny,:)

      IF(beltr%xperiodic)THEN
         CALL bicube_lsfit_xp(beltr%griderr)
      ELSEIF(beltr%yperiodic)THEN
         CALL bicube_lsfit_yp(beltr%griderr)
      ELSE
         CALL bicube_lsfit(beltr%griderr)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_griderr
c-----------------------------------------------------------------------
c     subprogram 9. beltrami_gmat1.
c     computes metric tensor.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_gmat1(x,y,gmat)

      REAL(r8), DIMENSION(0:beltr%np,0:beltr%np), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(4,0:beltr%np,0:beltr%np), INTENT(OUT) :: 
     $     gmat

      REAL(r8), DIMENSION(0:beltr%np,0:beltr%np) :: psi,psi_ksi,psi_eta
c-----------------------------------------------------------------------
c     compute gmat and gdet, align or both.
c-----------------------------------------------------------------------
      gmat = 0
      
      CALL beltrami_phi(x,y,psi,psi_ksi,psi_eta)
      gmat(1,:,:)=psi_eta/psi
      gmat(2,:,:)=psi_ksi/psi
      gmat(4,:,:)=psi_ksi*psi_eta
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_gmat1
c-----------------------------------------------------------------------
c     subprogram 10. beltrami_gmat2.
c     computes metric tensor.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_gmat2(x,y,gmat,gxmat,gymat)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: gmat,gxmat,gymat

      INTEGER :: ix,iy,jx,jy,inx,iny,np,nx,ny
      REAL(r8), DIMENSION(0:beltr%np) :: px,qx,py,qy
      REAL(r8), DIMENSION(:,:,:), POINTER :: gfit
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      np=beltr%np
      nx=beltr%nx
      ny=beltr%ny

      gfit => beltr%gfit
      gmat=0
      gxmat=0
      gymat=0
c-----------------------------------------------------------------------
c     compute output.
c-----------------------------------------------------------------------
      DO ix=1,SIZE(x,1)
         DO iy=1,SIZE(x,2)

            CALL jacobi_interp(nx,ny,x(ix,iy),y(ix,iy),jx,jy,
     $           beltr%basis,px,py,qx,qy)
            jx=(jx-1)*(np+1)
            jy=(jy-1)*(np+1)
c-----------------------------------------------------------------------
c     interpolate functions and derivates.
c-----------------------------------------------------------------------
            DO inx=0,np
               DO iny=0,np
                  gmat(:,ix,iy) = gmat(:,ix,iy) 
     $                 + gfit(:,jx+inx,jy+iny)*px(inx)*py(iny)
                  gxmat(:,ix,iy) = gxmat(:,ix,iy) 
     $                 + gfit(:,jx+inx,jy+iny)*qx(inx)*py(iny)
                  gymat(:,ix,iy) = gymat(:,ix,iy) 
     $                 + gfit(:,jx+inx,jy+iny)*px(inx)*qy(iny)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_gmat2
c-----------------------------------------------------------------------
c     subprogram 11. beltrami_gmat_fit.
c     initializes solution using p2_grid module subroutine.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_gmat_fit(gfit,alloc)
      
      LOGICAL, INTENT(IN) :: alloc
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: gfit
c-----------------------------------------------------------------------
c     pass gfit to local object "beltr".
c-----------------------------------------------------------------------
      IF(alloc)ALLOCATE(beltr%gfit(SIZE(gfit,1),0:SIZE(gfit,2)-1,
     $     0:SIZE(gfit,3)-1))
      beltr%gfit=0
      beltr%gfit=gfit 
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_gmat_fit
c-----------------------------------------------------------------------
c     subprogram 12. beltrami_getdata.
c     input grid data into object "beltr" from the main code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_getdata(xperiodic,yperiodic,nx,ny,np,
     $     mass,alloc)

      LOGICAL, INTENT(IN) :: xperiodic,yperiodic,alloc
      INTEGER, INTENT(IN) :: nx,ny,np
      REAL(r8), DIMENSION(0:np,0:np), INTENT(IN) :: mass
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      IF(alloc)THEN
         ALLOCATE(beltr%mass(0:np,0:np))
         CALL jacobi_alloc(beltr%basis,np,.FALSE.,.FALSE.,"gll")
      ENDIF
      beltr%mass=0
c-----------------------------------------------------------------------
c     pass variables to scalar fields of objects beltr.
c-----------------------------------------------------------------------
      beltr%nx=nx
      beltr%ny=ny
      beltr%np=np
      beltr%xperiodic=xperiodic
      beltr%yperiodic=yperiodic
c-----------------------------------------------------------------------
c     pass variables to array fields of objects beltr.
c-----------------------------------------------------------------------
      beltr%mass=mass
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_getdata
c-----------------------------------------------------------------------
c     subprogram 13. beltrami_dealloc.
c     deallocates pointers of "beltr" object.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_dealloc

      IF(ASSOCIATED(beltr%mass))THEN
         DEALLOCATE(beltr%mass,beltr%gfit)
         CALL bicube_dealloc(beltr%griderr)
         CALL jacobi_dealloc(beltr%basis)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_dealloc
c-----------------------------------------------------------------------
c     subprogram 14. beltrami_diagnose.
c     diagnoses flux function and its derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_diagnose

      INTEGER :: ix,iy,jx,jy,kx,ky,nxc,nyc,np
      REAL(r8) :: dx,dy
      REAL(r8), DIMENSION(0:beltr%nx) :: xx
      REAL(r8), DIMENSION(0:beltr%ny) :: yy
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: x,y,lambda,phi,phix,phiy
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: gmat,gxmat,gymat
c-----------------------------------------------------------------------
c     define sizes and allocate space.
c-----------------------------------------------------------------------
      np=beltr%np
      nxc=beltr%nx*np
      nyc=beltr%ny*np
      xx = (/(ix,ix=0,beltr%nx)/)/REAL(beltr%nx,r8)
      yy = (/(iy,iy=0,beltr%ny)/)/REAL(beltr%ny,r8)
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(x(0:nxc,0:nyc),y(0:nxc,0:nyc),
     $     gmat(4,0:nxc,0:nyc),gxmat(4,0:nxc,0:nyc),
     $     gymat(4,0:nxc,0:nyc),lambda(0:nxc,0:nyc))
c-----------------------------------------------------------------------
c     fill arrays.
c-----------------------------------------------------------------------
      jx=0
      DO ix=0,beltr%nx-1
         dx=one/REAL(nxc,8)
         jy=0
         DO iy=0,beltr%ny-1
            dy=one/REAL(nyc,8)
            DO ky=0,np
               x(jx:jx+np,jy+ky)=xx(ix)+(/(kx,kx=0,np)/)*dx
               y(jx:jx+np,jy+ky)=yy(iy)+ky*dy
            ENDDO
            jy=jy+np
         ENDDO
         jx=jx+np
      ENDDO
c-----------------------------------------------------------------------
c     write contour plots of phi.
c-----------------------------------------------------------------------
      ALLOCATE(phi(0:nxc,0:nyc),phix(0:nxc,0:nyc),phiy(0:nxc,0:nyc))
      CALL beltrami_phi(x,y,phi,phix,phiy)
      OPEN(UNIT=bin_unit,FILE="phi.bin",ACTION="WRITE",
     $     STATUS="REPLACE",FORM="UNFORMATTED")
      WRITE(bin_unit)1,0,1
      WRITE(bin_unit)nxc,nyc
      WRITE(bin_unit)REAL(x,4),REAL(y,4)
      WRITE(bin_unit)REAL(phi,4)
      CLOSE(UNIT=bin_unit)
      DEALLOCATE(phi,phix,phiy)
c-----------------------------------------------------------------------
c     compute gmat and its derivatives.
c-----------------------------------------------------------------------
      jx=0
      DO ix=0,beltr%nx-1
         jy=0
         DO iy=0,beltr%ny-1
            DO ky=0,np
               x(jx:jx+np,jy+ky)=(ix+(/(kx,kx=0,np)/)/REAL(np,r8))
     $              /REAL(beltr%nx,r8)
               y(jx:jx+np,jy+ky)=(iy+ky/REAL(np,r8))/REAL(beltr%ny,r8)
            ENDDO
            jy=jy+np
         ENDDO
         jx=jx+np
      ENDDO
      CALL beltrami_gmat2(x,y,gmat,gxmat,gymat)
c-----------------------------------------------------------------------
c     compute eigenvalue ratio.
c-----------------------------------------------------------------------
      lambda = (gmat(1,:,:) + gmat(2,:,:)
     $     + SQRT((gmat(1,:,:) - gmat(2,:,:))**2 + 4*gmat(3,:,:)**2))
     $     /(two*SQRT(ABS(gmat(4,:,:))))
c-----------------------------------------------------------------------
c     write contour plots of g11 and its derivatives.
c-----------------------------------------------------------------------
      OPEN(UNIT=bin_unit,FILE="gmat.bin",ACTION="WRITE",
     $     STATUS="REPLACE",FORM="UNFORMATTED")
      WRITE(bin_unit)1,0,5
      WRITE(bin_unit)nxc,nyc
      WRITE(bin_unit)REAL(x,4),REAL(y,4)
      WRITE(bin_unit)REAL(gmat(1,:,:),4)
      WRITE(bin_unit)REAL(gmat(2,:,:),4)
      WRITE(bin_unit)REAL(gmat(3,:,:),4)
      WRITE(bin_unit)REAL(gmat(4,:,:),4)
      WRITE(bin_unit)REAL(lambda,4)
      CLOSE(UNIT=bin_unit)
c-----------------------------------------------------------------------
c     deallocate arrays and terminate.
c-----------------------------------------------------------------------
      DEALLOCATE(x,y,gmat,gxmat,gymat,lambda)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_diagnose
      END MODULE beltrami_mod
