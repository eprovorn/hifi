c-----------------------------------------------------------------------
c     file p2_edge.F.
c     evaluation of edge boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_edge_mod.
c     1. p2_edge_rhs_lr.
c     2. p2_edge_rhs_tb.
c     3. p2_edge_drdu_lr.
c     4. p2_edge_drdu_tb.
c     5. p2_edge_mass_lr.
c     6. p2_edge_mass_tb.
c     7. p2_edge_interp.
c     8. p2_edge_ct.
c     9. p2_edge_utransf.
c-----------------------------------------------------------------------
c     subprogram 0. p2_edge_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_edge_mod
      USE p2_ct_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_edge_rhs_lr.
c     computes vertical boundary conditions for rhs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_rhs_lr(t_var,edge,uu,rr,crd_mat)

      REAL(r8), INTENT(IN) :: t_var
      TYPE(edge_type) :: edge
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%np,pt2%jymin:pt2%jymax), 
     $     INTENT(IN) :: uu
      REAL(r8), DIMENSION(ndim,0:pt2%np,pt2%jymin:pt2%jymax), 
     $     INTENT(IN) :: crd_mat
      REAL(r8), DIMENSION(pt2%nqty,pt2%jymin:pt2%jymax), 
     $     INTENT(INOUT) :: rr

      LOGICAL :: bc_flag
      INTEGER :: nqty,iqty,np,nx,ny,iymin,iymax,imy,jy,iy,iyp
      REAL(r8) :: xi
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,weight,ksi,eta,jac,
     $     x_ksi,y_ksi,x_eta,y_eta
      REAL(r8), DIMENSION(ndim,0:pt2%nq) :: nvec
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%nq) :: u,uk,ue,ux,uy,
     $     uxy,uxx,uyy,ukk,uee,uke
      REAL(r8), DIMENSION(0:pt2%nq,pt2%nqty) :: bc
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np) :: bq_amp
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np) :: crd_local
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%np,0:pt2%np) :: uu_local
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      nqty=pt2%nqty
      iymin=pt2%iymin
      iymax=pt2%iymax
      bq_amp=pt2%bq_amp(:,:,1)
c-----------------------------------------------------------------------
c     polar regularity conditions.
c-----------------------------------------------------------------------
      IF(edge%edgenum == 1 .AND. polar)THEN
         DO iqty=1,nqty
            IF(.NOT. pt2%polar_axis(iqty))CYCLE 
            DO jy=pt2%jymin,pt2%jymax-np,np+1
               rr(iqty,jy+1:jy+np-1)=0
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if necessary, zero interior contribution
c-----------------------------------------------------------------------
      bc_flag = .FALSE.
      DO iqty=1,nqty
         SELECT CASE(edge%bc_type(iqty))
         CASE("robin")
            rr(iqty,:)=0
            bc_flag=.TRUE.
         CASE("natural","normflux","robin+")
            bc_flag=.TRUE.
         END SELECT
      ENDDO
      IF(.NOT. bc_flag)RETURN
c-----------------------------------------------------------------------
c     define quadrature points.
c-----------------------------------------------------------------------
      weight=0.5*pt2%quad%weight/REAL(ny,r8)
      IF(pt2%quad%quadr=="gl0")THEN
         quad=pt2%quad%pzero
      ELSEIF(pt2%quad%quadr=="gll")THEN
         quad=pt2%quad%qzero
      ENDIF
      SELECT CASE(edge%edgenum)
      CASE(1)
         xi=-one
      CASE(3)
         xi=one
      END SELECT
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      DO iy=iymin,iymax-1
c-----------------------------------------------------------------------
c     interpolate positions and compute fluxes and sources.
c-----------------------------------------------------------------------
         iyp=iy*(np+1)
         DO iqty=1,nqty
            uu_local(iqty,:,:) = TRANSPOSE(uu(iqty,:,iyp:iyp+np))
         ENDDO
         crd_local(1,:,:) = TRANSPOSE(crd_mat(2,:,iyp:iyp+np))
         crd_local(2,:,:) = TRANSPOSE(crd_mat(1,:,iyp:iyp+np))
         CALL p2_edge_interp(ny,nx,xi,uu_local,u,uy,ux,uxy,uyy,uxx)
         CALL p2_edge_utransf(ny,nx,iy,quad,xi,crd_local,
     $        eta,ksi,y_eta,y_ksi,x_eta,x_ksi,jac,
     $        uy,ux,uxy,uyy,uxx,ue,uk,uee,ukk,uke)

         nvec(1,:)=xi*x_ksi
         nvec(2,:)=xi*x_eta

         CALL job2_edge_rhs(edge,t_var,ksi,eta,nvec,u,uk,ue,ukk,uee,uke,
     $        bc)
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
         DO imy=0,np
            jy=iy*(np+1)+imy
c-----------------------------------------------------------------------
c     loop over equations and sum over quadrature points.
c-----------------------------------------------------------------------
            DO iqty=1,nqty
               SELECT CASE(edge%bc_type(iqty))
               CASE("robin","robin+","natural","normflux")
                  rr(iqty,jy)=rr(iqty,jy) + SUM(weight*jac
     $                 *bq_amp(:,imy)*bc(:,iqty))
               END SELECT
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_rhs_lr
c-----------------------------------------------------------------------
c     subprogram 2. p2_edge_rhs_tb.
c     computes horizontal boundary conditions for rhs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_rhs_tb(t_var,edge,uu,rr,crd_mat)

      REAL(r8), INTENT(IN) :: t_var
      TYPE(edge_type) :: edge
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,0:pt2%np), 
     $     INTENT(IN) :: uu
      REAL(r8), DIMENSION(ndim,pt2%jxmin:pt2%jxmax,0:pt2%np), 
     $     INTENT(IN) :: crd_mat
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax), 
     $     INTENT(INOUT) :: rr

      LOGICAL :: bc_flag
      INTEGER :: nqty,iqty,np,nx,ny,ixmin,ixmax,imx,jx,ix,ixp
      REAL(r8) :: yi
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,weight,ksi,eta,jac,
     $     x_ksi,y_ksi,x_eta,y_eta
      REAL(r8), DIMENSION(ndim,0:pt2%nq) :: nvec
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%nq) :: u,uk,ue,ux,uy,
     $     uxy,uxx,uyy,ukk,uee,uke
      REAL(r8), DIMENSION(0:pt2%nq,pt2%nqty) :: bc
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np) :: bq_amp
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np) :: crd_local
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%np,0:pt2%np) :: uu_local
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      nqty=pt2%nqty
      ixmin=pt2%ixmin
      ixmax=pt2%ixmax
      bq_amp=pt2%bq_amp(:,:,1)
c-----------------------------------------------------------------------
c     if necessary, zero interior contribution
c-----------------------------------------------------------------------
      bc_flag = .FALSE.
      DO iqty=1,nqty
         SELECT CASE(edge%bc_type(iqty))
         CASE("robin")
            rr(iqty,:)=0
            bc_flag=.TRUE.
         CASE("natural","normflux","robin+")
            bc_flag=.TRUE.
         END SELECT
      ENDDO
      IF(.NOT. bc_flag)RETURN
c-----------------------------------------------------------------------
c     define quadrature points.
c-----------------------------------------------------------------------
      weight=0.5*pt2%quad%weight/REAL(nx,r8)
      IF(pt2%quad%quadr=="gl0")THEN
         quad=pt2%quad%pzero
      ELSEIF(pt2%quad%quadr=="gll")THEN
         quad=pt2%quad%qzero
      ENDIF
      SELECT CASE(edge%edgenum)
      CASE(4)
         yi=-one
      CASE(2)
         yi=one
      END SELECT
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      DO ix=ixmin,ixmax-1
c-----------------------------------------------------------------------
c     interpolate positions and compute fluxes and sources.
c-----------------------------------------------------------------------
         ixp=ix*(np+1)
         uu_local=uu(:,ixp:ixp+np,:)
         crd_local=crd_mat(:,ixp:ixp+np,:)
         CALL p2_edge_interp(nx,ny,yi,uu_local,u,ux,uy,uxy,uxx,uyy)
         CALL p2_edge_utransf(nx,ny,ix,quad,yi,crd_local,
     $        ksi,eta,x_ksi,x_eta,y_ksi,y_eta,jac,
     $        ux,uy,uxy,uxx,uyy,uk,ue,ukk,uee,uke)

         nvec(1,:)=yi*y_ksi
         nvec(2,:)=yi*y_eta

         CALL job2_edge_rhs(edge,t_var,ksi,eta,nvec,u,uk,ue,ukk,uee,uke,
     $        bc)
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
         DO imx=0,np
            jx=ix*(np+1)+imx
c-----------------------------------------------------------------------
c     loop over equations and sum over quadrature points.
c-----------------------------------------------------------------------
            DO iqty=1,nqty
               SELECT CASE(edge%bc_type(iqty))
               CASE("robin","robin+","natural","normflux")
                  rr(iqty,jx)=rr(iqty,jx) + SUM(weight*jac
     $                 *bq_amp(:,imx)*bc(:,iqty))
               END SELECT
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_rhs_tb
c-----------------------------------------------------------------------
c     subprogram 3. p2_edge_drdu_lr.
c     computes vertical boundary conditions for drdu.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_drdu_lr(t_var,edge,uu,drdu,crd_mat)

      REAL(r8), INTENT(IN) :: t_var
      TYPE(edge_type) :: edge
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%np,pt2%jymin:pt2%jymax), 
     $     INTENT(IN) :: uu
      REAL(r8), DIMENSION(ndim,0:pt2%np,pt2%jymin:pt2%jymax), 
     $     INTENT(IN) :: crd_mat
      REAL(r8), DIMENSION(pt2%elsize,pt2%nqty,pt2%jymin:pt2%jymax), 
     $     INTENT(INOUT) :: drdu

      LOGICAL :: bc_flag
      INTEGER :: nqty,iqty,jqty,np,nx,ny,iymin,iymax,imy,inx,iny,jy,iy,
     $     iyp,ii,inn
      REAL(r8) :: xi
      REAL(r8), DIMENSION(0:pt2%np) :: pbx
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,weight,ksi,eta,jac,
     $     x_ksi,y_ksi,x_eta,y_eta
      REAL(r8), DIMENSION(ndim,0:pt2%nq) :: nvec
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%nq) :: u,uk,ue,ux,uy,
     $     uxy,uxx,uyy,ukk,uee,uke
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np) :: bq_amp
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np) :: crd_local
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%np,0:pt2%np) :: uu_local
      REAL(r8), DIMENSION(0:pt2%nq,pt2%nqty,pt2%nqty) :: 
     $     bc_u,bc_uk,bc_ue,bc_ukk,bc_uee,bc_uke
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np,0:pt2%np,5) :: kappa
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      nqty=pt2%nqty
      iymin=pt2%iymin
      iymax=pt2%iymax
      bq_amp=pt2%bq_amp(:,:,1)
c-----------------------------------------------------------------------
c     polar regularity conditions.
c-----------------------------------------------------------------------
      IF(edge%edgenum == 1 .AND. polar)THEN
         DO iqty=1,nqty
            IF(.NOT. pt2%polar_axis(iqty))CYCLE 
            DO jy=pt2%jymin,pt2%jymax-np,np+1
               drdu(:,iqty,jy+1:jy+np-1) = 0
               DO iy=1,np-1
                  ii=pt2%elmap(iqty,0,iy)
                  drdu(ii,iqty,jy+iy) = 1.
               ENDDO
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if necessary, zero interior contribution
c-----------------------------------------------------------------------
      bc_flag = .FALSE.
      DO iqty=1,nqty
         SELECT CASE(edge%bc_type(iqty))
         CASE("robin")
            drdu(:,iqty,:)=0
            bc_flag=.TRUE.
         CASE("natural","normflux","robin+")
            bc_flag=.TRUE.
         END SELECT
      ENDDO
      IF(.NOT. bc_flag)RETURN
c-----------------------------------------------------------------------
c     define quadrature points.
c-----------------------------------------------------------------------
      weight=0.5*pt2%quad%weight/REAL(ny,r8)
      IF(pt2%quad%quadr=="gl0")THEN
         quad=pt2%quad%pzero
      ELSEIF(pt2%quad%quadr=="gll")THEN
         quad=pt2%quad%qzero
      ENDIF
      SELECT CASE(edge%edgenum)
      CASE(1)
         xi=-one
      CASE(3)
         xi=one
      END SELECT
      CALL jacobi_basis(xi,pt2%basis)
      pbx=pt2%basis%pb
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      DO iy=iymin,iymax-1
c-----------------------------------------------------------------------
c     interpolate positions and compute fluxes and sources.
c-----------------------------------------------------------------------
         iyp=iy*(np+1)
         DO iqty=1,nqty
            uu_local(iqty,:,:) = TRANSPOSE(uu(iqty,:,iyp:iyp+np))
         ENDDO
         crd_local(1,:,:) = TRANSPOSE(crd_mat(2,:,iyp:iyp+np))
         crd_local(2,:,:) = TRANSPOSE(crd_mat(1,:,iyp:iyp+np))      
         CALL p2_edge_interp(ny,nx,xi,uu_local,u,uy,ux,uxy,uyy,uxx)
         CALL p2_edge_utransf(ny,nx,iy,quad,xi,crd_local,
     $        eta,ksi,y_eta,y_ksi,x_eta,x_ksi,jac,
     $        uy,ux,uxy,uyy,uxx,ue,uk,uee,ukk,uke,kappa)
         DO jy=0,pt2%nq
            kappa(jy,:,:,1)=TRANSPOSE(kappa(jy,:,:,1))
            kappa(jy,:,:,2)=TRANSPOSE(kappa(jy,:,:,2))
            kappa(jy,:,:,3)=TRANSPOSE(kappa(jy,:,:,3))
            kappa(jy,:,:,4)=TRANSPOSE(kappa(jy,:,:,4))
            kappa(jy,:,:,5)=TRANSPOSE(kappa(jy,:,:,5))
         ENDDO
         kappa(:,:,:,1:2)=kappa(:,:,:,2:1:-1)
         kappa(:,:,:,3:4)=kappa(:,:,:,4:3:-1)

         nvec(1,:)=xi*x_ksi
         nvec(2,:)=xi*x_eta

         CALL job2_edge_drdu(edge,t_var,ksi,eta,nvec,u,uk,ue,ukk,uee,
     $        uke,bc_u,bc_uk,bc_ue,bc_ukk,bc_uee,bc_uke)
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
         DO imy=0,np
            jy=iy*(np+1)+imy
c-----------------------------------------------------------------------
c     loop over equations
c-----------------------------------------------------------------------
            DO iqty=1,nqty
               SELECT CASE(edge%bc_type(iqty))
               CASE("robin","robin+","natural","normflux")
                  DO inn=0,pt2%elsize-1
                     iny=inn/(nqty*(np+1))
                     inx=MOD(inn,nqty*(np+1))/nqty
                     jqty=MOD(inn,nqty)+1

                     drdu(inn+1,iqty,jy) = drdu(inn+1,iqty,jy)
     $                    + SUM(jac*weight*bq_amp(:,imy)
     $                    *(pbx(inx)*bq_amp(:,iny)*bc_u(:,jqty,iqty)
     $                    + kappa(:,inx,iny,1)*bc_uk(:,jqty,iqty)
     $                    + kappa(:,inx,iny,2)*bc_ue(:,jqty,iqty)
     $                    + kappa(:,inx,iny,3)*bc_ukk(:,jqty,iqty)
     $                    + kappa(:,inx,iny,4)*bc_uee(:,jqty,iqty)
     $                    + kappa(:,inx,iny,5)*bc_uke(:,jqty,iqty)))
                  ENDDO
               END SELECT
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_drdu_lr
c-----------------------------------------------------------------------
c     subprogram 4. p2_edge_drdu_tb.
c     computes horizontal boundary conditions for drdu.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_drdu_tb(t_var,edge,uu,drdu,crd_mat)

      REAL(r8), INTENT(IN) :: t_var
      TYPE(edge_type) :: edge
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,0:pt2%np), 
     $     INTENT(IN) :: uu
      REAL(r8), DIMENSION(ndim,pt2%jxmin:pt2%jxmax,0:pt2%np), 
     $     INTENT(IN) :: crd_mat
      REAL(r8), DIMENSION(pt2%elsize,pt2%nqty,pt2%jxmin:pt2%jxmax), 
     $     INTENT(INOUT) :: drdu

      LOGICAL :: bc_flag
      INTEGER :: nqty,iqty,jqty,np,nx,ny,ixmin,ixmax,imx,inx,iny,jx,ix,
     $     ixp,inn
      REAL(r8) :: yi
      REAL(r8), DIMENSION(0:pt2%np) :: pby
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,weight,ksi,eta,jac,
     $     x_ksi,y_ksi,x_eta,y_eta
      REAL(r8), DIMENSION(ndim,0:pt2%nq) :: nvec
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%nq) :: u,uk,ue,ux,uy,
     $     uxy,uxx,uyy,ukk,uee,uke
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np) :: bq_amp
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np) :: crd_local
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%np,0:pt2%np) :: uu_local
      REAL(r8), DIMENSION(0:pt2%nq,pt2%nqty,pt2%nqty) :: 
     $     bc_u,bc_uk,bc_ue,bc_ukk,bc_uee,bc_uke
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np,0:pt2%np,5) :: kappa
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      nqty=pt2%nqty
      ixmin=pt2%ixmin
      ixmax=pt2%ixmax
      bq_amp=pt2%bq_amp(:,:,1)
c-----------------------------------------------------------------------
c     if necessary, zero interior contribution
c-----------------------------------------------------------------------
      bc_flag = .FALSE.
      DO iqty=1,nqty
         SELECT CASE(edge%bc_type(iqty))
         CASE("robin")
            drdu(:,iqty,:)=0
            bc_flag=.TRUE.
         CASE("natural","normflux","robin+")
            bc_flag=.TRUE.
         END SELECT
      ENDDO
      IF(.NOT. bc_flag)RETURN
c-----------------------------------------------------------------------
c     define quadrature points.
c-----------------------------------------------------------------------
      weight=0.5*pt2%quad%weight/REAL(nx,r8)
      IF(pt2%quad%quadr=="gl0")THEN
         quad=pt2%quad%pzero
      ELSEIF(pt2%quad%quadr=="gll")THEN
         quad=pt2%quad%qzero
      ENDIF
      SELECT CASE(edge%edgenum)
      CASE(4)
         yi=-one
      CASE(2)
         yi=one
      END SELECT
      CALL jacobi_basis(yi,pt2%basis)
      pby=pt2%basis%pb
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      DO ix=ixmin,ixmax-1
c-----------------------------------------------------------------------
c     interpolate positions and compute fluxes and sources.
c-----------------------------------------------------------------------
         ixp=ix*(np+1)
         uu_local=uu(:,ixp:ixp+np,:)
         crd_local=crd_mat(:,ixp:ixp+np,:)
         CALL p2_edge_interp(nx,ny,yi,uu_local,u,ux,uy,uxy,uxx,uyy)
         CALL p2_edge_utransf(nx,ny,ix,quad,yi,crd_local,
     $        ksi,eta,x_ksi,x_eta,y_ksi,y_eta,jac,
     $        ux,uy,uxy,uxx,uyy,uk,ue,ukk,uee,uke,kappa)

         nvec(1,:)=yi*y_ksi
         nvec(2,:)=yi*y_eta

         CALL job2_edge_drdu(edge,t_var,ksi,eta,nvec,u,uk,ue,ukk,uee,
     $        uke,bc_u,bc_uk,bc_ue,bc_ukk,bc_uee,bc_uke)
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
         DO imx=0,np
            jx=ix*(np+1)+imx
c-----------------------------------------------------------------------
c     loop over equations
c-----------------------------------------------------------------------
            DO iqty=1,nqty
               SELECT CASE(edge%bc_type(iqty))
               CASE("robin","robin+","natural","normflux")
                  DO inn=0,pt2%elsize-1
                     iny=inn/(nqty*(np+1))
                     inx=MOD(inn,nqty*(np+1))/nqty
                     jqty=MOD(inn,nqty)+1

                     drdu(inn+1,iqty,jx) = drdu(inn+1,iqty,jx)
     $                    + SUM(jac*weight*bq_amp(:,imx)
     $                    *(bq_amp(:,inx)*pby(iny)*bc_u(:,jqty,iqty)
     $                    + kappa(:,inx,iny,1)*bc_uk(:,jqty,iqty)
     $                    + kappa(:,inx,iny,2)*bc_ue(:,jqty,iqty)
     $                    + kappa(:,inx,iny,3)*bc_ukk(:,jqty,iqty)
     $                    + kappa(:,inx,iny,4)*bc_uee(:,jqty,iqty)
     $                    + kappa(:,inx,iny,5)*bc_uke(:,jqty,iqty)))
                  ENDDO
               END SELECT
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_drdu_tb
c-----------------------------------------------------------------------
c     subprogram 5. p2_edge_mass_lr.
c     computes mass matrix for vertical boundaries.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_mass_lr(edge,mmat,crd_mat)

      TYPE(edge_type) :: edge
      REAL(r8), DIMENSION(ndim,0:pt2%np,pt2%jymin:pt2%jymax), 
     $     INTENT(IN) :: crd_mat
      REAL(r8), DIMENSION((pt2%np+1)**2,pt2%couple_msize,
     $     pt2%jymin:pt2%jymax), INTENT(INOUT) :: mmat

      LOGICAL :: bc_flag
      INTEGER :: nqty,iqty,jqty,np,nx,ny,iymin,iymax,imy,iny,inx,jy,iy,
     $     iyp,ii,inn
      REAL(r8) :: xi
      REAL(r8), DIMENSION(0:pt2%np) :: pbx,qbx
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,weight,ksi,eta,jac,pby,qby,
     $     tby,x_ksi,y_ksi,x_eta,y_eta
      REAL(r8), DIMENSION(ndim,0:pt2%nq) :: nvec
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np,2) :: bq_amp
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np) :: crd_local
      REAL(r8), DIMENSION(0:pt2%nq,pt2%nqty,pt2%nqty) ::
     $     mass,mass_k,mass_e
c-----------------------------------------------------------------------
c     define local quantities.
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      nqty=pt2%nqty
      iymin=pt2%iymin
      iymax=pt2%iymax
      bq_amp=pt2%bq_amp(:,:,1:2)
c-----------------------------------------------------------------------
c     polar regularity conditions.
c-----------------------------------------------------------------------
      IF(edge%edgenum == 1 .AND. polar)THEN
         DO ii=1,pt2%couple_msize
            iqty=pt2%couple_mlist(2,ii)
            IF(.NOT. pt2%polar_axis(iqty))CYCLE 
            DO jy=pt2%jymin,pt2%jymax-np,np+1
               mmat(:,ii,jy+1:jy+np-1) = 0
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if necessary, zero interior contribution
c-----------------------------------------------------------------------
      bc_flag = .FALSE.
      DO ii=1,pt2%couple_msize
         iqty=pt2%couple_mlist(2,ii)
         SELECT CASE(edge%bc_type(iqty))
         CASE("robin")
            mmat(:,ii,:)=0
            bc_flag=.TRUE.
         CASE("robin+")
            bc_flag=.TRUE.
         END SELECT
      ENDDO
      IF(.NOT. bc_flag)RETURN
c-----------------------------------------------------------------------
c     define quadrature points.
c-----------------------------------------------------------------------
      weight=0.5*pt2%quad%weight/REAL(ny,r8)
      IF(pt2%quad%quadr=="gl0")THEN
         quad=pt2%quad%pzero
      ELSEIF(pt2%quad%quadr=="gll")THEN
         quad=pt2%quad%qzero
      ENDIF
      SELECT CASE(edge%edgenum)
      CASE(1)
         xi=-one
      CASE(3)
         xi=one
      END SELECT
      CALL jacobi_basis(xi,pt2%basis)
      pbx=pt2%basis%pb
      qbx=pt2%basis%qb*(two*nx)
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      DO iy=iymin,iymax-1
c-----------------------------------------------------------------------
c     interpolate positions and compute mass matrices.
c-----------------------------------------------------------------------
         iyp=iy*(np+1)
         crd_local(1,:,:) = TRANSPOSE(crd_mat(2,:,iyp:iyp+np))
         crd_local(2,:,:) = TRANSPOSE(crd_mat(1,:,iyp:iyp+np))
         CALL p2_edge_ct(ny,nx,iy,quad,xi,crd_local,
     $        eta,ksi,y_eta,y_ksi,x_eta,x_ksi,jac)

         nvec(1,:)=xi*x_ksi
         nvec(2,:)=xi*x_eta

         CALL job2_edge_mass(edge,ksi,eta,nvec,mass,mass_k,mass_e)
c-----------------------------------------------------------------------
c     loop over equations
c-----------------------------------------------------------------------
         DO ii=1,pt2%couple_msize
            iqty=pt2%couple_mlist(2,ii)
            IF(.NOT. (edge%bc_type(iqty)=="robin"
     $           .OR. edge%bc_type(iqty)=="robin+"))CYCLE
            jqty=pt2%couple_mlist(1,ii)
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
            DO imy=0,np
               tby=bq_amp(:,imy,1)
               jy=iy*(np+1)+imy
               DO inn=0,(np+1)**2-1
                  inx=MOD(inn,np+1)
                  iny=inn/(np+1)

                  pby=bq_amp(:,iny,1)
                  qby=bq_amp(:,iny,2)*(two*ny)

                  mmat(inn+1,ii,jy) = mmat(inn+1,ii,jy)
     $                 + SUM(jac*weight*tby
     $                 *(mass(:,jqty,iqty)*pbx(inx)*pby
     $                 + mass_k(:,jqty,iqty)
     $                 *(x_ksi*qbx(inx)*pby+y_ksi*pbx(inx)*qby)
     $                 + mass_e(:,jqty,iqty)
     $                 *(x_eta*qbx(inx)*pby+y_eta*pbx(inx)*qby)))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_mass_lr
c-----------------------------------------------------------------------
c     subprogram 6. p2_edge_mass_tb.
c     computes mass matrix for horizontal boundaries.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_mass_tb(edge,mmat,crd_mat)

      TYPE(edge_type) :: edge
      REAL(r8), DIMENSION(ndim,pt2%jxmin:pt2%jxmax,0:pt2%np), 
     $     INTENT(IN) :: crd_mat
      REAL(r8), DIMENSION((pt2%np+1)**2,pt2%couple_msize,
     $     pt2%jxmin:pt2%jxmax), INTENT(INOUT) :: mmat

      LOGICAL :: bc_flag
      INTEGER :: nqty,iqty,jqty,np,nx,ny,ixmin,ixmax,imx,inx,iny,jx,ix,
     $     ixp,ii,inn
      REAL(r8) :: yi
      REAL(r8), DIMENSION(0:pt2%np) :: pby,qby
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,weight,ksi,eta,jac,pbx,qbx,
     $     tbx,x_ksi,y_ksi,x_eta,y_eta
      REAL(r8), DIMENSION(ndim,0:pt2%nq) :: nvec
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np,2) :: bq_amp
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np) :: crd_local
      REAL(r8), DIMENSION(0:pt2%nq,pt2%nqty,pt2%nqty) ::
     $     mass,mass_k,mass_e
c-----------------------------------------------------------------------
c     define local quantities.
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      nqty=pt2%nqty
      ixmin=pt2%ixmin
      ixmax=pt2%ixmax
      bq_amp=pt2%bq_amp(:,:,1:2)
c-----------------------------------------------------------------------
c     if necessary, zero interior contribution
c-----------------------------------------------------------------------
      bc_flag = .FALSE.
      DO ii=1,pt2%couple_msize
         iqty=pt2%couple_mlist(2,ii)
         SELECT CASE(edge%bc_type(iqty))
         CASE("robin")
            mmat(:,ii,:)=0
            bc_flag=.TRUE.
         CASE("robin+")
            bc_flag=.TRUE.
         END SELECT
      ENDDO
      IF(.NOT. bc_flag)RETURN
c-----------------------------------------------------------------------
c     define quadrature points.
c-----------------------------------------------------------------------
      weight=0.5*pt2%quad%weight/REAL(nx,r8)
      IF(pt2%quad%quadr=="gl0")THEN
         quad=pt2%quad%pzero
      ELSEIF(pt2%quad%quadr=="gll")THEN
         quad=pt2%quad%qzero
      ENDIF
      SELECT CASE(edge%edgenum)
      CASE(4)
         yi=-one
      CASE(2)
         yi=one
      END SELECT
      CALL jacobi_basis(yi,pt2%basis)
      pby=pt2%basis%pb
      qby=pt2%basis%qb*(two*ny)
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      DO ix=ixmin,ixmax-1
c-----------------------------------------------------------------------
c     map positions and compute jacobian of coordinate system.
c-----------------------------------------------------------------------
         ixp=ix*(np+1)
         crd_local=crd_mat(:,ixp:ixp+np,:)
         CALL p2_edge_ct(nx,ny,ix,quad,yi,crd_local,
     $        ksi,eta,x_ksi,x_eta,y_ksi,y_eta,jac)

         nvec(1,:)=yi*y_ksi
         nvec(2,:)=yi*y_eta

         CALL job2_edge_mass(edge,ksi,eta,nvec,mass,mass_k,mass_e)
c-----------------------------------------------------------------------
c     loop over equations
c-----------------------------------------------------------------------
         DO ii=1,pt2%couple_msize
            iqty=pt2%couple_mlist(2,ii)
            IF(.NOT. (edge%bc_type(iqty)=="robin"
     $           .OR. edge%bc_type(iqty)=="robin+"))CYCLE
            jqty=pt2%couple_mlist(1,ii)
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
            DO imx=0,np
               tbx=bq_amp(:,imx,1)
               jx=ix*(np+1)+imx
               DO inn=0,(np+1)**2-1
                  inx=MOD(inn,np+1)
                  iny=inn/(np+1)

                  pbx=bq_amp(:,inx,1)
                  qbx=bq_amp(:,inx,2)*(two*nx)

                  mmat(inn+1,ii,jx) = mmat(inn+1,ii,jx)
     $                 + SUM(jac*weight*tbx
     $                 *(mass(:,jqty,iqty)*pbx*pby(iny)
     $                 + mass_k(:,jqty,iqty)
     $                 *(x_ksi*qbx*pby(iny)+y_ksi*pbx*qby(iny))
     $                 + mass_e(:,jqty,iqty)
     $                 *(x_eta*qbx*pby(iny)+y_eta*pbx*qby(iny))))
               ENDDO
             ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_mass_tb
c-----------------------------------------------------------------------
c     subprogram 7. p2_edge_interp.
c     interpolate solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_interp(nx,ny,yi,uu,u,ux,uy,uxy,uxx,uyy)

      INTEGER, INTENT(IN) :: nx,ny
      REAL(r8), INTENT(IN) :: yi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(r8), DIMENSION(:,0:), INTENT(OUT) :: u,ux,uy,uxy,uxx,uyy

      INTEGER :: imx,jx,imy,np
      REAL(r8) :: dxf,dyf
      REAL(r8), DIMENSION(0:pt2%np) :: pby,qby,q2by,pbx,qbx,q2bx
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      np=pt2%np
c-----------------------------------------------------------------------
c     compute positions and increments.
c-----------------------------------------------------------------------
      dxf=two*nx
      dyf=two*ny

      CALL jacobi_basis(yi,pt2%basis)
      pby=pt2%basis%pb
      qby=pt2%basis%qb*dyf
      q2by=pt2%basis%qb2*dyf**2
c-----------------------------------------------------------------------
c     interpolate solutions and derivatives..
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      uxy=0
      uxx=0
      uyy=0
      DO jx=0,pt2%nq
         pbx=pt2%bq_amp(jx,:,1)
         qbx=pt2%bq_amp(jx,:,2)*dxf
         q2bx=pt2%bq_amp(jx,:,3)*dxf**2
         DO imy=0,np
            DO imx=0,np
               u(:,jx)=u(:,jx) + uu(:,imx,imy)*pbx(imx)*pby(imy)
               ux(:,jx)=ux(:,jx) + uu(:,imx,imy)*qbx(imx)*pby(imy)
               uy(:,jx)=uy(:,jx) + uu(:,imx,imy)*pbx(imx)*qby(imy)
               uxy(:,jx)=uxy(:,jx) + uu(:,imx,imy)*qbx(imx)*qby(imy)
               uxx(:,jx)=uxx(:,jx) + uu(:,imx,imy)*q2bx(imx)*pby(imy)
               uyy(:,jx)=uyy(:,jx) + uu(:,imx,imy)*pbx(imx)*q2by(imy)
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_interp
c-----------------------------------------------------------------------
c     subprogram 8. p2_edge_ct.
c     calculate coordinate transformation matrices for boundaries
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_ct(nx,ny,ix,quad,yi,crd_mat,ksi,eta,
     $     x_k,x_e,y_k,y_e,jac)
      
      INTEGER, INTENT(IN) :: nx,ny,ix
      REAL(r8), INTENT(IN) :: yi
      REAL(r8), DIMENSION(0:pt2%nq), INTENT(IN) :: quad
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np), INTENT(IN) :: crd_mat
      REAL(r8), DIMENSION(0:pt2%nq), INTENT(OUT) :: ksi,eta,jac,
     $     x_k,y_k,x_e,y_e
      
      REAL(r8), DIMENSION(0:pt2%nq) :: x
      REAL(r8), DIMENSION(ndim,0:pt2%nq) :: kt,kt_x,kt_y,kt_xy,kt_xx,
     $     kt_yy
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      ksi=0
      eta=0
      x_k=0
      y_k=0
      x_e=0
      y_e=0
      jac=0
c-----------------------------------------------------------------------
c     if uniform [0,1]x[0,1] grid -- assign default values.
c-----------------------------------------------------------------------
      IF(.NOT. pt2%curve_grid)THEN
         x_k = one
         x_e = zero
         y_k = zero
         y_e = one
         jac = one
         ksi = (REAL(ix,r8) + 0.5*(one+quad))/REAL(nx,r8)
         eta = 0.5*(yi+one)
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     interpolate physical positions and coordinates.
c-----------------------------------------------------------------------
      CALL p2_edge_interp(nx,ny,yi,crd_mat,kt,kt_x,kt_y,
     $     kt_xy,kt_xx,kt_yy)
c-----------------------------------------------------------------------
c     calculate jacobian of the coordinate transformation
c-----------------------------------------------------------------------
      jac = kt_x(1,:)*kt_y(2,:) - kt_x(2,:)*kt_y(1,:)
c-----------------------------------------------------------------------
c     calculate derivatives of current w.r.t. fixed coord. system
c-----------------------------------------------------------------------
      ksi = kt(1,:)
      eta = kt(2,:)
      x_k = kt_y(2,:)/jac
      y_k = -kt_x(2,:)/jac
      x_e = -kt_y(1,:)/jac
      y_e = kt_x(1,:)/jac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_ct
c-----------------------------------------------------------------------
c     subprogram 9. p2_edge_utransf.
c     transform physical variables from logical to physical space
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_edge_utransf(nx,ny,ix,quad,yi,crd_mat,ksi,eta,
     $     x_k,x_e,y_k,y_e,jac,ux,uy,uxy,uxx,uyy,uk,ue,ukk,uee,uke,
     $     kappa)
      
      INTEGER, INTENT(IN) :: nx,ny,ix
      REAL(r8), INTENT(IN) :: yi
      REAL(r8), DIMENSION(0:pt2%nq), INTENT(IN) :: quad
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np), INTENT(IN) :: crd_mat
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%nq), INTENT(IN) :: 
     $     ux,uy,uxy,uxx,uyy
      REAL(r8), DIMENSION(0:pt2%nq), INTENT(OUT) :: ksi,eta,jac,
     $     x_k,y_k,x_e,y_e
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%nq), INTENT(OUT) :: 
     $     uk,ue,ukk,uee,uke
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np,0:pt2%np,5), 
     $     INTENT(OUT), OPTIONAL :: kappa
      
      INTEGER :: iqty,np,imx,imy
      REAL(r8) :: dxf,dyf
      REAL(r8), DIMENSION(0:pt2%np) :: pby,qby,q2by
      REAL(r8), DIMENSION(0:pt2%nq) :: x,pbx,qbx,q2bx,jac_x,jac_y,
     $     x_kk,x_ee,x_ke,y_kk,y_ee,y_ke
      REAL(r8), DIMENSION(ndim,0:pt2%nq) :: kt,kt_x,kt_y,kt_xy,kt_xx,
     $     kt_yy
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np,3) :: bq_amp
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      uk=0
      ue=0
      ukk=0
      uee=0
      uke=0
      IF(PRESENT(kappa))kappa=0
      np=pt2%np
      dxf=two*nx
      dyf=two*ny
      bq_amp=pt2%bq_amp
c-----------------------------------------------------------------------
c     map positions and compute jacobian of coordinate system.
c-----------------------------------------------------------------------
      CALL p2_edge_ct(nx,ny,ix,quad,yi,crd_mat,ksi,eta,x_k,x_e,y_k,y_e,
     $     jac)

      CALL jacobi_basis(yi,pt2%basis)
      pby=pt2%basis%pb
      qby=pt2%basis%qb*dyf
      q2by=pt2%basis%qb2*dyf**2
c-----------------------------------------------------------------------
c     if uniform [0,1]x[0,1] grid -- assign default values.
c-----------------------------------------------------------------------
      IF(.NOT. pt2%curve_grid)THEN
         DO imy=0,np
            DO imx=0,np
               pbx=pt2%bq_amp(:,imx,1)
               qbx=pt2%bq_amp(:,imx,2)*dxf
               q2bx=pt2%bq_amp(:,imx,3)*dxf**2
               IF(PRESENT(kappa))THEN
                  kappa(:,imx,imy,1)=qbx*pby(imy)
                  kappa(:,imx,imy,2)=pbx*qby(imy)
                  kappa(:,imx,imy,3)=q2bx*pby(imy)
                  kappa(:,imy,imy,4)=pbx*q2by(imy)
                  kappa(:,imx,imy,5)=qbx*qby(imy)
               ENDIF
            ENDDO
         ENDDO
         uk=ux
         ue=uy
         ukk=uxx
         uee=uyy
         uke=uxy
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     interpolate physical positions and coordinates.
c-----------------------------------------------------------------------
      CALL p2_edge_interp(nx,ny,yi,crd_mat,kt,kt_x,kt_y,
     $     kt_xy,kt_xx,kt_yy)
c-----------------------------------------------------------------------
c     calculate jacobian of the coordinate systems
c-----------------------------------------------------------------------
      jac = kt_x(1,:)*kt_y(2,:) - kt_x(2,:)*kt_y(1,:)
      jac_x = kt_xx(1,:)*kt_y(2,:) + kt_x(1,:)*kt_xy(2,:)
     $     - kt_xx(2,:)*kt_y(1,:) - kt_x(2,:)*kt_xy(1,:)
      jac_y = kt_xy(1,:)*kt_y(2,:) + kt_x(1,:)*kt_yy(2,:)
     $     - kt_xy(2,:)*kt_y(1,:) - kt_x(2,:)*kt_yy(1,:)
c-----------------------------------------------------------------------
c     calculate second order transformation derivatives
c-----------------------------------------------------------------------
      x_kk=(kt_xy(2,:)*x_k + kt_yy(2,:)*y_k)/jac
     $     - kt_y(2,:)*(jac_x*x_k + jac_y*y_k)/jac**2
      x_ee=-(kt_xy(1,:)*x_e + kt_yy(1,:)*y_e)/jac
     $     + kt_y(1,:)*(jac_x*x_e + jac_y*y_e)/jac**2
      x_ke=(kt_xy(2,:)*x_e + kt_yy(2,:)*y_e)/jac
     $     - kt_y(2,:)*(jac_x*x_e + jac_y*y_e)/jac**2

      y_kk=-(kt_xx(2,:)*x_k + kt_xy(2,:)*y_k)/jac
     $     + kt_x(2,:)*(jac_x*x_k + jac_y*y_k)/jac**2
      y_ee=(kt_xx(1,:)*x_e + kt_xy(1,:)*y_e)/jac
     $     - kt_x(1,:)*(jac_x*x_e + jac_y*y_e)/jac**2
      y_ke=-(kt_xx(2,:)*x_e + kt_xy(2,:)*y_e)/jac
     $     + kt_x(2,:)*(jac_x*x_e + jac_y*y_e)/jac**2
c-----------------------------------------------------------------------
c     calculate kappa transformation matrix.
c-----------------------------------------------------------------------
      IF(PRESENT(kappa))THEN
         DO imy=0,np
            DO imx=0,np
               pbx=bq_amp(:,imx,1)
               qbx=bq_amp(:,imx,2)*dxf
               q2bx=bq_amp(:,imx,3)*dxf**2
               
               kappa(:,imx,imy,1)=x_k*qbx*pby(imy) + y_k*pbx*qby(imy)
               kappa(:,imx,imy,2)=x_e*qbx*pby(imy) + y_e*pbx*qby(imy)
               kappa(:,imx,imy,3)=q2bx*pby(imy)*x_k**2 
     $              + two*qbx*qby(imy)*x_k*y_k + pbx*q2by(imy)*y_k**2 
     $              + qbx*pby(imy)*x_kk + pbx*qby(imy)*y_kk
               kappa(:,imx,imy,4)=q2bx*pby(imy)*x_e**2 
     $              + two*qbx*qby(imy)*x_e*y_e + pbx*q2by(imy)*y_e**2 
     $              + qbx*pby(imy)*x_ee + pbx*qby(imy)*y_ee
               kappa(:,imx,imy,5)=q2bx*pby(imy)*x_k*x_e 
     $              + qbx*qby(imy)*(x_k*y_e + x_e*y_k) 
     $              + pbx*q2by(imy)*y_k*y_e 
     $              + qbx*pby(imy)*x_ke + pbx*qby(imy)*y_ke
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     transform physical variables from logical to physical space
c-----------------------------------------------------------------------
      DO iqty=1,pt2%nqty
         uk(iqty,:)=ux(iqty,:)*x_k + uy(iqty,:)*y_k
         ue(iqty,:)=ux(iqty,:)*x_e + uy(iqty,:)*y_e
         ukk(iqty,:)=uxx(iqty,:)*x_k**2 + two*uxy(iqty,:)*x_k*y_k
     $        + uyy(iqty,:)*y_k**2 + ux(iqty,:)*x_kk + uy(iqty,:)*y_kk
         uee(iqty,:)=uxx(iqty,:)*x_e**2 + two*uxy(iqty,:)*x_e*y_e
     $        + uyy(iqty,:)*y_e**2 + ux(iqty,:)*x_ee + uy(iqty,:)*y_ee
         uke(iqty,:)=uxx(iqty,:)*x_k*x_e + uyy(iqty,:)*y_k*y_e
     $        + uxy(iqty,:)*(x_k*y_e + x_e*y_k)
     $        + ux(iqty,:)*x_ke + uy(iqty,:)*y_ke
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_edge_utransf
      END MODULE p2_edge_mod
