c-----------------------------------------------------------------------
c     file p2_interior.F.
c     integration of rhs, jacobian and mass matrix over 2D interior.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_interior_mod.
c     1. p2_interior_rhs.
c     2. p2_interior_jac.
c     3. p2_interior_mass.
c-----------------------------------------------------------------------
c     subprogram 0. p2_interior_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_interior_mod
      USE p2_ct_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_interior_rhs.
c     computes right-hand side of differential equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_interior_rhs(t_var,crd,uu,rr)
      
      REAL(r8), INTENT(IN) :: t_var
      TYPE(coord_type) :: crd
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(IN) :: uu
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(OUT) :: rr

      INTEGER :: imm,mm,iqty,ix,iy,np,nqty,ii
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,x,y
      REAL(r8), DIMENSION((pt2%nq+1)**2) :: ksi,eta,
     $     x_ksi,y_ksi,x_eta,y_eta,jac_xy
      REAL(r8), DIMENSION(pt2%elsize) :: rr_1d
      REAL(r8), DIMENSION(pt2%nqty,(pt2%nq+1)**2) :: u,ux,uy,uk,ue
      REAL(r8), DIMENSION((pt2%nq+1)**2,pt2%nqty) :: fk,fe,s
      REAL(r8), DIMENSION((pt2%nq+1)**2,(pt2%np+1)**2) :: ppbq,pqbq,qpbq
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      np=pt2%np
      nqty=pt2%nqty
      ppbq=pt2%ppbq
      pqbq=pt2%pqbq
      qpbq=pt2%qpbq
c-----------------------------------------------------------------------
c     select quadrature points.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%quad%quadr)
      CASE("gl0")
         quad=pt2%quad%pzero
      CASE("gll")
         quad=pt2%quad%qzero
      END SELECT
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      ii=0
      rr=0
      DO iy=pt2%iymin,pt2%iymax-1
         DO ix=pt2%ixmin,pt2%ixmax-1
c-----------------------------------------------------------------------
c     interpolate positions and compute fluxes and sources.
c-----------------------------------------------------------------------
            CALL p2_sel_interp_all(ix,iy,quad,
     $           uu(:,ix*(np+1):ix*(np+1)+np,iy*(np+1):iy*(np+1)+np),
     $           x,y,u,ux,uy)
            IF(.NOT. pt2%curve_grid)THEN
               CALL p2_ct_default(x,y,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $              jac_xy)
            ELSE
               CALL p2_ct_transf(ii,crd,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $              jac_xy)
            ENDIF
            DO iqty=1,nqty
               uk(iqty,:)=ux(iqty,:)*x_ksi+uy(iqty,:)*y_ksi
               ue(iqty,:)=ux(iqty,:)*x_eta+uy(iqty,:)*y_eta
            ENDDO
            CALL job2_rhs(t_var,ksi,eta,u,uk,ue,fk,fe,s)
            jac_xy = jac_xy*pt2%weight
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
            DO imm=0,pt2%elsize-1
               iqty=MOD(imm,nqty)+1
               mm=imm/nqty+1
               rr_1d(imm+1) = SUM(jac_xy*(s(:,iqty)*ppbq(:,mm)
     $              + fk(:,iqty)*(qpbq(:,mm)*x_ksi + pqbq(:,mm)*y_ksi)
     $              + fe(:,iqty)*(qpbq(:,mm)*x_eta + pqbq(:,mm)*y_eta)))
            ENDDO
            rr(:,ix*(np+1):ix*(np+1)+np,iy*(np+1):iy*(np+1)+np)
     $           = RESHAPE(rr_1d,(/nqty,np+1,np+1/))
            ii = ii + 1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_interior_rhs
c-----------------------------------------------------------------------
c     subprogram 2. p2_interior_jac.
c     computes derivatives of rhs with respect to dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_interior_jac(t_var,crd,uu,drdu)

      REAL(r8), INTENT(IN) :: t_var
      TYPE(coord_type) :: crd
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(IN) :: uu
      REAL(r8), DIMENSION(pt2%elsize,pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(OUT) :: drdu

      INTEGER :: imm,mm,inn,nn,iqty,jqty,ix,iy,np,nqty,ii
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,x,y
      REAL(r8), DIMENSION((pt2%nq+1)**2) :: ksi,eta,
     $     x_ksi,y_ksi,x_eta,y_eta,jac_xy
      REAL(r8), DIMENSION(pt2%nqty,(pt2%nq+1)**2) :: u,ux,uy,uk,ue
      REAL(r8), DIMENSION(pt2%elsize,pt2%elsize) :: drdu_1d
      REAL(r8), DIMENSION((pt2%nq+1)**2,pt2%nqty,pt2%nqty) :: 
     $     fk_u,fk_uk,fk_ue,fe_u,fe_uk,fe_ue,s_u,s_uk,s_ue
      REAL(r8), DIMENSION((pt2%nq+1)**2,(pt2%np+1)**2) :: 
     $     kap1,kap2,ppbq,pqbq,qpbq
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      np=pt2%np
      nqty=pt2%nqty
      ppbq=pt2%ppbq
      pqbq=pt2%pqbq
      qpbq=pt2%qpbq
c-----------------------------------------------------------------------
c     select quadrature points.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%quad%quadr)
      CASE("gl0")
         quad=pt2%quad%pzero
      CASE("gll")
         quad=pt2%quad%qzero
      END SELECT
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      ii=0
      drdu=0
      DO iy=pt2%iymin,pt2%iymax-1
         DO ix=pt2%ixmin,pt2%ixmax-1
c-----------------------------------------------------------------------
c     interpolate positions and compute fluxes and sources.
c-----------------------------------------------------------------------
            CALL p2_sel_interp_all(ix,iy,quad,
     $           uu(:,ix*(np+1):ix*(np+1)+np,iy*(np+1):iy*(np+1)+np),
     $           x,y,u,ux,uy)
            IF(.NOT. pt2%curve_grid)THEN
               CALL p2_ct_default(x,y,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $              jac_xy)
            ELSE
               CALL p2_ct_transf(ii,crd,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $              jac_xy)
            ENDIF
            DO iqty=1,nqty
               uk(iqty,:)=ux(iqty,:)*x_ksi+uy(iqty,:)*y_ksi
               ue(iqty,:)=ux(iqty,:)*x_eta+uy(iqty,:)*y_eta
            ENDDO
            CALL job2_drdu(t_var,ksi,eta,u,uk,ue,fk_u,fk_uk,
     $           fk_ue,fe_u,fe_uk,fe_ue,s_u,s_uk,s_ue)
            jac_xy = jac_xy*pt2%weight
c-----------------------------------------------------------------------
c     calculate kappa transformation matrix.
c-----------------------------------------------------------------------
            kap1=0
            kap2=0
            DO imm=1,(np+1)**2
               kap1(:,imm) = x_ksi*qpbq(:,imm) + y_ksi*pqbq(:,imm)
               kap2(:,imm) = x_eta*qpbq(:,imm) + y_eta*pqbq(:,imm)
            ENDDO
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
            DO imm=0,pt2%elsize-1
               iqty=MOD(imm,nqty)+1
               mm=imm/nqty+1
c-----------------------------------------------------------------------
c     loop over basis polynomials.
c-----------------------------------------------------------------------
               DO inn=0,pt2%elsize-1
                  jqty=MOD(inn,nqty)+1
                  nn=inn/nqty+1
c-----------------------------------------------------------------------
c     compute derivatives of rhs with respect to dependent variables.
c-----------------------------------------------------------------------
                  drdu_1d(inn+1,imm+1) = SUM(jac_xy
     $                 *(ppbq(:,mm)*(ppbq(:,nn)*s_u(:,jqty,iqty)
     $                 + kap1(:,nn)*s_uk(:,jqty,iqty)
     $                 + kap2(:,nn)*s_ue(:,jqty,iqty))
     $                 + kap1(:,mm)*(ppbq(:,nn)*fk_u(:,jqty,iqty)
     $                 + kap1(:,nn)*fk_uk(:,jqty,iqty)
     $                 + kap2(:,nn)*fk_ue(:,jqty,iqty))
     $                 + kap2(:,mm)*(ppbq(:,nn)*fe_u(:,jqty,iqty)
     $                 + kap1(:,nn)*fe_uk(:,jqty,iqty)
     $                 + kap2(:,nn)*fe_ue(:,jqty,iqty))))
               ENDDO
            ENDDO
            drdu(:,:,ix*(np+1):ix*(np+1)+np,iy*(np+1):iy*(np+1)+np)
     $           = RESHAPE(drdu_1d,(/pt2%elsize,nqty,np+1,np+1/))
            ii = ii + 1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_interior_jac
c-----------------------------------------------------------------------
c     subprogram 3. p2_interior_mass.
c     computes main mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_interior_mass(crd,mmat)

      TYPE(coord_type) :: crd
      REAL(r8), DIMENSION((pt2%np+1)**2,pt2%couple_msize,
     $     pt2%jxmin:pt2%jxmax,pt2%jymin:pt2%jymax), INTENT(OUT) :: mmat

      INTEGER :: imm,inn,ix,iy,iqty,jqty,np,np1,nqty,ii,jj
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,x,y
      REAL(r8), DIMENSION((pt2%nq+1)**2) :: ksi,eta,x_ksi,y_ksi,x_eta,
     $     y_eta,jac_xy
      REAL(r8), DIMENSION((pt2%np+1)**2,pt2%couple_msize,(pt2%np+1)**2) 
     $     :: mmat_1d
      REAL(r8), DIMENSION((pt2%nq+1)**2,(pt2%np+1)**2) :: 
     $     kap1,kap2,ppbq,pqbq,qpbq
      REAL(r8), DIMENSION((pt2%nq+1)**2,pt2%nqty,pt2%nqty) :: 
     $     mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      np=pt2%np
      np1=np+1
      nqty=pt2%nqty
      ppbq=pt2%ppbq
      pqbq=pt2%pqbq
      qpbq=pt2%qpbq
c-----------------------------------------------------------------------
c     select quadrature points.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%quad%quadr)
      CASE("gl0")
         quad=pt2%quad%pzero
      CASE("gll")
         quad=pt2%quad%qzero
      END SELECT
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      ii=0
      mmat=0
      DO iy=pt2%iymin,pt2%iymax-1
         y = (REAL(iy,r8) + 0.5*(one+quad))/REAL(pt2%ny,r8)
         DO ix=pt2%ixmin,pt2%ixmax-1
            x = (REAL(ix,r8) + 0.5*(one+quad))/REAL(pt2%nx,r8)
c-----------------------------------------------------------------------
c     remap positions and compute mass matrix.
c-----------------------------------------------------------------------
            IF(.NOT. pt2%curve_grid)THEN
               CALL p2_ct_default(x,y,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $              jac_xy)
            ELSE
               CALL p2_ct_transf(ii,crd,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $              jac_xy)
            ENDIF
            CALL job2_mass(ksi,eta,mass,mass_x,mass_y)
            jac_xy = jac_xy*pt2%weight
c-----------------------------------------------------------------------
c     calculate kappa transformation matrix.
c-----------------------------------------------------------------------
            kap1=0
            kap2=0
            DO imm=1,np1**2
               kap1(:,imm) = x_ksi*qpbq(:,imm) + y_ksi*pqbq(:,imm)
               kap2(:,imm) = x_eta*qpbq(:,imm) + y_eta*pqbq(:,imm)
            ENDDO
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
            DO imm=1,np1**2
c-----------------------------------------------------------------------
c     loop over dependent variable couplings.
c-----------------------------------------------------------------------
               DO jj=1,pt2%couple_msize
                  jqty=pt2%couple_mlist(1,jj)
                  iqty=pt2%couple_mlist(2,jj)
c-----------------------------------------------------------------------
c     loop over basis polynomials.
c-----------------------------------------------------------------------
                  DO inn=1,np1**2
c-----------------------------------------------------------------------
c     accumulate mass matrix.
c-----------------------------------------------------------------------
                     mmat_1d(inn,jj,imm) = SUM(jac_xy*ppbq(:,imm)
     $                    *(mass(:,jqty,iqty)*ppbq(:,inn)
     $                    + mass_x(:,jqty,iqty)*kap1(:,inn)
     $                    + mass_y(:,jqty,iqty)*kap2(:,inn)))
                  ENDDO
               ENDDO
            ENDDO
            mmat(:,:,ix*np1:ix*np1+np,iy*np1:iy*np1+np)
     $           = RESHAPE(mmat_1d,(/np1**2,pt2%couple_msize,np1,np1/))
            ii = ii+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_interior_mass
      END MODULE p2_interior_mod
