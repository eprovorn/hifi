c-----------------------------------------------------------------------
c     file transport.f.
c     calculates transport coefficients.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. transport_mod.
c     1. transport_seteta.
c     2. transport_seteta_u.
c     3. transport_setkaniso.
c     4. transport_setkaniso_u.
c     5. transport_kbrag.
c     6. transport_kbrag_u.
c     7. transport_BdotT.
c     8. transport_BdotT_u.
c     9. transport_mubrag.
c     10. transport_mubrag_u.
c     11. transport_mun.
c     12. transport_mun_u.
c     13. transport_hexch.
c     14. transport_hexch_u.
c     15. transport_radloss.
c     16. transport_radloss_u.
c-----------------------------------------------------------------------
c     subprogram 0. transport_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE transport_mod
      USE local_mod
      IMPLICIT NONE

      REAL(r8), PARAMETER, PRIVATE :: fc=3._r8,
     $     gamma=5._r8/3._r8,Bmin=1.e-8,eps=1.e-12

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. transport_seteta.
c     sets resistivity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c p_e - electron pressure, p_i - ion pressure
c-----------------------------------------------------------------------
      SUBROUTINE transport_seteta(eta_case,x,rho,p_e,p_i,j,jc_norm,
     $     eta_anm_norm, chod_const,etas_norm,etac_norm,
     $     v_chod_norm,r_eta,eta,etavac,coeff,eta_local)
      
      CHARACTER(*), INTENT(IN) :: eta_case
      REAL(r8), INTENT(IN) :: etas_norm,etac_norm,v_chod_norm,r_eta,eta,
     $     etavac,chod_const, jc_norm, eta_anm_norm
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,rho,p_e,p_i,j
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: eta_local
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: coeff

      REAL(r8), DIMENSION(SIZE(rho,1),SIZE(rho,2)) :: Te_local,Ti_local, 
     $     dj, j_c, eta1_anm, deta_dj
c-----------------------------------------------------------------------
c     initialize values.
c-----------------------------------------------------------------------
      eta_local=eta
      Te_local = MAX(p_e/rho,eps)
      Ti_local = MAX(p_i/rho,eps)
      dj = 0.
      j_c = 0.
      eta1_anm=0.
      deta_dj=0. 
      coeff = 0.
c-----------------------------------------------------------------------
c     set resistivity.
c-----------------------------------------------------------------------
      SELECT CASE(eta_case)
c-----------------------------------------------------------------------
c     set eta = (spitzer resistivity + chodura resistivity)
c-----------------------------------------------------------------------
      CASE("spitzer-chodura","sc+r-layer")
         eta_local=etas_norm/Te_local**1.5
     $        + chod_const*etac_norm/SQRT(rho)
     $        *(one - EXP(-v_chod_norm*ABS(j)
     $        /(fc*rho*SQRT(gamma*Te_local))))
         eta_local = MIN(eta_local,etavac)
         IF(eta_case=="sc+r-layer")THEN
            WHERE(x < r_eta) eta_local = MAX(eta_local,
     $           etavac*half*(one-COS(pi*(one-x/r_eta))))
         ENDIF
      CASE("IAT-anomalous")
         eta_local = etas_norm/Te_local**1.5
         eta_local = MIN(eta_local,etavac)

         j_c = jc_norm*rho*SQRT(Te_local)
         dj = 0.05*j_c
         eta1_anm = eta_anm_norm*(Te_local/Ti_local)/SQRT(rho) 
         deta_dj =  eta1_anm/j_c

         WHERE (j > j_c .AND. j < (j_c + dj))
            coeff(1,:,:) = -(deta_dj*(-dj) - 2.*eta_local + 2.*eta1_anm)
     $           /dj**3 

            coeff(2,:,:) = -(-2.*j_c**2*deta_dj + j_c*(3.*eta_local 
     $           + (j_c+dj)*deta_dj - 3.*eta1_anm) 
     $           + (j_c+dj)*(3.*eta_local + (j_c+dj)*deta_dj 
     $           - 3.*eta1_anm))/dj**3

            coeff(3,:,:) = - j_c*(deta_dj*j_c**2 + j_c*(j_c+dj)*deta_dj 
     $           - 6.*eta_local*(j_c+dj) - 2.*deta_dj*(j_c+dj)**2 
     $           + 6.*(j_c+dj)*eta1_anm)/dj**3

            coeff(4,:,:) = -(j_c**3*(eta1_anm - (j_c+dj)*deta_dj) 
     $           + j_c**2*(j_c+dj)*((j_c+dj)*deta_dj - 3.*eta1_anm)
     $           + 3.*j_c*eta_local*(j_c+dj)**2 - eta_local*(j_c+dj)**3)
     $           /dj**3

            eta_local = coeff(1,:,:)*j**3 + coeff(2,:,:)*j**2 
     $           + coeff(3,:,:)*j + coeff(4,:,:)
            
         ELSEWHERE (j > j_c + dj)
            eta_local = deta_dj*j
         END WHERE
c-----------------------------------------------------------------------
c     define distance r_eta s.t. resistivity rises from "interior" 
c     value eta to "wall" value etavac for x < r_eta
c-----------------------------------------------------------------------
      CASE("x-layer")
         WHERE(x < r_eta)eta_local = eta + 
     $        (etavac - eta)*half*(one-COS(pi*(one-x/r_eta)))
      CASE("uniform")
         CONTINUE
      CASE DEFAULT
         CALL PROGRAM_STOP("transport_seteta: cannot recognize "
     $        //"eta case = "//TRIM(eta_case))
      END SELECT
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rho <= 0)eta_local=zero
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_seteta
c-----------------------------------------------------------------------
c     subprogram 2. transport_seteta_u.
c     sets resistivity and d/du(resistivity).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_seteta_u(eta_case,x,rho,p_e,p_i,j,jc_norm,
     $     eta_anm_norm,chod_const,
     $     etas_norm,etac_norm,v_chod_norm,r_eta,eta,etavac,
     $     eta_local,eta_rho,eta_p,eta_pi,eta_j)
      
      CHARACTER(*), INTENT(IN) :: eta_case
      REAL(r8), INTENT(IN) :: etas_norm,etac_norm,v_chod_norm,r_eta,eta,
     $     etavac,chod_const, jc_norm, eta_anm_norm
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,rho,p_e,p_i,j
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: eta_local,eta_rho,eta_p,
     $     eta_j, eta_pi

      REAL(r8), DIMENSION(4,SIZE(rho,1),SIZE(rho,2)) :: coeff
      REAL(r8), DIMENSION(SIZE(rho,1),SIZE(rho,2)) :: exp_chod,Te_local,
     $     Ti_local, j_c, dj, deta_dj
c-----------------------------------------------------------------------
c     initialize values.
c-----------------------------------------------------------------------
      eta_rho=zero
      eta_p=zero
      eta_j=zero
      eta_pi = zero
      Te_local = MAX(p_e/rho,eps)
      Ti_local = MAX(p_i/rho,eps)
c-----------------------------------------------------------------------
c     set resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta(eta_case,x,rho,p_e,p_i,j,jc_norm,
     $     eta_anm_norm,chod_const,etas_norm,etac_norm,v_chod_norm,
     $     r_eta,eta,etavac,coeff,eta_local)
c-----------------------------------------------------------------------
c     set d/du(resistivity).
c-----------------------------------------------------------------------
      SELECT CASE(eta_case)
      CASE("spitzer-chodura","sc+r-layer")
         exp_chod=EXP(-v_chod_norm*ABS(j)/(fc*rho*SQRT(gamma*Te_local)))

         eta_rho = 1.5_r8*etas_norm/(rho*Te_local**1.5)
     $        + chod_const*etac_norm*half/rho**1.5
     $        *(-one + exp_chod - v_chod_norm*ABS(j)*exp_chod
     $        /(fc*rho*SQRT(gamma*Te_local)))
         eta_p = -1.5_r8*etas_norm/(rho*Te_local**2.5)
     $        - chod_const*etac_norm*v_chod_norm*ABS(j)
     $        /(two*fc*rho**2.5*SQRT(gamma)*Te_local**1.5)*exp_chod
         eta_j = chod_const*etac_norm*v_chod_norm*SIGN(exp_chod,j)
     $        /(fc*rho**1.5*SQRT(gamma*Te_local))
c-----------------------------------------------------------------------
c     where eta_local is overwritten, set derivatives to zero.
c-----------------------------------------------------------------------
         WHERE(eta_local == etavac)
            eta_rho=zero
            eta_p=zero
            eta_j=zero
         ELSEWHERE(p_e/rho < eps)
            eta_rho = chod_const*etac_norm*half/rho**1.5
     $           *(-one + exp_chod - two*v_chod_norm*ABS(j)*exp_chod
     $           /(fc*rho*SQRT(gamma*Te_local)))
            eta_p = zero
         END WHERE

         IF(eta_case=="sc+r-layer")THEN
            WHERE(x < r_eta .AND.
     $           eta_local == etavac*half*(one-COS(pi*(one-x/r_eta))))
               eta_rho=zero
               eta_p=zero
               eta_j=zero
            END WHERE
         ENDIF
      CASE("IAT-anomalous")
         j_c = jc_norm*rho*SQRT(Te_local)
         dj = 0.05*j_c

         eta_rho = 1.5_r8*etas_norm/(rho*Te_local**1.5)
         eta_p = -1.5_r8*etas_norm/(rho*Te_local**2.5)
         eta_j = zero
         eta_pi = zero

         WHERE (j > j_c .AND. j < (j_c + dj))
            eta_j = 3.*coeff(1,:,:)*j**2 + 2.*coeff(2,:,:)*j 
     $           + coeff(3,:,:)
            eta_rho = zero
            eta_p = zero
            eta_pi = zero
         ELSEWHERE (j > j_c + dj)
            deta_dj =  eta_anm_norm*(Te_local/Ti_local)/SQRT(rho)/j_c
            eta_j = deta_dj
            eta_rho = - deta_dj*j/rho
            eta_p = 0.5*deta_dj*j/p_e
            eta_pi = - deta_dj*j/p_i
         END WHERE
      END SELECT
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         eta_local=zero
         eta_rho=zero
         eta_p=zero
         eta_j=zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_seteta_u
c-----------------------------------------------------------------------
c     subprogram 3. transport_setkaniso.
c     sets anisotropic thermal conductivity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_setkaniso(kappa_par,kappa_perp,Bsq,
     $     kperp,kfac)
      
      REAL(r8), INTENT(IN) :: kappa_perp,kappa_par
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: Bsq
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: kperp,kfac
c-----------------------------------------------------------------------
c     set thermal conductivity.
c-----------------------------------------------------------------------
      kperp=kappa_par*kappa_perp/(Bsq*kappa_par + kappa_perp)
      kfac=kappa_par-kperp
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_setkaniso
c-----------------------------------------------------------------------
c     subprogram 4. transport_setkaniso_u.
c     sets anisotropic thermal conductivity and 
c     d/du(thermal conductivity).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_setkaniso_u(kappa_par,kappa_perp,Bsq,
     $     kperp,kfac,kperp_bsq)

      REAL(r8), INTENT(IN) :: kappa_par,kappa_perp
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: Bsq
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: kperp,kfac,kperp_bsq
c-----------------------------------------------------------------------
c     set thermal conductivity and set d/du(thermal conductivity).
c-----------------------------------------------------------------------
      kperp=kappa_par*kappa_perp/(Bsq*kappa_par + kappa_perp)
      kfac=kappa_par-kperp
      
      kperp_bsq = -kappa_par**2*kappa_perp
     $     /(Bsq*kappa_par + kappa_perp)**2
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_setkaniso_u
c-----------------------------------------------------------------------
c     subprogram 5. transport_kbrag
c     sets perpendicular and parallel thermal conductivity per 
c     Braginskii.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_kbrag(rho,p,te_frac,Bsq,ke_norm,ki_norm,
     $     xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac)
      
      REAL(r8), INTENT(IN) :: ke_norm,ki_norm,xe_norm,xi_norm,te_frac,
     $     kappa_min,kappa_max
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,p,Bsq
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: kperp,kfac

      REAL(r8), PARAMETER :: ge0=11.92,ge1=4.664,de0=3.770,de1=14.79,
     $     gi0=2.645,gi1=2.0,di0=0.677,di1=2.7
      REAL(r8), DIMENSION(SIZE(rho,1),SIZE(rho,2)) :: Ti,Te,kpare,
     $     kperpe,kpari,kperpi,delta,x
c-----------------------------------------------------------------------
c     set thermal conductivity.
c-----------------------------------------------------------------------
      Ti = (one - te_frac)*MAX(p/rho,eps)
      Te = te_frac*MAX(p/rho,eps)

      kpare = ke_norm*Te**2.5*ge0/de0

      x = xe_norm*Te**1.5*SQRT(Bsq)/rho
      delta = x**4 + de1*x**2 + de0
      kperpe = ke_norm*Te**2.5*(ge1*x**2 + ge0)/delta

      kpari = ki_norm*Ti**2.5*gi0/di0

      x = xi_norm*Ti**1.5*SQRT(Bsq)/rho
      delta = x**4 + di1*x**2 + di0
      kperpi = ki_norm*Ti**2.5*(gi1*x**2 + gi0)/delta

      kperp = kperpe*te_frac + kperpi*(one-te_frac)
      kfac = kpare*te_frac + kpari*(one-te_frac)
c-----------------------------------------------------------------------
c     set minimum and maximum allowable perp & parallel heat conduction 
c-----------------------------------------------------------------------
      kperp = MIN(MAX(kperp, kappa_min), kappa_max)
      kfac = MIN(MAX(kfac, kappa_min), kappa_max) - kperp
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         kperp=zero
         kfac=zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_kbrag
c-----------------------------------------------------------------------
c     subprogram 6. transport_setkbrag_u.
c     sets perpendicular and parallel thermal conductivity and 
c     d/du(thermal conductivity) per Braginskii.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_kbrag_u(rho,p,te_frac,Bsq,ke_norm,ki_norm,
     $     xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac,kpar_rho,
     $     kpar_p,kperp_rho,kperp_p,kperp_bsq)

      REAL(r8), INTENT(IN) :: ke_norm,ki_norm,xe_norm,xi_norm,te_frac,
     $     kappa_min,kappa_max
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,p,Bsq
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: kperp,kfac,kpar_rho,
     $     kpar_p,kperp_rho,kperp_p,kperp_bsq
      
      REAL(r8), PARAMETER :: ge1=4.664,ge0=11.92,de0=3.770,de1=14.79,
     $     gi0=2.645,gi1=2.0,di0=0.677,di1=2.7
      REAL(r8), DIMENSION(SIZE(rho,1),SIZE(rho,2)) :: Ti,Te,invBsq,
     $     kpare,kperpe,kpari,kperpi,delta,x,fac_x,efac_rho,
     $     efac_p,efac_bsq,ifac_rho,ifac_p,ifac_bsq,p_inv
c-----------------------------------------------------------------------
c     set thermal conductivity and set d/du(thermal conductivity).
c-----------------------------------------------------------------------
      invBsq = one/Bsq

      Ti = (one - te_frac)*MAX(p/rho,eps)
      Te = te_frac*MAX(p/rho,eps)
      p_inv = one/(rho*MAX(p/rho,eps))

      kpare = ke_norm*Te**2.5*ge0/de0

      x = xe_norm*Te**1.5*SQRT(Bsq)/rho
      delta = x**4 + de1*x**2 + de0
      kperpe = ke_norm*Te**2.5*(ge1*x**2 + ge0)/delta

      fac_x = -two*x*(de1*ge0 - de0*ge1 + two*ge0*x**2 + ge1*x**4)
     $     /delta**2
      efac_rho = -2.5_r8*x*fac_x/rho
      WHERE(p/rho < eps)efac_rho = -x*fac_x/rho
      efac_p = 1.5_r8*x*p_inv*fac_x
      efac_bsq = half*fac_x*x*invBsq

      kpari = ki_norm*Ti**2.5*gi0/di0

      x = xi_norm*Ti**1.5*SQRT(Bsq)/rho
      delta = x**4 + di1*x**2 + di0
      kperpi = ki_norm*Ti**2.5*(gi1*x**2 + gi0)/delta

      fac_x = -two*x*(di1*gi0 - di0*gi1 + two*gi0*x**2 + gi1*x**4)
     $     /delta**2
      ifac_rho = -2.5_r8*x*fac_x/rho
      WHERE(p/rho < eps)ifac_rho = -x*fac_x/rho
      ifac_p = 1.5_r8*x*p_inv*fac_x
      ifac_bsq = half*fac_x*x*invBsq
      
      kperp = kperpe*te_frac + kperpi*(one-te_frac)
      kperp_rho = -2.5_r8*kperp/rho
     $     + te_frac*Te**2.5*ke_norm*efac_rho 
     $     + (one-te_frac)*Ti**2.5*ki_norm*ifac_rho
      kperp_p = 2.5_r8*p_inv*(te_frac*kperpe + (one-te_frac)*kperpi)
     $     + te_frac*Te**2.5*ke_norm*efac_p 
     $     + (one-te_frac)*Ti**2.5*ki_norm*ifac_p
      WHERE(p/rho < eps)
         kperp_rho = te_frac*Te**2.5*ke_norm*efac_rho 
     $        + (one-te_frac)*Ti**2.5*ki_norm*ifac_rho
         kperp_p = zero
      END WHERE
      kperp_bsq = te_frac*Te**2.5*ke_norm*efac_bsq 
     $     + (one-te_frac)*Ti**2.5*ki_norm*ifac_bsq
      WHERE(Bsq <= Bmin**2)kperp_bsq=zero

      kfac = kpare*te_frac + kpari*(one-te_frac)
      kpar_rho = -2.5_r8*kfac/rho
      kpar_p = 2.5_r8*p_inv*(te_frac*kpare + (one-te_frac)*kpari)
      WHERE(p/rho < eps)
         kpar_rho = zero
         kpar_p = zero
      END WHERE
c-----------------------------------------------------------------------
c     set minimum and maximum allowable perp & parallel heat conduction 
c-----------------------------------------------------------------------
      WHERE(kperp < kappa_min .OR. kperp > kappa_max)
         kperp_rho=zero
         kperp_p=zero
         kperp_bsq=zero
      END WHERE
      kperp = MIN(MAX(kperp, kappa_min), kappa_max)

      WHERE(kfac < kappa_min .OR. kfac > kappa_max)
         kpar_rho=zero
         kpar_p=zero
      END WHERE
      kfac = MIN(MAX(kfac, kappa_min), kappa_max) - kperp
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         kperp=zero
         kfac=zero
         kpar_rho=zero
         kpar_p=zero
         kperp_rho=zero
         kperp_p=zero
         kperp_bsq=zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_kbrag_u
c-----------------------------------------------------------------------
c     subprogram 7. transport_BdotT
c     sets b(b.(gradT)).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_BdotT(b1,b2,b3,Tx,Ty,BdotT)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: b1,b2,b3,Tx,Ty
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: BdotT

      REAL(r8), DIMENSION(SIZE(b1,1),SIZE(b1,2)) :: invBsq,Bsq
c-----------------------------------------------------------------------
c     set thermal conductivity.
c-----------------------------------------------------------------------
      Bsq = b1**2 + b2**2 + b3**2
      invBsq = one/Bsq

      BdotT(1,:,:) = b1*invBsq*(b1*Tx + b2*Ty)
      BdotT(2,:,:) = b2*invBsq*(b1*Tx + b2*Ty)
      BdotT(3,:,:) = b3*invBsq*(b1*Tx + b2*Ty)

      WHERE(Bsq <= Bmin**2)
         BdotT(1,:,:) = zero
         BdotT(2,:,:) = zero
         BdotT(3,:,:) = zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_BdotT
c-----------------------------------------------------------------------
c     subprogram 8. transport_BdotT_u.
c     sets d(b(b.(gradT)))/du
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_BdotT_u(b1,b2,b3,Tx,Ty,
     $     BdotT,BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: b1,b2,b3,Tx,Ty
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: BdotT,
     $     BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(b1,1),SIZE(b1,2)) :: invBsq,Bsq
c-----------------------------------------------------------------------
c     set derivatives of BdotT.
c-----------------------------------------------------------------------
      CALL transport_BdotT(b1,b2,b3,Tx,Ty,BdotT)
      Bsq = b1**2 + b2**2 + b3**2
      invBsq = one/Bsq
      
      BdotT_b1(1,:,:) = invBsq**2*(two*b1*(b2**2 + b3**2)*Tx 
     $     + b2*(-b1**2 + b2**2 + b3**2)*Ty)
      BdotT_b2(1,:,:) = -b1*invBsq**2
     $     *(two*b2*b1*Tx + (b2**2  - b1**2 - b3**2)*Ty)
      BdotT_b3(1,:,:) = -two*b3*b1*invBsq**2*(b1*Tx + b2*Ty)
      BdotT_Tx(1,:,:) = b1**2*invBsq
      BdotT_Ty(1,:,:) = b2*b1*invBsq

      BdotT_b1(2,:,:) = -b2*invBsq**2
     $     *(two*b2*b1*Ty + (b1**2 - b2**2 - b3**2)*Tx)
      BdotT_b2(2,:,:) = invBsq**2*(two*b2*(b1**2 + b3**2)*Ty
     $     + b1*(b1**2 - b2**2 + b3**2)*Tx)
      BdotT_b3(2,:,:) = -two*b3*b2*invBsq**2*(b2*Ty + b1*Tx)
      BdotT_Tx(2,:,:) = b2*b1*invBsq
      BdotT_Ty(2,:,:) = b2**2*invBsq
      
      BdotT_b1(3,:,:) = -b3*invBsq**2
     $     *(two*b2*b1*Ty + (b1**2 - b2**2 - b3**2)*Tx)
      BdotT_b2(3,:,:) = -b3*invBsq**2
     $     *(two*b2*b1*Tx + (b2**2  - b1**2 - b3**2)*Ty)
      BdotT_b3(3,:,:) = invBsq**2*(b1*Tx + b2*Ty)
     $     *(b1**2 + b2**2 - b3**2)
      BdotT_Tx(3,:,:) = b3*b1*invBsq
      BdotT_Ty(3,:,:) = b3*b2*invBsq

      DO i=1,SIZE(BdotT,1)
         WHERE(Bsq <= Bmin**2)
            BdotT(i,:,:)=zero
            BdotT_b1(i,:,:)=zero
            BdotT_b2(i,:,:)=zero
            BdotT_b3(i,:,:)=zero
            BdotT_Tx(i,:,:)=zero
            BdotT_Ty(i,:,:)=zero
         END WHERE
      ENDDO
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_BdotT_u
c-----------------------------------------------------------------------
c     subprogram 9. transport_mubrag.
c     sets isotropic viscosity per Braginskii assuming that ion 
c     viscosity is dominant and ions are not magnetized.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_mubrag(rho,p,ti_frac,mu_norm,mu_min,r_visc,r,
     $     visc)
      
      REAL(r8), INTENT(IN) :: mu_norm,ti_frac,mu_min,r_visc
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,p,r
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: visc

      REAL(r8), DIMENSION(SIZE(rho,1),SIZE(rho,2)) :: Ti
c-----------------------------------------------------------------------
c     set thermal viscosity.
c-----------------------------------------------------------------------
      Ti = ti_frac*MAX(p/rho,eps)
      visc = mu_norm*Ti**2.5
      visc = MAX(visc,mu_min)

      WHERE(r > r_visc)visc = visc 
     $     + 10*mu_min*half*(one - COS((r-r_visc)*pi/(one-r_visc)))
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         visc = zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_mubrag
c-----------------------------------------------------------------------
c     subprogram 10. transport_mubrag_u.
c     sets viscosity and d/du(viscoscity) per Braginskii.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_mubrag_u(rho,p,ti_frac,mu_norm,mu_min,r_visc,
     $     r,visc,visc_rho,visc_p)

      REAL(r8), INTENT(IN) :: mu_norm,ti_frac,mu_min,r_visc
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,p,r
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: visc,visc_rho,visc_p

      REAL(r8), DIMENSION(SIZE(rho,1),SIZE(rho,2)) :: Ti
c-----------------------------------------------------------------------
c     set viscosity and d/du(viscosity).
c-----------------------------------------------------------------------
      Ti = ti_frac*MAX(p/rho,eps)
      visc = mu_norm*Ti**2.5
      
      visc_rho = -2.5*visc/rho
      visc_p = 2.5*ti_frac*mu_norm*Ti**1.5/rho

      WHERE(visc < mu_min .OR. (p/rho) < eps)
         visc_rho = zero
         visc_p = zero
      END WHERE
      visc = MAX(visc,mu_min)

      WHERE(r > r_visc)visc = visc
     $     + 10*mu_min*half*(one - COS((r-r_visc)*pi/(one-r_visc)))
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         visc = zero
         visc_rho = zero
         visc_p = zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_mubrag_u
c-----------------------------------------------------------------------
c     subprogram 11. transport_mun.
c     sets isotropic viscosity using hard sphere model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_mun(rhon,pn,mun_norm,viscn)
      
      REAL(r8), INTENT(IN) :: mun_norm
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rhon,pn
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: viscn

      REAL(r8), DIMENSION(SIZE(rhon,1),SIZE(rhon,2)) :: Tn
c-----------------------------------------------------------------------
c     set viscosity.
c-----------------------------------------------------------------------
      Tn = MAX(pn/rhon,eps)
      viscn = mun_norm*SQRT(Tn)
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rhon <= 0)
         viscn = zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_mun
c-----------------------------------------------------------------------
c     subprogram 12. transport_mun_u.
c     sets viscosity and d/du(viscoscity) for hard sphere model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_mun_u(rhon,pn,mun_norm,viscn,viscn_rhon,
     $     viscn_pn)

      REAL(r8), INTENT(IN) :: mun_norm
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rhon,pn
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: viscn,viscn_rhon,viscn_pn

      REAL(r8), DIMENSION(SIZE(rhon,1),SIZE(rhon,2)) :: Tn
c-----------------------------------------------------------------------
c     set viscosity and d/du(viscosity).
c-----------------------------------------------------------------------
      Tn = MAX(pn/rhon,eps)
      viscn = mun_norm*SQRT(Tn)

      viscn_rhon = -half*viscn/rhon
      viscn_pn = half*mun_norm/(SQRT(Tn)*rhon)

      WHERE(pn/rhon < eps)
         viscn_rhon = zero
         viscn_pn = zero
      END WHERE
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rhon <= 0)
         viscn = zero
         viscn_rhon = zero
         viscn_pn = zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_mun_u
c-----------------------------------------------------------------------
c     subprogram 13. transport_hexch.
c     sets electron-ion heat exchange.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_hexch(ieheat,kappa_e,rho,pi,pe,heat_exch)
      
      REAL(r8), INTENT(IN) :: ieheat,kappa_e
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,pi,pe
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: heat_exch

      REAL(r8), PARAMETER :: max_exch=1._r8
      REAL(r8), DIMENSION(SIZE(rho,1),SIZE(rho,2)) :: Ti,Te
c-----------------------------------------------------------------------
c     initialize values.
c-----------------------------------------------------------------------
      Ti = MAX(pi/rho,eps)
      Te = MAX(pe/rho,eps)
c-----------------------------------------------------------------------
c     set heat exchange.
c-----------------------------------------------------------------------
      heat_exch = ieheat*rho**2*(Te - Ti)/(kappa_e*Te**1.5)
      heat_exch = SIGN(MIN(heat_exch,max_exch),(Te-Ti))
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         heat_exch = zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_hexch
c-----------------------------------------------------------------------
c     subprogram 14. transport_hexch_u.
c     sets electron-ion heat exchange jacobian
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_hexch_u(ieheat,kappa_e,rho,pi,pe,heat_exch,
     $     hexch_rho,hexch_pi,hexch_pe)
      
      REAL(r8), INTENT(IN) :: ieheat,kappa_e
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,pi,pe
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: heat_exch,hexch_rho,
     $     hexch_pi,hexch_pe

      REAL(r8), PARAMETER :: max_exch=1._r8
      REAL(r8), DIMENSION(SIZE(rho,1),SIZE(rho,2)) :: Ti,Te
c-----------------------------------------------------------------------
c     initialize values.
c-----------------------------------------------------------------------
      Ti = MAX(pi/rho,eps)
      Te = MAX(pe/rho,eps)
c-----------------------------------------------------------------------
c     set heat exchange.
c-----------------------------------------------------------------------
      CALL transport_hexch(ieheat,kappa_e,rho,pi,pe,heat_exch)
c-----------------------------------------------------------------------
c     set d/du.
c-----------------------------------------------------------------------
      hexch_rho=2.5_r8*heat_exch/rho
      hexch_pi=-ieheat*rho/(kappa_e*Te**1.5)
      hexch_pe=half*ieheat*rho*(3._r8*Ti - Te)/(kappa_e*Te**2.5)
c-----------------------------------------------------------------------
c     where necessary, set derivatives to zero.
c-----------------------------------------------------------------------
      WHERE(pi/rho < eps .AND. pe/rho < eps)
         hexch_rho = zero
         hexch_pi = zero
         hexch_pe = zero
      ELSEWHERE(pi/rho < eps)
         hexch_rho = ieheat*rho*(2.5_r8 - 3.5_r8*Ti/Te)
     $        /(kappa_e*SQRT(Te))
         hexch_pi = zero
      ELSEWHERE(pe/rho < eps)
         hexch_rho = ieheat*rho*(two - Ti/Te)/(kappa_e*SQRT(Te))
         hexch_pe = zero
      END WHERE

      WHERE(ABS(heat_exch) == max_exch)
         hexch_rho = zero
         hexch_pi = zero
         hexch_pe = zero
      END WHERE
c-----------------------------------------------------------------------
c     if density is non-positive, set transport coefficients to 0
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         heat_exch = zero
         hexch_rho = zero
         hexch_pi = zero
         hexch_pe = zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_hexch_u
c-----------------------------------------------------------------------
c     subprogram 15. transport_radloss.
c     sets energy loss term due to radiation 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_radloss(rad_fac,T0,rho,pe,alpha,rad_loss)

      REAL(r8), INTENT(IN) :: rad_fac, T0
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho, pe
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: alpha, rad_loss 

      REAL(r8), PARAMETER :: mp=1.673e-27 
      REAL(r8), DIMENSION(6) :: T_range
      REAL(r8), DIMENSION(7) :: hi, power
      REAL(r8), DIMENSION(SIZE(rho,1), SIZE(rho,2)) :: Te
c-----------------------------------------------------------------------
c     initialize values.
c-----------------------------------------------------------------------
      alpha = 0.0
      rad_loss = 0.0
      Te = MAX(pe/rho, eps)
c-----------------------------------------------------------------------
c     initialize prametrization from Klimchuk et al. 2008
c-----------------------------------------------------------------------      
       T_range = (/9.332e4_r8, 4.677e5_r8, 1.513e6_r8, 3.548e6_r8, 
     $            7.943e6_r8, 4.265e7_r8/)

       power = (/2.0, -1.0, 0.0, -1.5, 0.33, -1.0, 0.5 /)

       hi = (/1.09e-31_r8, 8.87e-17_r8, 1.9e-22_r8, 3.53e-13_r8, 
     $        3.46e-25_r8, 5.49e-16_r8, 1.96e-27_r8/)

      T_range = T_range/T0
       hi  = rad_fac * (T0**power) * hi * 1.e-13_r8

c-----------------------------------------------------------------------
c     set radiative loss function
c-----------------------------------------------------------------------
      WHERE(Te <= T_range(1))
         rad_loss = hi(1)*(Te**power(1))
         alpha = power(1)
      ELSEWHERE(T_range(1)<Te .AND. Te<= T_range(2))
         rad_loss = hi(2)*(Te**power(2))
         alpha = power(2)
      ELSEWHERE(T_range(2)<Te .AND. Te<=T_range(3))
         rad_loss = hi(3)*(Te**power(3))
         alpha = power(3)
      ELSEWHERE(T_range(3)<Te .AND. Te<=T_range(4))
         rad_loss = hi(4)*(Te**power(4))
         alpha = power(4)
      ELSEWHERE(T_range(4)<Te .AND. Te<=T_range(5))
         rad_loss = hi(5)*(Te**power(5))
         alpha = power(5)
      ELSEWHERE(T_range(5)<Te .AND. Te<=T_range(6))
         rad_loss = hi(6)*(Te**power(6))
         alpha = power(6)
      ELSEWHERE(T_range(6)<Te )
         rad_loss = hi(7)*(Te**power(7))
         alpha = power(7)
      END WHERE
 
      rad_loss = rad_loss * rho * rho
c-----------------------------------------------------------------------
c     if density is non-positive, set rad_loss to zero
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         rad_loss = zero
      END WHERE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE transport_radloss
c-----------------------------------------------------------------------
c     subprogram 16. transport_radloss_u.
c     sets radiative loss term jacobian
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transport_radloss_u(rad_fac,T0,rho,pe,rad_loss,
     $           rloss_rho,rloss_pe)

      REAL(r8), INTENT(IN) :: rad_fac, T0
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho, pe
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: rad_loss,
     $                          rloss_rho,rloss_pe
      REAL(r8), DIMENSION(SIZE(rho,1), SIZE(rho,2)) :: Te, alpha 
c-----------------------------------------------------------------------
c     initialize values.
c-----------------------------------------------------------------------
      rad_loss = 0.0
      rloss_rho = 0.0
      rloss_pe = 0.0

      Te = MAX(pe/rho, eps)

      CALL transport_radloss(rad_fac,T0,rho,pe,alpha,rad_loss) 
c-----------------------------------------------------------------------
c     set drdu
c-----------------------------------------------------------------------
      rloss_rho = (2.-alpha) * rad_loss/rho  
      rloss_pe = alpha * rad_loss/pe  
 
c-----------------------------------------------------------------------
c     if density is non-positive, set to zero
c-----------------------------------------------------------------------
      WHERE(rho <= 0)
         rad_loss = zero
         rloss_rho = zero
         rloss_pe = zero
      END WHERE
      
      RETURN
      END SUBROUTINE transport_radloss_u
      END MODULE transport_mod     


