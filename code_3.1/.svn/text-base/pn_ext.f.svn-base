c-----------------------------------------------------------------------
c     file pn_ext.f.
c     contains specifications for a reacting plasma-neutral model
c     with extended MHD effects in Ohm's Law (Hall and diamag. terms).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. pn_ext_mod.
c     1. pn_ext_equil.
c     2. pn_ext_grq.
c     3. pn_ext_grq_u.
c     4. pn_ext_bc_jac.
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
c     n. physics_schur
c     o. physics_dealloc
c     p. physics_main.
c-----------------------------------------------------------------------
c     subprogram 0. pn_ext_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE pn_ext_mod
      USE extra_mod
      USE transport_mod
      IMPLICIT NONE

      LOGICAL, PARAMETER :: ion_momentum=.TRUE.
      LOGICAL :: source=.FALSE.,cylinder=.FALSE.
      CHARACTER(16) :: init_type=".",kappa_case=".",eta_case=".",
     $     atom="hydrogen"
      INTEGER :: cyl_fac=0
      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,qe=1.602e-19,lnlam=10.,
     $     me=9.109e-31,ep0=8.854e-12,mu0=4.e-7*pi,mp=1.673e-27,
     $     gsun=2.74e2,gearth=9.81
      REAL(r8), PARAMETER :: eps_loc=1.e-12_r8
      REAL(r8) :: di=0.,nu=0.,eta=0.,r_eta=1.e10,etavac=1.,
     $     v_chod_norm=1.,etac_norm=1.,etas_norm=1.,mu=0.,kappa_par=0.,
     $     kappa_perp=0.,Dn=0.,beta0=0.,beta_e=0.,rhomin=0.,pmin=0.,
     $     pmax=1.,Tmin=1.,n0=1.,b0=1.,L0=1.,lx=0.,ly=0.,te_frac=0.,
     $     ti_frac=0.,lambda_psi=0.,lambda_phi=0.,kx=0.,ky=0.,delta=0.,
     $     deltatau=1.,gamma_fac=1.,gr_curve=0.,cc=.1,di_fac=1.,
     $     ion_fac=1.,recomb_fac=1.,cx_fac=1.,rin_fac=1.,kapn_fac=1.,
     $     civ_fac=0.,mfp_fac=0.,ion_norm=1.,recomb_norm=1.,cx_norm=1.,
     $     phi_ion=0.,phi_eff=0.,Top_norm=1.,initrho=1.,initrhon=1.,
     $     initpn=1.,v_civ=0.,ke_norm=1.,xe_norm=1.,ki_norm=1.,
     $     xi_norm=1.,kappa_min=0.,kappa_max=1.e10,cx_c1,cx_c2,cx_c3,
     $     a_vor,p_vor,x_vor,k_vor,mu_sv=0.,pert_frac=1.,equil_flow=1.,
     $     r_in_norm=1.,r_nn_norm=1.,r_en_norm=1.,gravity=0.,bplane=1.,
     $     bguide=0.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. pn_ext_equil.
c     sets equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pn_ext_equil(x,y,u)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      SELECT CASE(init_type)
      END SELECT 
cc-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pn_ext_equil
c-----------------------------------------------------------------------
c     subprogram 2. pn_ext_grq.
c     set reaction source rates (g), frictional momentum contributions
c     (r) and thermal energy transfer (q).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pn_ext_grq(p0,rho0,pn0,rhon0,vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: p0,pn0,rho0,rhon0
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: vi,vn
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: g_i,g_r,g_x,q_inx,q_nix
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: r_inx,r_nix

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(p0,1),SIZE(p0,2)) :: p,rho,pn,rhon,Top,
     $     poT,g_civ,dv,v_cx,vrep,vrepa,vrepb,s_cx
      REAL(r8), DIMENSION(SIZE(vi,1),SIZE(p0,1),SIZE(p0,2)) :: vin
c-----------------------------------------------------------------------
c     set local variables.
c-----------------------------------------------------------------------
      p=MAX(p0,eps_loc)
      rho=MAX(rho0,eps_loc)
      pn=MAX(pn0,eps_loc)
      rhon=MAX(rhon0,eps_loc)
c-----------------------------------------------------------------------
c     set reaction rates.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     note regarding critical ionization velocity (CIV):
c     ionization enhancement due to CIV occurs as if Te = 100 eV when 
c     the relative ion-neutral speed exceeds the critical ionization 
c     velocity,
c     v_civ = SQRT(2*phi_ion/mn).
c-----------------------------------------------------------------------
      dv = SQRT(SUM((vi - vn)**2,1))
      vin = vi - vn

      IF(civ_fac == 0)THEN
         g_civ = 0
      ELSE
         WHERE(dv > .75*v_civ .AND. dv < 1.25*v_civ)
            Top = 100./phi_ion*half
     $           *(1 - COS(pi*(dv - .75*v_civ)/(0.5*v_civ)))
            poT = one/MAX(Top,eps_loc)
            g_civ = civ_fac*ion_norm*rho*rhon*EXP(-poT)*
     $           (one + p_vor*SQRT(poT))/(x_vor + poT)*poT**k_vor
            WHERE(Top <= 0) g_civ = 0
         ELSEWHERE(dv >= 1.25*v_civ)
            poT = phi_ion/100.
            g_civ = civ_fac*ion_norm*rho*rhon*EXP(-poT)*
     $           (one + p_vor*SQRT(poT))/(x_vor + poT)*poT**k_vor
         ELSEWHERE
            g_civ = 0
         END WHERE
      ENDIF

      poT = rho/(Top_norm*te_frac*p)
      g_i = g_civ + ion_norm*rho*rhon*EXP(-poT)*poT**k_vor
     $     *(one + p_vor*SQRT(poT))/(x_vor + poT)
      g_r = recomb_norm*rho**2*SQRT(poT)

      v_cx = SQRT(MAX(8.*(ti_frac*p/rho + pn/rhon)/pi + dv**2,eps_loc))
      s_cx = MAX(-LOG(v_cx) + cx_c3,eps_loc)
      g_x = cx_norm*s_cx*v_cx*rhon*rho
c-----------------------------------------------------------------------
c     set frictional transfer for CX.
c-----------------------------------------------------------------------
      vrepa = two*pn/rhon
      vrepb = SQRT(MAX(4._r8*(8._r8*ti_frac*p/(rho*pi) + dv**2) 
     $     + 9._r8*pi*pn/(two*rhon),eps_loc))
      vrep = vrepa/vrepb
      DO i=1,3
         r_inx(i,:,:) = -cx_norm*s_cx*rho*rhon*vin(i,:,:)*vrep
      ENDDO

      vrepa = two*ti_frac*p/rho
      vrepb = SQRT(MAX(4._r8*(8._r8*pn/(rhon*pi) + dv**2) 
     $     + 9._r8*pi*ti_frac*p/(two*rho),eps_loc))
      vrep = vrepa/vrepb
      DO i=1,3
         r_nix(i,:,:) = cx_norm*s_cx*rho*rhon*vin(i,:,:)*vrep
      ENDDO
c-----------------------------------------------------------------------
c     set thermal energy transfer for CX.
c-----------------------------------------------------------------------
      vrepa = 3._r8*pn/(two*rhon)
      vrepb = SQRT(MAX(8._r8*ti_frac*p/(pi*rho) 
     $     + 128._r8*pn/(9._r8*pi*rhon) + dv**2,eps_loc))
      q_inx = cx_norm*s_cx*rho*rhon*vrepa*vrepb

      vrepa = 3._r8*ti_frac*p/(two*rho)
      vrepb = SQRT(MAX(8._r8*pn/(pi*rhon)
     $     + 128._r8*ti_frac*p/(9._r8*pi*rho) + dv**2,eps_loc))
      q_nix = cx_norm*s_cx*rho*rhon*vrepa*vrepb
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pn_ext_grq
c-----------------------------------------------------------------------
c     subprogram 3. pn_ext_grq_u.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pn_ext_grq_u(p0,rho0,pn0,rhon0,vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix,
     $     g_i_rho,g_i_rhon,g_i_p,g_r_rho,g_r_p,
     $     g_x_rho,g_x_rhon,g_x_p,g_x_pn,g_x_mi1,g_x_mi2,g_x_mi3,
     $     g_x_mn1,g_x_mn2,g_x_mn3,
     $     rinx_rho,rinx_p,rinx_rhon,rinx_pn,rinx_mi1,rinx_mi2,rinx_mi3,
     $     rinx_mn1,rinx_mn2,rinx_mn3,
     $     rnix_rho,rnix_p,rnix_rhon,rnix_pn,rnix_mi1,rnix_mi2,rnix_mi3,
     $     rnix_mn1,rnix_mn2,rnix_mn3,
     $     qinx_rho,qinx_p,qinx_rhon,qinx_pn,qinx_mi,qinx_mn,
     $     qnix_rho,qnix_p,qnix_rhon,qnix_pn,qnix_mi,qnix_mn)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: p0,rho0,pn0,rhon0
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: vi,vn
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: g_i,g_r,g_x,q_inx,q_nix,
     $     g_i_rho,g_i_rhon,g_i_p,g_r_rho,g_r_p,
     $     g_x_rho,g_x_rhon,g_x_p,g_x_pn,g_x_mi1,g_x_mi2,g_x_mi3,
     $     g_x_mn1,g_x_mn2,g_x_mn3,qinx_rho,qinx_p,qinx_rhon,qinx_pn,
     $     qnix_rho,qnix_p,qnix_rhon,qnix_pn
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: r_inx,r_nix,
     $     rinx_rho,rinx_p,rinx_rhon,rinx_pn,rinx_mi1,rinx_mi2,rinx_mi3,
     $     rinx_mn1,rinx_mn2,rinx_mn3,qinx_mi,qinx_mn,
     $     rnix_rho,rnix_p,rnix_rhon,rnix_pn,rnix_mi1,rnix_mi2,rnix_mi3,
     $     rnix_mn1,rnix_mn2,rnix_mn3,qnix_mi,qnix_mn

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(p0,1),SIZE(p0,2)) :: p,rho,pn,rhon,Top,
     $     poT,g_civ,dv,vcx,s_cx,g_x_vi1,g_x_vi2,g_x_vi3,
     $     v_cx,vrep,vrepa,vrepb,vrep_rho,vrep_p,vrep_rhon,vrep_pn,
     $     scx_p,scx_rho,scx_pn,scx_rhon,g_poT
      
      REAL(r8), DIMENSION(3,SIZE(p0,1),SIZE(p0,2)) :: vin,vrep_vi,scx_vi
c-----------------------------------------------------------------------
c     set local variables.
c-----------------------------------------------------------------------
      p=MAX(p0,eps_loc)
      rho=MAX(rho0,eps_loc)
      pn=MAX(pn0,eps_loc)
      rhon=MAX(rhon0,eps_loc)
c-----------------------------------------------------------------------
c     set reaction rates and derivatives of rates.
c-----------------------------------------------------------------------
      dv = SQRT(SUM((vi - vn)**2,1))
      vin = vi - vn

      IF(civ_fac == 0)THEN
         g_civ = 0
      ELSE
         WHERE(dv > .75*v_civ .AND. dv < 1.25*v_civ)
            Top = 100./phi_ion*half
     $           *(1 - COS(pi*(dv - .75*v_civ)/(0.5*v_civ)))
            poT = one/MAX(Top,eps_loc)
            g_civ = civ_fac*ion_norm*rho*rhon*EXP(-poT)*
     $           (one + p_vor*SQRT(poT))/(x_vor + poT)*poT**k_vor
            WHERE(Top <= 0) g_civ = 0
         ELSEWHERE(dv >= 1.25*v_civ)
            poT = phi_ion/100.
            g_civ = civ_fac*ion_norm*rho*rhon*EXP(-poT)*
     $           (one + p_vor*SQRT(poT))/(x_vor + poT)*poT**k_vor
         ELSEWHERE
            g_civ = 0
         END WHERE
      ENDIF

      poT = rho/(Top_norm*te_frac*p)

      CALL pn_ext_grq(p,rho,pn,rhon,vi,vn,g_i,g_r,g_x,r_inx,r_nix,
     $     q_inx,q_nix)
      
      g_poT = (g_i-g_civ)*(-one/(poT + x_vor) + k_vor/poT - one
     $     + p_vor*half/(SQRT(poT) + p_vor*poT))
      g_i_rho = (g_i + g_poT*poT)/rho
      g_i_rhon = g_i/rhon
      g_i_p = -g_poT*poT/p

      g_r_rho = 2.5_r8*g_r/rho
      g_r_p = -half*g_r/p

      v_cx = SQRT(MAX(8.*(ti_frac*p/rho + pn/rhon)/pi + dv**2,eps_loc))
      s_cx = MAX(-LOG(v_cx) + cx_c3,eps_loc)

      scx_p = -4._r8*ti_frac/(pi*rho*v_cx**2)
      scx_rho = (4._r8*ti_frac*p/(pi*rho) + vin(1,:,:)*vi(1,:,:) 
     $     + vin(2,:,:)*vi(2,:,:) + vin(3,:,:)*vi(3,:,:))/(v_cx**2*rho)

      scx_pn = -4._r8/(pi*rhon*v_cx**2)
      scx_rhon = (4._r8*pn/(pi*rhon) - vin(1,:,:)*vn(1,:,:) 
     $     - vin(2,:,:)*vn(2,:,:) - vin(3,:,:)*vn(3,:,:))/(v_cx**2*rhon)

      DO i=1,3
         scx_vi(i,:,:) = -vin(i,:,:)/v_cx**2
      ENDDO

      g_x_vi1 = g_x*(vin(1,:,:)/v_cx**2 + scx_vi(1,:,:)/s_cx)
      g_x_vi2 = g_x*(vin(2,:,:)/v_cx**2 + scx_vi(2,:,:)/s_cx)
      g_x_vi3 = g_x*(vin(3,:,:)/v_cx**2 + scx_vi(3,:,:)/s_cx)

      g_x_p = g_x*scx_p*(one - s_cx)/s_cx
      g_x_rho = g_x/rho + g_x*scx_rho*(one - s_cx)/s_cx
      g_x_pn = g_x*scx_pn*(one - s_cx)/s_cx
      g_x_rhon = g_x/rhon + g_x*scx_rhon*(one - s_cx)/s_cx
      g_x_mi1 = g_x_vi1/rho
      g_x_mi2 = g_x_vi2/rho
      g_x_mi3 = g_x_vi3/rho
      g_x_mn1 = -g_x_vi1/rhon
      g_x_mn2 = -g_x_vi2/rhon
      g_x_mn3 = -g_x_vi3/rhon
c-----------------------------------------------------------------------
c     set derivatives for CX frictional transfer, R_inx.
c-----------------------------------------------------------------------
      vrepa = two*pn/rhon
      vrepb = SQRT(MAX(4._r8*(8._r8*ti_frac*p/(rho*pi) + dv**2) 
     $     + 9._r8*pi*pn/(two*rhon),eps_loc))
      vrep = vrepa/vrepb

      DO i=1,3
         vrep_vi(i,:,:) = -4._r8*vrepa*vin(i,:,:)/vrepb**3
      ENDDO

      vrep_rho = 16._r8*vrepa*ti_frac*p/(pi*rho**2*vrepb**3)
     $     - SUM(vrep_vi*vi,1)/rho
      vrep_p = -16._r8*ti_frac*vrepa/(pi*vrepb**3*rho)

      vrep_rhon = 9._r8*pi*vrepa*pn/(4._r8*vrepb**3*rhon**2)
     $     + SUM(vrep_vi*vn,1)/rhon
     $     - two*pn/(rhon**2*vrepb)
      vrep_pn = -9._r8*pi*vrepa/(4._r8*vrepb**3*rhon) 
     $     + two/(rhon*vrepb)

      DO i=1,3
         rinx_rho(i,:,:) = r_inx(i,:,:)*(one/rho + scx_rho/s_cx
     $        + vrep_rho/vrep) + cx_norm*s_cx*rhon*vrep*vi(i,:,:)
         rinx_p(i,:,:) = r_inx(i,:,:)*(vrep_p/vrep + scx_p/s_cx)
         rinx_rhon(i,:,:) = r_inx(i,:,:)*(one/rhon + scx_rhon/s_cx
     $        + vrep_rhon/vrep) - cx_norm*s_cx*rho*vrep*vn(i,:,:)
         rinx_pn(i,:,:) = r_inx(i,:,:)*(vrep_pn/vrep + scx_pn/s_cx)

         rinx_mi1(i,:,:) = r_inx(i,:,:)/rho
     $        *(vrep_vi(1,:,:)/vrep + scx_vi(1,:,:)/s_cx)         
         rinx_mi2(i,:,:) = r_inx(i,:,:)/rho
     $        *(vrep_vi(2,:,:)/vrep + scx_vi(2,:,:)/s_cx)         
         rinx_mi3(i,:,:) = r_inx(i,:,:)/rho
     $        *(vrep_vi(3,:,:)/vrep + scx_vi(3,:,:)/s_cx)         

         rinx_mn1(i,:,:) = -r_inx(i,:,:)/rhon
     $        *(vrep_vi(1,:,:)/vrep + scx_vi(1,:,:)/s_cx)         
         rinx_mn2(i,:,:) = -r_inx(i,:,:)/rhon
     $        *(vrep_vi(2,:,:)/vrep + scx_vi(2,:,:)/s_cx)         
         rinx_mn3(i,:,:) = -r_inx(i,:,:)/rhon
     $        *(vrep_vi(3,:,:)/vrep + scx_vi(3,:,:)/s_cx)         
      ENDDO
      
      rinx_mi1(1,:,:) = rinx_mi1(1,:,:) - cx_norm*s_cx*rhon*vrep
      rinx_mi2(2,:,:) = rinx_mi2(2,:,:) - cx_norm*s_cx*rhon*vrep
      rinx_mi3(3,:,:) = rinx_mi3(3,:,:) - cx_norm*s_cx*rhon*vrep

      rinx_mn1(1,:,:) = rinx_mn1(1,:,:) + cx_norm*s_cx*rho*vrep
      rinx_mn2(2,:,:) = rinx_mn2(2,:,:) + cx_norm*s_cx*rho*vrep
      rinx_mn3(3,:,:) = rinx_mn3(3,:,:) + cx_norm*s_cx*rho*vrep
c-----------------------------------------------------------------------
c     set derivatives for CX frictional momentum transfer, R_nix.
c-----------------------------------------------------------------------
      vrepa = two*ti_frac*p/rho
      vrepb = SQRT(MAX(4._r8*(8._r8*pn/(rhon*pi) + dv**2) 
     $     + 9._r8*pi*ti_frac*p/(two*rho),eps_loc))
      vrep = vrepa/vrepb

      DO i=1,3
         vrep_vi(i,:,:) = -4._r8*vrepa*vin(i,:,:)/vrepb**3
      ENDDO

      vrep_rho = 9._r8*pi*vrepa*ti_frac*p/(4._r8*vrepb**3*rho**2)
     $     - SUM(vrep_vi*vi,1)/rho - two*ti_frac*p/(rho**2*vrepb)
      vrep_p = -9._r8*pi*vrepa*ti_frac/(4._r8*vrepb**3*rho)
     $     + two*ti_frac/(rho*vrepb)

      vrep_rhon = 16._r8*vrepa*pn/(pi*rhon**2*vrepb**3)
     $     + SUM(vrep_vi*vn,1)/rhon
      vrep_pn = -16._r8*vrepa/(pi*vrepb**3*rhon)

      DO i=1,3
         rnix_rho(i,:,:) = r_nix(i,:,:)*(one/rho + scx_rho/s_cx
     $        + vrep_rho/vrep) - cx_norm*s_cx*rhon*vrep*vi(i,:,:)
         rnix_p(i,:,:) = r_nix(i,:,:)*(vrep_p/vrep + scx_p/s_cx)
         rnix_rhon(i,:,:) = r_nix(i,:,:)*(one/rhon + scx_rhon/s_cx
     $        + vrep_rhon/vrep) + cx_norm*s_cx*rho*vrep*vn(i,:,:)
         rnix_pn(i,:,:) = r_nix(i,:,:)*(vrep_pn/vrep + scx_pn/s_cx)

         rnix_mi1(i,:,:) = r_nix(i,:,:)/rho
     $        *(vrep_vi(1,:,:)/vrep + scx_vi(1,:,:)/s_cx)         
         rnix_mi2(i,:,:) = r_nix(i,:,:)/rho
     $        *(vrep_vi(2,:,:)/vrep + scx_vi(2,:,:)/s_cx)         
         rnix_mi3(i,:,:) = r_nix(i,:,:)/rho
     $        *(vrep_vi(3,:,:)/vrep + scx_vi(3,:,:)/s_cx)         

         rnix_mn1(i,:,:) = -r_nix(i,:,:)/rhon
     $        *(vrep_vi(1,:,:)/vrep + scx_vi(1,:,:)/s_cx)         
         rnix_mn2(i,:,:) = -r_nix(i,:,:)/rhon
     $        *(vrep_vi(2,:,:)/vrep + scx_vi(2,:,:)/s_cx)         
         rnix_mn3(i,:,:) = -r_nix(i,:,:)/rhon
     $        *(vrep_vi(3,:,:)/vrep + scx_vi(3,:,:)/s_cx)         
      ENDDO
      
      rnix_mi1(1,:,:) = rnix_mi1(1,:,:) + cx_norm*s_cx*rhon*vrep
      rnix_mi2(2,:,:) = rnix_mi2(2,:,:) + cx_norm*s_cx*rhon*vrep
      rnix_mi3(3,:,:) = rnix_mi3(3,:,:) + cx_norm*s_cx*rhon*vrep

      rnix_mn1(1,:,:) = rnix_mn1(1,:,:) - cx_norm*s_cx*rho*vrep
      rnix_mn2(2,:,:) = rnix_mn2(2,:,:) - cx_norm*s_cx*rho*vrep
      rnix_mn3(3,:,:) = rnix_mn3(3,:,:) - cx_norm*s_cx*rho*vrep
c-----------------------------------------------------------------------
c     set derivatives for thermal energy  transfer, Q_inx.
c-----------------------------------------------------------------------
      vrepa = 3._r8*pn/(two*rhon)
      vrepb = SQRT(MAX(8._r8*ti_frac*p/(pi*rho) 
     $     + 128._r8*pn/(9._r8*pi*rhon) + dv**2,eps_loc))
      vrep = vrepa*vrepb

      DO i=1,3
         vrep_vi(i,:,:) = vrepa*vin(i,:,:)/vrepb
      ENDDO

      vrep_rho = -vrepa*4._r8*ti_frac*p/(pi*vrepb*rho**2)
     $     - SUM(vrep_vi*vi,1)/rho
      vrep_p = vrepa*4._r8*ti_frac/(pi*rho*vrepb)

      vrep_rhon = -vrepb*3._r8*pn/(two*rhon**2) 
     $     - vrepa*64._r8*pn/(9._r8*pi*rhon**2*vrepb)
     $     + SUM(vrep_vi*vn,1)/rhon
      vrep_pn = vrepb*3._r8/(two*rhon)
     $     + vrepa*64._r8/(9._r8*pi*rhon*vrepb)

      qinx_rho = q_inx*(one/rho + scx_rho/s_cx + vrep_rho/vrep)
      qinx_p = q_inx*(scx_p/s_cx + vrep_p/vrep)
      qinx_rhon = q_inx*(one/rhon + scx_rhon/s_cx + vrep_rhon/vrep)
      qinx_pn = q_inx*(scx_pn/s_cx + vrep_pn/vrep)
      
      DO i=1,3
         qinx_mi(i,:,:) = q_inx/rho
     $        *(scx_vi(i,:,:)/s_cx + vrep_vi(i,:,:)/vrep)
         qinx_mn(i,:,:) = -q_inx/rhon
     $        *(scx_vi(i,:,:)/s_cx + vrep_vi(i,:,:)/vrep)
      ENDDO
c-----------------------------------------------------------------------
c     set derivatives for thermal energy  transfer, Q_nix.
c-----------------------------------------------------------------------
      vrepa = 3._r8*ti_frac*p/(two*rho)
      vrepb = SQRT(MAX(8._r8*pn/(pi*rhon)
     $     + 128._r8*ti_frac*p/(9._r8*pi*rho) + dv**2,eps_loc))
      vrep = vrepa*vrepb

      DO i=1,3
         vrep_vi(i,:,:) = vrepa*vin(i,:,:)/vrepb
      ENDDO

      vrep_rho = -vrepb*3._r8*ti_frac*p/(two*rho**2)
     $     - vrepa*64._r8*ti_frac*p/(9._r8*pi*rho**2*vrepb)
     $     - SUM(vrep_vi*vi,1)/rho
      vrep_p = vrepb*3._r8*ti_frac/(two*rho)
     $     + vrepa*64._r8*ti_frac/(9._r8*pi*rho*vrepb)

      vrep_rhon = -vrepa*4._r8*pn/(pi*vrepb*rhon**2)
     $     + SUM(vrep_vi*vn,1)/rhon
      vrep_pn = vrepa*4._r8/(pi*rhon*vrepb)

      qnix_rho = q_nix*(one/rho + scx_rho/s_cx + vrep_rho/vrep)
      qnix_p = q_nix*(scx_p/s_cx + vrep_p/vrep)
      qnix_rhon = q_nix*(one/rhon + scx_rhon/s_cx + vrep_rhon/vrep)
      qnix_pn = q_nix*(scx_pn/s_cx + vrep_pn/vrep)
      
      DO i=1,3
         qnix_mi(i,:,:) = q_nix/rho
     $        *(scx_vi(i,:,:)/s_cx + vrep_vi(i,:,:)/vrep)
         qnix_mn(i,:,:) = -q_nix/rhon
     $        *(scx_vi(i,:,:)/s_cx + vrep_vi(i,:,:)/vrep)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pn_ext_grq_u
c-----------------------------------------------------------------------
c     subprogram 4. pn_ext_bc_jac.
c     computes jacobian for a boundary condition.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pn_ext_bc_jac(edge_type,t,x,y,nhat,
     $     u,ux,uy,uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy,uxx,uyy,uxy,nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy

      INTEGER :: i
      REAL(r8), PARAMETER :: du=1.e-6_r8
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: 
     $     u2,ux2,uy2,uxx2,uyy2,f,f2
c-----------------------------------------------------------------------
c     interface block
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_edge_rhs(lrtb,t,x,y,nhat,u,ux,uy,uxx,uyy,
     $        uxy,c)
         USE local_mod
         IMPLICIT NONE
         CHARACTER(*), INTENT(IN) :: lrtb
         REAL(r8), INTENT(IN) :: t
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
         REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,
     $        uxy
         REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c
         END SUBROUTINE physics_edge_rhs
      END INTERFACE
c-----------------------------------------------------------------------
c     compute the approximate boundary jacobian.
c-----------------------------------------------------------------------
      u2=u
      ux2=ux
      uy2=uy
      uxx2=uxx
      uyy2=uyy
      f = 0
      f2 = 0

      CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy,uxx,uyy,uxy,f)

      DO i=1,SIZE(u,1)
         u2(i,:,:) = u(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u2,ux,uy,uxx,uyy,
     $        uxy,f2)
         u2(i,:,:) = u(i,:,:)
         c_u(:,i,:,:) = (f2 - f)/du
      ENDDO

      DO i=1,SIZE(ux,1)
         ux2(i,:,:) = ux(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux2,uy,uxx,uyy,
     $        uxy,f2)
         ux2(i,:,:) = ux(i,:,:)
         c_ux(:,i,:,:) = (f2 - f)/du
      ENDDO

      DO i=1,SIZE(uy,1)
         uy2(i,:,:) = uy(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy2,uxx,uyy,
     $        uxy,f2)
         uy2(i,:,:) = uy(i,:,:)
         c_uy(:,i,:,:) = (f2 - f)/du
      ENDDO

      DO i=1,SIZE(uxx,1)
         uxx2(i,:,:) = uxx(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy,uxx2,uyy,
     $        uxy,f2)
         uxx2(i,:,:) = uxx(i,:,:)
         c_uxx(:,i,:,:) = (f2 - f)/du
      ENDDO

      DO i=1,SIZE(uyy,1)
         uyy2(i,:,:) = uyy(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy,uxx,uyy2,
     $        uxy,f2)
         uyy2(i,:,:) = uyy(i,:,:)
         c_uyy(:,i,:,:) = (f2 - f)/du
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pn_ext_bc_jac
      END MODULE pn_ext_mod
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
      USE pn_ext_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: physics_type
      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty,nqty_schur
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep,exit_flag
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax

      INTEGER :: myios
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/pn_ext_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,init_type,source,cylinder,eta,
     $     eta_case,r_eta,etavac,cc,mu,mu_sv,nu,kappa_par,kappa_perp,
     $     kappa_case,kappa_min,kappa_max,Dn,beta0,beta_e,rhomin,pmin,
     $     pmax,n0,b0,L0,lx,ly,lambda_psi,lambda_phi,delta,deltatau,
     $     atom,di_fac,ion_fac,recomb_fac,cx_fac,rin_fac,civ_fac,
     $     kapn_fac,mfp_fac,initrho,initrhon,initpn,pert_frac,
     $     equil_flow,bplane,bguide
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=pn_ext_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="pn_ext"
      nqty=13
      IF(nu /= 0)nqty=14
      nqty_schur=0
c-----------------------------------------------------------------------
c     set cylinder to true.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("GEM")
         cylinder=.FALSE.
         xperiodic=.FALSE.
         yperiodic=.FALSE.
      CASE("MRX","MRXpush")
         cylinder=.FALSE.
         xperiodic=.FALSE.
         yperiodic=.FALSE.
      CASE("forcefree")
         cylinder=.FALSE.
         xperiodic=.TRUE.
         yperiodic=.FALSE.
      CASE("RTsun")
         cylinder=.FALSE.
         xperiodic=.TRUE.
         yperiodic=.FALSE.
         di_fac=0
      CASE("RTearth")
         cylinder=.FALSE.
         xperiodic=.TRUE.
         yperiodic=.FALSE.
      END SELECT
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
      USE pn_ext_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty

      REAL(r8) :: tnorm,T0,mi
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(atom,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(eta_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(kappa_case,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(source,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(cylinder,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(r_eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(etavac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(cc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu_sv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_par,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_perp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Dn,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta_e,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rhomin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_min,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_max,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(L0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_psi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lambda_phi,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(delta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(deltatau,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ion_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(recomb_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(cx_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rin_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kapn_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mfp_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(civ_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initrho,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initrhon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initpn,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pert_frac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(equil_flow,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bplane,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bguide,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(7)=.TRUE.
      IF(nu /= 0)static(14)=.TRUE.
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      IF(cylinder)cyl_fac = 1
      SELECT CASE(atom)
      CASE("deuterium","hydrogen")
         IF(atom=="deuterium")mi = two*mp
         IF(atom=="hydrogen")mi = mp
         cx_c1 = 7.148e-20
         cx_c2 = 1.092e-18
         phi_ion = 13.6
         phi_eff = 33.*qe*n0/(b0**2/mu0)
         a_vor = 2.91e-14
         p_vor = 0.
         x_vor = .232
         k_vor = .39
      CASE("neon")
         mi = 3.351e-26
         cx_c1 = 5.646e-20
         cx_c2 = 7.954e-19
         phi_ion = 21.6
         phi_eff = 25.*qe*n0/(b0**2/mu0)
         a_vor = 1.5e-14
         p_vor = 1.
         x_vor = .0329
         k_vor = .43
      CASE DEFAULT
         CALL program_stop("physics_init_parameters cannot recognize "
     $        //"atom = "//TRIM(atom))
c-----------------------------------------------------------------------
c     Refer to G. S. Voronov, Atomic Data and Nuclear Tables, 65 (1997)
c     1-35 to determine constants (a_vor, x_vor, etc.) for your atom 
c     case.
c-----------------------------------------------------------------------
      END SELECT
c-----------------------------------------------------------------------
c     compute normalizations in MKS units.
c     for resistivity:
c     - v_chod_norm normalizes the ve/vs term in Chodura resistivity.
c     - etac_norm absorbs the constants in front of Chodura resistivity.
c     - etas_norm absorbs the constants in front of Spitzer resistivity.
c     for Braginskii thermal conduction:
c     - ke/ki_norm absorb constants in front of Braginskii
c     perpendicular and parallel thermal conduction.
c     - xe/xi_norm absorb constants in front of omega*tau. 
c     for viscosity:
c     for neutrals:
c     - ion_norm is tau_alfven/rho0 with additional factors as
c     required for ionization temperature dependence.
c     - recomb_norm is tau_alfven/rho0 with additional factors as
c     required for recombination temperature dependence.
c     - cx_norm is tau_alfven/rho0 with additional factors as
c     required for charge exchange temperature dependence.
c     - Top_norm is multiplied by normalized T to get Tev/phi_ion.
c     - cx_c3 is a constant required for normalizing sig_cx
c      (sig_cx = -c1*ln(v_cx)+ c2 = -c1*(ln(~v_cx) + ln(v0)) + c2
c       --> sig_cx = c1*~sig_cx, where ~sig_cx = ln(~v_cx) + c3,
c       and c3 = (c2 - c1*ln(v0))/c1,
c       and "~" refers to normalized quantities)
c     - r_in_norm normalizes the ion-neutral collision frequency
c         r_in_norm = 4/SQRT{pi}*n0*L0*sigma_in. sigma_in = 5.e-19 m^2
c     - r_nn_norm normalizes the neutral-neutral collision frequency
c         r_nn_norm = 4/SQRT{pi}*n0*L0*sigma_nn. sigma_nn = 5.e-19 m^2
c     - r_en_norm normalizes the electron-neutral collisional friction
c         with collisional cross-section sigma_en = 1.e-19 m^2
c-----------------------------------------------------------------------
      tnorm=L0*SQRT(mu0*n0*mi)/b0
      v_civ=SQRT(two*phi_ion*qe/mi)/(L0/tnorm)

      di = di_fac*2.28e8*SQRT(mi/(mp*n0))/L0
      T0 = b0**2/(n0*mu0*qe)
      v_chod_norm=SQRT(mi/(mu0*n0))/(qe*L0)
      etac_norm=me/(qe*L0*b0*SQRT(ep0*mu0))
      etas_norm=1.e-4*lnlam*SQRT(two*mi)*mu0*qe**1.5
     $     *n0**2/(L0*b0**4)
      ke_norm = 6.05e22/lnlam*T0**2.5*tnorm/(n0*L0**2)
      xe_norm = 6.05e22/lnlam*T0**1.5*b0/n0
      ki_norm = 3.35e-6/(lnlam*SQRT(mi*mp))*T0**2.5*tnorm/(n0*L0**2)
      xi_norm = 3.35e-6/(lnlam*SQRT(mi*mp))*T0**1.5*b0/n0
      ion_norm = ion_fac*a_vor*n0*tnorm
      recomb_norm = recomb_fac*2.6e-19*n0*tnorm/SQRT(phi_ion)
      cx_norm = cx_fac*L0*n0*7.15e-20
      cx_c3 = (cx_c2 - cx_c1*LOG(L0/tnorm))/cx_c1
      Top_norm = T0/phi_ion
      r_in_norm = rin_fac*5.e-19*4._r8*n0*L0/SQRT(pi)
      r_nn_norm = 5.e-19*4._r8*n0*L0/SQRT(pi)
      r_en_norm = 1.e-19*two*di*n0*L0*SQRT(me/(pi*mi))

      IF(1==2)THEN
         IF(mpi_rank==0)
     $        WRITE(*,*) "WARNING!! Overriding norms for FD test!"
         v_chod_norm = 1.1
         etac_norm = 1.2
         etas_norm = 1.3
         ke_norm = 1.4
         xe_norm = 1.5
         ki_norm = 1.6
         xi_norm = 1.7
         ion_norm = 2.1
         recomb_norm = 2.2
         cx_norm = 2.3
         r_in_norm = 2.4
         r_nn_norm = 2.5
         r_en_norm = 2.6
         di=2.7
      ENDIF

      gravity=0
      SELECT CASE(init_type)
      CASE("GEM","forcefree","MRX","MRXpush")
         mu = mu*ki_norm*(beta0 - beta_e)**2.5
         nu = nu*(di/initrho)**2*ke_norm*me/mp*beta_e**2.5
         IF(mpi_rank == 0)THEN
            WRITE(*,*)"di = ",di
            WRITE(*,*)"time_norm = ",tnorm
            WRITE(*,*)"Top_norm = ",Top_norm
            WRITE(*,*)"etas_norm = ",etas_norm
            WRITE(*,*)"ion_norm = ",ion_norm
            WRITE(*,*)"recomb_norm = ",recomb_norm
            WRITE(*,*)"cx_norm = ",cx_norm
            WRITE(*,*)"r_in_norm = ",r_in_norm
            WRITE(*,*)"r_nn_norm = ",r_nn_norm
            WRITE(*,*)"r_en_norm = ",r_en_norm
            WRITE(*,*)"mu = ",mu
            WRITE(*,*)"nu = ",nu
         ENDIF
      CASE("RTsun")
         beta_e=half*beta0
         gravity=gsun*tnorm**2/L0
         mu = mu*ki_norm*(beta0 - beta_e)**2.5
         IF(mpi_rank == 0)THEN
            WRITE(*,*)"time_norm = ",tnorm
            WRITE(*,*)"T0 = ",T0
            WRITE(*,*)"etas_norm = ",etas_norm
            WRITE(*,*)"ion_norm = ",ion_norm
            WRITE(*,*)"recomb_norm = ",recomb_norm
            WRITE(*,*)"cx_norm = ",cx_norm
            WRITE(*,*)"r_in_norm = ",r_in_norm
            WRITE(*,*)"r_nn_norm = ",r_nn_norm
            WRITE(*,*)"mu = ",mu
         ENDIF
      CASE("RTearth")
         beta_e=half*beta0
         gravity=gearth*tnorm**2/L0
         mu = mu*ki_norm*(beta0 - beta_e)**2.5
         nu = nu*di**2/(rhomin*initrho)*ke_norm*me/mp*beta_e**2.5
         IF(mpi_rank == 0)THEN
            WRITE(*,*)"di = ",di
            WRITE(*,*)"time_norm = ",tnorm
            WRITE(*,*)"T0 = ",T0
            WRITE(*,*)"etas_norm = ",etas_norm
            WRITE(*,*)"ion_norm = ",ion_norm
            WRITE(*,*)"recomb_norm = ",recomb_norm
            WRITE(*,*)"cx_norm = ",cx_norm
            WRITE(*,*)"r_in_norm = ",r_in_norm
            WRITE(*,*)"r_nn_norm = ",r_nn_norm
            WRITE(*,*)"r_en_norm = ",r_en_norm
            WRITE(*,*)"mu = ",mu
            WRITE(*,*)"nu = ",nu
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     set plasma compressibility parameters 
c     and adjust Spitzer resistivity constant
c-----------------------------------------------------------------------
      te_frac=beta_e/beta0
      ti_frac=one-te_frac
      IF(te_frac > 0)THEN
         etas_norm = etas_norm/te_frac**1.5
      ELSE
         etas_norm=0.
      ENDIF
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
      USE pn_ext_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u

      INTEGER :: i
      REAL(r8) :: lfac,ldelta
      REAL(r8), DIMENSION(SIZE(y,1),SIZE(y,2)) :: gauss,fchap,fchap2
c-----------------------------------------------------------------------
c     pn_ext initial conditions.
c-----------------------------------------------------------------------
      u = 0
      CALL pn_ext_equil(x,y,u)

      SELECT CASE(init_type)
      CASE("GEM")
         lfac=5
         ldelta=lfac*lambda_psi
         gauss = EXP(-((x/lfac)**2 + y**2)/lambda_psi**2)

         u(1,:,:) = initrho
     $        + half*pert_frac*bplane**2/(beta0*COSH(y/lambda_psi)**2)
         u(2,:,:) = -bplane*lambda_psi*LOG(COSH(y/lambda_psi))
ccc     $        + delta*bplane*gauss
ccc     $        *(one - half*EXP(-3._r8*y**2/lambda_psi**2))
         u(4,:,:) = delta*u(1,:,:)*SIN(x/ldelta)*gauss
     $        *(COS(y/lambda_psi) - two*y/lambda_psi*SIN(y/lambda_psi))
         u(5,:,:) = -delta*u(1,:,:)*SIN(y/lambda_psi)*gauss/lfac
     $        *(COS(x/ldelta) - two*x/ldelta*SIN(x/ldelta))
         u(7,:,:) = -bplane/(lambda_psi*COSH(y/lambda_psi)**2)
         u(8,:,:) = u(1,:,:)*beta0
         u(9,:,:) = initrhon*(one + half*(one-pert_frac)*bplane**2
     $        /(initpn*COSH(y/lambda_psi)**2))
         u(10,:,:) = delta*u(9,:,:)*SIN(x/ldelta)*gauss
     $        *(COS(y/lambda_psi) - two*y/lambda_psi*SIN(y/lambda_psi))
         u(11,:,:) = -delta*u(9,:,:)*SIN(y/lambda_psi)*gauss/lfac
     $        *(COS(x/ldelta) - two*x/ldelta*SIN(x/ldelta))
         u(13,:,:) = initpn
     $        + half*(one-pert_frac)*bplane**2/COSH(y/lambda_psi)**2

         IF(r_in_norm /= 0.)
     $        u(5,:,:) = u(5,:,:) + bplane**2*equil_flow*(pert_frac-one)
     $        *TANH(y/lambda_psi)/(lambda_psi*COSH(y/lambda_psi)**2)
     $        /(half*r_in_norm*u(9,:,:)
     $        *SQRT(half*(ti_frac*u(8,:,:)/u(1,:,:)
     $        + u(13,:,:)/u(9,:,:))))
      CASE("MRX")
         u(1,:,:) = initrho
     $        + half*bplane**2/(beta0*COSH(y/lambda_psi)**2)
         u(2,:,:) = bplane*lambda_psi*LOG(COSH(y/lambda_psi))
         u(7,:,:) = bplane/(lambda_psi*COSH(y/lambda_psi)**2)
         u(8,:,:) = u(1,:,:)*beta0
         u(9,:,:) = initrhon
         u(13,:,:) = initpn
      CASE("MRXpush")
         u(1,:,:) = initrho + COS(pi*x)**2*bplane**2/beta0
     $        *(0.625_r8*SIN(twopi*y)**2
     $        + half*(COS(twopi*y)**2
     $        - (SIN(twopi*y)/(twopi*lambda_psi*COSH(y/lambda_psi)**2)
     $        + COS(twopi*y)*TANH(y/lambda_psi))**2))
         u(2,:,:) = bplane/twopi
     $        *TANH(y/lambda_psi)*COS(pi*x)*SIN(twopi*y)
         u(7,:,:) = bplane*COS(pi*x)
     $        *(two*COS(twopi*y)/(lambda_psi*COSH(y/lambda_psi)**2)
     $        - pi*SIN(twopi*y)*TANH(y/lambda_psi)
     $        *(2.5_r8 + one/(pi*lambda_psi*COSH(y/lambda_psi))**2))
         u(8,:,:) = u(1,:,:)*beta0
         u(9,:,:) = initrhon
         u(13,:,:) = initpn
      CASE("forcefree")
         u(1,:,:) = initrho
         u(2,:,:) = bplane*lambda_psi*LOG(COSH(y/lambda_psi))
         u(3,:,:) = bplane/COSH(y/lambda_psi)
         u(7,:,:) = bplane/(lambda_psi*COSH(y/lambda_psi)**2)
         u(8,:,:) = u(1,:,:)*beta0
         u(9,:,:) = initrhon
         u(13,:,:) = initpn
      CASE("RTsun")
         lfac = one/20._r8
         fchap = initrho/COSH(two*y - one)**2
         fchap2 = lfac*COSH(y-half)**2/(lfac*COSH(y-half)**2 + one) 

         u(1,:,:) = EXP(-y*gravity/beta0)
         u(8,:,:) = beta0*u(1,:,:)*fchap2
         u(9,:,:) = rhomin + fchap
         u(13,:,:) = half*beta0*u(9,:,:)*fchap2

         u(3,:,:) = SQRT(one + beta0*((one-fchap2)*two*u(1,:,:)
     $        - u(9,:,:)*fchap2) 
     $        - gravity*initrho*(TANH(two*y-one) - one))

         fchap = 0
         DO i=1,5
            fchap =  fchap + 0.2*SIN(i*pi*x/lx)
         ENDDO
         u(9,:,:) = u(9,:,:)*(one + delta*EXP(-(two*y)**2)*fchap)
      CASE("RTearth")
         fchap = initrho/COSH(two*y - one)**2

         u(1,:,:) = rhomin + fchap
         u(8,:,:) = u(1,:,:)*beta0
         u(9,:,:) = EXP(-y*gravity*two/beta0)
         u(13,:,:) = u(9,:,:)*half*beta0

         u(3,:,:) = SQRT(one - two*beta0*fchap
     $        - gravity*initrho*(TANH(two*y-one) - one))

         fchap = 0
         DO i=1,5
            fchap =  fchap + 0.2*SIN(i*pi*x/lx)
         ENDDO
         u(1,:,:) = u(1,:,:)*(one + delta*EXP(-(two*y)**2)*fchap)
      END SELECT

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
      USE pn_ext_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     magnetic reconnection, boundary conditions.
c-----------------------------------------------------------------------
      top%bc_type="robin"
      top%static=.FALSE.
      bottom%bc_type="robin"
      bottom%static=.FALSE.
      left%bc_type="robin"
      left%static=.FALSE.
      right%bc_type="robin"
      right%static=.FALSE.

      SELECT CASE(init_type)
      CASE("GEM")
         edge_order=(/1,3,2,4/)
         bottom%bc_type="zeroflux"
         bottom%bc_type(3)="robin"
         bottom%static(3)=.TRUE.
         bottom%bc_type(5)="robin"
         bottom%static(5)=.TRUE.
         bottom%bc_type(11)="robin"
         bottom%static(11)=.TRUE.
         left%bc_type="zeroflux"
         left%bc_type(3)="robin"
         left%static(3)=.TRUE.
         left%bc_type(4)="robin"
         left%static(4)=.TRUE.
         left%bc_type(10)="robin"
         left%static(10)=.TRUE.
         right%bc_type="zeroflux"
         right%bc_type(3)="robin"
         right%static(3)=.TRUE.
         right%bc_type(4)="robin"
         right%static(4)=.TRUE.
         right%bc_type(10)="robin"
         right%static(10)=.TRUE.
c-----------------------------------------------------------------------
c     plasma density
c-----------------------------------------------------------------------
         top%bc_type(1)="zeroflux"
c-----------------------------------------------------------------------
c     poloidal magnetic flux
c-----------------------------------------------------------------------
         top%bc_type(2)="robin"
         top%static(2)=.FALSE.
c-----------------------------------------------------------------------
c     out of plane field
c-----------------------------------------------------------------------
         top%bc_type(3)="zeroflux"
c-----------------------------------------------------------------------
c     plasma momentum
c-----------------------------------------------------------------------
         top%bc_type(4:6)="robin"
         top%static(4:6)=.TRUE.
c-----------------------------------------------------------------------
c     current density
c-----------------------------------------------------------------------
         top%bc_type(7)="robin"
         top%static(7)=.TRUE.
c-----------------------------------------------------------------------
c     plasma pressure
c-----------------------------------------------------------------------
         top%bc_type(8)="zeroflux"
c-----------------------------------------------------------------------
c     neutral density
c-----------------------------------------------------------------------
         top%bc_type(9)="zeroflux"
c-----------------------------------------------------------------------
c     neutral momentum
c-----------------------------------------------------------------------
         top%bc_type(10:12)="robin"
         top%static(10:12)=.TRUE.
c-----------------------------------------------------------------------
c     neutral pressure
c-----------------------------------------------------------------------
         top%bc_type(13)="zeroflux"
c-----------------------------------------------------------------------
c     del^2(Bz) auxiliary variable
c-----------------------------------------------------------------------
         IF(nu /= 0)THEN
            top%bc_type(14)="robin"
            top%static(14)=.TRUE.
            bottom%bc_type(14)="robin"
            bottom%static(14)=.TRUE.
            left%bc_type(14)="robin"
            left%static(14)=.TRUE.
            right%bc_type(14)="robin"
            right%static(14)=.TRUE.
         ENDIF
      CASE("MRX")
         edge_order=(/1,3,2,4/)
         bottom%bc_type="zeroflux"
         bottom%bc_type(3)="robin"
         bottom%static(3)=.TRUE.
         bottom%bc_type(5)="robin"
         bottom%static(5)=.TRUE.
         bottom%bc_type(11)="robin"
         bottom%static(11)=.TRUE.
         left%bc_type="zeroflux"
         left%bc_type(3)="robin"
         left%static(3)=.TRUE.
         left%bc_type(4)="robin"
         left%static(4)=.TRUE.
         left%bc_type(10)="robin"
         left%static(10)=.TRUE.
         right%bc_type="zeroflux"
         right%bc_type(3)="robin"
         right%static(3)=.TRUE.
         right%bc_type(4)="robin"
         right%static(4)=.TRUE.
         right%bc_type(10)="robin"
         right%static(10)=.TRUE.
c-----------------------------------------------------------------------
c     plasma density
c-----------------------------------------------------------------------
         top%bc_type(1)="zeroflux"
c-----------------------------------------------------------------------
c     poloidal magnetic flux
c-----------------------------------------------------------------------
         top%bc_type(2)="robin"
         top%static(2)=.FALSE.
c-----------------------------------------------------------------------
c     out of plane field
c-----------------------------------------------------------------------
         top%bc_type(3)="zeroflux"
c-----------------------------------------------------------------------
c     plasma momentum
c-----------------------------------------------------------------------
         top%bc_type(4:6)="robin"
         top%static(4:6)=.TRUE.
c-----------------------------------------------------------------------
c     current density
c-----------------------------------------------------------------------
         top%bc_type(7)="robin"
         top%static(7)=.TRUE.
c-----------------------------------------------------------------------
c     plasma pressure
c-----------------------------------------------------------------------
         top%bc_type(8)="zeroflux"
c-----------------------------------------------------------------------
c     neutral density
c-----------------------------------------------------------------------
         top%bc_type(9)="zeroflux"
c-----------------------------------------------------------------------
c     neutral momentum
c-----------------------------------------------------------------------
         top%bc_type(10:12)="robin"
         top%static(10:12)=.TRUE.
c-----------------------------------------------------------------------
c     neutral pressure
c-----------------------------------------------------------------------
         top%bc_type(13)="zeroflux"
c-----------------------------------------------------------------------
c     del^2(Bz) auxiliary variable
c-----------------------------------------------------------------------
         IF(nu /= 0)THEN
            top%bc_type(14)="robin"
            top%static(14)=.TRUE.
            bottom%bc_type(14)="robin"
            bottom%static(14)=.TRUE.
            left%bc_type(14)="robin"
            left%static(14)=.TRUE.
            right%bc_type(14)="robin"
            right%static(14)=.TRUE.
         ENDIF
      CASE("MRXpush")
c-----------------------------------------------------------------------
c     plasma density
c-----------------------------------------------------------------------
         top%bc_type(1)="zeroflux"
         bottom%bc_type(1)="zeroflux"
         left%bc_type(1)="zeroflux"
         right%bc_type(1)="zeroflux"
c-----------------------------------------------------------------------
c     poloidal magnetic flux
c-----------------------------------------------------------------------
         top%bc_type(2)="robin"
         top%static(2)=.FALSE.
         bottom%bc_type(2)="robin"
         bottom%static(2)=.FALSE.
         left%bc_type(2)="robin"
         left%static(2)=.FALSE.
         right%bc_type(2)="robin"
         right%static(2)=.FALSE.
c-----------------------------------------------------------------------
c     out of plane field
c-----------------------------------------------------------------------
         top%bc_type(3)="robin"
         top%static(3)=.TRUE.
         bottom%bc_type(3)="robin"
         bottom%static(3)=.TRUE.
         left%bc_type(3)="robin"
         left%static(3)=.TRUE.
         right%bc_type(3)="robin"
         right%static(3)=.TRUE.
c-----------------------------------------------------------------------
c     plasma momentum
c-----------------------------------------------------------------------
         top%bc_type(4:6)="robin"
         top%static(4:6)=.TRUE.
         bottom%bc_type(4:6)="robin"
         bottom%static(4:6)=.TRUE.
         left%bc_type(4:6)="robin"
         left%static(4:6)=.TRUE.
         right%bc_type(4:6)="robin"
         right%static(4:6)=.TRUE.
c-----------------------------------------------------------------------
c     current density
c-----------------------------------------------------------------------
         top%bc_type(7)="robin"
         top%static(7)=.TRUE.
         bottom%bc_type(7)="robin"
         bottom%static(7)=.TRUE.
         left%bc_type(7)="robin"
         left%static(7)=.TRUE.
         right%bc_type(7)="robin"
         right%static(7)=.TRUE.
c-----------------------------------------------------------------------
c     plasma pressure
c-----------------------------------------------------------------------
         top%bc_type(8)="zeroflux"
         bottom%bc_type(8)="zeroflux"
         left%bc_type(8)="zeroflux"
         right%bc_type(8)="zeroflux"
c-----------------------------------------------------------------------
c     neutral density
c-----------------------------------------------------------------------
         top%bc_type(9)="zeroflux"
         bottom%bc_type(9)="zeroflux"
         left%bc_type(9)="zeroflux"
         right%bc_type(9)="zeroflux"
c-----------------------------------------------------------------------
c     neutral momentum
c-----------------------------------------------------------------------
         top%bc_type(10:12)="robin"
         top%static(10:12)=.TRUE.
         bottom%bc_type(10:12)="robin"
         bottom%static(10:12)=.TRUE.
         left%bc_type(10:12)="robin"
         left%static(10:12)=.TRUE.
         right%bc_type(10:12)="robin"
         right%static(10:12)=.TRUE.
c-----------------------------------------------------------------------
c     neutral pressure
c-----------------------------------------------------------------------
         top%bc_type(13)="zeroflux"
         bottom%bc_type(13)="zeroflux"
         left%bc_type(13)="zeroflux"
         right%bc_type(13)="zeroflux"
c-----------------------------------------------------------------------
c     del^2(Bz) auxiliary variable
c-----------------------------------------------------------------------
         IF(nu /= 0)THEN
            top%bc_type(14)="robin"
            top%static(14)=.TRUE.
            bottom%bc_type(14)="robin"
            bottom%static(14)=.TRUE.
            left%bc_type(14)="robin"
            left%static(14)=.TRUE.
            right%bc_type(14)="robin"
            right%static(14)=.TRUE.
         ENDIF
      CASE("forcefree")
c-----------------------------------------------------------------------
c     plasma density
c-----------------------------------------------------------------------
         top%bc_type(1)="zeroflux"
         bottom%bc_type(1)="zeroflux"
c-----------------------------------------------------------------------
c     poloidal magnetic flux
c-----------------------------------------------------------------------
         top%bc_type(2)="robin"
         top%static(2)=.FALSE.
         bottom%bc_type(2)="robin"
         bottom%static(2)=.FALSE.
c-----------------------------------------------------------------------
c     out of plane field
c-----------------------------------------------------------------------
         top%bc_type(3)="zeroflux"
         bottom%bc_type(3)="zeroflux"
c-----------------------------------------------------------------------
c     plasma momentum
c-----------------------------------------------------------------------
         top%bc_type(4:6)="robin"
         top%static(4:6)=.TRUE.
         bottom%bc_type(4:6)="robin"
         bottom%static(4:6)=.TRUE.
c-----------------------------------------------------------------------
c     current density
c-----------------------------------------------------------------------
         top%bc_type(7)="robin"
         top%static(7)=.TRUE.
         bottom%bc_type(7)="robin"
         bottom%static(7)=.TRUE.
c-----------------------------------------------------------------------
c     plasma pressure
c-----------------------------------------------------------------------
         top%bc_type(8)="zeroflux"
         bottom%bc_type(8)="zeroflux"
c-----------------------------------------------------------------------
c     neutral density
c-----------------------------------------------------------------------
         top%bc_type(9)="zeroflux"
         bottom%bc_type(9)="zeroflux"
c-----------------------------------------------------------------------
c     neutral momentum
c-----------------------------------------------------------------------
         top%bc_type(10:12)="robin"
         top%static(10:12)=.TRUE.
         bottom%bc_type(10:12)="robin"
         bottom%static(10:12)=.TRUE.
c-----------------------------------------------------------------------
c     neutral pressure
c-----------------------------------------------------------------------
         top%bc_type(13)="zeroflux"
         bottom%bc_type(13)="zeroflux"
c-----------------------------------------------------------------------
c     del^2(Bz) auxiliary variable
c-----------------------------------------------------------------------
         IF(nu /= 0)THEN
            top%bc_type(14)="robin"
            top%static(14)=.TRUE.
            bottom%bc_type(14)="robin"
            bottom%static(14)=.TRUE.
         ENDIF
      CASE("RTsun","RTearth")
c-----------------------------------------------------------------------
c     plasma density
c-----------------------------------------------------------------------
         top%bc_type(1)="zeroflux"
         bottom%bc_type(1)="zeroflux"
c-----------------------------------------------------------------------
c     poloidal magnetic flux
c-----------------------------------------------------------------------
         top%bc_type(2)="robin"
         top%static(2)=.FALSE.
         bottom%bc_type(2)="robin"
         bottom%static(2)=.FALSE.
c-----------------------------------------------------------------------
c     out of plane field
c-----------------------------------------------------------------------
         top%bc_type(3)="robin"
         top%static(3)=.TRUE.
         bottom%bc_type(3)="robin"
         bottom%static(3)=.TRUE.
c-----------------------------------------------------------------------
c     plasma momentum
c-----------------------------------------------------------------------
         top%bc_type(4:6)="robin"
         top%static(4:6)=.TRUE.
         bottom%bc_type(4:6)="robin"
         bottom%static(4:6)=.TRUE.
c-----------------------------------------------------------------------
c     current density
c-----------------------------------------------------------------------
         top%bc_type(7)="robin"
         top%static(7)=.TRUE.
         bottom%bc_type(7)="robin"
         bottom%static(7)=.TRUE.
c-----------------------------------------------------------------------
c     plasma pressure
c-----------------------------------------------------------------------
         top%bc_type(8)="zeroflux"
         bottom%bc_type(8)="zeroflux"
c-----------------------------------------------------------------------
c     neutral density
c-----------------------------------------------------------------------
         top%bc_type(9)="zeroflux"
         bottom%bc_type(9)="zeroflux"
c-----------------------------------------------------------------------
c     neutral momentum
c-----------------------------------------------------------------------
         top%bc_type(10:12)="robin"
         top%static(10:12)=.TRUE.
         bottom%bc_type(10:12)="robin"
         bottom%static(10:12)=.TRUE.
c-----------------------------------------------------------------------
c     neutral pressure
c-----------------------------------------------------------------------
         top%bc_type(13)="zeroflux"
         bottom%bc_type(13)="zeroflux"
c-----------------------------------------------------------------------
c     del^2(Bz) auxiliary variable
c-----------------------------------------------------------------------
         IF(nu /= 0)THEN
            top%bc_type(14)="robin"
            top%static(14)=.TRUE.
            bottom%bc_type(14)="robin"
            bottom%static(14)=.TRUE.
         ENDIF
      END SELECT
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
      USE pn_ext_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0

      SELECT CASE(init_type)
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top")
            c(4,:,:) = u(4,:,:)*uy(1,:,:) - u(1,:,:)*uy(4,:,:)
            c(5,:,:) = u(5,:,:)
            c(6,:,:) = u(6,:,:)*uy(1,:,:) - u(1,:,:)*uy(6,:,:)
            c(7,:,:) = uy(7,:,:)
            c(10,:,:) = u(10,:,:)*uy(9,:,:) - u(9,:,:)*uy(10,:,:)
            c(11,:,:) = u(11,:,:)
            c(12,:,:) = u(12,:,:)*uy(9,:,:) - u(9,:,:)*uy(12,:,:)
         CASE("left","right")
            c(3,:,:) = u(3,:,:)
            c(4,:,:) = u(4,:,:)
            c(10,:,:) = u(10,:,:)
         CASE("bottom")
            c(3,:,:) = u(3,:,:)
            c(5,:,:) = u(5,:,:)
            c(11,:,:) = u(11,:,:)
         END SELECT
         IF(nu /= 0)c(14,:,:) = u(14,:,:)
      CASE("MRX")
         SELECT CASE(lrtb)
         CASE("top")
            c(4,:,:) = u(4,:,:)*uy(1,:,:) - u(1,:,:)*uy(4,:,:)
            c(5,:,:) = u(5,:,:)
            c(6,:,:) = u(6,:,:)*uy(1,:,:) - u(1,:,:)*uy(6,:,:)
            c(7,:,:) = uy(7,:,:)
            c(10,:,:) = u(10,:,:)*uy(9,:,:) - u(9,:,:)*uy(10,:,:)
            c(11,:,:) = u(11,:,:)
            c(12,:,:) = u(12,:,:)*uy(9,:,:) - u(9,:,:)*uy(12,:,:)
         CASE("left","right")
            c(3,:,:) = u(3,:,:)
            c(4,:,:) = u(4,:,:)
            c(10,:,:) = u(10,:,:)
         CASE("bottom")
            c(3,:,:) = u(3,:,:)
            c(5,:,:) = u(5,:,:)
            c(11,:,:) = u(11,:,:)
         END SELECT
         IF(nu /= 0)c(14,:,:) = u(14,:,:)
      CASE("MRXpush")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(3,:,:) = uy(3,:,:)
            c(4,:,:) = u(4,:,:)*uy(1,:,:) - u(1,:,:)*uy(4,:,:)
            c(5,:,:) = u(5,:,:)
            c(6,:,:) = u(6,:,:)*uy(1,:,:) - u(1,:,:)*uy(6,:,:)
            c(7,:,:) = u(7,:,:)
            c(10,:,:) = u(10,:,:)*uy(9,:,:) - u(9,:,:)*uy(10,:,:)
            c(11,:,:) = u(11,:,:)
            c(12,:,:) = u(12,:,:)*uy(9,:,:) - u(9,:,:)*uy(12,:,:)
            IF(nu /= 0)c(14,:,:) = uy(14,:,:)
         CASE("left","right")
            c(3,:,:) = ux(3,:,:)
            c(4,:,:) = u(4,:,:)
            c(5,:,:) = u(5,:,:)*ux(1,:,:) - u(1,:,:)*ux(5,:,:)
            c(6,:,:) = u(6,:,:)*ux(1,:,:) - u(1,:,:)*ux(6,:,:)
            c(7,:,:) = u(7,:,:)
            c(10,:,:) = u(10,:,:)
            c(11,:,:) = u(11,:,:)*ux(9,:,:) - u(9,:,:)*ux(11,:,:)
            c(12,:,:) = u(12,:,:)*ux(9,:,:) - u(9,:,:)*ux(12,:,:)
            IF(nu /= 0)c(14,:,:) = ux(14,:,:)
         END SELECT
      CASE("forcefree")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(4,:,:) = u(4,:,:)*uy(1,:,:) - u(1,:,:)*uy(4,:,:)
            c(5,:,:) = u(5,:,:)
            c(6,:,:) = u(6,:,:)*uy(1,:,:) - u(1,:,:)*uy(6,:,:)
            c(7,:,:) = uy(7,:,:)
            c(10,:,:) = u(10,:,:)*uy(9,:,:) - u(9,:,:)*uy(10,:,:)
            c(11,:,:) = u(11,:,:)
            c(12,:,:) = u(12,:,:)*uy(9,:,:) - u(9,:,:)*uy(12,:,:)
         END SELECT
         IF(nu /= 0)c(14,:,:) = u(14,:,:)
      CASE("RTearth","RTsun")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c(3,:,:) = uy(3,:,:)
            c(4,:,:) = u(4,:,:)*uy(1,:,:) - u(1,:,:)*uy(4,:,:)
            c(5,:,:) = u(5,:,:)
            c(6,:,:) = u(6,:,:)
            c(7,:,:) = u(7,:,:)
            c(10,:,:) = u(10,:,:)*uy(9,:,:) - u(9,:,:)*uy(10,:,:)
            c(11,:,:) = u(11,:,:)
            c(12,:,:) = u(12,:,:)
            IF(nu /= 0)c(14,:,:) = uy(14,:,:)
         END SELECT
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
      USE pn_ext_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy
c-----------------------------------------------------------------------
c     zero arrays.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
      c_uxy=0

      SELECT CASE(init_type)
      CASE("GEM") 
         SELECT CASE(lrtb)
         CASE("top")
            c_u(4,1,:,:) = -uy(4,:,:)
            c_u(4,4,:,:) = uy(1,:,:)
            c_uy(4,1,:,:) = u(4,:,:)
            c_uy(4,4,:,:) = -u(1,:,:)

            c_u(5,5,:,:) = one

            c_u(6,1,:,:) = -uy(6,:,:)
            c_u(6,6,:,:) = uy(1,:,:)
            c_uy(6,1,:,:) = u(6,:,:)
            c_uy(6,6,:,:) = -u(1,:,:)

            c_uy(7,7,:,:) = one

            c_u(10,9,:,:) = -uy(10,:,:)
            c_u(10,10,:,:) = uy(9,:,:)
            c_uy(10,9,:,:) = u(10,:,:)
            c_uy(10,10,:,:) = -u(9,:,:)

            c_u(11,11,:,:) = one

            c_u(12,9,:,:) = -uy(12,:,:)
            c_u(12,12,:,:) = uy(9,:,:)
            c_uy(12,9,:,:) = u(12,:,:)
            c_uy(12,12,:,:) = -u(9,:,:)
         CASE("left","right")
            c_u(3,3,:,:) = one
            c_u(4,4,:,:) = one
            c_u(10,10,:,:) = one
         CASE("bottom")
            c_u(3,3,:,:) = one
            c_u(5,5,:,:) = one
            c_u(11,11,:,:) = one
         END SELECT
         IF(nu /= 0)c_u(14,14,:,:) = one
      CASE("MRX") 
         SELECT CASE(lrtb)
         CASE("top")
            c_u(4,1,:,:) = -uy(4,:,:)
            c_u(4,4,:,:) = uy(1,:,:)
            c_uy(4,1,:,:) = u(4,:,:)
            c_uy(4,4,:,:) = -u(1,:,:)

            c_u(5,5,:,:) = one

            c_u(6,1,:,:) = -uy(6,:,:)
            c_u(6,6,:,:) = uy(1,:,:)
            c_uy(6,1,:,:) = u(6,:,:)
            c_uy(6,6,:,:) = -u(1,:,:)

            c_uy(7,7,:,:) = one

            c_u(10,9,:,:) = -uy(10,:,:)
            c_u(10,10,:,:) = uy(9,:,:)
            c_uy(10,9,:,:) = u(10,:,:)
            c_uy(10,10,:,:) = -u(9,:,:)

            c_u(11,11,:,:) = one

            c_u(12,9,:,:) = -uy(12,:,:)
            c_u(12,12,:,:) = uy(9,:,:)
            c_uy(12,9,:,:) = u(12,:,:)
            c_uy(12,12,:,:) = -u(9,:,:)
         CASE("left","right")
            c_u(3,3,:,:) = one
            c_u(4,4,:,:) = one
            c_u(10,10,:,:) = one
         CASE("bottom")
            c_u(3,3,:,:) = one
            c_u(5,5,:,:) = one
            c_u(11,11,:,:) = one
         END SELECT
         IF(nu /= 0)c_u(14,14,:,:) = one
      CASE("MRXpush") 
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_uy(3,3,:,:) = one

            c_u(4,1,:,:) = -uy(4,:,:)
            c_u(4,4,:,:) = uy(1,:,:)
            c_uy(4,1,:,:) = u(4,:,:)
            c_uy(4,4,:,:) = -u(1,:,:)

            c_u(5,5,:,:) = one

            c_u(6,1,:,:) = -uy(6,:,:)
            c_u(6,6,:,:) = uy(1,:,:)
            c_uy(6,1,:,:) = u(6,:,:)
            c_uy(6,6,:,:) = -u(1,:,:)

            c_u(7,7,:,:) = one

            c_u(10,9,:,:) = -uy(10,:,:)
            c_u(10,10,:,:) = uy(9,:,:)
            c_uy(10,9,:,:) = u(10,:,:)
            c_uy(10,10,:,:) = -u(9,:,:)

            c_u(11,11,:,:) = one

            c_u(12,9,:,:) = -uy(12,:,:)
            c_u(12,12,:,:) = uy(9,:,:)
            c_uy(12,9,:,:) = u(12,:,:)
            c_uy(12,12,:,:) = -u(9,:,:)
            IF(nu /= 0)c_uy(14,14,:,:) = one
         CASE("left","right")
            c_ux(3,3,:,:) = one

            c_u(4,4,:,:) = one

            c_u(5,1,:,:) = -ux(5,:,:)
            c_u(5,5,:,:) = ux(1,:,:)
            c_ux(5,1,:,:) = u(5,:,:)
            c_ux(5,5,:,:) = -u(1,:,:)

            c_u(6,1,:,:) = -ux(6,:,:)
            c_u(6,6,:,:) = ux(1,:,:)
            c_ux(6,1,:,:) = u(6,:,:)
            c_ux(6,6,:,:) = -u(1,:,:)

            c_u(7,7,:,:) = one

            c_u(10,10,:,:) = one

            c_u(11,9,:,:) = -ux(11,:,:)
            c_u(11,11,:,:) = ux(9,:,:)
            c_ux(11,9,:,:) = u(11,:,:)
            c_ux(11,11,:,:) = -u(9,:,:)

            c_u(12,9,:,:) = -ux(12,:,:)
            c_u(12,12,:,:) = ux(9,:,:)
            c_ux(12,9,:,:) = u(12,:,:)
            c_ux(12,12,:,:) = -u(9,:,:)
            IF(nu /= 0)c_ux(14,14,:,:) = one
         END SELECT
      CASE("forcefree") 
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_u(4,1,:,:) = -uy(4,:,:)
            c_u(4,4,:,:) = uy(1,:,:)
            c_uy(4,1,:,:) = u(4,:,:)
            c_uy(4,4,:,:) = -u(1,:,:)

            c_u(5,5,:,:) = one

            c_u(6,1,:,:) = -uy(6,:,:)
            c_u(6,6,:,:) = uy(1,:,:)
            c_uy(6,1,:,:) = u(6,:,:)
            c_uy(6,6,:,:) = -u(1,:,:)

            c_uy(7,7,:,:) = one

            c_u(10,9,:,:) = -uy(10,:,:)
            c_u(10,10,:,:) = uy(9,:,:)
            c_uy(10,9,:,:) = u(10,:,:)
            c_uy(10,10,:,:) = -u(9,:,:)

            c_u(11,11,:,:) = one

            c_u(12,9,:,:) = -uy(12,:,:)
            c_u(12,12,:,:) = uy(9,:,:)
            c_uy(12,9,:,:) = u(12,:,:)
            c_uy(12,12,:,:) = -u(9,:,:)
         END SELECT
         IF(nu /= 0)c_u(14,14,:,:) = one
      CASE("RTearth","RTsun") 
         SELECT CASE(lrtb)
         CASE("top","bottom")
            c_uy(3,3,:,:) = one

            c_u(4,1,:,:) = -uy(4,:,:)
            c_u(4,4,:,:) = uy(1,:,:)
            c_uy(4,1,:,:) = u(4,:,:)
            c_uy(4,4,:,:) = -u(1,:,:)

            c_u(5,5,:,:) = one

            c_u(6,6,:,:) = one

            c_u(7,7,:,:) = one

            c_u(10,9,:,:) = -uy(10,:,:)
            c_u(10,10,:,:) = uy(9,:,:)
            c_uy(10,9,:,:) = u(10,:,:)
            c_uy(10,10,:,:) = -u(9,:,:)

            c_u(11,11,:,:) = one

            c_u(12,12,:,:) = one

            IF(nu /= 0)c_uy(14,14,:,:) = one
         END SELECT
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
      USE pn_ext_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     set mass,mass_x,mass_y(iqty,jqty,:,:), ... to coupling 
c     mass matrices for du/dt, d(du/dx)/dt, and d(du/dy)/dt
c-----------------------------------------------------------------------
      mass = zero
      mass_x = zero
      mass_y = zero

      SELECT CASE(init_type)
      CASE("GEM")
         SELECT CASE(lrtb)
         CASE("top")
            mass(2,2,:,:) = one
         END SELECT
      CASE("MRX")
         SELECT CASE(lrtb)
         CASE("top")
            mass(2,2,:,:) = one
         END SELECT
      CASE("MRXpush")
         SELECT CASE(lrtb)
         CASE("top","bottom","left","right")
            mass(2,2,:,:) = one
         END SELECT
      CASE("forcefree")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            mass(2,2,:,:) = one
         END SELECT
      CASE("RTsun","RTearth")
         SELECT CASE(lrtb)
         CASE("top","bottom")
            mass(2,2,:,:) = one
         END SELECT
      END SELECT
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
      USE pn_ext_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Tix,Tiy,kperp,kfac,
     $     ve,r_fac,r_faci,j1,j2,j3,jtot,eta_local,b1,b2,Bsq,
     $     n_inv,nx_inv,ny_inv,g_i,g_r,g_x,nn_inv,nnx_inv,qt_in,Tin, 
     $     nny_inv,Tnx,Tny,kediff,q_inx,q_nix,q_in,Ti,Tn,kapn,mun,r_en
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: vi,vix,viy,vn,vnx,
     $     vny,BdotT,r_inx,r_nix,r_in
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
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
      ny_inv=-uy(1,:,:)*n_inv**2
      nn_inv=one/u(9,:,:)
      nnx_inv=-ux(9,:,:)*nn_inv**2
      nny_inv=-uy(9,:,:)*nn_inv**2
c-----------------------------------------------------------------------
c     temperature gradients.
c-----------------------------------------------------------------------
      Tix = ux(8,:,:)*n_inv + u(8,:,:)*nx_inv
      Tiy = uy(8,:,:)*n_inv + u(8,:,:)*ny_inv
      Tnx = ux(13,:,:)*nn_inv + u(13,:,:)*nnx_inv
      Tny = uy(13,:,:)*nn_inv + u(13,:,:)*nny_inv
c-----------------------------------------------------------------------
c     ion & neutral temperatures.
c-----------------------------------------------------------------------
      Ti = ti_frac*u(8,:,:)*n_inv
      Tn = u(13,:,:)*nn_inv
      Tin= half*(Ti+Tn)
c-----------------------------------------------------------------------
c     magnetic fields and currents.
c-----------------------------------------------------------------------
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2 + u(3,:,:)**2
      j1 = uy(3,:,:) + cyl_fac*r_faci*u(3,:,:)
      j2 = -ux(3,:,:)
      j3 = u(7,:,:)
      jtot = SQRT(j1**2 + j2**2 + j3**2)
c-----------------------------------------------------------------------
c     resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta(eta_case,ly-y,u(1,:,:),u(8,:,:),jtot,
     $     cc,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,
     $     etavac,eta_local)
c-----------------------------------------------------------------------
c     ionized fluid heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT(b1,b2,u(3,:,:),Tix,Tiy,BdotT)
         CALL transport_kbrag(u(1,:,:),u(8,:,:),te_frac,Bsq,ke_norm,
     $        ki_norm,xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac)
      CASE("anisotropic")
         CALL transport_BdotT(b1,b2,u(3,:,:),Tix,Tiy,BdotT)
         CALL transport_setkaniso(kappa_par,kappa_perp,Bsq,kperp,kfac)
      CASE("isotropic")
         BdotT=0
         kperp=kappa_par
         kfac=0
      CASE DEFAULT
         CALL program_stop("pn_ext_rhs: cannot recognize "
     $        //"kappa_case = "//TRIM(kappa_case))
      END SELECT
c-----------------------------------------------------------------------
c     neutral fluid viscosity and heat conduction.
c-----------------------------------------------------------------------
      kapn = kapn_fac*4._r8/r_nn_norm*SQRT(Tn)
      mun = one/r_nn_norm*SQRT(Tn)
      IF(init_type == "forcefree")THEN
         WHERE(ABS(y) >= 0.8*ly)
            mun = mun + half*(one + COS(5._r8*pi*(ly-ABS(y))/ly))
         END WHERE
      ENDIF
c-----------------------------------------------------------------------
c     velocities and their gradients.
c     also, ion-neutral & electron-neutral friction terms
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:)=u(i+3,:,:)*n_inv
         vix(i,:,:)=ux(i+3,:,:)*n_inv + u(i+3,:,:)*nx_inv
         viy(i,:,:)=uy(i+3,:,:)*n_inv + u(i+3,:,:)*ny_inv

         vn(i,:,:)=u(i+9,:,:)*nn_inv
         vnx(i,:,:)=ux(i+9,:,:)*nn_inv + u(i+9,:,:)*nnx_inv
         vny(i,:,:)=uy(i+9,:,:)*nn_inv + u(i+9,:,:)*nny_inv

         r_in(i,:,:) = half*r_in_norm*SQRT(Tin)*u(1,:,:)*u(9,:,:)*
     $        (vn(i,:,:)-vi(i,:,:))
      ENDDO

      ve=(u(6,:,:)-di*u(7,:,:))*n_inv
      r_en = r_en_norm*SQRT(te_frac*u(8,:,:)*n_inv + Tn)*u(9,:,:)*n_inv
c-----------------------------------------------------------------------
c     reaction collision terms.
c-----------------------------------------------------------------------
      CALL pn_ext_grq(u(8,:,:),u(1,:,:),u(13,:,:),u(9,:,:),vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix)
c-----------------------------------------------------------------------
c     ion-neutral collision heating terms
c-----------------------------------------------------------------------
      kediff = half*SUM((vi - vn)**2,1)
      q_in = r_in_norm*half*SQRT(Tin)*u(1,:,:)*u(9,:,:)*kediff
      qt_in= 0.75_r8*r_in_norm*SQRT(Tin)*u(1,:,:)*u(9,:,:)*(Tn-Ti)
c-----------------------------------------------------------------------
c     density equation.
c-----------------------------------------------------------------------
      fx(1,:,:) = r_fac*(u(4,:,:) - Dn*ux(1,:,:))
      fy(1,:,:) = r_fac*(u(5,:,:) - Dn*uy(1,:,:))

      s(1,:,:) = r_fac*(g_i - g_r)
c-----------------------------------------------------------------------
c     -Az (or -Aphi)
c-----------------------------------------------------------------------
      fx(2,:,:) = r_fac*nu*ux(7,:,:)
      fy(2,:,:) = r_fac*nu*uy(7,:,:)
         
      s(2,:,:) = r_fac*(vi(2,:,:)*b1 - vi(1,:,:)*b2
     $     - di*(j2*b1 - j1*b2)*n_inv + eta_local*j3 
     $     + r_en*(vn(3,:,:)*u(1,:,:) - u(6,:,:) + di*j3)) 
     $     + cyl_fac*r_faci*nu*u(7,:,:)


      IF((init_type == "forcefree" .OR. init_type == "MRX") 
     $     .AND. t <= deltatau)s(2,:,:) = s(2,:,:) 
     $     + delta*bplane*EXP(-((0.25_r8*x)**2 + y**2)/lambda_psi**2)
     $     *(one - half*EXP(-3._r8*y**2/lambda_psi**2))
     $     *half*(one - COS(twopi*t/deltatau))
c-----------------------------------------------------------------------
c     out-of-plane Bz (or Bphi) magnetic field equation.
c-----------------------------------------------------------------------
      fx(3,:,:) = u(3,:,:)*vi(1,:,:) - ve*b1
     $     - di*te_frac*uy(8,:,:)*n_inv + eta_local*j2 
     $     + r_en*(vn(2,:,:)*u(1,:,:) - u(5,:,:) + di*j2)
      
      fy(3,:,:) = u(3,:,:)*vi(2,:,:) - ve*b2
     $     + di*te_frac*ux(8,:,:)*n_inv - eta_local*j1
     $     - r_en*(vn(1,:,:)*u(1,:,:) - u(4,:,:) + di*j1)
      
      s(3,:,:) = di*u(3,:,:)*(j2*ny_inv + j1*nx_inv) 
     $     + cyl_fac*two*r_faci*di*u(3,:,:)*ux(3,:,:)*n_inv

      IF(nu /= 0)THEN
         fx(3,:,:)=fx(3,:,:) + nu*ux(14,:,:)
         fy(3,:,:)=fy(3,:,:) 
     $        + nu*(cyl_fac*r_faci*u(14,:,:) + uy(14,:,:))
      ENDIF

      IF(ion_momentum)THEN
c-----------------------------------------------------------------------
c     x-component of ion momentum equation.
c-----------------------------------------------------------------------
         fx(4,:,:) = r_fac*(u(4,:,:)*vi(1,:,:) - two*mu*vix(1,:,:) 
     $        + half*u(3,:,:)**2 + u(8,:,:)
     $        - two*mu_sv*vix(1,:,:)*ABS(vix(1,:,:)))
         
         fy(4,:,:) = r_fac*(u(4,:,:)*vi(2,:,:) 
     $        - mu*(viy(1,:,:) + vix(2,:,:)))
         
         s(4,:,:) = r_fac*(-j3*b2 + g_i*vn(1,:,:) - g_r*vi(1,:,:) 
     $        - g_x*(vi(1,:,:) - vn(1,:,:))
     $        - r_nix(1,:,:) + r_inx(1,:,:) + r_in(1,:,:))
c-----------------------------------------------------------------------
c     y-component of ion momentum equation.
c-----------------------------------------------------------------------
         fx(5,:,:) = r_fac*(u(5,:,:)*vi(1,:,:)
     $        - mu*(vix(2,:,:) + viy(1,:,:)))
         
         fy(5,:,:) = r_fac*(u(5,:,:)*vi(2,:,:) - two*mu*viy(2,:,:)
     $        + half*u(3,:,:)**2 - two*mu_sv*viy(2,:,:)*ABS(viy(2,:,:)))
         
         s(5,:,:) = r_fac*(j3*b1 - uy(8,:,:) - gravity*u(1,:,:)
     $        + g_i*vn(2,:,:) - g_r*vi(2,:,:) 
     $        - g_x*(vi(2,:,:) - vn(2,:,:))
     $        - r_nix(2,:,:) + r_inx(2,:,:) + r_in(2,:,:))
     $        + cyl_fac*(u(6,:,:)*vi(3,:,:) - half*u(3,:,:)**2
     $        - two*r_faci*(mu*vi(2,:,:)
     $        + mu_sv*r_faci*vi(2,:,:)*ABS(vi(2,:,:))))
c-----------------------------------------------------------------------
c     z-component of ion momentum equation.
c-----------------------------------------------------------------------
         fx(6,:,:) = r_fac*(u(6,:,:)*vi(1,:,:) - mu*vix(3,:,:))
         
         fy(6,:,:) = r_fac*(u(6,:,:)*vi(2,:,:) - mu*viy(3,:,:))
         
         s(6,:,:) = r_fac*(j1*b2 - j2*b1
     $        + g_i*vn(3,:,:) - g_r*vi(3,:,:)
     $        - g_x*(vi(3,:,:) - vn(3,:,:))
     $        - r_nix(3,:,:) + r_inx(3,:,:) + r_in(3,:,:))
     $        - cyl_fac*(u(6,:,:)*vi(2,:,:) + r_faci*mu*vi(3,:,:))
      ELSE
c-----------------------------------------------------------------------
c     x-component of total fluid momentum equation.
c-----------------------------------------------------------------------
         fx(4,:,:) = r_fac*(u(4,:,:)*vi(1,:,:) + u(10,:,:)*vn(1,:,:)
     $        - two*(mu*vix(1,:,:) + mun*vnx(1,:,:))
     $        + half*u(3,:,:)**2 + u(8,:,:) + u(13,:,:)
     $        - two*mu_sv*vix(1,:,:)*ABS(vix(1,:,:)))
         
         fy(4,:,:) = r_fac*(u(4,:,:)*vi(2,:,:) + u(10,:,:)*vn(2,:,:)
     $        - mu*(viy(1,:,:) + vix(2,:,:))
     $        - mun*(vny(1,:,:) + vnx(2,:,:)))
         
         s(4,:,:) = -r_fac*(j3*b2 
     $        + mfp_fac*(g_x*vi(1,:,:) + r_nix(1,:,:)))
c-----------------------------------------------------------------------
c     y-component of total fluid momentum equation.
c-----------------------------------------------------------------------
         fx(5,:,:) = r_fac*(u(5,:,:)*vi(1,:,:) + u(11,:,:)*vn(1,:,:)
     $        - mu*(vix(2,:,:) + viy(1,:,:))
     $        - mun*(vnx(2,:,:) + vny(1,:,:)))
         
         fy(5,:,:) = r_fac*(u(5,:,:)*vi(2,:,:) + u(11,:,:)*vn(2,:,:)
     $        + half*u(3,:,:)**2 - two*(mu*viy(2,:,:) + mun*vny(2,:,:)) 
     $        - two*mu_sv*viy(2,:,:)*ABS(viy(2,:,:)))
         
         s(5,:,:) = r_fac*(j3*b1 - uy(8,:,:) - uy(13,:,:)
     $        - mfp_fac*(g_x*vi(2,:,:) + r_nix(2,:,:))
     $        - gravity*(u(1,:,:) + u(9,:,:)))
     $        + cyl_fac*(u(6,:,:)*vi(3,:,:) + u(12,:,:)*vn(3,:,:)
     $        - half*u(3,:,:)**2
     $        - two*r_faci*(mu*vi(2,:,:) + mun*vn(2,:,:)
     $        + mu_sv*r_faci*vi(2,:,:)*ABS(vi(2,:,:))))
c-----------------------------------------------------------------------
c     z-component of total fluid momentum equation.
c-----------------------------------------------------------------------
         fx(6,:,:) = r_fac*(u(6,:,:)*vi(1,:,:) + u(12,:,:)*vn(1,:,:)
     $        - mu*vix(3,:,:) - mun*vnx(3,:,:))
         
         fy(6,:,:) = r_fac*(u(6,:,:)*vi(2,:,:) + u(12,:,:)*vn(2,:,:)
     $        - mu*viy(3,:,:) - mun*vny(3,:,:))
         
         s(6,:,:) = r_fac*(j1*b2 - j2*b1
     $        - mfp_fac*(g_x*vi(3,:,:) + r_nix(3,:,:)))
     $        - cyl_fac*(u(6,:,:)*vi(2,:,:) + u(12,:,:)*vn(2,:,:)
     $        + r_faci*(mu*vi(3,:,:) + mun*vn(3,:,:)))
      ENDIF
c-----------------------------------------------------------------------
c     out-of-plane electron momentum (current) equation.
c-----------------------------------------------------------------------
      fx(7,:,:)=b2
      fy(7,:,:)=-b1
      s(7,:,:)=j3
c-----------------------------------------------------------------------
c     plasma pressure equation.
c-----------------------------------------------------------------------
      fx(8,:,:)=r_fac*(gamma_fac*u(8,:,:)
     $     *(vi(1,:,:) - te_frac*di*j1*n_inv)
     $     - kfac*BdotT(1,:,:) - kperp*Tix)
      fy(8,:,:)=r_fac*(gamma_fac*u(8,:,:)
     $     *(vi(2,:,:) - te_frac*di*j2*n_inv)
     $     - kfac*BdotT(2,:,:) - kperp*Tiy)
      s(8,:,:)=r_fac*(ux(8,:,:)*(vi(1,:,:) - te_frac*di*j1*n_inv) 
     $     + uy(8,:,:)*(vi(2,:,:) - te_frac*di*j2*n_inv)
     $     + jtot**2*(eta_local + di*r_en)
     $     + mu*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:)
     $     + vix(3,:,:)**2 + viy(3,:,:)**2)
     $     + cyl_fac*r_faci**2*mu*(two*vi(2,:,:)**2 + vi(3,:,:)**2) 
     $     + mu_sv*two*(vix(1,:,:)**2*ABS(vix(1,:,:))
     $     + viy(2,:,:)**2*ABS(viy(2,:,:)))
     $     + cyl_fac*mu_sv*two*r_faci**3*vi(2,:,:)**2*ABS(vi(2,:,:))
     $     + (g_i + g_x)*kediff - g_i*phi_eff
     $     + 1.5_r8*(g_i*u(13,:,:)*nn_inv - g_r*ti_frac*u(8,:,:)*n_inv)
     $     + r_inx(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + r_inx(2,:,:)*(vn(2,:,:) - vi(2,:,:))
     $     + r_inx(3,:,:)*(vn(3,:,:) - vi(3,:,:))
     $     + q_inx - q_nix + q_in + qt_in)
c-----------------------------------------------------------------------
c     total fluid density.
c-----------------------------------------------------------------------
      fx(9,:,:) = r_fac*(u(4,:,:) + u(10,:,:) 
     $     - Dn*(ux(1,:,:) + ux(9,:,:)))
      fy(9,:,:) = r_fac*(u(5,:,:) + u(11,:,:) 
     $     - Dn*(uy(1,:,:) + uy(9,:,:)))

      s(9,:,:) = -r_fac*mfp_fac*g_x
c-----------------------------------------------------------------------
c     x-component of neutral momentum equation.
c-----------------------------------------------------------------------
      fx(10,:,:) = r_fac*(u(10,:,:)*vn(1,:,:)
     $     - two*mun*vnx(1,:,:) + u(13,:,:))
      fy(10,:,:) = r_fac*(u(10,:,:)*vn(2,:,:) 
     $     - mun*(vny(1,:,:) + vnx(2,:,:)))
      s(10,:,:) = r_fac*(-g_i*vn(1,:,:) + g_r*vi(1,:,:) 
     $     + g_x*((one-mfp_fac)*vi(1,:,:) - vn(1,:,:))
     $     + (one-mfp_fac)*r_nix(1,:,:) - r_inx(1,:,:) - r_in(1,:,:))
c-----------------------------------------------------------------------
c     y-component of neutral momentum equation.
c-----------------------------------------------------------------------
      fx(11,:,:) = r_fac*(u(11,:,:)*vn(1,:,:)
     $     - mun*(vnx(2,:,:) + vny(1,:,:)))
      fy(11,:,:) = r_fac*(u(11,:,:)*vn(2,:,:) - two*mun*vny(2,:,:))
      s(11,:,:) = r_fac*(-uy(13,:,:) - g_i*vn(2,:,:) 
     $     + g_r*vi(2,:,:) + g_x*((one-mfp_fac)*vi(2,:,:) 
     $     - vn(2,:,:)) + (one-mfp_fac)*r_nix(2,:,:) - r_inx(2,:,:) 
     $     - r_in(2,:,:) - gravity*u(9,:,:))
     $     + cyl_fac*(u(12,:,:)*vn(3,:,:) - two*mun*vn(2,:,:)*r_faci)
c-----------------------------------------------------------------------
c     z-component of neutral momentum equation.
c-----------------------------------------------------------------------
      fx(12,:,:) = r_fac*(u(12,:,:)*vn(1,:,:) - mun*vnx(3,:,:))
      fy(12,:,:) = r_fac*(u(12,:,:)*vn(2,:,:) - mun*vny(3,:,:))
      s(12,:,:) = r_fac*(-g_i*vn(3,:,:) + g_r*vi(3,:,:)
     $     + g_x*((one-mfp_fac)*vi(3,:,:) - vn(3,:,:))
     $     + (one-mfp_fac)*r_nix(3,:,:) - r_inx(3,:,:) - r_in(3,:,:))
     $     - cyl_fac*(u(12,:,:)*vn(2,:,:) + r_faci*mun*vn(3,:,:))
c-----------------------------------------------------------------------
c     neutral pressure.
c-----------------------------------------------------------------------
      fx(13,:,:) = r_fac*(gamma_fac*u(13,:,:)*vn(1,:,:) - kapn*Tnx)
      fy(13,:,:) = r_fac*(gamma_fac*u(13,:,:)*vn(2,:,:) - kapn*Tny)
      s(13,:,:) = r_fac*(vn(1,:,:)*ux(13,:,:) + vn(2,:,:)*uy(13,:,:)
     $     + mun*(two*vnx(1,:,:)**2 + two*vny(2,:,:)**2 + vny(1,:,:)**2
     $     + vnx(2,:,:)**2 + two*vny(1,:,:)*vnx(2,:,:)
     $     + vnx(3,:,:)**2 + vny(3,:,:)**2)
     $     + cyl_fac*r_faci**2*mun*(two*vn(2,:,:)**2 + vn(3,:,:)**2)
     $     + (g_r + g_x*(one-mfp_fac))*kediff 
     $     + 1.5_r8*(-g_i*u(13,:,:)*nn_inv + g_r*ti_frac*u(8,:,:)*n_inv) 
     $     + (one-mfp_fac)*(r_nix(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + r_nix(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + r_nix(3,:,:)*(vi(3,:,:) - vn(3,:,:)) + q_nix)
     $     - q_inx + q_in - qt_in)
c-----------------------------------------------------------------------
c     del^2(out-of-plane B) auxiliary equation.
c-----------------------------------------------------------------------
      IF(nu /= 0)THEN
         fx(14,:,:)=ux(3,:,:)
         fy(14,:,:)=cyl_fac*r_faci*u(3,:,:) + uy(3,:,:)
         s(14,:,:)=u(14,:,:)
      ENDIF
c-----------------------------------------------------------------------
c     zero out-of-plane equations when init_type="RTearth" or "RTsun"
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("RTearth","RTsun")
         fx(2,:,:)=0
         fy(2,:,:)=0
         s(2,:,:)=0
         fx(6,:,:)=0
         fy(6,:,:)=0
         s(6,:,:)=0
         fx(12,:,:)=0
         fy(12,:,:)=0
         s(12,:,:)=0
         fx(7,:,:)=0
         fy(7,:,:)=0
         s(7,:,:)=u(7,:,:)
      END SELECT
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
      USE pn_ext_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac,r_faci,Tix,Tiy,
     $     Tix_un,Tiy_un,Ti_un,kperp,kfac,kperp_bsq,kpar_rho,kpar_p,
     $     kperp_rho,kperp_p,ve,j1,j2,j3,jtot,eta_local,b1,b2,
     $     Bsq,eta_rho,eta_p,eta_j,j_u3,j_ux3,j_uy3,j_u7,
     $     g_i,g_i_rho,g_i_p,g_i_rhon,
     $     g_r,g_r_rho,g_r_p,g_x,g_x_p,g_x_pn,g_x_rho,g_x_rhon,
     $     q_inx,q_nix,n_inv,nx_inv,ny_inv,nn_inv,nnx_inv,nny_inv,
     $     kediff,kediff_rho,kediff_rhon,Tnx,Tny,Tn_un,Tnx_un,Tny_un,
     $     g_x_mi1,g_x_mi2,g_x_mi3,g_x_mn1,g_x_mn2,g_x_mn3,
     $     qinx_rho,qinx_p,qinx_rhon,qinx_pn,
     $     qnix_rho,qnix_p,qnix_rhon,qnix_pn, 
     $     r_in_u1,r_in_u9,q_in_8_1,q_in_8_8,q_in_8_13,q_in_8_9,
     $     qt_in_8_1,qt_in_8_8,qt_in_8_9,qt_in_8_13,Ti,Tn,Tin,Te,
     $     kapn,mun,mun_un,mun_pn,r_en,r_en_ut

      REAL(r8), DIMENSION(3,SIZE(x,1),SIZE(x,2)) :: vi,vix,viy,r_inx,
     $     r_nix,vi_un,vix_un,viy_un,BdotT,BdotT_b1,BdotT_b2,BdotT_b3,
     $     BdotT_Tx,BdotT_Ty,vn,vnx,vny,vn_un,vnx_un,vny_un,
     $     rinx_rho,rinx_p,rinx_rhon,rinx_pn,rinx_mi1,rinx_mi2,
     $     rinx_mi3,rinx_mn1,rinx_mn2,rinx_mn3,qinx_mi,qinx_mn,
     $     rnix_rho,rnix_p,rnix_rhon,rnix_pn,rnix_mi1,rnix_mi2,rnix_mi3,
     $     rnix_mn1,rnix_mn2,rnix_mn3,qnix_mi,qnix_mn,kediff_mi,
     $     kediff_mn,q_in_vi,q_in_vn,r_in_rho,r_in_rhon,r_in_p,r_in_pn,
     $     r_en_rho,r_en_rhon,r_en_p,r_en_pn
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
      ny_inv=-uy(1,:,:)*n_inv**2

      nn_inv=one/u(9,:,:)
      nnx_inv=-ux(9,:,:)*nn_inv**2
      nny_inv=-uy(9,:,:)*nn_inv**2
c-----------------------------------------------------------------------
c     temperature gradients and derivatives.
c-----------------------------------------------------------------------
      Tix = ux(8,:,:)*n_inv + u(8,:,:)*nx_inv
      Tiy = uy(8,:,:)*n_inv + u(8,:,:)*ny_inv
      Ti_un=-u(8,:,:)*n_inv**2
      Tix_un=-(two*u(8,:,:)*nx_inv + ux(8,:,:)*n_inv)*n_inv
      Tiy_un=-(two*u(8,:,:)*ny_inv + uy(8,:,:)*n_inv)*n_inv

      Tnx = ux(13,:,:)*nn_inv + u(13,:,:)*nnx_inv
      Tny = uy(13,:,:)*nn_inv + u(13,:,:)*nny_inv
      Tn_un=-u(13,:,:)*nn_inv**2
      Tnx_un=-(two*u(13,:,:)*nnx_inv + ux(13,:,:)*nn_inv)*nn_inv
      Tny_un=-(two*u(13,:,:)*nny_inv + uy(13,:,:)*nn_inv)*nn_inv
c-----------------------------------------------------------------------
c     ion & neutral temperatures.
c-----------------------------------------------------------------------
      Ti = ti_frac*u(8,:,:)*n_inv
      Tn = u(13,:,:)*nn_inv
      Tin= half*(Ti+Tn)
c-----------------------------------------------------------------------
c     magnetic fields and currents.
c-----------------------------------------------------------------------
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2 + u(3,:,:)**2
      j1 = uy(3,:,:) + cyl_fac*r_faci*u(3,:,:)
      j2 = -ux(3,:,:)
      j3 = u(7,:,:)
      jtot = SQRT(j1**2 + j2**2 + j3**2)

      j_u3 = cyl_fac*r_faci*j1/jtot
      j_ux3 = -j2/jtot
      j_uy3 = j1/jtot
      j_u7 = j3/jtot
      WHERE(jtot == 0)
         j_u3 = cyl_fac*r_faci
         j_ux3 = -one
         j_uy3 = one
         j_u7 = one
      END WHERE
c-----------------------------------------------------------------------
c     resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta_u(eta_case,ly-y,u(1,:,:),u(8,:,:),jtot,cc,
     $     etas_norm,etac_norm,v_chod_norm,r_eta,eta,etavac,eta_local,
     $     eta_rho,eta_p,eta_j)
c-----------------------------------------------------------------------
c     ionized fluid heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT_u(b1,b2,u(3,:,:),Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty)
         CALL transport_kbrag_u(u(1,:,:),u(8,:,:),te_frac,Bsq,ke_norm,
     $        ki_norm,xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac,
     $        kpar_rho,kpar_p,kperp_rho,kperp_p,kperp_bsq)
      CASE("anisotropic")
         kpar_rho=0
         kpar_p=0
         kperp_rho=0
         kperp_p=0
         CALL transport_BdotT_u(b1,b2,u(3,:,:),Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,BdotT_b3,BdotT_Tx,BdotT_Ty)
         CALL transport_setkaniso_u(kappa_par,kappa_perp,Bsq,kperp,kfac,
     $        kperp_bsq)
      CASE("isotropic")
         BdotT=0
         BdotT_b1=0
         BdotT_b2=0
         BdotT_b3=0
         BdotT_Tx=0
         BdotT_Ty=0
         kperp=kappa_par
         kfac=0
         kpar_rho=0
         kpar_p=0
         kperp_rho=0
         kperp_p=0
         kperp_bsq=0
      END SELECT
c-----------------------------------------------------------------------
c     neutral fluid viscosity and heat conduction.
c-----------------------------------------------------------------------
      kapn = kapn_fac*4._r8/r_nn_norm*SQRT(Tn)
      mun = one/r_nn_norm*SQRT(Tn)
      IF(init_type == "forcefree")THEN
         WHERE(ABS(y) >= 0.8*ly)
            mun = mun + half*(one + COS(5._r8*pi*(ly-ABS(y))/ly))
         END WHERE
      ENDIF
      mun_un = -half*mun*nn_inv
      mun_pn = half*mun/u(13,:,:)
c-----------------------------------------------------------------------
c     velocities and their derivatives.
c-----------------------------------------------------------------------
      DO i=1,3
         vi(i,:,:)=u(i+3,:,:)*n_inv
         vix(i,:,:)=ux(i+3,:,:)*n_inv + u(i+3,:,:)*nx_inv
         viy(i,:,:)=uy(i+3,:,:)*n_inv + u(i+3,:,:)*ny_inv

         vi_un(i,:,:)=-vi(i,:,:)*n_inv
         vix_un(i,:,:) = -(vi(i,:,:)*nx_inv + vix(i,:,:)*n_inv)
         viy_un(i,:,:) = -(vi(i,:,:)*ny_inv + viy(i,:,:)*n_inv)

         vn(i,:,:)=u(i+9,:,:)*nn_inv
         vnx(i,:,:)=ux(i+9,:,:)*nn_inv + u(i+9,:,:)*nnx_inv
         vny(i,:,:)=uy(i+9,:,:)*nn_inv + u(i+9,:,:)*nny_inv

         vn_un(i,:,:)=-vn(i,:,:)*nn_inv
         vnx_un(i,:,:) = -(vn(i,:,:)*nnx_inv + vnx(i,:,:)*nn_inv)
         vny_un(i,:,:) = -(vn(i,:,:)*nny_inv + vny(i,:,:)*nn_inv)

         kediff_mi(i,:,:) = (vi(i,:,:) - vn(i,:,:))*n_inv
         kediff_mn(i,:,:) = (vn(i,:,:) - vi(i,:,:))*nn_inv
      ENDDO

      ve=(u(6,:,:)-di*u(7,:,:))*n_inv
c-----------------------------------------------------------------------
c     kinetic energies and derivatives.
c-----------------------------------------------------------------------
      kediff = half*SUM((vi - vn)**2,1)

      kediff_rho = SUM((vi - vn)*vi_un,1)
      kediff_rhon = SUM((vn - vi)*vn_un,1)
c-----------------------------------------------------------------------
c     atomic reaction rates and derivatives.
c-----------------------------------------------------------------------
      CALL pn_ext_grq_u(u(8,:,:),u(1,:,:),u(13,:,:),u(9,:,:),vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix,
     $     g_i_rho,g_i_rhon,g_i_p,g_r_rho,g_r_p,
     $     g_x_rho,g_x_rhon,g_x_p,g_x_pn,g_x_mi1,g_x_mi2,g_x_mi3,
     $     g_x_mn1,g_x_mn2,g_x_mn3,
     $     rinx_rho,rinx_p,rinx_rhon,rinx_pn,rinx_mi1,rinx_mi2,rinx_mi3,
     $     rinx_mn1,rinx_mn2,rinx_mn3,
     $     rnix_rho,rnix_p,rnix_rhon,rnix_pn,rnix_mi1,rnix_mi2,rnix_mi3,
     $     rnix_mn1,rnix_mn2,rnix_mn3,
     $     qinx_rho,qinx_p,qinx_rhon,qinx_pn,qinx_mi,qinx_mn,
     $     qnix_rho,qnix_p,qnix_rhon,qnix_pn,qnix_mi,qnix_mn)
c-----------------------------------------------------------------------
c     compute i-n collision terms.
c-----------------------------------------------------------------------
      r_in_u1 = half*r_in_norm*u(1,:,:)*SQRT(Tin)
      r_in_u9 = -half*r_in_norm*u(9,:,:)*SQRT(Tin)

      q_in_8_8 = -r_in_u9*0.25_r8*ti_frac*kediff/Tin
      q_in_8_13 = r_in_u1*0.25_r8*kediff/Tin

      q_in_8_1 = 0._r8
      q_in_8_9 = 0._r8
      DO i=1,3
         r_in_rho(i,:,:) = -r_in_u9
     $        *(vn(i,:,:)*(one - 0.25_r8*Ti/Tin) 
     $        + vi(i,:,:)*0.25_r8*Ti/Tin)
         r_in_rhon(i,:,:) = -r_in_u1
     $        *(vi(i,:,:)*(one - 0.25_r8*Tn/Tin) 
     $        + vn(i,:,:)*0.25_r8*Tn/Tin)
         r_in_p(i,:,:) = -ti_frac*r_in_norm
     $        *0.125_r8*u(1,:,:)*u(9,:,:)*kediff_mi(i,:,:)/SQRT(Tin)
         r_in_pn(i,:,:) = r_in_norm
     $        *0.125_r8*u(1,:,:)*u(9,:,:)*kediff_mn(i,:,:)/SQRT(Tin)  
         q_in_vi(i,:,:) = r_in_u9*(vn(i,:,:) - vi(i,:,:))
         q_in_vn(i,:,:) = r_in_u1*(vn(i,:,:) - vi(i,:,:))

         q_in_8_1 = q_in_8_1 + r_in_u9*half*(vi(i,:,:) - vn(i,:,:))
     $        *(vn(i,:,:)*(one - 0.25_r8*Ti/Tin)
     $        + vi(i,:,:)*(one + 0.25_r8*Ti/Tin))

         q_in_8_9 = q_in_8_9 - r_in_u1*half*(vn(i,:,:) - vi(i,:,:))
     $        *(vn(i,:,:)*(one + 0.25_r8*Tn/Tin)
     $        + vi(i,:,:)*(one - 0.25_r8*Tn/Tin))
      ENDDO

c     NOTE - q_in_13's are same as q_in_8's as q_in is same for both 
c     (if m=min)
      qt_in_8_1 = -0.375_r8*r_in_u9*(4.0_r8*Tn - (Tn-Ti)*Ti/Tin)
      qt_in_8_9 = -0.375_r8*r_in_u1*(4.0_r8*Ti + (Tn-Ti)*Tn/Tin)
      qt_in_8_8 = 0.375_r8*ti_frac*r_in_u9*(4.0_r8 - (Tn-Ti)/Tin)
      qt_in_8_13 = 0.375_r8*r_in_u1*(4.0_r8 + (Tn-Ti)/Tin)
c     NOTE - qt_in_13's are inverse of  qt_in_8's (if m=min)
c-----------------------------------------------------------------------
c     compute e-n collision terms.
c-----------------------------------------------------------------------
      Te = te_frac*u(8,:,:)*n_inv

      r_en = r_en_norm*SQRT(Te + Tn)*u(9,:,:)*n_inv
      r_en_ut = r_en/(Te + Tn)

      r_en_rho(1,:,:) = r_en_ut*n_inv*(-half*te_frac*u(8,:,:)*vn(1,:,:)
     $     + (u(4,:,:) - di*j1)*(Tn + 1.5_r8*Te))
      r_en_rho(2,:,:) = r_en_ut*n_inv*(-half*te_frac*u(8,:,:)*vn(2,:,:)
     $     + (u(5,:,:) - di*j2)*(Tn + 1.5_r8*Te))
      r_en_rho(3,:,:) = r_en_ut*n_inv*(-half*te_frac*u(8,:,:)*vn(3,:,:)
     $     + (u(6,:,:) - di*j3)*(Tn + 1.5_r8*Te))

      r_en_rhon(1,:,:) = -r_en_ut*nn_inv*(half*Tn*vn(1,:,:)*u(1,:,:)
     $     + (u(4,:,:) - di*j1)*(half*Tn + Te))
      r_en_rhon(2,:,:) = -r_en_ut*nn_inv*(half*Tn*vn(2,:,:)*u(1,:,:)
     $     + (u(5,:,:) - di*j2)*(half*Tn + Te))
      r_en_rhon(3,:,:) = -r_en_ut*nn_inv*(half*Tn*vn(3,:,:)*u(1,:,:)
     $     + (u(6,:,:) - di*j3)*(half*Tn + Te))

      r_en_p(1,:,:) = r_en_ut*half*te_frac
     $     *(vn(1,:,:) - vi(1,:,:) + di*j1*n_inv)
      r_en_p(2,:,:) = r_en_ut*half*te_frac
     $     *(vn(2,:,:) - vi(2,:,:) + di*j2*n_inv)
      r_en_p(3,:,:) = r_en_ut*half*te_frac
     $     *(vn(3,:,:) - vi(3,:,:) + di*j3*n_inv)

      r_en_pn(1,:,:) = r_en_ut*half*nn_inv
     $     *(vn(1,:,:)*u(1,:,:) - u(4,:,:) + di*j1)
      r_en_pn(2,:,:) = r_en_ut*half*nn_inv
     $     *(vn(2,:,:)*u(1,:,:) - u(5,:,:) + di*j2)
      r_en_pn(3,:,:) = r_en_ut*half*nn_inv
     $     *(vn(3,:,:)*u(1,:,:) - u(6,:,:) + di*j3)
c-----------------------------------------------------------------------
c     plasma density equation.
c-----------------------------------------------------------------------
      fx_u(1,4,:,:) = r_fac
      fx_ux(1,1,:,:) = -r_fac*Dn
      fy_u(1,5,:,:) = r_fac
      fy_uy(1,1,:,:) = -r_fac*Dn

      s_u(1,1,:,:)=r_fac*(g_i_rho-g_r_rho)
      s_u(1,8,:,:)=r_fac*(g_i_p-g_r_p)
      s_u(1,9,:,:)=r_fac*g_i_rhon
c-----------------------------------------------------------------------
c     -Az (or -Aphi)
c-----------------------------------------------------------------------
      fx_ux(2,7,:,:) = r_fac*nu
      fy_uy(2,7,:,:) = r_fac*nu

      s_u(2,1,:,:)=r_fac*(vi_un(2,:,:)*b1 - vi_un(1,:,:)*b2
     $     + di*(j2*b1 - j1*b2)*n_inv**2 + eta_rho*j3 + r_en_rho(3,:,:))
      s_u(2,2,:,:)=-cyl_fac*(vi(2,:,:)-di*j2*n_inv)
      s_u(2,3,:,:)=di*cyl_fac*b2*n_inv + r_fac*eta_j*j_u3*j3
      s_u(2,4,:,:)=-r_fac*b2*n_inv
      s_u(2,5,:,:)=r_fac*b1*n_inv
      s_u(2,6,:,:)=-r_fac*r_en 
      s_u(2,7,:,:)=r_fac*(eta_j*j_u7*j3 + eta_local + di*r_en)
     $     + cyl_fac*r_faci*nu
      s_u(2,8,:,:)=r_fac*(eta_p*j3 + r_en_p(3,:,:))
      s_u(2,9,:,:)=r_fac*r_en_rhon(3,:,:)
      s_u(2,12,:,:)=r_fac*r_en
      s_u(2,13,:,:)=r_fac*r_en_pn(3,:,:)
      s_ux(2,2,:,:)=r_fac*(di*j1*n_inv - vi(1,:,:))
      s_ux(2,3,:,:)=r_fac*(di*b1*n_inv + eta_j*j_ux3*j3)
      s_uy(2,2,:,:)=-r_fac*(vi(2,:,:)-di*j2*n_inv)
      s_uy(2,3,:,:)=r_fac*(di*b2*n_inv + eta_j*j_uy3*j3)
c-----------------------------------------------------------------------
c     out-of-plane magnetic field equation.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=u(3,:,:)*vi_un(1,:,:) + ve*n_inv*b1 
     $     + di*te_frac*uy(8,:,:)*n_inv**2 +  eta_rho*j2
     $     + r_en_rho(2,:,:)
      fx_u(3,2,:,:)=cyl_fac*r_faci*ve 
      fx_u(3,3,:,:)=vi(1,:,:) + eta_j*j_u3*j2
      fx_u(3,4,:,:)=u(3,:,:)*n_inv
      fx_u(3,5,:,:)=-r_en
      fx_u(3,6,:,:)=-n_inv*b1
      fx_u(3,7,:,:)=di*n_inv*b1 + eta_j*j_u7*j2
      fx_u(3,8,:,:)= eta_p*j2 + r_en_p(2,:,:)
      fx_u(3,9,:,:)=r_en_rhon(2,:,:)
      fx_u(3,11,:,:)=r_en
      fx_u(3,13,:,:)=r_en_pn(2,:,:)
      fx_ux(3,3,:,:)=-eta_local + eta_j*j_ux3*j2 - di*r_en
      fx_uy(3,2,:,:)=ve
      fx_uy(3,3,:,:)=eta_j*j_uy3*j2
      fx_uy(3,8,:,:)=-di*te_frac*n_inv
      
      fy_u(3,1,:,:)=u(3,:,:)*vi_un(2,:,:) + ve*n_inv*b2 
     $     - di*te_frac*ux(8,:,:)*n_inv**2 - eta_rho*j1
     $     - r_en_rho(1,:,:)
      fy_u(3,3,:,:)=vi(2,:,:) - eta_j*j_u3*j1 
     $     - cyl_fac*r_faci*(eta_local + di*r_en)
      fy_u(3,4,:,:)=r_en
      fy_u(3,5,:,:)=u(3,:,:)*n_inv
      fy_u(3,6,:,:)=-n_inv*b2
      fy_u(3,7,:,:)=di*n_inv*b2 - eta_j*j_u7*j1
      fy_u(3,8,:,:)= -eta_p*j1 - r_en_p(1,:,:)
      fy_u(3,9,:,:)=-r_en_rhon(1,:,:)
      fy_u(3,10,:,:)=-r_en
      fy_u(3,13,:,:)=-r_en_pn(1,:,:)
      fy_ux(3,2,:,:)=-ve
      fy_ux(3,3,:,:)=-eta_j*j_ux3*j1
      fy_ux(3,8,:,:)=di*te_frac*n_inv
      fy_uy(3,3,:,:)=-eta_local - eta_j*j_uy3*j1 - di*r_en
      
      s_u(3,1,:,:)=-two*di*u(3,:,:)*(j2*ny_inv + j1*nx_inv)*n_inv
     $     - cyl_fac*two*r_faci*di*u(3,:,:)*ux(3,:,:)*n_inv**2
      s_u(3,3,:,:)=di*(j2*ny_inv + j1*nx_inv)
     $     + cyl_fac*r_faci*di*(u(3,:,:)*nx_inv 
     $     + two*ux(3,:,:)*n_inv)
      s_ux(3,1,:,:)=-di*u(3,:,:)*j1*n_inv**2
      s_ux(3,3,:,:)=-di*u(3,:,:)*ny_inv
     $     + cyl_fac*two*r_faci*di*u(3,:,:)*n_inv
      s_uy(3,1,:,:)=-di*u(3,:,:)*j2*n_inv**2
      s_uy(3,3,:,:)=di*u(3,:,:)*nx_inv

      IF(nu /= 0)THEN
         fx_ux(3,14,:,:) = nu
         fy_u(3,14,:,:) = cyl_fac*r_faci*nu
         fy_uy(3,14,:,:) = nu
      ENDIF

      IF(ion_momentum)THEN
c-----------------------------------------------------------------------
c     x-component of ion fluid momentum equation.
c-----------------------------------------------------------------------
         fx_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(1,:,:) 
     $        - mu*two*vix_un(1,:,:)
     $        - 4.*mu_sv*vix_un(1,:,:)*ABS(vix(1,:,:)))
         fx_u(4,3,:,:)=r_fac*u(3,:,:)
         fx_u(4,4,:,:)=r_fac*two*(vi(1,:,:) - mu*nx_inv
     $        - two*mu_sv*nx_inv*ABS(vix(1,:,:)))
         fx_u(4,8,:,:)=r_fac
         fx_ux(4,1,:,:)=-r_fac*two*(mu*vi_un(1,:,:)
     $        + mu_sv*two*vi_un(1,:,:)*ABS(vix(1,:,:)))
         fx_ux(4,4,:,:)=-r_fac*two*(mu*n_inv
     $        + mu_sv*two*n_inv*ABS(vix(1,:,:)))
         
         fy_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(2,:,:)
     $        - mu*(viy_un(1,:,:) + vix_un(2,:,:)))
         fy_u(4,4,:,:)=r_fac*(vi(2,:,:) - mu*ny_inv)
         fy_u(4,5,:,:)=r_fac*(vi(1,:,:) - mu*nx_inv)
         fy_ux(4,1,:,:)=-r_fac*mu*vi_un(2,:,:)
         fy_ux(4,5,:,:)=-r_fac*mu*n_inv
         fy_uy(4,1,:,:)=-r_fac*mu*vi_un(1,:,:)
         fy_uy(4,4,:,:)=-r_fac*mu*n_inv
         
         s_u(4,1,:,:) = -r_fac*(-g_i_rho*vn(1,:,:)
     $        + g_r*vi_un(1,:,:) + g_r_rho*vi(1,:,:)
     $        + g_x*vi_un(1,:,:) + g_x_rho*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_rho(1,:,:) - rinx_rho(1,:,:) 
     $        - r_in_rho(1,:,:))
         s_u(4,4,:,:) = -r_fac*(n_inv*(g_r + g_x)
     $        + g_x_mi1*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_mi1(1,:,:) - rinx_mi1(1,:,:) - r_in_u9)
         s_u(4,5,:,:) = -r_fac*(g_x_mi2*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_mi2(1,:,:) - rinx_mi2(1,:,:))
         s_u(4,6,:,:) = -r_fac*(g_x_mi3*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_mi3(1,:,:) - rinx_mi3(1,:,:))
         s_u(4,7,:,:) = -r_fac*b2
         s_u(4,8,:,:) = -r_fac*(-g_i_p*vn(1,:,:) + g_r_p*vi(1,:,:)
     $        + g_x_p*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_p(1,:,:) - rinx_p(1,:,:) - r_in_p(1,:,:))
         s_u(4,9,:,:)=-r_fac*(-g_i*vn_un(1,:,:) - g_i_rhon*vn(1,:,:)
     $        - g_x*vn_un(1,:,:) + g_x_rhon*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_rhon(1,:,:) - rinx_rhon(1,:,:) - r_in_rhon(1,:,:))
         s_u(4,10,:,:)=r_fac*(nn_inv*(g_i + g_x)
     $        + g_x_mn1*(vn(1,:,:) - vi(1,:,:))
     $        - rnix_mn1(1,:,:) + rinx_mn1(1,:,:) + r_in_u1)
         s_u(4,11,:,:)=-r_fac*(g_x_mn2*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_mn2(1,:,:) - rinx_mn2(1,:,:))
         s_u(4,12,:,:)=-r_fac*(g_x_mn3*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_mn3(1,:,:) - rinx_mn3(1,:,:))
         s_u(4,13,:,:)=-r_fac*(g_x_pn*(vi(1,:,:) - vn(1,:,:))
     $        + rnix_pn(1,:,:) - rinx_pn(1,:,:) - r_in_pn(1,:,:))         
         s_ux(4,2,:,:)=-r_fac*j3
c-----------------------------------------------------------------------
c     y-component of ion fluid momentum equation.
c-----------------------------------------------------------------------
         fx_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(1,:,:)
     $        - mu*(vix_un(2,:,:) + viy_un(1,:,:)))
         fx_u(5,4,:,:)=r_fac*(u(5,:,:)*n_inv - mu*ny_inv)
         fx_u(5,5,:,:)=r_fac*(vi(1,:,:)-mu*nx_inv)
         fx_ux(5,1,:,:)=-r_fac*mu*vi_un(2,:,:)
         fx_ux(5,5,:,:)=-r_fac*mu*n_inv
         fx_uy(5,1,:,:)=-r_fac*mu*vi_un(1,:,:)
         fx_uy(5,4,:,:)=-r_fac*mu*n_inv
         
         fy_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(2,:,:)
     $        - mu*two*viy_un(2,:,:)
     $        - 4.*mu_sv*viy_un(2,:,:)*ABS(viy(2,:,:)))
         fy_u(5,3,:,:)=r_fac*u(3,:,:)
         fy_u(5,5,:,:)=r_fac*two*(vi(2,:,:) - mu*ny_inv
     $        - two*mu_sv*ny_inv*ABS(viy(2,:,:)))
         fy_uy(5,1,:,:)=-r_fac*two*(mu*vi_un(2,:,:)
     $        + two*mu_sv*vi_un(2,:,:)*ABS(viy(2,:,:)))
         fy_uy(5,5,:,:)=-r_fac*two*(mu*n_inv
     $        + two*mu_sv*n_inv*ABS(viy(2,:,:)))
         
         s_u(5,1,:,:) = -r_fac*(gravity 
     $        - g_i_rho*vn(2,:,:) + g_r*vi_un(2,:,:) + g_r_rho*vi(2,:,:)
     $        + g_x*vi_un(2,:,:) + g_x_rho*(vi(2,:,:) - vn(2,:,:))
     $        + rnix_rho(2,:,:) - rinx_rho(2,:,:) - r_in_rho(2,:,:))
     $        + cyl_fac*(u(6,:,:)*vi_un(3,:,:) 
     $        - r_faci*two*(mu*vi_un(2,:,:)
     $        + mu_sv*two*r_faci*vi_un(2,:,:)*ABS(vi(2,:,:))))
         s_u(5,2,:,:) = -cyl_fac*j3
         s_u(5,3,:,:) = -cyl_fac*u(3,:,:)
         s_u(5,4,:,:) = -r_fac*(g_x_mi1*(vi(2,:,:) - vn(2,:,:))
     $        + rnix_mi1(2,:,:) - rinx_mi1(2,:,:))
         s_u(5,5,:,:) = -cyl_fac*two*r_faci*n_inv
     $        *(mu + mu_sv*r_faci*two*ABS(vi(2,:,:)))
     $        - r_fac*(n_inv*(g_r + g_x)
     $        + g_x_mi2*(vi(2,:,:) - vn(2,:,:))
     $        + rnix_mi2(2,:,:) - rinx_mi2(2,:,:) - r_in_u9)
         s_u(5,6,:,:) = cyl_fac*two*vi(3,:,:) 
     $        - r_fac*(g_x_mi3*(vi(2,:,:)-vn(2,:,:))
     $        + rnix_mi3(2,:,:) - rinx_mi3(2,:,:))
         s_u(5,7,:,:) = r_fac*b1
         s_u(5,8,:,:) = -r_fac*(-g_i_p*vn(2,:,:) + g_r_p*vi(2,:,:)
     $        + g_x_p*(vi(2,:,:) - vn(2,:,:))
     $        + rnix_p(2,:,:) - rinx_p(2,:,:) - r_in_p(2,:,:))
         s_u(5,9,:,:) = -r_fac*(-g_i*vn_un(2,:,:) - g_i_rhon*vn(2,:,:)
     $        - g_x*vn_un(2,:,:) + g_x_rhon*(vi(2,:,:) - vn(2,:,:))
     $        + rnix_rhon(2,:,:) - rinx_rhon(2,:,:) - r_in_rhon(2,:,:))
         s_u(5,10,:,:) = -r_fac*(g_x_mn1*(vi(2,:,:) - vn(2,:,:))
     $        + rnix_mn1(2,:,:) - rinx_mn1(2,:,:))
         s_u(5,11,:,:) = r_fac*(nn_inv*(g_i + g_x)
     $        + g_x_mn2*(vn(2,:,:) - vi(2,:,:))
     $        - rnix_mn2(2,:,:) + rinx_mn2(2,:,:) + r_in_u1)
         s_u(5,12,:,:) = -r_fac*(g_x_mn3*(vi(2,:,:) - vn(2,:,:))
     $        + rnix_mn3(2,:,:) - rinx_mn3(2,:,:))
         s_u(5,13,:,:) = -r_fac*(g_x_pn*(vi(2,:,:) - vn(2,:,:))
     $        + rnix_pn(2,:,:) - rinx_pn(2,:,:) - r_in_pn(2,:,:))
         
         s_uy(5,2,:,:)=-r_fac*j3
         s_uy(5,8,:,:)=-r_fac
         s_uy(5,13,:,:)=-r_fac
c-----------------------------------------------------------------------
c     z-component of ion fluid momentum equation.
c-----------------------------------------------------------------------
         fx_u(6,1,:,:)=r_fac*(u(6,:,:)*vi_un(1,:,:) - mu*vix_un(3,:,:))
         fx_u(6,4,:,:)=r_fac*u(6,:,:)*n_inv
         fx_u(6,6,:,:)=r_fac*(vi(1,:,:) - mu*nx_inv)
         fx_ux(6,1,:,:)=-r_fac*mu*vi_un(3,:,:)
         fx_ux(6,6,:,:)=-r_fac*mu*n_inv
         
         fy_u(6,1,:,:)=r_fac*(u(6,:,:)*vi_un(2,:,:) - mu*viy_un(3,:,:))
         fy_u(6,5,:,:)=r_fac*u(6,:,:)*n_inv
         fy_u(6,6,:,:)=r_fac*(vi(2,:,:) - mu*ny_inv)
         fy_uy(6,1,:,:)=-r_fac*mu*vi_un(3,:,:)
         fy_uy(6,6,:,:)=-r_fac*mu*n_inv
         
         s_u(6,1,:,:) = cyl_fac*(vi(2,:,:)*vi(3,:,:) 
     $        - r_faci*mu*vi_un(3,:,:)) - r_fac*(g_r*vi_un(3,:,:) 
     $        + g_r_rho*vi(3,:,:) - g_i_rho*vn(3,:,:) 
     $        + g_x*vi_un(3,:,:) + g_x_rho*(vi(3,:,:) - vn(3,:,:))
     $        + rnix_rho(3,:,:) - rinx_rho(3,:,:) - r_in_rho(3,:,:))
         s_u(6,2,:,:) = cyl_fac*j2
         s_u(6,3,:,:) = cyl_fac*b2
         s_u(6,4,:,:) = -r_fac*(g_x_mi1*(vi(3,:,:) - vn(3,:,:))
     $        + rnix_mi1(3,:,:) - rinx_mi1(3,:,:))
         s_u(6,5,:,:) = -cyl_fac*vi(3,:,:)
     $        - r_fac*(g_x_mi2*(vi(3,:,:) - vn(3,:,:))
     $        + rnix_mi2(3,:,:) - rinx_mi2(3,:,:))
         s_u(6,6,:,:) = -cyl_fac*(vi(2,:,:) + r_faci*mu*n_inv) 
     $        - r_fac*(n_inv*(g_r + g_x)
     $        + g_x_mi3*(vi(3,:,:) - vn(3,:,:))
     $        + rnix_mi3(3,:,:) - rinx_mi3(3,:,:) - r_in_u9)
         s_u(6,8,:,:) = -r_fac*(g_r_p*vi(3,:,:) - g_i_p*vn(3,:,:)
     $        + g_x_p*(vi(3,:,:) - vn(3,:,:))
     $        + rnix_p(3,:,:) - rinx_p(3,:,:) - r_in_p(3,:,:))
         s_u(6,9,:,:) = r_fac*(g_i*vn_un(3,:,:) + g_i_rhon*vn(3,:,:)
     $        + g_x*vn_un(3,:,:) + g_x_rhon*(vn(3,:,:) - vi(3,:,:))
     $        - rnix_rhon(3,:,:) + rinx_rhon(3,:,:) + r_in_rhon(3,:,:))
         s_u(6,10,:,:) = -r_fac*(g_x_mn1*(vi(3,:,:) - vn(3,:,:))
     $        + rnix_mn1(3,:,:) - rinx_mn1(3,:,:))
         s_u(6,11,:,:) = -r_fac*(g_x_mn2*(vi(3,:,:) - vn(3,:,:))
     $        + rnix_mn2(3,:,:) - rinx_mn2(3,:,:))
         s_u(6,12,:,:) = -r_fac*(g_x_mn3*(vi(3,:,:) - vn(3,:,:))
     $        - nn_inv*(g_i + g_x)
     $        + rnix_mn3(3,:,:) - rinx_mn3(3,:,:) - r_in_u1)
         s_u(6,13,:,:) = -r_fac*(g_x_pn*(vi(3,:,:) - vn(3,:,:))
     $        + rnix_pn(3,:,:) - rinx_pn(3,:,:) - r_in_pn(3,:,:))
         
         s_ux(6,2,:,:)=r_fac*j1
         s_ux(6,3,:,:)=r_fac*b1
         s_uy(6,2,:,:)=r_fac*j2
         s_uy(6,3,:,:)=r_fac*b2
      ELSE
c-----------------------------------------------------------------------
c     x-component of total fluid momentum equation.
c-----------------------------------------------------------------------
         fx_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(1,:,:) 
     $        - mu*two*vix_un(1,:,:)
     $        - 4.*mu_sv*vix_un(1,:,:)*ABS(vix(1,:,:)))
         fx_u(4,3,:,:)=r_fac*u(3,:,:)
         fx_u(4,4,:,:)=r_fac*two*(vi(1,:,:) - mu*nx_inv
     $        - two*mu_sv*nx_inv*ABS(vix(1,:,:)))
         fx_u(4,8,:,:)=r_fac
         fx_u(4,9,:,:)=r_fac*(u(10,:,:)*vn_un(1,:,:)
     $        - two*(mun*vnx_un(1,:,:) + mun_un*vnx(1,:,:)))
         fx_u(4,10,:,:)=r_fac*two*(vn(1,:,:) - mun*nnx_inv)
         fx_u(4,13,:,:)=r_fac*(one - two*mun_pn*vnx(1,:,:))
         fx_ux(4,1,:,:)=-r_fac*two*(mu*vi_un(1,:,:)
     $        + mu_sv*two*vi_un(1,:,:)*ABS(vix(1,:,:)))
         fx_ux(4,4,:,:)=-r_fac*two*(mu*n_inv
     $        + mu_sv*two*n_inv*ABS(vix(1,:,:)))
         fx_ux(4,9,:,:)=-r_fac*mun*two*vn_un(1,:,:)
         fx_ux(4,10,:,:)=-r_fac*mun*two*nn_inv
         
         fy_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(2,:,:)
     $        - mu*(viy_un(1,:,:) + vix_un(2,:,:)))
         fy_u(4,4,:,:)=r_fac*(vi(2,:,:) - mu*ny_inv)
         fy_u(4,5,:,:)=r_fac*(vi(1,:,:) - mu*nx_inv)
         fy_u(4,9,:,:)=r_fac*(u(10,:,:)*vn_un(2,:,:) 
     $        - mun*(vny_un(1,:,:) + vnx_un(2,:,:)) 
     $        - mun_un*(vny(1,:,:) + vnx(2,:,:)))
         fy_u(4,10,:,:)=r_fac*(vn(2,:,:) - mun*nny_inv)
         fy_u(4,11,:,:)=r_fac*(vn(1,:,:) - mun*nnx_inv)
         fy_u(4,13,:,:)=-r_fac*mun_pn*(vny(1,:,:) + vnx(2,:,:))
         fy_ux(4,1,:,:)=-r_fac*mu*vi_un(2,:,:)
         fy_ux(4,5,:,:)=-r_fac*mu*n_inv
         fy_ux(4,9,:,:)=-r_fac*mun*vn_un(2,:,:)
         fy_ux(4,11,:,:)=-r_fac*mun*nn_inv
         fy_uy(4,1,:,:)=-r_fac*mu*vi_un(1,:,:)
         fy_uy(4,4,:,:)=-r_fac*mu*n_inv
         fy_uy(4,9,:,:)=-r_fac*mun*vn_un(1,:,:)
         fy_uy(4,10,:,:)=-r_fac*mun*nn_inv
         
         s_u(4,1,:,:) = -r_fac*mfp_fac
     $        *(g_x*vi_un(1,:,:) + g_x_rho*vi(1,:,:) + rnix_rho(1,:,:))
         s_u(4,4,:,:) = -r_fac*mfp_fac
     $        *(n_inv*g_x + g_x_mi1*vi(1,:,:) + rnix_mi1(1,:,:))
         s_u(4,5,:,:) = -r_fac*mfp_fac
     $        *(g_x_mi2*vi(1,:,:) + rnix_mi2(1,:,:))
         s_u(4,6,:,:) = -r_fac*mfp_fac
     $        *(g_x_mi3*vi(1,:,:) + rnix_mi3(1,:,:))
         s_u(4,7,:,:) = -r_fac*b2
         s_u(4,8,:,:) = -r_fac*mfp_fac*(g_x_p*vi(1,:,:) + rnix_p(1,:,:))
         s_u(4,9,:,:) = -r_fac*mfp_fac
     $        *(g_x_rhon*vi(1,:,:) + rnix_rhon(1,:,:))
         s_u(4,10,:,:) = -r_fac*mfp_fac
     $        *(g_x_mn1*vi(1,:,:) + rnix_mn1(1,:,:))
         s_u(4,11,:,:) = -r_fac*mfp_fac
     $        *(g_x_mn2*vi(1,:,:) + rnix_mn2(1,:,:))
         s_u(4,12,:,:) = -r_fac*mfp_fac
     $        *(g_x_mn3*vi(1,:,:) + rnix_mn3(1,:,:))
         s_u(4,13,:,:) = -r_fac*mfp_fac
     $        *(g_x_pn*vi(1,:,:) + rnix_pn(1,:,:))
         
         s_ux(4,2,:,:)=-r_fac*j3
c-----------------------------------------------------------------------
c     y-component of total fluid momentum equation.
c-----------------------------------------------------------------------
         fx_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(1,:,:)
     $        - mu*(vix_un(2,:,:) + viy_un(1,:,:)))
         fx_u(5,4,:,:)=r_fac*(u(5,:,:)*n_inv - mu*ny_inv)
         fx_u(5,5,:,:)=r_fac*(vi(1,:,:)-mu*nx_inv)
         fx_u(5,9,:,:)=r_fac*(u(11,:,:)*vn_un(1,:,:) 
     $        - mun*(vnx_un(2,:,:) + vny_un(1,:,:))
     $        - mun_un*(vnx(2,:,:) + vny(1,:,:)))
         fx_u(5,10,:,:)=r_fac*(vn(2,:,:) - mun*nny_inv)
         fx_u(5,11,:,:)=r_fac*(vn(1,:,:) - mun*nnx_inv)
         fx_u(5,13,:,:)=-r_fac*mun_pn*(vnx(2,:,:) + vny(1,:,:))
         fx_ux(5,1,:,:)=-r_fac*mu*vi_un(2,:,:)
         fx_ux(5,5,:,:)=-r_fac*mu*n_inv
         fx_ux(5,9,:,:)=-r_fac*mun*vn_un(2,:,:)
         fx_ux(5,11,:,:)=-r_fac*mun*nn_inv
         fx_uy(5,1,:,:)=-r_fac*mu*vi_un(1,:,:)
         fx_uy(5,4,:,:)=-r_fac*mu*n_inv
         fx_uy(5,9,:,:)=-r_fac*mun*vn_un(1,:,:)
         fx_uy(5,10,:,:)=-r_fac*mun*nn_inv
         
         fy_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(2,:,:)
     $        - mu*two*viy_un(2,:,:)
     $        - 4.*mu_sv*viy_un(2,:,:)*ABS(viy(2,:,:)))
         fy_u(5,3,:,:)=r_fac*u(3,:,:)
         fy_u(5,5,:,:)=r_fac*two*(vi(2,:,:) - mu*ny_inv
     $        - two*mu_sv*ny_inv*ABS(viy(2,:,:)))
         fy_u(5,9,:,:)=r_fac*(u(11,:,:)*vn_un(2,:,:) 
     $        - two*(mun*vny_un(2,:,:) + mun_un*vny(2,:,:)))
         fy_u(5,11,:,:)=r_fac*two*(vn(2,:,:) - mun*nny_inv)
         fy_u(5,13,:,:)=-r_fac*two*mun_pn*vny(2,:,:)
         fy_uy(5,1,:,:)=-r_fac*two*(mu*vi_un(2,:,:)
     $        + two*mu_sv*vi_un(2,:,:)*ABS(viy(2,:,:)))
         fy_uy(5,5,:,:)=-r_fac*two*(mu*n_inv
     $        + two*mu_sv*n_inv*ABS(viy(2,:,:)))
         fy_uy(5,9,:,:)=-r_fac*mun*two*vn_un(2,:,:)
         fy_uy(5,11,:,:)=-r_fac*mun*two*nn_inv
         
         s_u(5,1,:,:) = -r_fac*gravity + cyl_fac*(u(6,:,:)*vi_un(3,:,:) 
     $        - r_faci*two*(mu*vi_un(2,:,:)
     $        + mu_sv*two*r_faci*vi_un(2,:,:)*ABS(vi(2,:,:))))
     $        - r_fac*mfp_fac*(g_x*vi_un(2,:,:) + g_x_rho*vi(2,:,:)
     $        + rnix_rho(2,:,:))
         s_u(5,2,:,:) = -cyl_fac*j3
         s_u(5,3,:,:) = -cyl_fac*u(3,:,:)
         s_u(5,4,:,:) = -r_fac*mfp_fac
     $        *(g_x_mi1*vi(2,:,:) + rnix_mi1(2,:,:))
         s_u(5,5,:,:) = -cyl_fac*two*r_faci*n_inv
     $        *(mu + mu_sv*r_faci*two*ABS(vi(2,:,:)))
     $        - r_fac*mfp_fac
     $        *(g_x*n_inv + g_x_mi2*vi(2,:,:) + rnix_mi2(2,:,:))
         s_u(5,6,:,:) = cyl_fac*two*vi(3,:,:) 
     $        - r_fac*mfp_fac*(g_x_mi3*vi(2,:,:) + rnix_mi3(2,:,:))
         s_u(5,7,:,:) = r_fac*b1
         s_u(5,8,:,:) = -r_fac*mfp_fac*(g_x_p*vi(2,:,:) + rnix_p(2,:,:))
         s_u(5,9,:,:) = -r_fac*gravity + cyl_fac*(u(12,:,:)*vn_un(3,:,:) 
     $        - r_faci*two*(mun*vn_un(2,:,:) + mun_un*vn(2,:,:)))
     $        - r_fac*mfp_fac*(g_x_rhon*vi(2,:,:) + rnix_rhon(2,:,:))
         s_u(5,10,:,:) = -r_fac*mfp_fac
     $        *(g_x_mn1*vi(2,:,:) + rnix_mn1(2,:,:))
         s_u(5,11,:,:) = -cyl_fac*two*mun*nn_inv*r_faci
     $        - r_fac*mfp_fac*(g_x_mn2*vi(2,:,:) + rnix_mn2(2,:,:))
         s_u(5,12,:,:) = cyl_fac*two*vn(3,:,:)
     $        - r_fac*mfp_fac*(g_x_mn3*vi(2,:,:) + rnix_mn3(2,:,:))
         s_u(5,13,:,:) = -r_fac*mfp_fac*(g_x_pn*vi(2,:,:)
     $        + rnix_pn(2,:,:)) - cyl_fac*two*r_faci*mun_pn*vn(2,:,:)
         
         s_uy(5,2,:,:)=-r_fac*j3
         s_uy(5,8,:,:)=-r_fac
         s_uy(5,13,:,:)=-r_fac
c-----------------------------------------------------------------------
c     z-component of total fluid momentum equation.
c-----------------------------------------------------------------------
         fx_u(6,1,:,:)=r_fac*(u(6,:,:)*vi_un(1,:,:) - mu*vix_un(3,:,:))
         fx_u(6,4,:,:)=r_fac*u(6,:,:)*n_inv
         fx_u(6,6,:,:)=r_fac*(vi(1,:,:) - mu*nx_inv)
         fx_u(6,9,:,:)=r_fac*(u(12,:,:)*vn_un(1,:,:) 
     $        - mun*vnx_un(3,:,:) - mun_un*vnx(3,:,:))
         fx_u(6,10,:,:)=r_fac*u(12,:,:)*nn_inv      
         fx_u(6,12,:,:)=r_fac*(vn(1,:,:) - mun*nnx_inv)
         fx_u(6,13,:,:)=-r_fac*mun_pn*vnx(3,:,:)
         fx_ux(6,1,:,:)=-r_fac*mu*vi_un(3,:,:)
         fx_ux(6,6,:,:)=-r_fac*mu*n_inv
         fx_ux(6,9,:,:)=-r_fac*mun*vn_un(3,:,:)
         fx_ux(6,12,:,:)=-r_fac*mun*nn_inv
         
         fy_u(6,1,:,:)=r_fac*(u(6,:,:)*vi_un(2,:,:) - mu*viy_un(3,:,:))
         fy_u(6,5,:,:)=r_fac*u(6,:,:)*n_inv
         fy_u(6,6,:,:)=r_fac*(vi(2,:,:) - mu*ny_inv)
         fy_u(6,9,:,:)=r_fac*(u(12,:,:)*vn_un(2,:,:) 
     $        - mun*vny_un(3,:,:) - mun_un*vny(3,:,:))
         fy_u(6,11,:,:)=r_fac*u(12,:,:)*nn_inv
         fy_u(6,12,:,:)=r_fac*(vn(2,:,:) - mun*nny_inv)
         fy_u(6,13,:,:)=-r_fac*mun_pn*vny(3,:,:)
         fy_uy(6,1,:,:)=-r_fac*mu*vi_un(3,:,:)
         fy_uy(6,6,:,:)=-r_fac*mu*n_inv
         fy_uy(6,9,:,:)=-r_fac*mun*vn_un(3,:,:)
         fy_uy(6,12,:,:)=-r_fac*mun*nn_inv
         
         s_u(6,1,:,:) = cyl_fac*(vi(2,:,:)*vi(3,:,:) 
     $        - r_faci*mu*vi_un(3,:,:)) - r_fac*mfp_fac
     $        *(g_x*vi_un(3,:,:) + g_x_rho*vi(3,:,:) + rnix_rho(3,:,:))
         s_u(6,2,:,:) = cyl_fac*j2
         s_u(6,3,:,:) = cyl_fac*b2
         s_u(6,4,:,:) = -r_fac*mfp_fac
     $        *(g_x_mi1*vi(3,:,:) + rnix_mi1(3,:,:))
         s_u(6,5,:,:) = -cyl_fac*vi(3,:,:)
     $        - r_fac*mfp_fac*(g_x_mi2*vi(3,:,:) + rnix_mi2(3,:,:))
         s_u(6,6,:,:) = -cyl_fac*(vi(2,:,:) + r_faci*mu*n_inv) 
     $        - r_fac*mfp_fac
     $        *(n_inv*g_x + g_x_mi3*vi(3,:,:) + rnix_mi3(3,:,:))
         s_u(6,8,:,:) = -r_fac*mfp_fac
     $        *(g_x_p*vi(3,:,:) + rnix_p(3,:,:))
         s_u(6,9,:,:) = cyl_fac*(vn(2,:,:)*vn(3,:,:) 
     $        - r_faci*(mun*vn_un(3,:,:) + mun_un*vn(3,:,:)))
     $        - r_fac*mfp_fac*(g_x_rhon*vi(3,:,:) + rnix_rhon(3,:,:))
         s_u(6,10,:,:) = -r_fac*mfp_fac
     $        *(g_x_mn1*vi(3,:,:) + rnix_mn1(3,:,:))
         s_u(6,11,:,:) = -cyl_fac*vn(3,:,:)
     $        - r_fac*mfp_fac*(g_x_mn2*vi(3,:,:) + rnix_mn2(3,:,:))
         s_u(6,12,:,:) = -cyl_fac*(vn(2,:,:) + r_faci*mun*nn_inv)
     $        - r_fac*mfp_fac*(g_x_mn3*vi(3,:,:) + rnix_mn3(3,:,:))
         s_u(6,13,:,:) = -r_fac*mfp_fac*(g_x_pn*vi(3,:,:)
     $        + rnix_pn(3,:,:)) - cyl_fac*r_faci*mun_pn*vn(3,:,:)
         
         s_ux(6,2,:,:)=r_fac*j1
         s_ux(6,3,:,:)=r_fac*b1
         s_uy(6,2,:,:)=r_fac*j2
         s_uy(6,3,:,:)=r_fac*b2
      ENDIF
c-----------------------------------------------------------------------
c     out-of-plane electron momentum (current) equation.
c-----------------------------------------------------------------------
      fx_ux(7,2,:,:)=one
      fy_u(7,2,:,:)=cyl_fac*r_faci
      fy_uy(7,2,:,:)=one
      s_u(7,7,:,:)=one
c-----------------------------------------------------------------------
c     plasma pressure equation.
c-----------------------------------------------------------------------
      fx_u(8,1,:,:)=r_fac*(gamma_fac*u(8,:,:)
     $     *(vi_un(1,:,:) + te_frac*di*j1*n_inv**2)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*Tix_un 
     $     - kfac*BdotT_Ty(1,:,:)*Tiy_un
     $     - (kpar_rho - kperp_rho)*BdotT(1,:,:) - kperp_rho*Tix)
      fx_u(8,2,:,:)=-cyl_fac*(-kfac*BdotT_b1(1,:,:)
     $     + two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_u(8,3,:,:)=r_fac*(-kfac*BdotT_b3(1,:,:)
     $     + two*u(3,:,:)*kperp_bsq*(BdotT(1,:,:) - Tix))
     $     - cyl_fac*gamma_fac*di*Te
      fx_u(8,4,:,:)=r_fac*gamma_fac*u(8,:,:)*n_inv
      fx_u(8,8,:,:)=r_fac*(gamma_fac*(vi(1,:,:) - te_frac*di*j1*n_inv)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*nx_inv
     $     - kfac*BdotT_Ty(1,:,:)*ny_inv
     $     - (kpar_p - kperp_p)*BdotT(1,:,:) - kperp_p*Tix)

      fx_ux(8,1,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*Ti_un
      fx_ux(8,2,:,:)=r_fac*(-kfac*BdotT_b2(1,:,:)
     $     + two*b2*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_ux(8,8,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*n_inv

      fx_uy(8,1,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*Ti_un
      fx_uy(8,2,:,:)=-r_fac*(-kfac*BdotT_b1(1,:,:)
     $     + two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_uy(8,3,:,:)=-r_fac*gamma_fac*di*Te
      fx_uy(8,8,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*n_inv

      fy_u(8,1,:,:)=r_fac*(gamma_fac*u(8,:,:)
     $     *(vi_un(2,:,:) + te_frac*di*j2*n_inv**2)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*Tiy_un
     $     - kfac*BdotT_Tx(2,:,:)*Tix_un
     $     - (kpar_rho - kperp_rho)*BdotT(2,:,:) - kperp_rho*Tiy)
      fy_u(8,2,:,:)=-cyl_fac*(-kfac*BdotT_b1(2,:,:)
     $     + two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(8,3,:,:)=r_fac*(-kfac*BdotT_b3(2,:,:)
     $     + two*u(3,:,:)*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(8,5,:,:)=r_fac*gamma_fac*u(8,:,:)*n_inv
      fy_u(8,8,:,:)=r_fac*(gamma_fac*(vi(2,:,:) - te_frac*di*j2*n_inv)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*ny_inv
     $     - kfac*BdotT_Tx(2,:,:)*nx_inv
     $     - (kpar_p - kperp_p)*BdotT(2,:,:) - kperp_p*Tiy)

      fy_ux(8,1,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*Ti_un
      fy_ux(8,2,:,:)=r_fac*(-kfac*BdotT_b2(2,:,:)
     $     + two*b2*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_ux(8,3,:,:)=r_fac*gamma_fac*di*Te
      fy_ux(8,8,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*n_inv

      fy_uy(8,1,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*Ti_un
      fy_uy(8,2,:,:)=-r_fac*(-kfac*BdotT_b1(2,:,:)
     $     + two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_uy(8,8,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*n_inv

      s_u(8,1,:,:)=r_fac
     $     *(ux(8,:,:)*(vi_un(1,:,:) + te_frac*di*j1*n_inv**2)
     $     + uy(8,:,:)*(vi_un(2,:,:) + te_frac*di*j2*n_inv**2)
     $     + (eta_rho - di*r_en_ut*n_inv*(1.5_r8*Te+Tn))*jtot**2 
     $     + mu*two*(two*vix(1,:,:)*vix_un(1,:,:) 
     $     + two*viy(2,:,:)*viy_un(2,:,:) 
     $     + viy(1,:,:)*viy_un(1,:,:) + vix(2,:,:)*vix_un(2,:,:)
     $     + viy(1,:,:)*vix_un(2,:,:) + vix(2,:,:)*viy_un(1,:,:)
     $     + vix(3,:,:)*vix_un(3,:,:) + viy(3,:,:)*viy_un(3,:,:))
     $     + cyl_fac*two*r_faci**2*mu
     $     *(two*vi(2,:,:)*vi_un(2,:,:) + vi(3,:,:)*vi_un(3,:,:))
     $     + mu_sv*(6.*vix(1,:,:)*vix_un(1,:,:)*ABS(vix(1,:,:))
     $     + 6.*viy(2,:,:)*viy_un(2,:,:)*ABS(viy(2,:,:)))
     $     + cyl_fac*mu_sv*6.*r_faci**3*vi(2,:,:)*vi_un(2,:,:)
     $     *ABS(vi(2,:,:))
     $     + (g_i + g_x)*kediff_rho + (g_i_rho + g_x_rho)*kediff 
     $     - g_i_rho*phi_eff + 1.5_r8*(g_i_rho*u(13,:,:)*nn_inv
     $     - 1.5_r8*ti_frac*g_r*u(8,:,:)*n_inv**2)
     $     + n_inv*(r_inx(1,:,:)*vi(1,:,:) + r_inx(2,:,:)*vi(2,:,:)
     $     + r_inx(3,:,:)*vi(3,:,:))
     $     + rinx_rho(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_rho(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + rinx_rho(3,:,:)*(vn(3,:,:) - vi(3,:,:)) 
     $     + qinx_rho - qnix_rho + q_in_8_1 + qt_in_8_1)
      s_u(8,3,:,:)=r_fac*(eta_local*two + eta_j*jtot + di*r_en*two)
     $     *jtot*j_u3 - cyl_fac*te_frac*di*ux(8,:,:)*n_inv
      s_u(8,4,:,:)=r_fac*(ux(8,:,:)*n_inv
     $     + mu*(4._r8*vix(1,:,:)*nx_inv
     $     + two*viy(1,:,:)*ny_inv + two*vix(2,:,:)*ny_inv)
     $     + mu_sv*6.*vix(1,:,:)*nx_inv*ABS(vix(1,:,:))
     $     + (g_i + g_x)*kediff_mi(1,:,:)
     $     + g_x_mi1*kediff - r_inx(1,:,:)*n_inv
     $     + rinx_mi1(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mi1(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + rinx_mi1(3,:,:)*(vn(3,:,:) - vi(3,:,:)) 
     $     + qinx_mi(1,:,:) - qnix_mi(1,:,:) + q_in_vi(1,:,:))
      s_u(8,5,:,:)=r_fac*(uy(8,:,:)*n_inv
     $     + mu*(4._r8*viy(2,:,:)*ny_inv
     $     + two*vix(2,:,:)*nx_inv + two*viy(1,:,:)*nx_inv)
     $     + cyl_fac*4._r8*mu*r_faci**2*vi(2,:,:)*n_inv
     $     + mu_sv*6.*viy(2,:,:)*ny_inv*ABS(viy(2,:,:))
     $     + cyl_fac*6.*mu_sv*r_faci
     $     *SIGN((vi(2,:,:)*r_faci)**2*n_inv,vi(2,:,:))
     $     + (g_i + g_x)*kediff_mi(2,:,:)
     $     + g_x_mi2*kediff - r_inx(2,:,:)*n_inv
     $     + rinx_mi2(1,:,:)*(vn(1,:,:) - vi(1,:,:)) 
     $     + rinx_mi2(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + rinx_mi2(3,:,:)*(vn(3,:,:) - vi(3,:,:)) 
     $     + qinx_mi(2,:,:) - qnix_mi(2,:,:) + q_in_vi(2,:,:))
      s_u(8,6,:,:)=r_fac*(two*(mu*(vix(3,:,:)*nx_inv 
     $     + viy(3,:,:)*ny_inv) + cyl_fac*r_faci**2*mu*vi(3,:,:)*n_inv)
     $     + (g_i + g_x)*kediff_mi(3,:,:)
     $     + g_x_mi3*kediff - r_inx(3,:,:)*n_inv
     $     + rinx_mi3(1,:,:)*(vn(1,:,:) - vi(1,:,:)) 
     $     + rinx_mi3(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + rinx_mi3(3,:,:)*(vn(3,:,:) - vi(3,:,:))
     $     + qinx_mi(3,:,:) - qnix_mi(3,:,:) + q_in_vi(3,:,:))
      s_u(8,7,:,:)=r_fac*(eta_local*two + eta_j*jtot + di*r_en*two)*j3
      s_u(8,8,:,:)=r_fac*(
     $     (eta_p + half*te_frac*di*r_en_ut*n_inv)*jtot**2
     $     + (g_i_p + g_x_p)*kediff - g_i_p*phi_eff
     $     + 1.5_r8*(g_i_p*u(13,:,:)*nn_inv 
     $     - ti_frac*(g_r_p*u(8,:,:)*n_inv + g_r*n_inv))
     $     + rinx_p(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_p(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + rinx_p(3,:,:)*(vn(3,:,:) - vi(3,:,:)) 
     $     + qinx_p - qnix_p + q_in_8_8 + qt_in_8_8)
      s_u(8,9,:,:)=r_fac*(di*r_en_ut*nn_inv*(half*Tn + Te)*jtot**2 
     $     - g_i_rhon*phi_eff
     $     + (g_i + g_x)*kediff_rhon + (g_i_rhon + g_x_rhon)*kediff
     $     + 1.5_r8*u(13,:,:)*nn_inv*(g_i_rhon - g_i*nn_inv)
     $     - nn_inv*(r_inx(1,:,:)*vn(1,:,:) + r_inx(2,:,:)*vn(2,:,:)
     $     + r_inx(3,:,:)*vn(3,:,:))
     $     + rinx_rhon(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_rhon(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + rinx_rhon(3,:,:)*(vn(3,:,:) - vi(3,:,:)) 
     $     + qinx_rhon - qnix_rhon + q_in_8_9 + qt_in_8_9)
      s_u(8,10,:,:)=r_fac*((g_i + g_x)*kediff_mn(1,:,:)
     $     + g_x_mn1*kediff + r_inx(1,:,:)*nn_inv
     $     + rinx_mn1(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mn1(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + rinx_mn1(3,:,:)*(vn(3,:,:) - vi(3,:,:)) 
     $     + qinx_mn(1,:,:) - qnix_mn(1,:,:) + q_in_vn(1,:,:))
      s_u(8,11,:,:)=r_fac*((g_i + g_x)*kediff_mn(2,:,:)
     $     + g_x_mn2*kediff + r_inx(2,:,:)*nn_inv
     $     + rinx_mn2(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mn2(2,:,:)*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_mn2(3,:,:)*(vn(3,:,:) - vi(3,:,:))
     $     + qinx_mn(2,:,:) - qnix_mn(2,:,:) + q_in_vn(2,:,:))
      s_u(8,12,:,:)=r_fac*((g_i + g_x)*kediff_mn(3,:,:)
     $     + g_x_mn3*kediff + r_inx(3,:,:)*nn_inv
     $     + rinx_mn3(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mn3(2,:,:)*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_mn3(3,:,:)*(vn(3,:,:) - vi(3,:,:))
     $     + qinx_mn(3,:,:) - qnix_mn(3,:,:) + q_in_vn(3,:,:))
      s_u(8,13,:,:)=r_fac*(half*di*r_en_ut*nn_inv*jtot**2
     $     + g_x_pn*kediff + 1.5_r8*g_i*nn_inv
     $     + rinx_pn(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_pn(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + rinx_pn(3,:,:)*(vn(3,:,:) - vi(3,:,:)) 
     $     + qinx_pn - qnix_pn + q_in_8_13 + qt_in_8_13)
      
      s_ux(8,1,:,:)=r_fac*(mu*two*(two*vix(1,:,:)*vi_un(1,:,:)  
     $     + vi_un(2,:,:)*(vix(2,:,:) + viy(1,:,:))
     $     + vix(3,:,:)*vi_un(3,:,:))
     $     + mu_sv*6.*vix(1,:,:)*vi_un(1,:,:)*ABS(vix(1,:,:)))
      s_ux(8,3,:,:)=-r_fac*((eta_local*two + eta_j*jtot + di*r_en*two)
     $     *j2 - te_frac*di*uy(8,:,:)*n_inv)
      s_ux(8,4,:,:)=r_fac*(mu*4._r8*vix(1,:,:)*n_inv
     $     + mu_sv*6.*vix(1,:,:)*n_inv*ABS(vix(1,:,:)))
      s_ux(8,5,:,:)=r_fac*mu*two*n_inv*(vix(2,:,:) + viy(1,:,:))
      s_ux(8,6,:,:)=r_fac*two*mu*vix(3,:,:)*n_inv
      s_ux(8,8,:,:)=r_fac*(vi(1,:,:) - te_frac*di*j1*n_inv)

      s_uy(8,1,:,:)=r_fac*(mu*two*(two*viy(2,:,:)*vi_un(2,:,:)  
     $     + vi_un(1,:,:)*(viy(1,:,:) + vix(2,:,:))
     $     + viy(3,:,:)*vi_un(3,:,:))
     $     + mu_sv*6.*viy(2,:,:)*vi_un(2,:,:)*ABS(viy(2,:,:)))
      s_uy(8,3,:,:)=r_fac*((eta_local*two + eta_j*jtot + di*r_en*two)
     $     *j1 - te_frac*di*ux(8,:,:)*n_inv)
      s_uy(8,4,:,:)=r_fac*mu*two*n_inv*(viy(1,:,:) + vix(2,:,:))
      s_uy(8,5,:,:)=r_fac*(mu*4._r8*viy(2,:,:)*n_inv
     $     + mu_sv*6.*viy(2,:,:)*n_inv*ABS(viy(2,:,:)))
      s_uy(8,6,:,:)=r_fac*two*mu*viy(3,:,:)*n_inv
      s_uy(8,8,:,:)=r_fac*(vi(2,:,:) - te_frac*di*j2*n_inv)
c-----------------------------------------------------------------------
c     total fluid density.
c-----------------------------------------------------------------------
      fx_u(9,4,:,:) = r_fac
      fx_u(9,10,:,:) = r_fac
      fx_ux(9,1,:,:) = -r_fac*Dn
      fx_ux(9,9,:,:) = -r_fac*Dn
      fy_u(9,5,:,:) = r_fac
      fy_u(9,11,:,:) = r_fac
      fy_uy(9,1,:,:) = -r_fac*Dn
      fy_uy(9,9,:,:) = -r_fac*Dn
      
      s_u(9,1,:,:) = -r_fac*mfp_fac*g_x_rho
      s_u(9,4,:,:) = -r_fac*mfp_fac*g_x_mi1
      s_u(9,5,:,:) = -r_fac*mfp_fac*g_x_mi2
      s_u(9,6,:,:) = -r_fac*mfp_fac*g_x_mi3
      s_u(9,8,:,:) = -r_fac*mfp_fac*g_x_p
      s_u(9,9,:,:) = -r_fac*mfp_fac*g_x_rhon
      s_u(9,10,:,:) = -r_fac*mfp_fac*g_x_mn1
      s_u(9,11,:,:) = -r_fac*mfp_fac*g_x_mn2
      s_u(9,12,:,:) = -r_fac*mfp_fac*g_x_mn3
      s_u(9,13,:,:) = -r_fac*mfp_fac*g_x_pn
c-----------------------------------------------------------------------
c     x-component of neutral momentum equation.
c-----------------------------------------------------------------------
      fx_u(10,9,:,:)=r_fac*(u(10,:,:)*vn_un(1,:,:)
     $     - two*(mun*vnx_un(1,:,:) + mun_un*vnx(1,:,:)))
      fx_u(10,10,:,:)=r_fac*two*(vn(1,:,:) - mun*nnx_inv)
      fx_u(10,13,:,:)=r_fac*(one - two*mun_pn*vnx(1,:,:))
      fx_ux(10,9,:,:)=-r_fac*mun*two*vn_un(1,:,:)
      fx_ux(10,10,:,:)=-r_fac*mun*two*nn_inv
      
      fy_u(10,9,:,:)=r_fac*(u(10,:,:)*vn_un(2,:,:) 
     $     - mun*(vny_un(1,:,:) + vnx_un(2,:,:))
     $     - mun_un*(vny(1,:,:) + vnx(2,:,:)))
      fy_u(10,10,:,:)=r_fac*(vn(2,:,:) - mun*nny_inv)
      fy_u(10,11,:,:)=r_fac*(vn(1,:,:) - mun*nnx_inv)
      fy_u(10,13,:,:)=-r_fac*mun_pn*(vny(1,:,:) + vnx(2,:,:))
      fy_ux(10,9,:,:)=-r_fac*mun*vn_un(2,:,:)
      fy_ux(10,11,:,:)=-r_fac*mun*nn_inv
      fy_uy(10,9,:,:)=-r_fac*mun*vn_un(1,:,:)
      fy_uy(10,10,:,:)=-r_fac*mun*nn_inv
      
      s_u(10,1,:,:)=r_fac*(-g_i_rho*vn(1,:,:)
     $     + g_r*vi_un(1,:,:) + g_r_rho*vi(1,:,:)
     $     + g_x*(one-mfp_fac)*vi_un(1,:,:) 
     $     + g_x_rho*((one-mfp_fac)*vi(1,:,:) - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_rho(1,:,:) - rinx_rho(1,:,:) 
     $     - r_in_rho(1,:,:))
      s_u(10,4,:,:)=r_fac*(n_inv*(g_r + g_x*(one-mfp_fac))
     $     + g_x_mi1*((one-mfp_fac)*vi(1,:,:) - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_mi1(1,:,:) - rinx_mi1(1,:,:)
     $     - r_in_u9)
      s_u(10,5,:,:)=r_fac*(g_x_mi2*((one-mfp_fac)*vi(1,:,:) - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_mi2(1,:,:) - rinx_mi2(1,:,:))
      s_u(10,6,:,:)=r_fac*(g_x_mi3*((one-mfp_fac)*vi(1,:,:) - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_mi3(1,:,:) - rinx_mi3(1,:,:))
      s_u(10,8,:,:)=r_fac*(-g_i_p*vn(1,:,:) + g_r_p*vi(1,:,:)
     $     + g_x_p*((one-mfp_fac)*vi(1,:,:) - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_p(1,:,:) - rinx_p(1,:,:) 
     $     - r_in_p(1,:,:))
      s_u(10,9,:,:)=r_fac*(-g_i*vn_un(1,:,:) - g_i_rhon*vn(1,:,:)
     $     - g_x*vn_un(1,:,:) 
     $     + g_x_rhon*((one-mfp_fac)*vi(1,:,:) - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_rhon(1,:,:) - rinx_rhon(1,:,:) 
     $     - r_in_rhon(1,:,:))
      s_u(10,10,:,:)=-r_fac*(nn_inv*(g_i + g_x)
     $     + g_x_mn1*(vn(1,:,:) - (one-mfp_fac)*vi(1,:,:))
     $     - (one-mfp_fac)*rnix_mn1(1,:,:) + rinx_mn1(1,:,:)
     $     + r_in_u1)
      s_u(10,11,:,:)=r_fac*(g_x_mn2*((one-mfp_fac)*vi(1,:,:) 
     $     - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_mn2(1,:,:) - rinx_mn2(1,:,:))
      s_u(10,12,:,:)=r_fac*(g_x_mn3*((one-mfp_fac)*vi(1,:,:) 
     $     - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_mn3(1,:,:) - rinx_mn3(1,:,:))
      s_u(10,13,:,:)=r_fac*(g_x_pn*((one-mfp_fac)*vi(1,:,:) - vn(1,:,:))
     $     + (one-mfp_fac)*rnix_pn(1,:,:) - rinx_pn(1,:,:) 
     $     - r_in_pn(1,:,:))
c-----------------------------------------------------------------------
c     y-component of neutral momentum equation.
c-----------------------------------------------------------------------
      fx_u(11,9,:,:)=r_fac*(u(11,:,:)*vn_un(1,:,:) 
     $     - mun*(vnx_un(2,:,:) + vny_un(1,:,:))
     $     - mun_un*(vnx(2,:,:) + vny(1,:,:)))
      fx_u(11,10,:,:)=r_fac*(vn(2,:,:) - mun*nny_inv)
      fx_u(11,11,:,:)=r_fac*(vn(1,:,:) - mun*nnx_inv)
      fx_u(11,13,:,:)=-r_fac*mun_pn*(vnx(2,:,:) + vny(1,:,:))
      fx_ux(11,9,:,:)=-r_fac*mun*vn_un(2,:,:)
      fx_ux(11,11,:,:)=-r_fac*mun*nn_inv
      fx_uy(11,9,:,:)=-r_fac*mun*vn_un(1,:,:)
      fx_uy(11,10,:,:)=-r_fac*mun*nn_inv
      
      fy_u(11,9,:,:)=r_fac*(u(11,:,:)*vn_un(2,:,:) 
     $     - two*(mun*vny_un(2,:,:) + mun_un*vny(2,:,:)))
      fy_u(11,11,:,:)=r_fac*two*(vn(2,:,:) - mun*nny_inv)
      fy_u(11,13,:,:)=-r_fac*two*mun_pn*vny(2,:,:)
      fy_uy(11,9,:,:)=-r_fac*mun*two*vn_un(2,:,:)
      fy_uy(11,11,:,:)=-r_fac*mun*two*nn_inv
      
      s_u(11,1,:,:)=r_fac*(-g_i_rho*vn(2,:,:)
     $     + g_r*vi_un(2,:,:) + g_r_rho*vi(2,:,:)
     $     + g_x*(one-mfp_fac)*vi_un(2,:,:) 
     $     + g_x_rho*((one-mfp_fac)*vi(2,:,:) - vn(2,:,:))
     $     + (one-mfp_fac)*rnix_rho(2,:,:) - rinx_rho(2,:,:) 
     $     - r_in_rho(2,:,:))
      s_u(11,4,:,:)=r_fac*(g_x_mi1*((one-mfp_fac)*vi(2,:,:)-vn(2,:,:))
     $     + (one-mfp_fac)*rnix_mi1(2,:,:) - rinx_mi1(2,:,:))
      s_u(11,5,:,:)=r_fac*(n_inv*(g_r + (one-mfp_fac)*g_x)
     $     + g_x_mi2*((one-mfp_fac)*vi(2,:,:) - vn(2,:,:))
     $     + (one-mfp_fac)*rnix_mi2(2,:,:) - rinx_mi2(2,:,:)
     $     - r_in_u9)
      s_u(11,6,:,:)=r_fac*(g_x_mi3*((one-mfp_fac)*vi(2,:,:)-vn(2,:,:))
     $     + (one-mfp_fac)*rnix_mi3(2,:,:) - rinx_mi3(2,:,:))
      s_u(11,8,:,:)=r_fac*(-g_i_p*vn(2,:,:) + g_r_p*vi(2,:,:)
     $     + g_x_p*((one-mfp_fac)*vi(2,:,:) - vn(2,:,:))
     $     + (one-mfp_fac)*rnix_p(2,:,:) - rinx_p(2,:,:) 
     $     - r_in_p(2,:,:))
      s_u(11,9,:,:)=cyl_fac*(u(12,:,:)*vn_un(3,:,:) 
     $     - r_faci*two*(mun*vn_un(2,:,:) + mun_un*vn(2,:,:)))
     $     + r_fac*(-g_i*vn_un(2,:,:) - g_i_rhon*vn(2,:,:)
     $     - g_x*vn_un(2,:,:) 
     $     + g_x_rhon*((one-mfp_fac)*vi(2,:,:) - vn(2,:,:))
     $     + (one-mfp_fac)*rnix_rhon(2,:,:) - rinx_rhon(2,:,:) 
     $     - r_in_rhon(2,:,:) - gravity)
      s_u(11,10,:,:)=r_fac*(
     $     g_x_mn1*((one-mfp_fac)*vi(2,:,:) - vn(2,:,:))
     $     + (one-mfp_fac)*rnix_mn1(2,:,:) - rinx_mn1(2,:,:))
      s_u(11,11,:,:)=-r_fac*(nn_inv*(g_i + g_x)
     $     + g_x_mn2*(vn(2,:,:) - (one-mfp_fac)*vi(2,:,:))
     $     - (one-mfp_fac)*rnix_mn2(2,:,:) + rinx_mn2(2,:,:)
     $     + r_in_u1) - cyl_fac*two*mun*nn_inv*r_faci
      s_u(11,12,:,:)=cyl_fac*two*vn(3,:,:) 
     $     + r_fac*(g_x_mn3*((one-mfp_fac)*vi(2,:,:) - vn(2,:,:))
     $     + (one-mfp_fac)*rnix_mn3(2,:,:) - rinx_mn3(2,:,:))
      s_u(11,13,:,:)=r_fac*(g_x_pn*((one-mfp_fac)*vi(2,:,:) - vn(2,:,:))
     $     + (one-mfp_fac)*rnix_pn(2,:,:) - rinx_pn(2,:,:) 
     $     - r_in_pn(2,:,:)) - cyl_fac*two*r_faci*mun_pn*vn(2,:,:)
      
      s_uy(11,13,:,:)=-r_fac
c-----------------------------------------------------------------------
c     z-component of neutral momentum equation.
c-----------------------------------------------------------------------
      fx_u(12,9,:,:)=r_fac*(u(12,:,:)*vn_un(1,:,:) 
     $     - mun*vnx_un(3,:,:) - mun_un*vnx(3,:,:))
      fx_u(12,10,:,:)=r_fac*u(12,:,:)*nn_inv      
      fx_u(12,12,:,:)=r_fac*(vn(1,:,:) - mun*nnx_inv)
      fx_u(12,13,:,:)=-r_fac*mun_pn*vnx(3,:,:)
      fx_ux(12,9,:,:)=-r_fac*mun*vn_un(3,:,:)
      fx_ux(12,12,:,:)=-r_fac*mun*nn_inv
      
      fy_u(12,9,:,:)=r_fac*(u(12,:,:)*vn_un(2,:,:) 
     $     - mun*vny_un(3,:,:) - mun_un*vny(3,:,:))
      fy_u(12,11,:,:)=r_fac*u(12,:,:)*nn_inv
      fy_u(12,12,:,:)=r_fac*(vn(2,:,:) - mun*nny_inv)
      fy_u(12,13,:,:)=-r_fac*mun_pn*vny(3,:,:)
      fy_uy(12,9,:,:)=-r_fac*mun*vn_un(3,:,:)
      fy_uy(12,12,:,:)=-r_fac*mun*nn_inv
      
      s_u(12,1,:,:)=r_fac*(g_r*vi_un(3,:,:) 
     $     + g_r_rho*vi(3,:,:) - g_i_rho*vn(3,:,:) 
     $     + g_x*(one-mfp_fac)*vi_un(3,:,:) 
     $     + g_x_rho*((one-mfp_fac)*vi(3,:,:) - vn(3,:,:))
     $     + (one-mfp_fac)*rnix_rho(3,:,:) - rinx_rho(3,:,:) 
     $     - r_in_rho(3,:,:))
      s_u(12,4,:,:)=r_fac*(g_x_mi1*((one-mfp_fac)*vi(3,:,:) 
     $     - vn(3,:,:))
     $     + (one-mfp_fac)*rnix_mi1(3,:,:) - rinx_mi1(3,:,:))
      s_u(12,5,:,:)=r_fac*(g_x_mi2*((one-mfp_fac)*vi(3,:,:) 
     $     - vn(3,:,:))
     $     + (one-mfp_fac)*rnix_mi2(3,:,:) - rinx_mi2(3,:,:))
      s_u(12,6,:,:)=r_fac*(n_inv*(g_r + (one-mfp_fac)*g_x)
     $     + g_x_mi3*((one-mfp_fac)*vi(3,:,:) - vn(3,:,:))
     $     + (one-mfp_fac)*rnix_mi3(3,:,:) - rinx_mi3(3,:,:) 
     $     - r_in_u9)
      s_u(12,8,:,:)=r_fac*(g_r_p*vi(3,:,:) - g_i_p*vn(3,:,:)
     $     + g_x_p*((one-mfp_fac)*vi(3,:,:) - vn(3,:,:))
     $     + (one-mfp_fac)*rnix_p(3,:,:) - rinx_p(3,:,:) 
     $     - r_in_p(3,:,:))
      s_u(12,9,:,:)=cyl_fac*(vn(2,:,:)*vn(3,:,:) 
     $     - r_faci*(mun*vn_un(3,:,:) + mun_un*vn(3,:,:)))
     $     - r_fac*(g_i*vn_un(3,:,:) + g_i_rhon*vn(3,:,:)
     $     + g_x*vn_un(3,:,:) + g_x_rhon*(vn(3,:,:) 
     $     - (one-mfp_fac)*vi(3,:,:))
     $     - (one-mfp_fac)*rnix_rhon(3,:,:) + rinx_rhon(3,:,:) 
     $     + r_in_rhon(3,:,:))
      s_u(12,10,:,:)=r_fac*(g_x_mn1*((one-mfp_fac)*vi(3,:,:) 
     $     - vn(3,:,:))
     $     + (one-mfp_fac)*rnix_mn1(3,:,:) - rinx_mn1(3,:,:))
      s_u(12,11,:,:)=-cyl_fac*vn(3,:,:)
     $     + r_fac*(g_x_mn2*((one-mfp_fac)*vi(3,:,:) - vn(3,:,:))
     $     + (one-mfp_fac)*rnix_mn2(3,:,:) - rinx_mn2(3,:,:))
      s_u(12,12,:,:)=-cyl_fac*(vn(2,:,:) + r_faci*mun*nn_inv)
     $     + r_fac*(g_x_mn3*((one-mfp_fac)*vi(3,:,:) - vn(3,:,:))
     $     - nn_inv*(g_i + g_x)
     $     + (one-mfp_fac)*rnix_mn3(3,:,:) - rinx_mn3(3,:,:) 
     $     - r_in_u1)
      s_u(12,13,:,:)=r_fac*(g_x_pn*((one-mfp_fac)*vi(3,:,:) 
     $     - vn(3,:,:))
     $     + (one-mfp_fac)*rnix_pn(3,:,:) - rinx_pn(3,:,:) 
     $     - r_in_pn(3,:,:))
     $     - cyl_fac*r_faci*mun_pn*vn(3,:,:)
c-----------------------------------------------------------------------
c     neutral pressure.
c-----------------------------------------------------------------------
      fx_u(13,9,:,:)=r_fac*(gamma_fac*u(13,:,:)*vn_un(1,:,:)
     $     - kapn*(Tnx_un - half*Tnx*nn_inv))
      fx_u(13,10,:,:)=r_fac*gamma_fac*u(13,:,:)*nn_inv
      fx_u(13,13,:,:)=r_fac*(gamma_fac*vn(1,:,:)
     $     - kapn*(nnx_inv + half*Tnx/u(13,:,:)))
      fx_ux(13,9,:,:)=-r_fac*kapn*Tn_un
      fx_ux(13,13,:,:)=-r_fac*kapn*nn_inv

      fy_u(13,9,:,:)=r_fac*(gamma_fac*u(13,:,:)*vn_un(2,:,:)
     $     - kapn*(Tny_un - half*Tny*nn_inv))
      fy_u(13,11,:,:)=r_fac*gamma_fac*u(13,:,:)*nn_inv
      fy_u(13,13,:,:)=r_fac*(gamma_fac*vn(2,:,:)
     $     - kapn*(nny_inv + half*Tny/u(13,:,:)))
      fy_uy(13,9,:,:)=-r_fac*kapn*Tn_un
      fy_uy(13,13,:,:)=-r_fac*kapn*nn_inv

      s_u(13,1,:,:)=r_fac*(g_r*kediff_rho + g_r_rho*kediff
     $     + (one-mfp_fac)*(g_x*kediff_rho + g_x_rho*kediff)
     $     + 1.5_r8*(-g_i_rho*u(13,:,:)*nn_inv
     $     + 1.5_r8*g_r*ti_frac*u(8,:,:)*n_inv**2)
     $     + (one-mfp_fac)*(rnix_rho(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_rho(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_rho(3,:,:)*(vi(3,:,:) - vn(3,:,:))
     $     - n_inv*(r_nix(1,:,:)*vi(1,:,:) + r_nix(2,:,:)*vi(2,:,:)
     $     + r_nix(3,:,:)*vi(3,:,:)) + qnix_rho) - qinx_rho + q_in_8_1 
     $     - qt_in_8_1)
      s_u(13,4,:,:)=r_fac*((g_r + (one-mfp_fac)*g_x)*kediff_mi(1,:,:) 
     $     + g_x_mi1*(one-mfp_fac)*kediff
     $     + (one-mfp_fac)*(r_nix(1,:,:)*n_inv 
     $     + rnix_mi1(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mi1(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mi1(3,:,:)*(vi(3,:,:) - vn(3,:,:)) + qnix_mi(1,:,:))
     $     - qinx_mi(1,:,:) + q_in_vi(1,:,:))
      s_u(13,5,:,:)=r_fac*((g_r + (one-mfp_fac)*g_x)*kediff_mi(2,:,:) 
     $     + g_x_mi2*(one-mfp_fac)*kediff
     $     + (one-mfp_fac)*(r_nix(2,:,:)*n_inv 
     $     + rnix_mi2(1,:,:)*(vi(1,:,:) - vn(1,:,:)) 
     $     + rnix_mi2(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mi2(3,:,:)*(vi(3,:,:) - vn(3,:,:)) + qnix_mi(2,:,:))
     $     - qinx_mi(2,:,:) + q_in_vi(2,:,:))
      s_u(13,6,:,:)=r_fac*((g_r + (one-mfp_fac)*g_x)*kediff_mi(3,:,:) 
     $     + g_x_mi3*(one-mfp_fac)*kediff
     $     + (one-mfp_fac)*(r_nix(3,:,:)*n_inv 
     $     + rnix_mi3(1,:,:)*(vi(1,:,:) - vn(1,:,:)) 
     $     + rnix_mi3(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mi3(3,:,:)*(vi(3,:,:) - vn(3,:,:)) + qnix_mi(3,:,:))
     $     - qinx_mi(3,:,:) + q_in_vi(3,:,:))
      s_u(13,8,:,:)=r_fac*((g_r_p + g_x_p*(one-mfp_fac))*kediff
     $     + 1.5_r8*(-g_i_p*u(13,:,:)*nn_inv 
     $     + ti_frac*(g_r_p*u(8,:,:)*n_inv + g_r*n_inv)) 
     $     + (one-mfp_fac)*(rnix_p(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_p(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_p(3,:,:)*(vi(3,:,:) - vn(3,:,:)) + qnix_p) - qinx_p + 
     $       q_in_8_8 - qt_in_8_8)
      s_u(13,9,:,:)=r_fac*(vn_un(1,:,:)*ux(13,:,:)
     $     + vn_un(2,:,:)*uy(13,:,:)
     $     + two*mun*(two*vnx(1,:,:)*vnx_un(1,:,:) 
     $     + two*vny(2,:,:)*vny_un(2,:,:) 
     $     + vny(1,:,:)*vny_un(1,:,:) + vnx(2,:,:)*vnx_un(2,:,:)
     $     + vny(1,:,:)*vnx_un(2,:,:) + vnx(2,:,:)*vny_un(1,:,:)
     $     + vnx(3,:,:)*vnx_un(3,:,:) + vny(3,:,:)*vny_un(3,:,:))
     $     + cyl_fac*two*r_faci**2*mun*(two*vn(2,:,:)*vn_un(2,:,:)
     $     + vn(3,:,:)*vn_un(3,:,:))
     $     + mun_un*(two*vnx(1,:,:)**2 + two*vny(2,:,:)**2 
     $     + vny(1,:,:)**2 + vnx(2,:,:)**2 + two*vny(1,:,:)*vnx(2,:,:)
     $     + vnx(3,:,:)**2 + vny(3,:,:)**2)
     $     + cyl_fac*r_faci**2*mun_un*(two*vn(2,:,:)**2 + vn(3,:,:)**2)
     $     + (g_r + g_x*(one-mfp_fac))*kediff_rhon
     $     + g_x_rhon*(one-mfp_fac)*kediff
     $     + 1.5_r8*nn_inv*(-g_i_rhon*u(13,:,:) 
     $     + g_i*u(13,:,:)*nn_inv)
     $     + (one-mfp_fac)*(rnix_rhon(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_rhon(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_rhon(3,:,:)*(vi(3,:,:) - vn(3,:,:))
     $     + nn_inv*(r_nix(1,:,:)*vn(1,:,:) + r_nix(2,:,:)*vn(2,:,:)
     $     + r_nix(3,:,:)*vn(3,:,:)) + qnix_rhon) - qinx_rhon + q_in_8_9 
     $     - qt_in_8_9)
      s_u(13,10,:,:)=r_fac*(nn_inv*ux(13,:,:) + two*mun*
     $     (two*vnx(1,:,:)*nnx_inv + (vny(1,:,:) + vnx(2,:,:))*nny_inv)
     $     + (g_r + g_x*(one-mfp_fac))*kediff_mn(1,:,:)
     $     + g_x_mn1*(one-mfp_fac)*kediff
     $     + (one-mfp_fac)*(rnix_mn1(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mn1(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mn1(3,:,:)*(vi(3,:,:) - vn(3,:,:))
     $     - r_nix(1,:,:)*nn_inv + qnix_mn(1,:,:)) 
     $     - qinx_mn(1,:,:) + q_in_vn(1,:,:))
      s_u(13,11,:,:)=r_fac*(uy(13,:,:)*nn_inv + two*mun*
     $     (two*vny(2,:,:)*nny_inv + (vnx(2,:,:) + vny(1,:,:))*nnx_inv)
     $     + cyl_fac*4._r8*mun*r_faci**2*vn(2,:,:)*nn_inv
     $     + (g_r + g_x*(one-mfp_fac))*kediff_mn(2,:,:)
     $     + g_x_mn2*(one-mfp_fac)*kediff
     $     + (one-mfp_fac)*(rnix_mn2(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mn2(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mn2(3,:,:)*(vi(3,:,:) - vn(3,:,:))
     $     - r_nix(2,:,:)*nn_inv + qnix_mn(2,:,:))
     $     - qinx_mn(2,:,:) + q_in_vn(2,:,:))
      s_u(13,12,:,:)=r_fac*(two*(mun*(vnx(3,:,:)*nnx_inv 
     $     + vny(3,:,:)*nny_inv)
     $     + cyl_fac*r_faci**2*mun*vn(3,:,:)*nn_inv)
     $     + (g_r + g_x*(one-mfp_fac))*kediff_mn(3,:,:)
     $     + g_x_mn3*(one-mfp_fac)*kediff
     $     + (one-mfp_fac)*(rnix_mn3(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mn3(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mn3(3,:,:)*(vi(3,:,:) - vn(3,:,:))
     $     - r_nix(3,:,:)*nn_inv + qnix_mn(3,:,:)) 
     $     - qinx_mn(3,:,:) + q_in_vn(3,:,:))
      s_u(13,13,:,:)=r_fac*(g_x_pn*(one-mfp_fac)*kediff
     $     - 1.5_r8*g_i*nn_inv + (one-mfp_fac)*(qnix_pn
     $     + rnix_pn(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_pn(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_pn(3,:,:)*(vi(3,:,:) - vn(3,:,:))) - qinx_pn 
     $     + q_in_8_13 - qt_in_8_13
     $     + mun_pn*(two*vnx(1,:,:)**2 + two*vny(2,:,:)**2 
     $     + vny(1,:,:)**2 + vnx(2,:,:)**2 + two*vny(1,:,:)*vnx(2,:,:)
     $     + vnx(3,:,:)**2 + vny(3,:,:)**2)
     $     + cyl_fac*r_faci**2*mun_pn*(two*vn(2,:,:)**2 + vn(3,:,:)**2))

      s_ux(13,9,:,:)=r_fac*two*mun*(two*vnx(1,:,:)*vn_un(1,:,:)  
     $     + vn_un(2,:,:)*(vnx(2,:,:) + vny(1,:,:))
     $     + vnx(3,:,:)*vn_un(3,:,:))
      s_ux(13,10,:,:)=r_fac*mun*4._r8*vnx(1,:,:)*nn_inv
      s_ux(13,11,:,:)=r_fac*mun*two*nn_inv*(vnx(2,:,:) + vny(1,:,:))
      s_ux(13,12,:,:)=r_fac*two*mun*vnx(3,:,:)*nn_inv
      s_ux(13,13,:,:)=r_fac*vn(1,:,:)

      s_uy(13,9,:,:)=r_fac*two*mun*(two*vny(2,:,:)*vn_un(2,:,:)  
     $     + vn_un(1,:,:)*(vny(1,:,:) + vnx(2,:,:))
     $     + vny(3,:,:)*vn_un(3,:,:))
      s_uy(13,10,:,:)=r_fac*mun*two*nn_inv*(vny(1,:,:) + vnx(2,:,:))
      s_uy(13,11,:,:)=r_fac*mun*4._r8*vny(2,:,:)*nn_inv
      s_uy(13,12,:,:)=r_fac*two*mun*vny(3,:,:)*nn_inv
      s_uy(13,13,:,:)=r_fac*vn(2,:,:)
c-----------------------------------------------------------------------
c     del^2(out-of-plane B) auxiliary equation.
c-----------------------------------------------------------------------
      IF(nu /= 0)THEN
         fx_ux(14,3,:,:)=one
         fy_u(14,3,:,:)=cyl_fac*r_faci
         fy_uy(14,3,:,:)=one
         s_u(14,14,:,:)=one
      ENDIF
c-----------------------------------------------------------------------
c     zero out-of-plane equations when init_type="RTearth" or "RTsun"
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("RTearth","RTsun")
         fx_u(2,:,:,:)=0
         fy_u(2,:,:,:)=0
         s_u(2,:,:,:)=0
         fx_u(6,:,:,:)=0
         fy_u(6,:,:,:)=0
         s_u(6,:,:,:)=0
         fx_u(12,:,:,:)=0
         fy_u(12,:,:,:)=0
         s_u(12,:,:,:)=0
         fx_u(7,:,:,:)=0
         fy_u(7,:,:,:)=0
         s_u(7,:,:,:)=0
         s_u(7,7,:,:)=one

         fx_ux(2,:,:,:)=0
         fy_ux(2,:,:,:)=0
         s_ux(2,:,:,:)=0
         fx_ux(6,:,:,:)=0
         fy_ux(6,:,:,:)=0
         s_ux(6,:,:,:)=0
         fx_ux(12,:,:,:)=0
         fy_ux(12,:,:,:)=0
         s_ux(12,:,:,:)=0
         fx_ux(7,:,:,:)=0
         fy_ux(7,:,:,:)=0
         s_ux(7,:,:,:)=0

         fx_uy(2,:,:,:)=0
         fy_uy(2,:,:,:)=0
         s_uy(2,:,:,:)=0
         fx_uy(6,:,:,:)=0
         fy_uy(6,:,:,:)=0
         s_uy(6,:,:,:)=0
         fx_uy(12,:,:,:)=0
         fy_uy(12,:,:,:)=0
         s_uy(12,:,:,:)=0
         fx_uy(7,:,:,:)=0
         fy_uy(7,:,:,:)=0
         s_uy(7,:,:,:)=0
      END SELECT
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
      USE pn_ext_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y

      REAL(r8), DIMENSION(SIZE(y,1),SIZE(y,2)) :: r_fac
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
c        cylindrical to cartesian relationships:
c        1: z --> x
c        2: r --> y
c        3: phi --> z
c-----------------------------------------------------------------------
      r_fac=y**cyl_fac

      mass(1,1,:,:)=r_fac
      mass(2,2,:,:)=r_fac
      mass(3,3,:,:)=one
      mass(4,4,:,:)=r_fac
      mass(5,5,:,:)=r_fac
      mass(6,6,:,:)=r_fac
      mass(8,8,:,:)=r_fac/(gamma-one)
      mass(9,1,:,:)=r_fac
      mass(9,9,:,:)=r_fac
      mass(10,10,:,:)=r_fac
      mass(11,11,:,:)=r_fac
      mass(12,12,:,:)=r_fac
      mass(13,13,:,:)=r_fac/(gamma-one)
      IF(.NOT. ion_momentum)THEN
         mass(4,10,:,:)=r_fac
         mass(5,11,:,:)=r_fac
         mass(6,12,:,:)=r_fac
      ENDIF
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
      SUBROUTINE physics_grid(x,y,ksi,phi)
      USE pn_ext_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,phi
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("GEM")
         ksi=lx*x
         phi=ly*(gr_curve*y**4 + y)/(gr_curve+one)
      CASE("MRX")
         ksi=lx*(gr_curve*x**2 + x)/(gr_curve+one)
         phi=ly*(gr_curve*y**4 + y)/(gr_curve+one)
      CASE("MRXpush")
         ksi= x - half
         phi=(gr_curve*(y-half)**3 + (y-half))
     $        /(half**2*gr_curve + one)
      CASE("forcefree")
         ksi=two*lx*(x - half)
         phi=two*ly*(gr_curve*(y-half)**3 + (y-half))
     $        /(half**2*gr_curve + one)
      CASE("RTsun","RTearth")
         ksi=two*lx*(x - half)
         phi=two*ly*(y - half)
      END SELECT
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
      USE pn_ext_mod
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
      USE pn_ext_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     deallocate appropriate arrays.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      END SELECT
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
