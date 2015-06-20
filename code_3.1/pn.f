c-----------------------------------------------------------------------
c     file pn.f.
c     contains specifications for a reacting plasma-neutral model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. pn_mod.
c     1. pn_equil.
c     2. pn_grq.
c     3. pn_grq_u.
c     4. pn_bc_jac.
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
c     o. physics_dealloc.
c     p. physics_main.
c-----------------------------------------------------------------------
c     subprogram 0. physics_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE pn_mod
      USE extra_mod
      USE transport_mod
      IMPLICIT NONE

      LOGICAL :: cylinder=.FALSE.,rad=.FALSE.,ifcoils=.FALSE.
      CHARACTER(16) :: init_type=".",equilfile=".",atom="deuterium",
     $     coilfile=".",interp="bicube",eta_case=".",kappa_case=".",
     $     targ_type="block",visc_case="uniform",viscn_case="uniform",
     $     kappan_case="uniform"
      INTEGER :: cyl_fac=0
      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,qe=1.602e-19,lnlam=10.,
     $     me=9.109e-31,mp=1.673e-27,ep0=8.854e-12,mu0=4.e-7*pi,
     $     rho_rad=2.e-1_r8
      REAL(r8) :: eta=0.,r_eta=1.e10,etavac=1.,v_chod_norm=1.,mi=0.,
     $     etac_norm=1.,etas_norm=1.,mu=0.,kappa_par=0.,kappa_perp=0.,
     $     ddiff=0.,beta0=0.,rhomin=0.,pmin=0.,n0=1.e20,mu_min=0.,
     $     b0=1.,L0=1.,bx=0.,lx=1.,ly=1.,xmin=0.,mu_norm=1.,
     $     epsilon=0.,gamma_fac=1.,ke_norm=1.,xe_norm=1.,ki_norm=1.,
     $     xi_norm=1.,initrho=1.,initp=1.,ion_fac=1.,initrhon=1.,
     $     ion_norm=1.,Top_norm=1.,cx_norm=1.,recomb_norm=1.,hs_diam=0.,
     $     recomb_fac=1.,psource=0.,phi_ion=0.,phi_eff=0.,mun=0.,
     $     initpn=1.,cx_fac=1.,rhonmin=0.,te_frac=.5,cx_c1=0.,cx_c2=0.,
     $     cc=0.01,v_civ=0.,civ_fac=0.,ti_frac=.5,initv=0.,cx_c3,
     $     kappa_n=0.,a_vor=0.,p_vor=0.,x_vor=0.,k_vor=0.,mun_norm=1.,
     $     kn_norm,initTn=1.,kappa_min=0.,gr_curve=0.,r_visc=1.e10,
     $     Tfac=1._r8,eta_hr=0.,kappa_max=1.e10,mu_sv=0.,mun_sv=0.,
     $     rhofac=1._r8,etafac=1._r8

      TYPE(coil_type) :: coils
      TYPE(bicube_type) :: equil_bc
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: equil,equilxy

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. pn_equil.
c     sets equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pn_equil(x,y,u)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(INOUT) :: u

      REAL(r8) :: Tmin
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Ttot
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: eq
c-----------------------------------------------------------------------
c     set equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     eq(1) is flux, eq(2) is pressure, and eq(3) is jphi
c-----------------------------------------------------------------------
      CALL extra_interp(interp,x,y,equil_bc,equil,equilxy,eq)
c-----------------------------------------------------------------------
c     apply cosine smoothing to floor pressure (pmin).
c-----------------------------------------------------------------------
      WHERE(eq(2,:,:) < pmin)
         u(3,:,:)=pmin
      ELSEWHERE(eq(2,:,:) >= pmin .AND. eq(2,:,:) < two*pmin)
         u(3,:,:)=pmin 
     $        + half*pmin*(one - COS(pi*(eq(2,:,:)-pmin)/pmin))
      ELSEWHERE
         u(3,:,:)=eq(2,:,:)
      ENDWHERE

c-----------------------------------------------------------------------
c     assign magnetic vector potential and current.
c-----------------------------------------------------------------------
      WHERE(y == 0)
         u(2,:,:) = zero
      ELSEWHERE
         u(2,:,:) = -eq(1,:,:)/y
      ENDWHERE
      u(6,:,:) = eq(3,:,:)

c-----------------------------------------------------------------------
c     assign density based on sqaure root of pressure.
c-----------------------------------------------------------------------
      u(1,:,:) = SQRT(u(3,:,:))*rhofac

c-----------------------------------------------------------------------
c     modify pressure/temperature/density and maintain equilibrium by
c     appropriately changing magnetic components.
c-----------------------------------------------------------------------
      u(2,:,:) = u(2,:,:)*SQRT(Tfac*rhofac)
      u(3,:,:) = u(3,:,:)*Tfac*rhofac
      u(6,:,:) = u(6,:,:)*SQRT(Tfac*rhofac)

c-----------------------------------------------------------------------
c     where T doesn't approach desired boundary temperature, smoothly
c     correct it by modifying density.
c-----------------------------------------------------------------------
      Tmin = pmin/rhomin
      Ttot = u(3,:,:)/u(1,:,:)
      WHERE(y > .8) Ttot = 
     $     Ttot - (Ttot - Tmin)*half*(one - COS(pi*(y - 0.8)/0.2))
      u(1,:,:) = u(3,:,:)/Ttot

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pn_equil
c-----------------------------------------------------------------------
c     subprogram 2. pn_grq.
c     set reaction source rates (g), frictional momentum contributions
c     (r) and thermal energy transfer (q).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pn_grq(p0,rho0,pn0,rhon0,vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: p0,pn0,rho0,rhon0
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: vi,vn
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: g_i,g_r,g_x,q_inx,q_nix
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: r_inx,r_nix

      INTEGER :: i
      REAL(r8), PARAMETER :: eps_loc=1.e-12_r8
      REAL(r8), DIMENSION(SIZE(p0,1),SIZE(p0,2)) :: p,rho,pn,rhon,Top,
     $     poT,g_civ,dv,v_cx,vrep,vrepa,vrepb,s_cx
      REAL(r8), DIMENSION(2,SIZE(p0,1),SIZE(p0,2)) :: vin

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
      dv = 0
      DO i=1,SIZE(vi,1)
         dv = dv + (vi(i,:,:) - vn(i,:,:))**2
         vin(i,:,:) = vi(i,:,:) - vn(i,:,:)
      ENDDO
      WHERE(dv /= 0) dv = SQRT(dv)
      WHERE(ABS(vin) < 1.e-10_r8) vin = 1.e-10_r8

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
         ENDWHERE
      ENDIF

      toP = (Top_norm*te_frac*p)/rho
      poT = rho/(Top_norm*te_frac*p)
      g_i = g_civ + ion_fac*ion_norm*rho*rhon*EXP(-poT)*
     $     poT**k_vor*(one + p_vor*SQRT(poT))/(x_vor + poT)

      WHERE(rhon < rhonmin)g_i=0
      g_r=recomb_fac*recomb_norm*rho**2*SQRT(poT)

      v_cx = SQRT(MAX(8.*(ti_frac*p/rho + pn/rhon)/pi + dv**2,eps_loc))
      s_cx = MAX(-LOG(v_cx) + cx_c3,eps_loc)
      g_x = cx_fac*cx_norm*s_cx*v_cx*rhon*rho

c-----------------------------------------------------------------------
c     set frictional transfer for CX.
c-----------------------------------------------------------------------
      vrepa = two*pn/rhon
      vrepb = SQRT(MAX(4._r8*(8._r8*ti_frac*p/(rho*pi) + dv**2) 
     $     + 9._r8*pi*pn/(two*rhon),eps_loc))
      vrep = vrepa/vrepb
      DO i=1,2
         r_inx(i,:,:) = -cx_fac*cx_norm*s_cx*rho*rhon*vin(i,:,:)*vrep
      ENDDO

      vrepa = two*ti_frac*p/rho
      vrepb = SQRT(MAX(4._r8*(8._r8*pn/(rhon*pi) + dv**2) 
     $     + 9._r8*pi*ti_frac*p/(two*rho),eps_loc))
      vrep = vrepa/vrepb
      DO i=1,2
         r_nix(i,:,:) = cx_fac*cx_norm*s_cx*rho*rhon*vin(i,:,:)*vrep
      ENDDO

c-----------------------------------------------------------------------
c     set thermal energy transfer for CX.
c-----------------------------------------------------------------------
      vrepa = 3._r8*pn/(two*rhon)
      vrepb = SQRT(MAX(8._r8*ti_frac*p/(pi*rho) 
     $     + 128._r8*pn/(9._r8*pi*rhon) + dv**2,eps_loc))
      vrep = vrepa*vrepb
      q_inx = cx_fac*cx_norm*s_cx*rho*rhon*vrep

      vrepa = 3._r8*ti_frac*p/(two*rho)
      vrepb = SQRT(MAX(8._r8*pn/(pi*rhon)
     $     + 128._r8*ti_frac*p/(9._r8*pi*rho) + dv**2,eps_loc))
      vrep = vrepa*vrepb
      q_nix = cx_fac*cx_norm*s_cx*rho*rhon*vrep

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pn_grq
c-----------------------------------------------------------------------
c     subprogram 3. pn_grq_u.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pn_grq_u(p0,rho0,pn0,rhon0,vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix,
     $     g_i_rho,g_i_rhon,g_i_p,g_r_rho,g_r_p,
     $     g_x_rho,g_x_rhon,g_x_p,g_x_pn,g_x_mi1,g_x_mi2,g_x_mn1,
     $     g_x_mn2,
     $     rinx_rho,rinx_p,rinx_rhon,rinx_pn,rinx_mi1,
     $     rinx_mi2,rinx_mn1,rinx_mn2,
     $     rnix_rho,rnix_p,rnix_rhon,rnix_pn,rnix_mi1,
     $     rnix_mi2,rnix_mn1,rnix_mn2,
     $     qinx_rho,qinx_p,qinx_rhon,qinx_pn,qinx_mi1,
     $     qinx_mi2,qinx_mn1,qinx_mn2,
     $     qnix_rho,qnix_p,qnix_rhon,qnix_pn,qnix_mi1,
     $     qnix_mi2,qnix_mn1,qnix_mn2)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: p0,rho0,pn0,rhon0
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: vi,vn
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: g_i,g_r,g_x,q_inx,q_nix,
     $     g_i_rho,g_i_rhon,g_i_p,g_r_rho,g_r_p,
     $     g_x_rho,g_x_rhon,g_x_p,g_x_pn,g_x_mi1,g_x_mi2,g_x_mn1,
     $     g_x_mn2,
     $     qinx_rho,qinx_p,qinx_rhon,qinx_pn,qinx_mi1,
     $     qinx_mi2,qinx_mn1,qinx_mn2,
     $     qnix_rho,qnix_p,qnix_rhon,qnix_pn,qnix_mi1,
     $     qnix_mi2,qnix_mn1,qnix_mn2
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: r_inx,r_nix,
     $     rinx_rho,rinx_p,rinx_rhon,rinx_pn,rinx_mi1,
     $     rinx_mi2,rinx_mn1,rinx_mn2,
     $     rnix_rho,rnix_p,rnix_rhon,rnix_pn,rnix_mi1,
     $     rnix_mi2,rnix_mn1,rnix_mn2

      INTEGER :: i
      REAL(r8), PARAMETER :: eps_loc=1.e-12_r8
      REAL(r8), DIMENSION(SIZE(p0,1),SIZE(p0,2)) :: p,rho,pn,rhon,Top,
     $     poT,g_civ,dv,vcx,s_cx,g_x_vi1,g_x_vi2,g_x_vn1,g_x_vn2,v_cx,
     $     vrep,vrepa,vrepb,vrep_vi1,vrep_vi2,vrep_vn1,vrep_vn2,
     $     vrep_rho,vrep_p,vrep_rhon,vrep_pn,scx_p,scx_rho,scx_pn,
     $     scx_rhon,scx_vi1,scx_vi2,scx_vn1,scx_vn2,g_poT
      REAL(r8), DIMENSION(2,SIZE(p0,1),SIZE(p0,2)) :: vin

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
      dv = 0
      DO i=1,SIZE(vi,1)
         dv = dv + (vi(i,:,:) - vn(i,:,:))**2
         vin(i,:,:) = vi(i,:,:) - vn(i,:,:)
      ENDDO
      WHERE(dv /= 0) dv = SQRT(dv)
      WHERE(ABS(vin) < 1.e-10_r8) vin = 1.e-10_r8

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
         ENDWHERE
      ENDIF

      Top = (Top_norm*te_frac*p)/rho
      poT = rho/(Top_norm*te_frac*p)

      CALL pn_grq(p,rho,pn,rhon,vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix)
      
      g_poT = (g_i-g_civ)*(-one/(poT + x_vor) + k_vor/poT - one
     $     + p_vor*half/(SQRT(poT) + p_vor*poT))
      g_i_rho = (g_i-g_civ)/rho + g_poT/(Top_norm*te_frac*p) + g_civ/rho
      g_i_rhon = (g_i-g_civ)/rhon + g_civ/rhon
      g_i_p = -g_poT*rho/(Top_norm*te_frac*p**2)
      WHERE(rhon < rhonmin)
         g_i_rho=0
         g_i_rhon=0
         g_i_p=0
      ENDWHERE

      g_poT = half*g_r/poT
      g_r_rho = two*g_r/rho + g_poT/(Top_norm*te_frac*p)
      g_r_p = -g_poT*rho/(Top_norm*te_frac*p**2)

      v_cx = SQRT(MAX(8.*(ti_frac*p/rho + pn/rhon)/pi + dv**2,eps_loc))
      s_cx = MAX(-LOG(v_cx) + cx_c3,eps_loc)

      scx_p = -4._r8*ti_frac/(pi*rho*v_cx**2)
      scx_rho = (4._r8*ti_frac*p/(pi*rho) 
     $     + vin(1,:,:)*vi(1,:,:) + vin(2,:,:)*vi(2,:,:))/(v_cx**2*rho)

      scx_pn = -4._r8/(pi*rhon*v_cx**2)
      scx_rhon = (4._r8*pn/(pi*rhon) 
     $     - vin(1,:,:)*vn(1,:,:) - vin(2,:,:)*vn(2,:,:))/(v_cx**2*rhon)

      scx_vi1 = -vin(1,:,:)/v_cx**2
      scx_vi2 = -vin(2,:,:)/v_cx**2
      scx_vn1 = vin(1,:,:)/v_cx**2
      scx_vn2 = vin(2,:,:)/v_cx**2

      g_x_vi1 = g_x*(vin(1,:,:)/v_cx**2 + scx_vi1/s_cx)
      g_x_vi2 = g_x*(vin(2,:,:)/v_cx**2 + scx_vi2/s_cx)
      g_x_vn1 = g_x*(scx_vn1/s_cx - vin(1,:,:)/v_cx**2)
      g_x_vn2 = g_x*(scx_vn2/s_cx - vin(2,:,:)/v_cx**2)

      g_x_p = g_x*scx_p*(one - s_cx)/s_cx
      g_x_rho = g_x/rho + g_x*scx_rho*(one - s_cx)/s_cx
      g_x_pn = g_x*scx_pn*(one - s_cx)/s_cx
      g_x_rhon = g_x/rhon + g_x*scx_rhon*(one - s_cx)/s_cx
      g_x_mi1 = g_x_vi1/rho
      g_x_mi2 = g_x_vi2/rho
      g_x_mn1 = g_x_vn1/rhon
      g_x_mn2 = g_x_vn2/rhon

c-----------------------------------------------------------------------
c     set derivatives for CX frictional transfer, R_inx.
c-----------------------------------------------------------------------
      vrepa = two*pn/rhon
      vrepb = SQRT(MAX(4._r8*(8._r8*ti_frac*p/(rho*pi) + dv**2) 
     $     + 9._r8*pi*pn/(two*rhon),eps_loc))
      vrep = vrepa/vrepb

      vrep_vi1 = -4._r8*vrepa*vin(1,:,:)/vrepb**3
      vrep_vi2 = -4._r8*vrepa*vin(2,:,:)/vrepb**3
      vrep_vn1 = -vrep_vi1
      vrep_vn2 = -vrep_vi2

      vrep_rho = 16._r8*vrepa*ti_frac*p/(pi*rho**2*vrepb**3)
     $     - vrep_vi1*vi(1,:,:)/rho - vrep_vi2*vi(2,:,:)/rho
      vrep_p = -16._r8*ti_frac*vrepa/(pi*vrepb**3*rho)

      vrep_rhon = 9._r8*pi*vrepa*pn/(4._r8*vrepb**3*rhon**2)
     $     - vrep_vn1*vn(1,:,:)/rhon - vrep_vn2*vn(2,:,:)/rhon
     $     - two*pn/(rhon**2*vrepb)
      vrep_pn = -9._r8*pi*vrepa/(4._r8*vrepb**3*rhon) 
     $     + two/(rhon*vrepb)

      DO i=1,2
         rinx_rho(i,:,:) = r_inx(i,:,:)*(one/rho + scx_rho/s_cx
     $        - vi(i,:,:)/(rho*vin(i,:,:)) + vrep_rho/vrep)
     $        
         rinx_p(i,:,:) = r_inx(i,:,:)*(vrep_p/vrep + scx_p/s_cx)
         rinx_rhon(i,:,:) = r_inx(i,:,:)*(one/rhon + scx_rhon/s_cx
     $        + vn(i,:,:)/(rhon*vin(i,:,:)) + vrep_rhon/vrep)
         rinx_pn(i,:,:) = r_inx(i,:,:)*(vrep_pn/vrep + scx_pn/s_cx)
      ENDDO
      
      rinx_mi1(1,:,:) = r_inx(1,:,:)*(one/(vin(1,:,:)*rho)
     $     + vrep_vi1/(vrep*rho) + scx_vi1/(s_cx*rho))
      rinx_mi1(2,:,:) = r_inx(2,:,:)*(vrep_vi1/(vrep*rho)
     $     + scx_vi1/(s_cx*rho))

      rinx_mi2(1,:,:) = r_inx(1,:,:)*(vrep_vi2/(vrep*rho)
     $     + scx_vi2/(s_cx*rho))
      rinx_mi2(2,:,:) = r_inx(2,:,:)*(one/(vin(2,:,:)*rho)
     $     + vrep_vi2/(vrep*rho) + scx_vi2/(s_cx*rho))

      rinx_mn1(1,:,:) = r_inx(1,:,:)*(-one/(vin(1,:,:)*rhon)
     $     + vrep_vn1/(vrep*rhon) + scx_vn1/(s_cx*rhon))
      rinx_mn1(2,:,:) = r_inx(2,:,:)*(vrep_vn1/(vrep*rhon)
     $     + scx_vn1/(s_cx*rhon))

      rinx_mn2(1,:,:) = r_inx(1,:,:)*(vrep_vn2/(vrep*rhon)
     $     + scx_vn2/(s_cx*rhon))
      rinx_mn2(2,:,:) = r_inx(2,:,:)*(-one/(vin(2,:,:)*rhon)
     $     + vrep_vn2/(vrep*rhon) + scx_vn2/(s_cx*rhon))

c-----------------------------------------------------------------------
c     set derivatives for CX frictional momentum transfer, R_nix.
c-----------------------------------------------------------------------
      vrepa = two*ti_frac*p/rho
      vrepb = SQRT(MAX(4._r8*(8._r8*pn/(rhon*pi) + dv**2) 
     $     + 9._r8*pi*ti_frac*p/(two*rho),eps_loc))
      vrep = vrepa/vrepb

      vrep_vi1 = -4._r8*vrepa*vin(1,:,:)/vrepb**3
      vrep_vi2 = -4._r8*vrepa*vin(2,:,:)/vrepb**3
      vrep_vn1 = -vrep_vi1
      vrep_vn2 = -vrep_vi2

      vrep_rho = 9._r8*pi*vrepa*ti_frac*p/(4._r8*vrepb**3*rho**2)
     $     - vrep_vi1*vi(1,:,:)/rho - vrep_vi2*vi(2,:,:)/rho
     $     - two*ti_frac*p/(rho**2*vrepb)
      vrep_p = -9._r8*pi*vrepa*ti_frac/(4._r8*vrepb**3*rho)
     $     + two*ti_frac/(rho*vrepb)

      vrep_rhon = 16._r8*vrepa*pn/(pi*rhon**2*vrepb**3)
     $     - vrep_vn1*vn(1,:,:)/rhon - vrep_vn2*vn(2,:,:)/rhon
      vrep_pn = -16._r8*vrepa/(pi*vrepb**3*rhon)

      DO i=1,2
         rnix_rho(i,:,:) = r_nix(i,:,:)*(one/rho + scx_rho/s_cx
     $        - vi(i,:,:)/(rho*vin(i,:,:)) + vrep_rho/vrep)
         rnix_p(i,:,:) = r_nix(i,:,:)*(vrep_p/vrep + scx_p/s_cx)
         rnix_rhon(i,:,:) = r_nix(i,:,:)*(one/rhon + scx_rhon/s_cx
     $        + vn(i,:,:)/(rhon*vin(i,:,:)) + vrep_rhon/vrep)
         rnix_pn(i,:,:) = r_nix(i,:,:)*(vrep_pn/vrep + scx_pn/s_cx)
      ENDDO
      
      rnix_mi1(1,:,:) = r_nix(1,:,:)*(one/(vin(1,:,:)*rho)
     $     + vrep_vi1/(vrep*rho) + scx_vi1/(s_cx*rho))
      rnix_mi1(2,:,:) = r_nix(2,:,:)*(vrep_vi1/(vrep*rho)
     $     + scx_vi1/(s_cx*rho))

      rnix_mi2(1,:,:) = r_nix(1,:,:)*(vrep_vi2/(vrep*rho)
     $     + scx_vi2/(s_cx*rho))
      rnix_mi2(2,:,:) = r_nix(2,:,:)*(one/(vin(2,:,:)*rho)
     $     + vrep_vi2/(vrep*rho) + scx_vi2/(s_cx*rho))

      rnix_mn1(1,:,:) = r_nix(1,:,:)*(-one/(vin(1,:,:)*rhon)
     $     + vrep_vn1/(vrep*rhon) + scx_vn1/(s_cx*rhon))
      rnix_mn1(2,:,:) = r_nix(2,:,:)*(vrep_vn1/(vrep*rhon)
     $     + scx_vn1/(s_cx*rhon))

      rnix_mn2(1,:,:) = r_nix(1,:,:)*(vrep_vn2/(vrep*rhon)
     $     + scx_vn2/(s_cx*rhon))
      rnix_mn2(2,:,:) = r_nix(2,:,:)*(-one/(vin(2,:,:)*rhon)
     $     + vrep_vn2/(vrep*rhon) + scx_vn2/(s_cx*rhon))

c-----------------------------------------------------------------------
c     set derivatives for thermal energy  transfer, Q_inx.
c-----------------------------------------------------------------------
      vrepa = 3._r8*pn/(two*rhon)
      vrepb = SQRT(MAX(8._r8*ti_frac*p/(pi*rho) 
     $     + 128._r8*pn/(9._r8*pi*rhon) + dv**2,eps_loc))
      vrep = vrepa*vrepb

      vrep_vi1 = vrepa*vin(1,:,:)/vrepb
      vrep_vi2 = vrepa*vin(2,:,:)/vrepb
      vrep_vn1 = -vrep_vi1
      vrep_vn2 = -vrep_vi2

      vrep_rho = -vrepa*4._r8*ti_frac*p/(pi*vrepb*rho**2)
     $     - vrep_vi1*vi(1,:,:)/rho - vrep_vi2*vi(2,:,:)/rho
      vrep_p = vrepa*4._r8*ti_frac/(pi*rho*vrepb)

      vrep_rhon = -vrepb*3._r8*pn/(two*rhon**2) 
     $     - vrepa*64._r8*pn/(9._r8*pi*rhon**2*vrepb)
     $     - vrep_vn1*vn(1,:,:)/rhon - vrep_vn2*vn(2,:,:)/rhon
      vrep_pn = vrepb*3._r8/(two*rhon)
     $     + vrepa*64._r8/(9._r8*pi*rhon*vrepb)

      qinx_rho = q_inx*(one/rho + scx_rho/s_cx + vrep_rho/vrep)
      qinx_p = q_inx*(scx_p/s_cx + vrep_p/vrep)
      qinx_rhon = q_inx*(one/rhon + scx_rhon/s_cx + vrep_rhon/vrep)
      qinx_pn = q_inx*(scx_pn/s_cx + vrep_pn/vrep)
      
      qinx_mi1 = q_inx*(scx_vi1/(rho*s_cx) + vrep_vi1/(rho*vrep))
      qinx_mi2 = q_inx*(scx_vi2/(rho*s_cx) + vrep_vi2/(rho*vrep))
      qinx_mn1 = q_inx*(scx_vn1/(rhon*s_cx) + vrep_vn1/(rhon*vrep))
      qinx_mn2 = q_inx*(scx_vn2/(rhon*s_cx) + vrep_vn2/(rhon*vrep))

c-----------------------------------------------------------------------
c     set derivatives for thermal energy  transfer, Q_nix.
c-----------------------------------------------------------------------
      vrepa = 3._r8*ti_frac*p/(two*rho)
      vrepb = SQRT(MAX(8._r8*pn/(pi*rhon)
     $     + 128._r8*ti_frac*p/(9._r8*pi*rho) + dv**2,eps_loc))
      vrep = vrepa*vrepb

      vrep_vi1 = vrepa*vin(1,:,:)/vrepb
      vrep_vi2 = vrepa*vin(2,:,:)/vrepb
      vrep_vn1 = -vrep_vi1
      vrep_vn2 = -vrep_vi2

      vrep_rho = -vrepb*3._r8*ti_frac*p/(two*rho**2)
     $     - vrepa*64._r8*ti_frac*p/(9._r8*pi*rho**2*vrepb)
     $     - vrep_vi1*vi(1,:,:)/rho - vrep_vi2*vi(2,:,:)/rho
      vrep_p = vrepb*3._r8*ti_frac/(two*rho)
     $     + vrepa*64._r8*ti_frac/(9._r8*pi*rho*vrepb)

      vrep_rhon = -vrepa*4._r8*pn/(pi*vrepb*rhon**2)
     $     - vrep_vn1*vn(1,:,:)/rhon - vrep_vn2*vn(2,:,:)/rhon
      vrep_pn = vrepa*4._r8/(pi*rhon*vrepb)

      qnix_rho = q_nix*(one/rho + scx_rho/s_cx + vrep_rho/vrep)
      qnix_p = q_nix*(scx_p/s_cx + vrep_p/vrep)
      qnix_rhon = q_nix*(one/rhon + scx_rhon/s_cx + vrep_rhon/vrep)
      qnix_pn = q_nix*(scx_pn/s_cx + vrep_pn/vrep)
      
      qnix_mi1 = q_nix*(scx_vi1/(rho*s_cx) + vrep_vi1/(rho*vrep))
      qnix_mi2 = q_nix*(scx_vi2/(rho*s_cx) + vrep_vi2/(rho*vrep))
      qnix_mn1 = q_nix*(scx_vn1/(rhon*s_cx) + vrep_vn1/(rhon*vrep))
      qnix_mn2 = q_nix*(scx_vn2/(rhon*s_cx) + vrep_vn2/(rhon*vrep))

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pn_grq_u
c-----------------------------------------------------------------------
c     subprogram 4. pn_bc_jac.
c     computes jacobian for a boundary condition.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pn_bc_jac(edge_type,t,x,y,nhat,
     $     u,ux,uy,uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy,uxx,uyy,uxy,nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy

      INTEGER :: i
      REAL(r8), PARAMETER :: du=1.e-6_r8
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: 
     $     u2,ux2,uy2,uxx2,uyy2,uxy2,f,f2
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
      uxy2=uxy

      CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy,uxx,uyy,uxy,f)

      DO i=1,SIZE(u,1)
         u2(i,:,:) = u(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u2,ux,uy,uxx,uyy,
     $        uxy,f2)
         u2(i,:,:) = u(i,:,:)
         c_u(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

      DO i=1,SIZE(ux,1)
         ux2(i,:,:) = ux(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux2,uy,uxx,uyy,
     $        uxy,f2)
         ux2(i,:,:) = ux(i,:,:)
         c_ux(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

      DO i=1,SIZE(uy,1)
         uy2(i,:,:) = uy(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy2,uxx,uyy,
     $        uxy,f2)
         uy2(i,:,:) = uy(i,:,:)
         c_uy(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

      DO i=1,SIZE(uxx,1)
         uxx2(i,:,:) = uxx(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy,uxx2,uyy,
     $        uxy,f2)
         uxx2(i,:,:) = uxx(i,:,:)
         c_uxx(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

      DO i=1,SIZE(uyy,1)
         uyy2(i,:,:) = uyy(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy,uxx,uyy2,
     $        uxy,f2)
         uyy2(i,:,:) = uyy(i,:,:)
         c_uyy(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

      DO i=1,SIZE(uxy,1)
         uxy2(i,:,:) = uxy(i,:,:) + du
         CALL physics_edge_rhs(edge_type,t,x,y,nhat,u,ux,uy,uxx,uyy,
     $        uxy2,f2)
         uxy2(i,:,:) = uxy(i,:,:)
         c_uxy(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pn_bc_jac
      END MODULE pn_mod
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
      USE pn_mod
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
      NAMELIST/pn_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,init_type,cylinder,eta,eta_case,mu_min,kappa_min,
     $     r_eta,etavac,mu,mun,kappa_par,kappa_perp,visc_case,r_visc,
     $     ddiff,beta0,rhomin,pmin,b0,n0,L0,bx,lx,ly,cc,viscn_case,Tfac,
     $     xmin,epsilon,equilfile,interp,coilfile,atom,kappan_case,
     $     kappa_case,initrho,initp,ion_fac,recomb_fac,targ_type,eta_hr,
     $     initrhon,psource,kappa_n,initpn,cx_fac,rhonmin,te_frac,rad,
     $     civ_fac,initv,initTn,gr_curve,kappa_max,mu_sv,mun_sv,rhofac,
     $     etafac,ifcoils
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=pn_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="pn"
      nqty=10
      nqty_schur=0
c-----------------------------------------------------------------------
c     set cylinder to true.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("trans_test")
         cylinder=.TRUE.
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
      USE pn_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty

      REAL(r8) :: tnorm,T0
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      static(6)=.TRUE.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(atom,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(equilfile,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(interp,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(coilfile,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(eta_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(kappa_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(kappan_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(visc_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(viscn_case,16,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(targ_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(cylinder,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(rad,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(ifcoils,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(etafac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(cc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(r_eta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(r_visc,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Tfac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rhofac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(eta_hr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu_sv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mun_sv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(etavac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu_min,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_min,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_max,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mun,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_n,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_par,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kappa_perp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ddiff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rhomin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rhonmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(te_frac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(b0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(n0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(L0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(bx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(xmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(epsilon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initrho,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initp,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ion_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(civ_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(recomb_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(cx_fac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initrhon,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initpn,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(initTn,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(psource,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     define problem dependent parameters
c-----------------------------------------------------------------------
      IF(cylinder)cyl_fac=1
      SELECT CASE(atom)
      CASE("deuterium")
         mi = 3.345e-27
         hs_diam = .240e-9
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
         hs_diam = .308e-9
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

      tnorm=L0*SQRT(mu0*n0*mi)/b0
      v_civ=SQRT(two*phi_ion*qe/mi)/(L0/tnorm)
      ti_frac = one - te_frac
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
c     - mu_norm absorbs the constants in front of Braginskii viscosity.
c     - mun_norm absorbs the constants in front of hard-sphere
c     viscosity.
c     for neutral thermal conduction:
c     - kn_norm absorbs the constants in front of hard-sphere thermal
c     conduction.
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
c-----------------------------------------------------------------------
      T0 = b0**2/(n0*mu0*qe)
      v_chod_norm=SQRT(mi/(mu0*n0))/(qe*L0)
      etac_norm=me/(qe*L0*b0*SQRT(ep0*mu0))
      etas_norm=5.e-5*lnlam*tnorm
     $     *(two*n0*mu0*qe)**1.5/(mu0*L0**2*b0**3)
      IF(te_frac > 0)THEN
         etas_norm = etas_norm/te_frac**1.5
      ELSE
         etas_norm=0.
      ENDIF
      ke_norm = 6.05e22/lnlam*T0**2.5*tnorm/(n0*L0**2)
      xe_norm = 6.05e22/lnlam*T0**1.5*b0/n0
      ki_norm = 3.35e-6/(lnlam*SQRT(mi*mp))*T0**2.5*tnorm/(n0*L0**2)
      xi_norm = 3.35e-6/(lnlam*SQRT(mi*mp))*T0**1.5*b0/n0
      kn_norm = 1.25*SQRT(qe*T0)/(hs_diam**2*SQRT(twopi*mi))
     $     *tnorm/(n0*L0**2)
      mu_norm = 3.214e-6/lnlam*SQRT(mi/mp)*T0**2.5*tnorm/(n0*mi*L0**2)
      mun_norm = SQRT(qe*T0)/(two*hs_diam**2*SQRT(twopi/mi))
     $     *tnorm/(L0**2*n0*mi)
      ion_norm=a_vor*n0*tnorm
      recomb_norm = 2.6e-19*n0*tnorm/SQRT(phi_ion)
      cx_norm=L0*n0*7.15e-20
      cx_c3=(cx_c2 - cx_c1*LOG(L0/tnorm))/cx_c1
      Top_norm=T0/phi_ion

      IF(1==2)THEN
         IF(mpi_rank==0)
     $        write(*,*) "WARNING!! Overriding norms for FD test!"
         v_chod_norm = 1.1
         etac_norm = 1.2
         etas_norm = 1.3
         ke_norm = 1.4
         xe_norm = 1.5
         ki_norm = 1.6
         xi_norm = 1.7
         kn_norm = 1.8
         mu_norm = 1.9
         mun_norm = 2.0
         ion_norm = 2.1
         recomb_norm = 2.2
         cx_norm = 2.3
      ENDIF

      SELECT CASE(init_type)
      CASE("trans_test")
         CALL extra_read_marklin(equilfile,interp,1._r8,
     $        equil_bc,equil,equilxy)
c-----------------------------------------------------------------------
c     set rhomin and pmin to default values if not set by user.
c-----------------------------------------------------------------------
         IF(pmin==0)pmin=5.e-3
         IF(rhomin==0)rhomin=SQRT(pmin)
c-----------------------------------------------------------------------
c     read and normalize coil information.
c-----------------------------------------------------------------------
         IF(ifcoils)THEN
            CALL extra_coilalloc(coilfile,coils)
            coils%tfire=coils%tfire/tnorm
            coils%tqtr=coils%tqtr/tnorm
            coils%tcrow=coils%tcrow/tnorm
            coils%vc0=coils%vc0/(L0**2*b0/tnorm)
            coils%zb=coils%zb/L0
            coils%ze=coils%ze/L0
         ENDIF
      END SELECT
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
      USE pn_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u
      
      REAL(r8) :: psimax,x1,x2,x3,x4,dx,lam
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: psi,trans
c-----------------------------------------------------------------------
c     initialize u.
c-----------------------------------------------------------------------
      u=zero
      SELECT CASE(init_type)
      CASE("trans_test")
         CALL pn_equil(x,y,u)
c-----------------------------------------------------------------------
c     initialize uniform FRC motion within an axial region that
c     encompasses the FRC.  Speed smoothly drops to zero outside that
c     region.
c-----------------------------------------------------------------------
         psi = -y*u(2,:,:)
         psimax = .3064_r8*SQRT(Tfac*rhofac)
         u(4,:,:) = u(1,:,:)*initv

         WHERE(psi > 0)
     $        u(4,:,:) = u(1,:,:)*initv*half*(one + COS(pi*psi/psimax))

         x1 = 3.55_r8
         x2 = 4.43_r8
         dx = x2 - x1
         WHERE(ABS(x) >= x2)
            u(4,:,:) = 0
         ELSEWHERE(ABS(x) > x1)
            u(4,:,:) = u(4,:,:)*half*(one + COS(pi*ABS(x)/dx))
         ENDWHERE

         x1 = 4.43_r8
         x2 = 4.88_r8
         dx = x2 - x1
         x3 = 7.53_r8
         x4 = x3 + dx
         SELECT CASE(targ_type)
         CASE("block","gauss")
            u(7,:,:) = initrhon - rhomin
            IF(targ_type == "gauss") u(7,:,:) = u(7,:,:)*EXP(-(3.*y)**2)
            WHERE(x > x4)
               u(7,:,:) = 0
            ELSEWHERE(x > x3)
               trans = half*(1 + COS(pi*(x - x3)/dx))
               u(7,:,:) = u(7,:,:)*trans
            ELSEWHERE(x > x2)
               u(7,:,:) = u(7,:,:)
            ELSEWHERE(x > x1)
               trans = half*(1 - COS(pi*(x - x1)/dx))
               u(7,:,:) = u(7,:,:)*trans
            ELSEWHERE   
               u(7,:,:) = 0
            ENDWHERE
         CASE("gauss-rz")
            x1=6.28_r8
            lam=.709_r8
            u(7,:,:) = initrhon
            u(7,:,:) = u(7,:,:)*EXP(-(y/lam)**2)*EXP(-((x-x1)/lam)**2)
     $           + rhomin
         CASE DEFAULT
            u(7,:,:) = initrhon - rhomin
         END SELECT
         u(7,:,:) = u(7,:,:) + rhomin
         u(10,:,:) = initTn*u(7,:,:)
      CASE("uniform")
         u(1,:,:) = initrho
         u(2,:,:) = -half*y*b0
         u(3,:,:) = initp
         u(7,:,:) = initrhon
      CASE("civ")
         u(1,:,:) = initrho
         u(2,:,:) = x*1.e-1
         u(3,:,:) = initp
         u(7,:,:) = 1.e-4
         u(10,:,:) = 1.e-4
         WHERE(x > 2. .AND. x < 2.2)
            u(7,:,:) = u(7,:,:) + initrhon*half
     $           *(1-COS(pi*(x - lx/2)/.2))
            u(10,:,:) = u(10,:,:) + initpn*half
     $           *(1-COS(pi*(x - lx/2)/.2))
         ELSEWHERE(x >= 2.2)
            u(7,:,:) = u(7,:,:) + initrhon
            u(10,:,:) = u(10,:,:) + initpn
         ENDWHERE
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
      USE pn_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     set boundary conditions.
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
      CASE("trans_test")
         top%bc_type(1)="normflux"
         top%bc_type(4)="zeroflux"
         top%bc_type(7)="zeroflux"
         top%bc_type(8)="zeroflux"
         top%static=.TRUE.
         top%static(2)=.FALSE.

         bottom%static=.TRUE.

         left%bc_type(7)="zeroflux"
         left%static=.TRUE.
         left%static(2)=.FALSE.

         right%bc_type(7)="zeroflux"
         right%static=.TRUE.
         right%static(2)=.FALSE.
      CASE("uniform")
         top%bc_type(1)="natural"
         top%bc_type(7)="natural"
         top%static=.TRUE.
         top%static(2)=.FALSE.

         bottom%static=.TRUE.
      CASE("uniform2")
         top%bc_type(7)="natural"
         top%static=.TRUE.
         top%static(1:2)=.FALSE.

         bottom%static=.TRUE.
      CASE("civ")
c         left%bc_type(7)="zeroflux"
         left%static=.TRUE.
         left%static(2)=.FALSE.

c         right%bc_type(7)="zeroflux"
         right%static=.TRUE.
         right%static(2)=.FALSE.
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
      USE pn_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: volt,br,bz,bsq
      REAL(r8), DIMENSION(2,SIZE(x,1),SIZE(x,2)) :: vi
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
c-----------------------------------------------------------------------
c     boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("trans_test")
         SELECT CASE(lrtb)
         CASE("top")
            volt=0
            IF(ifcoils) CALL extra_coileval(t,x,coils,volt)
            br = ux(2,:,:)
            bz = -(uy(2,:,:) + u(2,:,:)/y)
            bsq = br**2 + bz**2
c-----------------------------------------------------------------------
c           density flux = rhomin*vr.
c-----------------------------------------------------------------------
            WHERE(u(5,:,:)/MAX(u(1,:,:),1.e-8_r8) < 0)
               c(1,:,:) = rhomin*u(5,:,:)/u(1,:,:)
            ELSEWHERE
               c(1,:,:) = u(5,:,:)
            ENDWHERE
c-----------------------------------------------------------------------
c           specify electric field based on Ephi = -V/(2 pi R).
c           Ephi = -d/dt(Aphi) = d/dt(u2), so assuming R=1, 
c           d/dt(u2) = -V/(2 pi).
c-----------------------------------------------------------------------
            c(2,:,:) = -volt/(twopi*y)

c-----------------------------------------------------------------------
c           dT/dn = 0.
c-----------------------------------------------------------------------
            c(3,:,:) = u(3,:,:)*uy(1,:,:) - u(1,:,:)*uy(3,:,:)

c-----------------------------------------------------------------------
c           T = T0.
c-----------------------------------------------------------------------
c            c(3,:,:) = u(3,:,:)/u(1,:,:) - pmin/rhomin

c-----------------------------------------------------------------------
c           zero normal flux for axial momentum.
c-----------------------------------------------------------------------
            c(4,:,:) = 0
c-----------------------------------------------------------------------
c           zero axial momentum.
c-----------------------------------------------------------------------
c            c(4,:,:) = u(4,:,:)

c-----------------------------------------------------------------------
c           flux convection BC for radial momentum.
c           set vr to ExB/B^2
c-----------------------------------------------------------------------
            c(5,:,:) = u(5,:,:)/u(1,:,:) + volt*bz/(twopi*y*bsq)
c-----------------------------------------------------------------------
c           jphi=0.
c-----------------------------------------------------------------------
            c(6,:,:) = u(6,:,:)
c-----------------------------------------------------------------------
c           neutral BC:
c           zero flux for density.
c           zero flux for axial mom.
c           zero normal flow.
c           dT/dn = 0.
c-----------------------------------------------------------------------
            c(7,:,:) = 0
            c(8,:,:) = 0
            c(9,:,:) = u(9,:,:)
            c(10,:,:) = u(10,:,:)*uy(7,:,:) - u(7,:,:)*uy(10,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=uy(7,:,:)
            c(8,:,:)=uy(8,:,:)
            c(9,:,:)=u(9,:,:)
            c(10,:,:)=uy(10,:,:)
         END SELECT
      CASE("civ")
         DO i=1,2
            vi(i,:,:)=u(i+3,:,:)/u(1,:,:)
         ENDDO
         
         SELECT CASE(lrtb)
         CASE("left")
            c(1,:,:) = ux(1,:,:)
            IF(t < 1.e-1)THEN
               volt = beta0*t
            ELSE
               volt = beta0*1.e-1
            ENDIF
            c(2,:,:) = -volt/twopi
            c(3,:,:) = ux(3,:,:)
            c(4,:,:) = volt/twopi - vi(1,:,:)*ux(2,:,:)
            c(5,:,:) = u(5,:,:)
            c(6,:,:) = u(6,:,:)

            c(7,:,:)=ux(7,:,:)
            c(8,:,:)=u(8,:,:)
            c(9,:,:)=ux(9,:,:)
            c(10,:,:)=ux(10,:,:)
         CASE("right")
            c(1,:,:) = u(1,:,:) - rhomin
            c(3,:,:) = u(3,:,:) - pmin
            c(4,:,:) = u(4,:,:)
            c(5,:,:) = u(5,:,:)
            c(6,:,:) = u(6,:,:)

            c(7,:,:)=ux(7,:,:)
            c(8,:,:)=u(8,:,:)
            c(9,:,:)=ux(9,:,:)
            c(10,:,:)=ux(10,:,:)
         END SELECT
      CASE("uniform")
         SELECT CASE(lrtb)
         CASE("top")
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=uy(7,:,:)
         END SELECT
      CASE("uniform2")
         SELECT CASE(lrtb)
         CASE("top")
            c(3,:,:)=u(3,:,:)-0.1022*u(1,:,:)
            c(4,:,:)=u(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
         CASE("bottom")
            c(1,:,:)=uy(1,:,:)
            c(2,:,:)=u(2,:,:)
            c(3,:,:)=uy(3,:,:)
            c(4,:,:)=uy(4,:,:)
            c(5,:,:)=u(5,:,:)
            c(6,:,:)=u(6,:,:)
            c(7,:,:)=uy(7,:,:)
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
      USE pn_mod
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

      CALL pn_bc_jac(lrtb,t,x,y,nhat,
     $     u,ux,uy,uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)
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
      USE pn_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y

      INTEGER :: iqty
c-----------------------------------------------------------------------
c     set mass,mass_x,mass_y(iqty,jqty,:,:), ... to coupling 
c     mass matrices for du/dt, d(du/dx)/dt, and d(du/dy)/dt
c-----------------------------------------------------------------------
      mass=0
      mass_x=0
      mass_y=0

      SELECT CASE(init_type)
      CASE("trans_test")
         SELECT CASE(lrtb)
         CASE("top")
            mass(2,2,:,:)=one
         END SELECT
      CASE("uniform")
         SELECT CASE(lrtb)
         CASE("top")
            mass(2,2,:,:)=one
         END SELECT
      CASE("uniform2")
         SELECT CASE(lrtb)
         CASE("top")
            mass(1,1,:,:)=one
            mass(2,2,:,:)=one
         END SELECT
      CASE("civ")
         SELECT CASE(lrtb)
         CASE("left","right")
            mass(2,2,:,:)=one
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
      USE pn_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      INTEGER :: i
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: Tix,Tiy,kperp,kfac,
     $     eta_local,r_fac,r_faci,visc,viscn,nil,b1,b2,Bsq,g_i,g_r,
     $     g_x,q_inx,q_nix,n_inv,nx_inv,ny_inv,nn_inv,nnx_inv,nny_inv,
     $     kei,ken,Tnx,Tny,vivn,kn_local,ones,rad_fac
      REAL(r8), DIMENSION(2,SIZE(u,2),SIZE(u,3)) :: vi,vix,viy,vn,vnx,
     $     vny,r_inx,r_nix
      REAL(r8), DIMENSION(3,SIZE(u,2),SIZE(u,3)) :: BdotT
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
      nn_inv=one/u(7,:,:)
      nnx_inv=-ux(7,:,:)*nn_inv**2
      nny_inv=-uy(7,:,:)*nn_inv**2
c-----------------------------------------------------------------------
c     viscosity.
c-----------------------------------------------------------------------
      SELECT CASE(visc_case)
      CASE("braginskii")
         CALL transport_mubrag(u(1,:,:),u(3,:,:),ti_frac,mu_norm,mu_min,
     $        r_visc,y,visc)
      CASE("uniform")
         visc = mu
      CASE DEFAULT
         CALL program_stop("physics_rhs cannot recognize "
     $        //"visc_case = "//TRIM(visc_case))
      END SELECT

      SELECT CASE(viscn_case)
      CASE("hard_sphere")
         CALL transport_mun(u(7,:,:),u(10,:,:),mun_norm,viscn)
      CASE("uniform")
         viscn = mun
      CASE DEFAULT
         CALL program_stop("physics_rhs cannot recognize "
     $        //"viscn_case = "//TRIM(viscn_case))
      END SELECT
c-----------------------------------------------------------------------
c     resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta(eta_case,ly-y,u(1,:,:),u(3,:,:),u(6,:,:),cc,
     $     etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,etavac,
     $     eta_local)
      eta_local=eta_local*etafac
c-----------------------------------------------------------------------
c     radiation.
c-----------------------------------------------------------------------
      rad_fac = zero
      IF(rad) WHERE(u(1,:,:) < rho_rad) rad_fac =
     $     half*(one + COS(pi*(rho_rad - u(1,:,:))/rho_rad))
c      rad_fac=one
c-----------------------------------------------------------------------
c     temperature gradients.
c-----------------------------------------------------------------------
      Tix = ux(3,:,:)*n_inv + u(3,:,:)*nx_inv
      Tiy = uy(3,:,:)*n_inv + u(3,:,:)*ny_inv
      Tnx = ux(10,:,:)*nn_inv + u(10,:,:)*nnx_inv
      Tny = uy(10,:,:)*nn_inv + u(10,:,:)*nny_inv
c-----------------------------------------------------------------------
c     magnetic fields.
c-----------------------------------------------------------------------
      nil=0
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT(b1,b2,nil,Tix,Tiy,BdotT)
         CALL transport_kbrag(u(1,:,:),u(3,:,:),te_frac,Bsq,ke_norm,
     $        ki_norm,xe_norm,xi_norm,kappa_min,kappa_max,kperp,kfac)
      CASE("anisotropic")
         CALL transport_BdotT(b1,b2,nil,Tix,Tiy,BdotT)
         CALL transport_setkaniso(kappa_par,kappa_perp,Bsq,kperp,kfac)
      CASE("isotropic")
         BdotT=0
         kperp=kappa_par
         kfac=0
      CASE DEFAULT
         CALL program_stop("physics_rhs: cannot recognize "
     $        //"kappa_case = "//TRIM(kappa_case))
      END SELECT

      SELECT CASE(kappan_case)
      CASE("hard_sphere")
         kn_local = viscn*kn_norm/mun_norm
      CASE("uniform")
         kn_local = kappa_n
      CASE DEFAULT
         CALL program_stop("physics_rhs cannot recognize "
     $        //"kappan_case = "//TRIM(kappan_case))
      END SELECT
c-----------------------------------------------------------------------
c     velocities and their gradients.
c-----------------------------------------------------------------------
      DO i=1,2
         vi(i,:,:)=u(i+3,:,:)*n_inv
         vix(i,:,:)=ux(i+3,:,:)*n_inv + u(i+3,:,:)*nx_inv
         viy(i,:,:)=uy(i+3,:,:)*n_inv + u(i+3,:,:)*ny_inv

         vn(i,:,:)=u(i+7,:,:)*nn_inv
         vnx(i,:,:)=ux(i+7,:,:)*nn_inv + u(i+7,:,:)*nnx_inv
         vny(i,:,:)=uy(i+7,:,:)*nn_inv + u(i+7,:,:)*nny_inv
      ENDDO
      kei = half*(vi(1,:,:)**2 + vi(2,:,:)**2)
      ken = half*(vn(1,:,:)**2 + vn(2,:,:)**2)
      vivn = vi(1,:,:)*vn(1,:,:) + vi(2,:,:)*vn(2,:,:)
c-----------------------------------------------------------------------
c     compute reaction collision terms.
c-----------------------------------------------------------------------
      CALL pn_grq(u(3,:,:),u(1,:,:),u(10,:,:),u(7,:,:),vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix)

c-----------------------------------------------------------------------
c     equations for density and -Az (cartesian) or -Aphi (cylindrical).
c-----------------------------------------------------------------------
      fx(1,:,:) = r_fac*(u(4,:,:) - ddiff*ux(1,:,:))
      fy(1,:,:) = r_fac*(u(5,:,:) - ddiff*uy(1,:,:))
      s(1,:,:) = r_fac*(g_i - g_r)

      fx(2,:,:) = eta_hr*ux(6,:,:)
      fy(2,:,:) = eta_hr*uy(6,:,:)
      s(2,:,:) = vi(2,:,:)*(b1+bx) - vi(1,:,:)*b2 + eta_local*u(6,:,:)
     $     + eta_hr*cyl_fac*r_faci*(-uy(6,:,:) + u(6,:,:)*r_faci)
c-----------------------------------------------------------------------
c     pressure equation.
c-----------------------------------------------------------------------
      fx(3,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi(1,:,:)
     $     - kfac*BdotT(1,:,:) - kperp*Tix)
      fy(3,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi(2,:,:)
     $     - kfac*BdotT(2,:,:) - kperp*Tiy)
      s(3,:,:)=r_fac*(vi(1,:,:)*ux(3,:,:) + vi(2,:,:)*uy(3,:,:)
     $     + rad_fac*
     $     (eta_local*u(6,:,:)**2
     $     + visc*(two*vix(1,:,:)**2 
     $     + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2
     $     + vix(2,:,:)**2 
     $     + two*viy(1,:,:)*vix(2,:,:))
     $     + cyl_fac*two*visc*(vi(2,:,:)*r_faci)**2
     $     + mu_sv*(two*vix(1,:,:)**2*ABS(vix(1,:,:))
     $     + two*viy(2,:,:)**2*ABS(viy(2,:,:)))
     $     + cyl_fac*mu_sv*two*r_faci**3*vi(2,:,:)**2*ABS(vi(2,:,:))
     $     + eta_hr*(ux(6,:,:)**2 + uy(6,:,:)**2)
     $     + cyl_fac*eta_hr*(u(6,:,:)*r_faci)**2)
     $     + (g_i + g_x)*(kei + ken - vivn) - g_i*phi_eff
     $     + 3._r8*half*(g_i*u(10,:,:)*nn_inv 
     $     - g_r*ti_frac*u(3,:,:)*n_inv)
     $     + r_inx(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + r_inx(2,:,:)*(vn(2,:,:) - vi(2,:,:))
     $     + q_inx - q_nix + psource)
c-----------------------------------------------------------------------
c     momentum equations.
c-----------------------------------------------------------------------
      fx(4,:,:)=r_fac*(u(4,:,:)*vi(1,:,:)
     $     - two*visc*vix(1,:,:) + u(3,:,:)
     $     - two*mu_sv*vix(1,:,:)*ABS(vix(1,:,:)))
      fy(4,:,:)=r_fac*(u(4,:,:)*vi(2,:,:) 
     $     - visc*(viy(1,:,:) + vix(2,:,:)))
      s(4,:,:)=r_fac*(-u(6,:,:)*b2 - g_r*vi(1,:,:)
     $     + g_i*vn(1,:,:) + g_x*(vn(1,:,:) - vi(1,:,:))
     $     + r_inx(1,:,:) - r_nix(1,:,:))

      fx(5,:,:)=r_fac*(u(5,:,:)*vi(1,:,:)
     $     - visc*(vix(2,:,:) + viy(1,:,:)))
      fy(5,:,:)=r_fac*(u(5,:,:)*vi(2,:,:) - two*visc*viy(2,:,:)
     $     - two*mu_sv*viy(2,:,:)*ABS(viy(2,:,:)))
      s(5,:,:)=r_fac*(u(6,:,:)*(b1+bx) - uy(3,:,:) - g_r*vi(2,:,:) 
     $     + g_i*vn(2,:,:) + g_x*(vn(2,:,:) - vi(2,:,:))
     $     + r_inx(2,:,:) - r_nix(2,:,:))
     $     - cyl_fac*two*r_faci*(visc*vi(2,:,:)
     $     + mu_sv*r_faci*vi(2,:,:)*ABS(vi(2,:,:)))
c-----------------------------------------------------------------------
c     current.
c-----------------------------------------------------------------------
      fx(6,:,:)=b2
      fy(6,:,:)=-b1
      s(6,:,:)=u(6,:,:)
c-----------------------------------------------------------------------
c     neutral density.
c-----------------------------------------------------------------------
      fx(7,:,:)=r_fac*(u(8,:,:) - ddiff*ux(7,:,:))
      fy(7,:,:)=r_fac*(u(9,:,:) - ddiff*uy(7,:,:))

      s(7,:,:)=r_fac*(-g_i+g_r)
c-----------------------------------------------------------------------
c     neutral momentum.
c-----------------------------------------------------------------------
      fx(8,:,:)=r_fac*(u(8,:,:)*vn(1,:,:)
     $     - two*viscn*vnx(1,:,:) + u(10,:,:)
     $     - two*mun_sv*vnx(1,:,:)*ABS(vnx(1,:,:)))
      fy(8,:,:)=r_fac*(u(8,:,:)*vn(2,:,:) 
     $     - viscn*(vny(1,:,:) + vnx(2,:,:)))
      s(8,:,:)=r_fac*(-g_i*vn(1,:,:) + g_r*vi(1,:,:) 
     $     + g_x*(vi(1,:,:) - vn(1,:,:))
     $     + r_nix(1,:,:) - r_inx(1,:,:))

      fx(9,:,:)=r_fac*(u(9,:,:)*vn(1,:,:)
     $     - viscn*(vnx(2,:,:) + vny(1,:,:)))
      fy(9,:,:)=r_fac*(u(9,:,:)*vn(2,:,:) - two*viscn*vny(2,:,:)
     $     - two*mun_sv*vny(2,:,:)*ABS(vny(2,:,:)))
      s(9,:,:)=r_fac*(-g_i*vn(2,:,:) + g_r*vi(2,:,:) - uy(10,:,:)
     $     + g_x*(vi(2,:,:) - vn(2,:,:))
     $     + r_nix(2,:,:) - r_inx(2,:,:))
     $     - cyl_fac*two*r_faci*(viscn*vn(2,:,:)
     $     + mun_sv*r_faci*vn(2,:,:)*ABS(vn(2,:,:)))
c-----------------------------------------------------------------------
c     neutral pressure.
c-----------------------------------------------------------------------
      fx(10,:,:)=r_fac*(gamma_fac*u(10,:,:)*vn(1,:,:)
     $     - kn_local*Tnx)
      fy(10,:,:)=r_fac*(gamma_fac*u(10,:,:)*vn(2,:,:)
     $     - kn_local*Tny)
      s(10,:,:)=r_fac*(vn(1,:,:)*ux(10,:,:) + vn(2,:,:)*uy(10,:,:)
     $     + viscn*(two*vnx(1,:,:)**2 + two*vny(2,:,:)**2 
     $     + vny(1,:,:)**2 + vnx(2,:,:)**2 + two*vny(1,:,:)*vnx(2,:,:))
     $     + cyl_fac*two*viscn*(vn(2,:,:)*r_faci)**2
     $     + (g_r + g_x)*(ken + kei - vivn)
     $     + mun_sv*(two*vnx(1,:,:)**2*ABS(vnx(1,:,:))
     $     + two*vny(2,:,:)**2*ABS(vny(2,:,:)))
     $     + cyl_fac*mun_sv*two*r_faci**3*vn(2,:,:)**2*ABS(vn(2,:,:))
     $     + 3._r8*half*(-g_i*u(10,:,:)*nn_inv 
     $     + g_r*ti_frac*u(3,:,:)*n_inv)
     $     + r_nix(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + r_nix(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + q_nix - q_inx + psource)

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
      USE pn_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      INTEGER :: i,j
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3)) :: r_fac,r_faci,Tix,Tiy,
     $     Tix_un,Tiy_un,Ti_un,kperp,kfac,kpar_u1,kperp_bsq,kperp_u1,
     $     kperp_u3,kpar_u3,eta_local,visc,visc_rho,visc_p,viscn,
     $     viscn_rhon,viscn_pn,eta_u1,eta_u3,eta_u6,nil,b1,b2,Bsq,
     $     kn_local,kn_rhon,kn_pn,
     $     g_i,g_i_rho,g_i_p,g_i_rhon,
     $     g_r,g_r_rho,g_r_p,
     $     g_x,g_x_p,g_x_pn,g_x_rho,g_x_rhon,q_inx,q_nix,
     $     n_inv,nx_inv,ny_inv,nn_inv,nnx_inv,nny_inv,kei,
     $     ken,kei_mix,kei_miy,ken_mnx,ken_mny,kei_rho,ken_rhon,Tnx,Tny,
     $     Tn_un,Tnx_un,Tny_un,vivn,vivn_rho,vivn_rhon,vivn_mix,
     $     vivn_miy,vivn_mnx,vivn_mny,g_x_mi1,g_x_mi2,g_x_mn1,
     $     g_x_mn2,
     $     qinx_rho,qinx_p,qinx_rhon,qinx_pn,qinx_mi1,
     $     qinx_mi2,qinx_mn1,qinx_mn2,
     $     qnix_rho,qnix_p,qnix_rhon,qnix_pn,qnix_mi1,
     $     qnix_mi2,qnix_mn1,qnix_mn2,
     $     rad_fac

      REAL(r8), DIMENSION(2,SIZE(u,2),SIZE(u,3)) :: vi,vix,viy,r_inx,
     $     r_nix,vi_un,vix_un,viy_un,vn,vnx,vny,vn_un,vnx_un,vny_un,
     $     rinx_rho,rinx_p,rinx_rhon,rinx_pn,rinx_mi1,
     $     rinx_mi2,rinx_mn1,rinx_mn2,
     $     rnix_rho,rnix_p,rnix_rhon,rnix_pn,rnix_mi1,
     $     rnix_mi2,rnix_mn1,rnix_mn2
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

      nn_inv=one/u(7,:,:)
      nnx_inv=-ux(7,:,:)*nn_inv**2
      nny_inv=-uy(7,:,:)*nn_inv**2
c-----------------------------------------------------------------------
c     viscosity.
c-----------------------------------------------------------------------
      SELECT CASE(visc_case)
      CASE("braginskii")
         CALL transport_mubrag_u(u(1,:,:),u(3,:,:),ti_frac,mu_norm,
     $        mu_min,r_visc,y,visc,visc_rho,visc_p)
      CASE("uniform")
         visc = mu
         visc_rho = 0
         visc_p = 0
      END SELECT

      SELECT CASE(viscn_case)
      CASE("hard_sphere")
         CALL transport_mun_u(u(7,:,:),u(10,:,:),mun_norm,viscn,
     $        viscn_rhon,viscn_pn)
      CASE("uniform")
         viscn = mun
         viscn_rhon = 0
         viscn_pn = 0
      END SELECT
c-----------------------------------------------------------------------
c     resistivity.
c-----------------------------------------------------------------------
      CALL transport_seteta_u(eta_case,ly-y,u(1,:,:),u(3,:,:),u(6,:,:),
     $     cc,etas_norm,etac_norm,v_chod_norm,ly-r_eta,eta,etavac,
     $     eta_local,eta_u1,eta_u3,eta_u6)
      eta_local=eta_local*etafac
      eta_u1=eta_u1*etafac
      eta_u3=eta_u3*etafac
      eta_u6=eta_u6*etafac
c-----------------------------------------------------------------------
c     radiation.
c-----------------------------------------------------------------------
      rad_fac = zero
      IF(rad) WHERE(u(1,:,:) < rho_rad) rad_fac =
     $     half*(one + COS(pi*(rho_rad - u(1,:,:))/rho_rad))
c      rad_fac=one
c-----------------------------------------------------------------------
c     temperature gradients and derivatives.
c-----------------------------------------------------------------------
      Tix = ux(3,:,:)*n_inv + u(3,:,:)*nx_inv
      Tiy = uy(3,:,:)*n_inv + u(3,:,:)*ny_inv
      Ti_un=-u(3,:,:)*n_inv**2
      Tix_un=-(two*u(3,:,:)*nx_inv + ux(3,:,:)*n_inv)*n_inv
      Tiy_un=-(two*u(3,:,:)*ny_inv + uy(3,:,:)*n_inv)*n_inv

      Tnx = ux(10,:,:)*nn_inv + u(10,:,:)*nnx_inv
      Tny = uy(10,:,:)*nn_inv + u(10,:,:)*nny_inv
      Tn_un=-u(10,:,:)*nn_inv**2
      Tnx_un=-(two*u(10,:,:)*nnx_inv + ux(10,:,:)*nn_inv)*nn_inv
      Tny_un=-(two*u(10,:,:)*nny_inv + uy(10,:,:)*nn_inv)*nn_inv
c-----------------------------------------------------------------------
c     magnetic fields.
c-----------------------------------------------------------------------
      nil=0
      nil3=0
      b1 = -(cyl_fac*r_faci*u(2,:,:) + uy(2,:,:))
      b2 = ux(2,:,:)
      Bsq = b1**2 + b2**2
c-----------------------------------------------------------------------
c     heat conduction.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         CALL transport_BdotT_u(b1,b2,nil,Tix,Tiy,BdotT,
     $        BdotT_b1,BdotT_b2,nil3,BdotT_Tx,BdotT_Ty)
         CALL transport_kbrag_u(u(1,:,:),u(3,:,:),te_frac,Bsq,ke_norm,
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

      SELECT CASE(kappan_case)
      CASE("hard_sphere")
         kn_local = viscn*kn_norm/mun_norm
         kn_rhon = viscn_rhon*kn_norm/mun_norm
         kn_pn = viscn_pn*kn_norm/mun_norm
      CASE("uniform")
         kn_local = kappa_n
         kn_rhon = 0
         kn_pn = 0
      END SELECT
c-----------------------------------------------------------------------
c     velocities and their derivatives.
c-----------------------------------------------------------------------
      DO i=1,2
         vi(i,:,:)=u(i+3,:,:)*n_inv
         vix(i,:,:)=ux(i+3,:,:)*n_inv + u(i+3,:,:)*nx_inv
         viy(i,:,:)=uy(i+3,:,:)*n_inv + u(i+3,:,:)*ny_inv

         vi_un(i,:,:)=-vi(i,:,:)*n_inv
         vix_un(i,:,:) = -(vi(i,:,:)*nx_inv + vix(i,:,:)*n_inv)
         viy_un(i,:,:) = -(vi(i,:,:)*ny_inv + viy(i,:,:)*n_inv)

         vn(i,:,:)=u(i+7,:,:)*nn_inv
         vnx(i,:,:)=ux(i+7,:,:)*nn_inv + u(i+7,:,:)*nnx_inv
         vny(i,:,:)=uy(i+7,:,:)*nn_inv + u(i+7,:,:)*nny_inv

         vn_un(i,:,:)=-vn(i,:,:)*nn_inv
         vnx_un(i,:,:) = -(vn(i,:,:)*nnx_inv + vnx(i,:,:)*nn_inv)
         vny_un(i,:,:) = -(vn(i,:,:)*nny_inv + vny(i,:,:)*nn_inv)
      ENDDO
c-----------------------------------------------------------------------
c     kinetic energies and derivatives.
c-----------------------------------------------------------------------
      kei = half*(vi(1,:,:)**2 + vi(2,:,:)**2)
      kei_mix = vi(1,:,:)/u(1,:,:)
      kei_miy = vi(2,:,:)/u(1,:,:)
      kei_rho = -two*kei/u(1,:,:)

      ken = half*(vn(1,:,:)**2 + vn(2,:,:)**2)
      ken_mnx = vn(1,:,:)/u(7,:,:)
      ken_mny = vn(2,:,:)/u(7,:,:)
      ken_rhon = -two*ken/u(7,:,:)

      vivn = vi(1,:,:)*vn(1,:,:) + vi(2,:,:)*vn(2,:,:)
      vivn_rho = vn(1,:,:)*vi_un(1,:,:) + vn(2,:,:)*vi_un(2,:,:)
      vivn_rhon = vn_un(1,:,:)*vi(1,:,:) + vn_un(2,:,:)*vi(2,:,:)
      vivn_mix = n_inv*vn(1,:,:)
      vivn_miy = n_inv*vn(2,:,:)
      vivn_mnx = nn_inv*vi(1,:,:)
      vivn_mny = nn_inv*vi(2,:,:)
c-----------------------------------------------------------------------
c     compute reaction collision terms.
c-----------------------------------------------------------------------
      CALL pn_grq_u(u(3,:,:),u(1,:,:),u(10,:,:),u(7,:,:),vi,vn,
     $     g_i,g_r,g_x,r_inx,r_nix,q_inx,q_nix,
     $     g_i_rho,g_i_rhon,g_i_p,g_r_rho,g_r_p,
     $     g_x_rho,g_x_rhon,g_x_p,g_x_pn,g_x_mi1,g_x_mi2,g_x_mn1,
     $     g_x_mn2,
     $     rinx_rho,rinx_p,rinx_rhon,rinx_pn,rinx_mi1,
     $     rinx_mi2,rinx_mn1,rinx_mn2,
     $     rnix_rho,rnix_p,rnix_rhon,rnix_pn,rnix_mi1,
     $     rnix_mi2,rnix_mn1,rnix_mn2,
     $     qinx_rho,qinx_p,qinx_rhon,qinx_pn,qinx_mi1,
     $     qinx_mi2,qinx_mn1,qinx_mn2,
     $     qnix_rho,qnix_p,qnix_rhon,qnix_pn,qnix_mi1,
     $     qnix_mi2,qnix_mn1,qnix_mn2)

c-----------------------------------------------------------------------
c     equations for density and -Az (cartesian) or -Aphi (cylindrical).
c-----------------------------------------------------------------------
      fx_u(1,4,:,:)=r_fac
      fy_u(1,5,:,:)=r_fac
      fx_ux(1,1,:,:)=-r_fac*ddiff
      fy_uy(1,1,:,:)=-r_fac*ddiff

      s_u(1,1,:,:)=r_fac*(g_i_rho-g_r_rho)
      s_u(1,3,:,:)=r_fac*(g_i_p-g_r_p)
      s_u(1,7,:,:)=r_fac*g_i_rhon

      fx_ux(2,6,:,:) = eta_hr
      fy_uy(2,6,:,:) = eta_hr

      s_u(2,1,:,:)=vi_un(2,:,:)*(b1+bx) - vi_un(1,:,:)*b2
     $     + eta_u1*u(6,:,:)
      s_u(2,2,:,:)=-vi(2,:,:)*cyl_fac*r_faci
      s_u(2,3,:,:)=eta_u3*u(6,:,:)
      s_u(2,4,:,:)=-n_inv*b2
      s_u(2,5,:,:)=n_inv*(b1+bx)
      s_u(2,6,:,:)=eta_local + u(6,:,:)*eta_u6 
     $     + eta_hr*cyl_fac*r_faci**2
      s_ux(2,2,:,:)=-vi(1,:,:)
      s_uy(2,2,:,:)=-vi(2,:,:)
      s_uy(2,6,:,:)=-eta_hr*cyl_fac*r_faci

c-----------------------------------------------------------------------
c     pressure equation.
c-----------------------------------------------------------------------
      fx_u(3,1,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi_un(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*Tix_un 
     $     - kfac*BdotT_Ty(1,:,:)*Tiy_un
     $     - (kpar_u1 - kperp_u1)*BdotT(1,:,:) - kperp_u1*Tix)
      fx_u(3,2,:,:)=cyl_fac*(kfac*BdotT_b1(1,:,:)
     $     - two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_u(3,3,:,:)=r_fac*(gamma_fac*vi(1,:,:)
     $     - (kfac*BdotT_Tx(1,:,:) + kperp)*nx_inv
     $     - kfac*BdotT_Ty(1,:,:)*ny_inv
     $     - (kpar_u3 - kperp_u3)*BdotT(1,:,:) - kperp_u3*Tix)
      fx_u(3,4,:,:)=r_fac*gamma_fac*u(3,:,:)*n_inv

      fx_ux(3,1,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*Ti_un
      fx_ux(3,2,:,:)=r_fac*(-kfac*BdotT_b2(1,:,:)
     $     + two*b2*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_ux(3,3,:,:)=-r_fac*(kfac*BdotT_Tx(1,:,:) + kperp)*n_inv

      fx_uy(3,1,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*Ti_un
      fx_uy(3,2,:,:)=r_fac*(kfac*BdotT_b1(1,:,:)
     $     - two*b1*kperp_bsq*(BdotT(1,:,:) - Tix))
      fx_uy(3,3,:,:)=-r_fac*kfac*BdotT_Ty(1,:,:)*n_inv

      fy_u(3,1,:,:)=r_fac*(gamma_fac*u(3,:,:)*vi_un(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*Tiy_un
     $     - kfac*BdotT_Tx(2,:,:)*Tix_un
     $     - (kpar_u1 - kperp_u1)*BdotT(2,:,:) - kperp_u1*Tiy)
      fy_u(3,2,:,:)=cyl_fac*(kfac*BdotT_b1(2,:,:)
     $     - two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_u(3,3,:,:)=r_fac*(gamma_fac*vi(2,:,:)
     $     - (kfac*BdotT_Ty(2,:,:) + kperp)*ny_inv
     $     - kfac*BdotT_Tx(2,:,:)*nx_inv
     $     - (kpar_u3 - kperp_u3)*BdotT(2,:,:) - kperp_u3*Tiy)
      fy_u(3,5,:,:)=r_fac*gamma_fac*u(3,:,:)*n_inv

      fy_ux(3,1,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*Ti_un
      fy_ux(3,2,:,:)=r_fac*(-kfac*BdotT_b2(2,:,:)
     $     + two*b2*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_ux(3,3,:,:)=-r_fac*kfac*BdotT_Tx(2,:,:)*n_inv

      fy_uy(3,1,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*Ti_un
      fy_uy(3,2,:,:)=r_fac*(kfac*BdotT_b1(2,:,:)
     $     - two*b1*kperp_bsq*(BdotT(2,:,:) - Tiy))
      fy_uy(3,3,:,:)=-r_fac*(kfac*BdotT_Ty(2,:,:) + kperp)*n_inv

      s_u(3,1,:,:)=r_fac*(vi_un(1,:,:)*ux(3,:,:)
     $     + vi_un(2,:,:)*uy(3,:,:)
     $     + rad_fac*
     $     (eta_u1*u(6,:,:)**2
     $     + two*visc*(two*vix(1,:,:)*vix_un(1,:,:) 
     $     + two*viy(2,:,:)*viy_un(2,:,:) 
     $     + viy(1,:,:)*viy_un(1,:,:) + vix(2,:,:)*vix_un(2,:,:)
     $     + viy(1,:,:)*vix_un(2,:,:) + vix(2,:,:)*viy_un(1,:,:))
     $     + visc_rho*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:))
     $     + cyl_fac*two*r_faci**2*(visc*two*vi(2,:,:)*vi_un(2,:,:)
     $     + visc_rho*vi(2,:,:)**2)
     $     + mu_sv*(6.*vix(1,:,:)*vix_un(1,:,:)*ABS(vix(1,:,:))
     $     + 6.*viy(2,:,:)*viy_un(2,:,:)*ABS(viy(2,:,:)))
     $     + cyl_fac*mu_sv*6.*r_faci**3*vi(2,:,:)*vi_un(2,:,:)
     $     *ABS(vi(2,:,:)))
     $     + (g_i + g_x)*(kei_rho - vivn_rho) 
     $     + (g_i_rho + g_x_rho)*(kei + ken - vivn) - g_i_rho*phi_eff
     $     + 3._r8*half*(g_i_rho*u(10,:,:)*nn_inv
     $     + ti_frac*(g_r*u(3,:,:)*n_inv**2 - g_r_rho*u(3,:,:)*n_inv))
     $     + n_inv*(r_inx(1,:,:)*vi(1,:,:) + r_inx(2,:,:)*vi(2,:,:))
     $     + rinx_rho(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_rho(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + qinx_rho - qnix_rho)
      s_u(3,3,:,:)=r_fac*(rad_fac*(eta_u3*u(6,:,:)**2
     $     + visc_p*(two*vix(1,:,:)**2 + two*viy(2,:,:)**2 
     $     + viy(1,:,:)**2 + vix(2,:,:)**2 + two*viy(1,:,:)*vix(2,:,:))
     $     + cyl_fac*two*visc_p*(vi(2,:,:)*r_faci)**2)
     $     + (g_i_p + g_x_p)*(kei + ken - vivn) - g_i_p*phi_eff
     $     + 3._r8*half*(g_i_p*u(10,:,:)*nn_inv 
     $     - ti_frac*(g_r_p*u(3,:,:)*n_inv + g_r*n_inv))
     $     + rinx_p(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_p(2,:,:)*(vn(2,:,:) - vi(2,:,:)) + qinx_p - qnix_p)
      s_u(3,4,:,:)=r_fac*(ux(3,:,:)*n_inv 
     $     + rad_fac*(two*visc*(two*vix(1,:,:)*nx_inv 
     $     + (viy(1,:,:) + vix(2,:,:))*ny_inv)
     $     + mu_sv*6.*vix(1,:,:)*nx_inv*ABS(vix(1,:,:)))
     $     + (g_i + g_x)*(kei_mix - vivn_mix)
     $     + g_x_mi1*(kei + ken - vivn) - r_inx(1,:,:)*n_inv
     $     + rinx_mi1(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mi1(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + qinx_mi1 - qnix_mi1)
      s_u(3,5,:,:)=r_fac*(uy(3,:,:)*n_inv 
     $     + rad_fac*(two*visc*
     $     (two*viy(2,:,:)*ny_inv + (vix(2,:,:) + viy(1,:,:))*nx_inv)
     $     + cyl_fac*4._r8*visc*r_faci**2*vi(2,:,:)*n_inv
     $     + mu_sv*6.*viy(2,:,:)*ny_inv*ABS(viy(2,:,:))
     $     + cyl_fac*6.*mu_sv*r_faci
     $     *SIGN((vi(2,:,:)*r_faci)**2*n_inv,vi(2,:,:)))
     $     + (g_i + g_x)*(kei_miy - vivn_miy)
     $     + g_x_mi2*(kei + ken - vivn) - r_inx(2,:,:)*n_inv
     $     + rinx_mi2(1,:,:)*(vn(1,:,:) - vi(1,:,:)) 
     $     + rinx_mi2(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + qinx_mi2 - qnix_mi2)
      s_u(3,6,:,:)=r_fac*rad_fac*(two*eta_local*u(6,:,:) 
     $     + eta_u6*u(6,:,:)**2 + cyl_fac*eta_hr*two*u(6,:,:)*r_faci**2)
      s_u(3,7,:,:)=r_fac*((g_i + g_x)*(ken_rhon - vivn_rhon)
     $     + (g_i_rhon + g_x_rhon)*(kei + ken - vivn) 
     $     - g_i_rhon*phi_eff
     $     + 3._r8*half*u(10,:,:)*nn_inv*(g_i_rhon - g_i*nn_inv)
     $     - nn_inv*(r_inx(1,:,:)*vn(1,:,:) + r_inx(2,:,:)*vn(2,:,:))
     $     + rinx_rhon(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_rhon(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + qinx_rhon - qnix_rhon)
      s_u(3,8,:,:)=r_fac*((g_i + g_x)*(ken_mnx - vivn_mnx) 
     $     + g_x_mn1*(kei + ken - vivn) + r_inx(1,:,:)*nn_inv
     $     + rinx_mn1(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mn1(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + qinx_mn1 - qnix_mn1)
      s_u(3,9,:,:)=r_fac*((g_i + g_x)*(ken_mny - vivn_mny) 
     $     + g_x_mn2*(kei + ken - vivn) + r_inx(2,:,:)*nn_inv
     $     + rinx_mn2(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mn2(2,:,:)*(vn(2,:,:) - vi(2,:,:)) 
     $     + qinx_mn2 - qnix_mn2)
      s_u(3,10,:,:)=r_fac*(g_x_pn*(kei + ken - vivn) 
     $     + 3._r8*half*g_i*nn_inv
     $     + rinx_pn(1,:,:)*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_pn(2,:,:)*(vn(2,:,:) - vi(2,:,:)) + qinx_pn - qnix_pn)
      
      s_ux(3,1,:,:)=r_fac*(two*visc*(two*vix(1,:,:)*vi_un(1,:,:)  
     $     + vi_un(2,:,:)*(vix(2,:,:) + viy(1,:,:)))
     $     + mu_sv*6.*vix(1,:,:)*vi_un(1,:,:)*ABS(vix(1,:,:)))
      s_ux(3,3,:,:)=r_fac*vi(1,:,:)
      s_ux(3,4,:,:)=r_fac*(visc*4._r8*vix(1,:,:)*n_inv
     $     + mu_sv*6.*vix(1,:,:)*n_inv*ABS(vix(1,:,:)))
      s_ux(3,5,:,:)=r_fac*(visc*two*n_inv*(vix(2,:,:) + viy(1,:,:)))
      s_ux(3,6,:,:)=r_fac*eta_hr*two*ux(6,:,:)

      s_uy(3,1,:,:)=r_fac*(two*visc*(two*viy(2,:,:)*vi_un(2,:,:)  
     $     + vi_un(1,:,:)*(viy(1,:,:) + vix(2,:,:)))
     $     + mu_sv*6.*viy(2,:,:)*vi_un(2,:,:)*ABS(viy(2,:,:)))
      s_uy(3,3,:,:)=r_fac*vi(2,:,:)
      s_uy(3,4,:,:)=r_fac*visc*two*n_inv*(viy(1,:,:) + vix(2,:,:))
      s_uy(3,5,:,:)=r_fac*(visc*4._r8*viy(2,:,:)*n_inv
     $     + mu_sv*6.*viy(2,:,:)*n_inv*ABS(viy(2,:,:)))
      s_uy(3,6,:,:)=r_fac*eta_hr*two*uy(6,:,:)
c-----------------------------------------------------------------------
c     momentum equations.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(1,:,:) 
     $     - two*visc*vix_un(1,:,:) - two*visc_rho*vix(1,:,:)
     $     - 4.*mu_sv*vix_un(1,:,:)*ABS(vix(1,:,:)))
      fx_u(4,3,:,:)=r_fac*(one - two*visc_p*vix(1,:,:))
      fx_u(4,4,:,:)=r_fac*two*(vi(1,:,:) - visc*nx_inv
     $     - two*mu_sv*nx_inv*ABS(vix(1,:,:)))
      fx_ux(4,1,:,:)=-r_fac*two*(visc*vi_un(1,:,:)
     $     + mu_sv*two*vi_un(1,:,:)*ABS(vix(1,:,:)))
      fx_ux(4,4,:,:)=-r_fac*two*(visc*n_inv
     $     + mu_sv*two*n_inv*ABS(vix(1,:,:)))

      fy_u(4,1,:,:)=r_fac*(u(4,:,:)*vi_un(2,:,:) 
     $     - visc*(viy_un(1,:,:) + vix_un(2,:,:))
     $     - visc_rho*(viy(1,:,:) + vix(2,:,:)))
      fy_u(4,3,:,:)=-r_fac*visc_p*(viy(1,:,:) + vix(2,:,:))
      fy_u(4,4,:,:)=r_fac*(vi(2,:,:) - visc*ny_inv)
      fy_u(4,5,:,:)=r_fac*(vi(1,:,:) - visc*nx_inv)

      fy_ux(4,1,:,:)=-r_fac*visc*vi_un(2,:,:)
      fy_ux(4,5,:,:)=-r_fac*visc*n_inv
      fy_uy(4,1,:,:)=-r_fac*visc*vi_un(1,:,:)
      fy_uy(4,4,:,:)=-r_fac*visc*n_inv

      s_u(4,1,:,:)=r_fac*(-g_r*vi_un(1,:,:) 
     $     - g_r_rho*vi(1,:,:)
     $     + g_i_rho*vn(1,:,:) - g_x*vi_un(1,:,:) 
     $     + g_x_rho*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_rho(1,:,:) - rnix_rho(1,:,:))
      s_u(4,3,:,:)=r_fac*(-g_r_p*vi(1,:,:) + g_i_p*vn(1,:,:)
     $     + g_x_p*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_p(1,:,:) - rnix_p(1,:,:))
      s_u(4,4,:,:)=-r_fac*(n_inv*(g_r + g_x)
     $     + g_x_mi1*(vi(1,:,:) - vn(1,:,:))
     $     - rinx_mi1(1,:,:) + rnix_mi1(1,:,:))
      s_u(4,5,:,:)=r_fac*(g_x_mi2*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mi2(1,:,:) - rnix_mi2(1,:,:))
      s_u(4,6,:,:)=-r_fac*b2
      s_u(4,7,:,:)=r_fac*(g_i*vn_un(1,:,:) + g_i_rhon*vn(1,:,:)
     $     + g_x*vn_un(1,:,:) + g_x_rhon*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_rhon(1,:,:) - rnix_rhon(1,:,:))
      s_u(4,8,:,:)=r_fac*(nn_inv*(g_i + g_x)
     $     + g_x_mn1*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mn1(1,:,:) - rnix_mn1(1,:,:))
      s_u(4,9,:,:)=r_fac*(g_x_mn2*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_mn2(1,:,:) - rnix_mn2(1,:,:))
      s_u(4,10,:,:)=r_fac*(g_x_pn*(vn(1,:,:) - vi(1,:,:))
     $     + rinx_pn(1,:,:) - rnix_pn(1,:,:))

      s_ux(4,2,:,:)=-r_fac*u(6,:,:)

      fx_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(1,:,:) 
     $     - visc*(vix_un(2,:,:) + viy_un(1,:,:))
     $     - visc_rho*(vix(2,:,:) + viy(1,:,:)))
      fx_u(5,3,:,:)=-r_fac*visc_p*(vix(2,:,:) + viy(1,:,:))
      fx_u(5,4,:,:)=r_fac*(vi(2,:,:) - visc*ny_inv)
      fx_u(5,5,:,:)=r_fac*(vi(1,:,:) - visc*nx_inv)

      fx_ux(5,1,:,:)=-r_fac*visc*vi_un(2,:,:)
      fx_ux(5,5,:,:)=-r_fac*visc*n_inv
      fx_uy(5,1,:,:)=-r_fac*visc*vi_un(1,:,:)
      fx_uy(5,4,:,:)=-r_fac*visc*n_inv

      fy_u(5,1,:,:)=r_fac*(u(5,:,:)*vi_un(2,:,:) 
     $     - two*visc*viy_un(2,:,:) - two*visc_rho*viy(2,:,:)
     $     - 4.*mu_sv*viy_un(2,:,:)*ABS(viy(2,:,:)))
      fy_u(5,3,:,:)=-r_fac*two*visc_p*viy(2,:,:)
      fy_u(5,5,:,:)=r_fac*two*(vi(2,:,:) - visc*ny_inv
     $     - two*mu_sv*ny_inv*ABS(viy(2,:,:)))

      fy_uy(5,1,:,:)=-r_fac*two*(visc*vi_un(2,:,:)
     $     + two*mu_sv*vi_un(2,:,:)*ABS(viy(2,:,:)))
      fy_uy(5,5,:,:)=-r_fac*two*(visc*n_inv
     $     + two*mu_sv*n_inv*ABS(viy(2,:,:)))

      s_u(5,1,:,:)=-cyl_fac*two*r_faci
     $     *(visc*vi_un(2,:,:) + visc_rho*vi(2,:,:)
     $     + mu_sv*two*r_faci*vi_un(2,:,:)*ABS(vi(2,:,:)))
     $     + r_fac*(-g_r*vi_un(2,:,:) - g_r_rho*vi(2,:,:)
     $     + g_i_rho*vn(2,:,:) - g_x*vi_un(2,:,:) 
     $     + g_x_rho*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_rho(2,:,:) - rnix_rho(2,:,:))
      s_u(5,2,:,:)=-cyl_fac*u(6,:,:)
      s_u(5,3,:,:)=r_fac*(-g_r_p*vi(2,:,:) + g_i_p*vn(2,:,:)
     $     + g_x_p*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_p(2,:,:) - rnix_p(2,:,:))
     $     - cyl_fac*two*visc_p*vi(2,:,:)*r_faci
      s_u(5,4,:,:)=r_fac*(g_x_mi1*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_mi1(2,:,:) - rnix_mi1(2,:,:))
      s_u(5,5,:,:)=-n_inv*(cyl_fac*two*r_faci*(visc
     $     + mu_sv*r_faci*two*ABS(vi(2,:,:)))
     $     + r_fac*(g_r + g_x))
     $     + r_fac*(g_x_mi2*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_mi2(2,:,:) - rnix_mi2(2,:,:))
      s_u(5,6,:,:)=r_fac*(b1+bx)
      s_u(5,7,:,:)=r_fac*(g_i*vn_un(2,:,:) + g_i_rhon*vn(2,:,:)
     $     + g_x*vn_un(2,:,:) + g_x_rhon*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_rhon(2,:,:) - rnix_rhon(2,:,:))
      s_u(5,8,:,:)=r_fac*(g_x_mn1*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_mn1(2,:,:) - rnix_mn1(2,:,:))
      s_u(5,9,:,:)=r_fac*(nn_inv*(g_i + g_x)
     $     + g_x_mn2*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_mn2(2,:,:) - rnix_mn2(2,:,:))
      s_u(5,10,:,:)=r_fac*(g_x_pn*(vn(2,:,:) - vi(2,:,:))
     $     + rinx_pn(2,:,:) - rnix_pn(2,:,:))
      s_uy(5,2,:,:)=-r_fac*u(6,:,:)
      s_uy(5,3,:,:)=-r_fac
c-----------------------------------------------------------------------
c     current.
c-----------------------------------------------------------------------
      fx_ux(6,2,:,:)=one
      fy_u(6,2,:,:)=cyl_fac*r_faci
      fy_uy(6,2,:,:)=one
      s_u(6,6,:,:)=one
c-----------------------------------------------------------------------
c     neutral density.
c-----------------------------------------------------------------------
      fx_u(7,8,:,:)=r_fac
      fy_u(7,9,:,:)=r_fac
      fx_ux(7,7,:,:)=-r_fac*ddiff
      fy_uy(7,7,:,:)=-r_fac*ddiff
      
      s_u(7,1,:,:)=r_fac*(-g_i_rho+g_r_rho)
      s_u(7,3,:,:)=r_fac*(-g_i_p+g_r_p)
      s_u(7,7,:,:)=-r_fac*g_i_rhon
c-----------------------------------------------------------------------
c     neutral momentum.
c-----------------------------------------------------------------------
      fx_u(8,7,:,:)=r_fac*(u(8,:,:)*vn_un(1,:,:)
     $     - two*(viscn*vnx_un(1,:,:) + viscn_rhon*vnx(1,:,:))
     $     - 4.*mun_sv*vnx_un(1,:,:)*ABS(vnx(1,:,:)))
      fx_u(8,8,:,:)=r_fac*two*(vn(1,:,:) - viscn*nnx_inv
     $     - two*mun_sv*nnx_inv*ABS(vnx(1,:,:)))
      fx_u(8,10,:,:)=r_fac*(one - two*viscn_pn*vnx(1,:,:))
      fx_ux(8,7,:,:)=-r_fac*two*(viscn*vn_un(1,:,:)
     $     + mun_sv*two*vn_un(1,:,:)*ABS(vnx(1,:,:)))
      fx_ux(8,8,:,:)=-r_fac*two*(viscn*nn_inv
     $     + mun_sv*two*nn_inv*ABS(vnx(1,:,:)))

      fy_u(8,7,:,:)=r_fac*(u(8,:,:)*vn_un(2,:,:) 
     $     - viscn*(vny_un(1,:,:) + vnx_un(2,:,:))
     $     - viscn_rhon*(vny(1,:,:) + vnx(2,:,:)))
      fy_u(8,8,:,:)=r_fac*(vn(2,:,:) - viscn*nny_inv)
      fy_u(8,9,:,:)=r_fac*(vn(1,:,:) - viscn*nnx_inv)
      fy_u(8,10,:,:) =-r_fac*viscn_pn*(vny(1,:,:) + vnx(2,:,:))
      fy_ux(8,7,:,:)=-r_fac*viscn*vn_un(2,:,:)
      fy_ux(8,9,:,:)=-r_fac*viscn*nn_inv
      fy_uy(8,7,:,:)=-r_fac*viscn*vn_un(1,:,:)
      fy_uy(8,8,:,:)=-r_fac*viscn*nn_inv

      s_u(8,1,:,:)=r_fac*(-g_i_rho*vn(1,:,:)
     $     + g_r*vi_un(1,:,:) + g_r_rho*vi(1,:,:)
     $     + g_x*vi_un(1,:,:) + g_x_rho*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_rho(1,:,:) - rinx_rho(1,:,:))
      s_u(8,3,:,:)=r_fac*(-g_i_p*vn(1,:,:) + g_r_p*vi(1,:,:)
     $     + g_x_p*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_p(1,:,:) - rinx_p(1,:,:))
      s_u(8,4,:,:)=r_fac*(n_inv*(g_r + g_x)
     $     + g_x_mi1*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mi1(1,:,:) - rinx_mi1(1,:,:))
      s_u(8,5,:,:)=r_fac*(g_x_mi2*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mi2(1,:,:) - rinx_mi2(1,:,:))
      s_u(8,7,:,:)=r_fac*(-g_i*vn_un(1,:,:) - g_i_rhon*vn(1,:,:)
     $     - g_x*vn_un(1,:,:) + g_x_rhon*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_rhon(1,:,:) - rinx_rhon(1,:,:))
      s_u(8,8,:,:)=-r_fac*(nn_inv*(g_i + g_x)
     $     + g_x_mn1*(vn(1,:,:) - vi(1,:,:))
     $     - rnix_mn1(1,:,:) + rinx_mn1(1,:,:))
      s_u(8,9,:,:)=r_fac*(g_x_mn2*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mn2(1,:,:) - rinx_mn2(1,:,:))
      s_u(8,10,:,:)=r_fac*(g_x_pn*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_pn(1,:,:) - rinx_pn(1,:,:))

      fx_u(9,7,:,:)=r_fac*(u(9,:,:)*vn_un(1,:,:) 
     $     - viscn*(vnx_un(2,:,:) + vny_un(1,:,:))
     $     - viscn_rhon*(vnx(2,:,:) + vny(1,:,:)))
      fx_u(9,8,:,:)=r_fac*(vn(2,:,:) - viscn*nny_inv)
      fx_u(9,9,:,:)=r_fac*(vn(1,:,:) - viscn*nnx_inv)
      fx_u(9,10,:,:)=-r_fac*viscn_pn*(vnx(2,:,:) + vny(1,:,:))
      fx_ux(9,7,:,:)=-r_fac*viscn*vn_un(2,:,:)
      fx_ux(9,9,:,:)=-r_fac*viscn*nn_inv
      fx_uy(9,7,:,:)=-r_fac*viscn*vn_un(1,:,:)
      fx_uy(9,8,:,:)=-r_fac*viscn*nn_inv
      
      fy_u(9,7,:,:)=r_fac*(u(9,:,:)*vn_un(2,:,:) 
     $     - two*(viscn*vny_un(2,:,:) + viscn_rhon*vny(2,:,:))
     $     - 4.*mun_sv*vny_un(2,:,:)*ABS(vny(2,:,:)))
      fy_u(9,9,:,:)=r_fac*two*(vn(2,:,:) - viscn*nny_inv
     $     - two*mun_sv*nny_inv*ABS(vny(2,:,:)))
      fy_u(9,10,:,:)=-r_fac*two*viscn_pn*vny(2,:,:)
      fy_uy(9,7,:,:)=-r_fac*two*(viscn*vn_un(2,:,:)
     $     + two*mun_sv*vn_un(2,:,:)*ABS(vny(2,:,:)))
      fy_uy(9,9,:,:)=-r_fac*two*(viscn*nn_inv
     $     + two*mun_sv*nn_inv*ABS(vny(2,:,:)))

      s_u(9,1,:,:)=r_fac*(-g_i_rho*vn(2,:,:)
     $     + g_r*vi_un(2,:,:) + g_r_rho*vi(2,:,:)
     $     + g_x*vi_un(2,:,:) + g_x_rho*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_rho(2,:,:) - rinx_rho(2,:,:))
      s_u(9,3,:,:)=r_fac*(-g_i_p*vn(2,:,:) + g_r_p*vi(2,:,:)
     $     + g_x_p*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_p(2,:,:) - rinx_p(2,:,:))
      s_u(9,4,:,:)=r_fac*(g_x_mi1*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mi1(2,:,:) - rinx_mi1(2,:,:))
      s_u(9,5,:,:)=r_fac*(n_inv*(g_r + g_x)
     $     + g_x_mi2*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mi2(2,:,:) - rinx_mi2(2,:,:))
      s_u(9,7,:,:)=r_fac*(-g_i*vn_un(2,:,:) - g_i_rhon*vn(2,:,:)
     $     - g_x*vn_un(2,:,:) + g_x_rhon*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_rhon(2,:,:) - rinx_rhon(2,:,:))
     $     - cyl_fac*two*r_faci*(viscn*vn_un(2,:,:) 
     $     + viscn_rhon*vn(2,:,:)
     $     + mun_sv*two*r_faci*vn_un(2,:,:)*ABS(vn(2,:,:)))
      s_u(9,8,:,:)=r_fac*(g_x_mn1*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_mn1(2,:,:) - rinx_mn1(2,:,:))
      s_u(9,9,:,:)= -nn_inv*cyl_fac*two*r_faci*(viscn
     $     + mun_sv*r_faci*two*ABS(vn(2,:,:)))
     $     - r_fac*(nn_inv*(g_i + g_x)
     $     + g_x_mn2*(vn(2,:,:) - vi(2,:,:))
     $     - rnix_mn2(2,:,:) + rinx_mn2(2,:,:))
      s_u(9,10,:,:)=r_fac*(g_x_pn*(vi(2,:,:) - vn(2,:,:))
     $     + rnix_pn(2,:,:) - rinx_pn(2,:,:))
     $     - cyl_fac*two*viscn_pn*vn(2,:,:)*r_faci
      s_uy(9,10,:,:)=-r_fac
c-----------------------------------------------------------------------
c     neutral pressure.
c-----------------------------------------------------------------------
      fx_u(10,7,:,:)=r_fac*(gamma_fac*u(10,:,:)*vn_un(1,:,:)
     $     - kn_local*Tnx_un - kn_rhon*Tnx)
      fx_u(10,8,:,:)=r_fac*gamma_fac*u(10,:,:)*nn_inv
      fx_u(10,10,:,:)=r_fac*(gamma_fac*vn(1,:,:)
     $     - kn_local*nnx_inv - kn_pn*Tnx)
      fx_ux(10,7,:,:)=-r_fac*kn_local*Tn_un
      fx_ux(10,10,:,:)=-r_fac*kn_local*nn_inv

      fy_u(10,7,:,:)=r_fac*(gamma_fac*u(10,:,:)*vn_un(2,:,:)
     $     - kn_local*Tny_un - kn_rhon*Tny)
      fy_u(10,9,:,:)=r_fac*gamma_fac*u(10,:,:)*nn_inv
      fy_u(10,10,:,:)=r_fac*(gamma_fac*vn(2,:,:)
     $     - kn_local*nny_inv - kn_pn*Tny)
      fy_uy(10,7,:,:)=-r_fac*kn_local*Tn_un
      fy_uy(10,10,:,:)=-r_fac*kn_local*nn_inv

      s_u(10,1,:,:)=r_fac*((g_r + g_x)*(kei_rho - vivn_rho) 
     $     + (g_r_rho + g_x_rho)*(ken + kei - vivn)
     $     + 3._r8*half*(-g_i_rho*u(10,:,:)*nn_inv 
     $     + ti_frac*(g_r_rho*u(3,:,:)*n_inv - g_r*u(3,:,:)*n_inv**2))
     $     - n_inv*(r_nix(1,:,:)*vi(1,:,:) + r_nix(2,:,:)*vi(2,:,:))
     $     + rnix_rho(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_rho(2,:,:)*(vi(2,:,:) - vn(2,:,:)) 
     $     + qnix_rho - qinx_rho)
      s_u(10,3,:,:)=r_fac*((g_r_p + g_x_p)*(ken + kei - vivn)
     $     + 3._r8*half*(-g_i_p*u(10,:,:)*nn_inv 
     $     + ti_frac*(g_r_p*u(3,:,:)*n_inv + g_r*n_inv)) 
     $     + qnix_p - qinx_p
     $     + rnix_p(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_p(2,:,:)*(vi(2,:,:) - vn(2,:,:)))
      s_u(10,4,:,:)=r_fac*((g_r + g_x)*(kei_mix - vivn_mix)
     $     + g_x_mi1*(ken + kei - vivn) + r_nix(1,:,:)*n_inv
     $     + rnix_mi1(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mi1(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + qnix_mi1 - qinx_mi1)
      s_u(10,5,:,:)=r_fac*((g_r + g_x)*(kei_miy - vivn_miy)
     $     + g_x_mi2*(ken + kei - vivn) + r_nix(2,:,:)*n_inv
     $     + rnix_mi2(1,:,:)*(vi(1,:,:) - vn(1,:,:)) 
     $     + rnix_mi2(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + qnix_mi2 - qinx_mi2)
      s_u(10,7,:,:)=r_fac*(vn_un(1,:,:)*ux(10,:,:)
     $     + vn_un(2,:,:)*uy(10,:,:)
     $     + two*viscn*(two*vnx(1,:,:)*vnx_un(1,:,:) 
     $     + two*vny(2,:,:)*vny_un(2,:,:) 
     $     + vny(1,:,:)*vny_un(1,:,:) + vnx(2,:,:)*vnx_un(2,:,:)
     $     + vny(1,:,:)*vnx_un(2,:,:) + vnx(2,:,:)*vny_un(1,:,:))
     $     + viscn_rhon*(two*vnx(1,:,:)**2 + two*vny(2,:,:)**2 
     $     + vny(1,:,:)**2 + vnx(2,:,:)**2 + two*vny(1,:,:)*vnx(2,:,:))
     $     + cyl_fac*two*r_faci**2*(viscn*two*vn(2,:,:)*vn_un(2,:,:)
     $     + viscn_rhon*vn(2,:,:)**2)
     $     + mun_sv*(6.*vnx(1,:,:)*vnx_un(1,:,:)*ABS(vnx(1,:,:))
     $     + 6.*vny(2,:,:)*vny_un(2,:,:)*ABS(vny(2,:,:)))
     $     + cyl_fac*mun_sv*6.*r_faci**3*vn(2,:,:)*vn_un(2,:,:)
     $     *ABS(vn(2,:,:))
     $     + (g_r + g_x)*(ken_rhon - vivn_rhon)
     $     + g_x_rhon*(kei + ken - vivn)
     $     + 3._r8*half*nn_inv*(-g_i_rhon*u(10,:,:) 
     $     + g_i*u(10,:,:)*nn_inv)
     $     + nn_inv*(r_nix(1,:,:)*vn(1,:,:) + r_nix(2,:,:)*vn(2,:,:))
     $     + rnix_rhon(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_rhon(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + qnix_rhon - qinx_rhon)
      s_u(10,8,:,:)=r_fac*(nn_inv*ux(10,:,:) + two*viscn*
     $     (two*vnx(1,:,:)*nnx_inv + (vny(1,:,:) + vnx(2,:,:))*nny_inv)
     $     + mun_sv*6.*vnx(1,:,:)*nnx_inv*ABS(vnx(1,:,:))
     $     + (g_r + g_x)*(ken_mnx - vivn_mnx)
     $     + g_x_mn1*(ken + kei - vivn) - r_nix(1,:,:)*nn_inv
     $     + rnix_mn1(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mn1(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + qnix_mn1 - qinx_mn1)
      s_u(10,9,:,:)=r_fac*(uy(10,:,:)*nn_inv + two*viscn*
     $     (two*vny(2,:,:)*nny_inv + (vnx(2,:,:) + vny(1,:,:))*nnx_inv)
     $     + cyl_fac*4._r8*viscn*r_faci**2*vn(2,:,:)*nn_inv
     $     + mun_sv*6.*vny(2,:,:)*nny_inv*ABS(vny(2,:,:))
     $     + cyl_fac*6.*mun_sv*r_faci
     $     *SIGN((vn(2,:,:)*r_faci)**2*nn_inv,vn(2,:,:))
     $     + (g_r + g_x)*(ken_mny - vivn_mny)
     $     + g_x_mn2*(ken + kei - vivn) - r_nix(2,:,:)*nn_inv
     $     + rnix_mn2(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_mn2(2,:,:)*(vi(2,:,:) - vn(2,:,:))
     $     + qnix_mn2 - qinx_mn2)
      s_u(10,10,:,:)=r_fac*(viscn_pn*(two*vnx(1,:,:)**2 
     $     + two*vny(2,:,:)**2 + vny(1,:,:)**2 + vnx(2,:,:)**2 
     $     + two*vny(1,:,:)*vnx(2,:,:))
     $     + cyl_fac*two*viscn_pn*(vn(2,:,:)*r_faci)**2
     $     + g_x_pn*(kei + ken - vivn)
     $     - 3._r8*half*g_i*nn_inv + qnix_pn - qinx_pn
     $     + rnix_pn(1,:,:)*(vi(1,:,:) - vn(1,:,:))
     $     + rnix_pn(2,:,:)*(vi(2,:,:) - vn(2,:,:)))

      s_ux(10,7,:,:)=r_fac*(two*viscn*(two*vnx(1,:,:)*vn_un(1,:,:)  
     $     + vn_un(2,:,:)*(vnx(2,:,:) + vny(1,:,:)))
     $     + mun_sv*6.*vnx(1,:,:)*vn_un(1,:,:)*ABS(vnx(1,:,:)))
      s_ux(10,8,:,:)=r_fac*(viscn*4._r8*vnx(1,:,:)*nn_inv
     $     + mun_sv*6.*vnx(1,:,:)*nn_inv*ABS(vnx(1,:,:)))
      s_ux(10,9,:,:)=r_fac*viscn*two*nn_inv*(vnx(2,:,:) + vny(1,:,:))
      s_ux(10,10,:,:)=r_fac*vn(1,:,:)

      s_uy(10,7,:,:)=r_fac*(two*viscn*(two*vny(2,:,:)*vn_un(2,:,:)  
     $     + vn_un(1,:,:)*(vny(1,:,:) + vnx(2,:,:)))
     $     + mun_sv*6.*vny(2,:,:)*vn_un(2,:,:)*ABS(vny(2,:,:)))
      s_uy(10,8,:,:)=r_fac*viscn*two*nn_inv*(vny(1,:,:) + vnx(2,:,:))
      s_uy(10,9,:,:)=r_fac*(viscn*4._r8*vny(2,:,:)*nn_inv
     $     + mun_sv*6.*vny(2,:,:)*nn_inv*ABS(vny(2,:,:)))
      s_uy(10,10,:,:)=r_fac*vn(2,:,:)

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
      USE pn_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: r_fac
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     cylindrical to cartesian relationships:
c     1: z --> x
c     2: r --> y
c     3: phi --> z
c-----------------------------------------------------------------------
      r_fac=y**cyl_fac

      mass(1,1,:,:)=r_fac
      mass(3,3,:,:)=r_fac/(gamma-one)
      mass(4,4,:,:)=r_fac
      mass(5,5,:,:)=r_fac

      mass(7,7,:,:)=r_fac
      mass(8,8,:,:)=r_fac
      mass(9,9,:,:)=r_fac
      mass(10,10,:,:)=r_fac/(gamma-one)
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
      SUBROUTINE physics_grid(ksi,etag,x,y)
      USE pn_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: ksi,etag
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: x,y
c-----------------------------------------------------------------------
c     set grid according to init_type.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("trans_test","uniform","uniform2","civ")
         x=lx*ksi+xmin
         IF(gr_curve == 0.)THEN
            y=ly*etag
         ELSE
            y=ly*(one + (etag-one)*(one+gr_curve-etag)/(one+gr_curve))
         ENDIF
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
      USE pn_mod
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
      USE pn_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     deallocate appropriate arrays.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("trans_test")
         CALL extra_equil_dealloc(interp,equil_bc,equil,equilxy)
         IF(ifcoils) CALL extra_coildealloc(coils)
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
