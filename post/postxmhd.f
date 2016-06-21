c-----------------------------------------------------------------------
c     file postxmhd.f.
c     post-processes output from sel code for xmhd problems.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. postxmhd_mod.
c     1. postxmhd_init.
c     2. postxmhd_geteta.
c     3. postxmhd_ploteta.
c     4. postxmhd_getkappa.
c     5. postxmhd_plotkappa.
c     6. postxmhd_screw.
c-----------------------------------------------------------------------
c     subprogram 0. postxmhd_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE postxmhd_mod
      USE plotter_mod
      USE transport_mod
      IMPLICIT NONE

      REAL(r8), PARAMETER, PRIVATE :: gamma=5._r8/3._r8,
     $     me=9.109e-31,ep0=8.854e-12,mp=1.673e-27,qe=1.602e-19

      LOGICAL, PRIVATE :: flag_eta,flag_kappa
      CHARACTER(20), PRIVATE :: init_type,eta_case,kappa_case
      INTEGER, PRIVATE :: neta,nkappa
      REAL(r8), PRIVATE :: b0,L0,n0,eta,r_eta=0.,etavac=1.,mi=0.,
     $     te_frac=.5,chod_const=0.1,kappa_par=0.,kappa_perp=0.,
     $     kappa_min=0.,kappa_max=1.e8

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. postxmhd_init.
c     initialize post-processing parameters from post.in 
c     and input parameters from sel.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postxmhd_init(indir,norm_plist,eta_name,eta_plist,
     $     kappa_name,kappa_plist)

      CHARACTER(*), INTENT(IN) :: indir,eta_name,kappa_name
      REAL(r8), DIMENSION(:), INTENT(IN) :: norm_plist,eta_plist,
     $     kappa_plist

      NAMELIST/xmhd_input/flag_eta,flag_kappa
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      CALL system("rm -f eta.bin kappa.bin")
      OPEN(UNIT=in_unit,FILE="post.in",STATUS="OLD")
      READ(in_unit,NML=xmhd_input)
      CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     initialize parameters.
c-----------------------------------------------------------------------
      b0 = norm_plist(1)
      L0 = norm_plist(2)
      n0 = norm_plist(3)
      te_frac = norm_plist(4)
      mi=mp*norm_plist(5)

      eta_case = eta_name
      eta = eta_plist(1)
      chod_const = eta_plist(2)
      etavac = eta_plist(3)
      r_eta = eta_plist(4)

      kappa_case = kappa_name
      kappa_par = kappa_plist(1)
      kappa_perp = kappa_plist(2)
      kappa_min = kappa_plist(3)
      kappa_max = kappa_plist(4)

      neta=6
      nkappa=9

      IF(flag2d .AND. flag_eta .AND. out_type/="hdf5")THEN
         OPEN(UNIT=Ucontour_unit,FILE="eta.bin",
     $        STATUS="UNKNOWN",FORM="UNFORMATTED")
         WRITE(Ucontour_unit)1,0,neta
         CLOSE(UNIT=Ucontour_unit)
      ENDIF
      IF(flag2d .AND. flag_kappa .AND. out_type/="hdf5")THEN
         OPEN(UNIT=Ucontour_unit,FILE="kappa.bin",
     $        STATUS="UNKNOWN",FORM="UNFORMATTED")
         WRITE(Ucontour_unit)1,0,nkappa
         CLOSE(UNIT=Ucontour_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postxmhd_init
c-----------------------------------------------------------------------
c     subprogram 2. postxmhd_geteta.
c     calculates quantities of interest for plotting resistivity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postxmhd_geteta(xyw,rho,p,j,uu)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,p,j
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: uu

      INTEGER :: i
      REAL(r8) :: etas_norm,etac_norm,v_chod_norm,tnorm
c-----------------------------------------------------------------------
c     zero uu.
c-----------------------------------------------------------------------
      IF(.NOT.flag_eta) RETURN
      uu=0._r8
c-----------------------------------------------------------------------
c     compute resistivity quantities:
c     eta(x,y)          -->  uu(1,:,:)
c     eta_spitzer(x,y)  -->  uu(2,:,:)
c     eta_chodura(x,y)  -->  uu(3,:,:)
c     j                 -->  uu(4,:,:)
c     T                 -->  uu(5,:,:)
c     eta*j^2           -->  uu(6,:,:)
c-----------------------------------------------------------------------
      tnorm=L0*SQRT(mu0*n0*mi)/b0
      v_chod_norm=SQRT(mi/(mu0*n0))/(qe*L0)
      etac_norm=me/(qe*L0*b0*SQRT(ep0*mu0))
      etas_norm=5.e-5*17.*tnorm*(two*n0*mu0*qe)**1.5/(mu0*L0**2*b0**3)
      IF(eta_case == "spitzer") etac_norm = zero
      IF(eta_case == "chodura") etas_norm = zero
      IF(te_frac > 0)THEN
         etas_norm = etas_norm/te_frac**1.5
      ELSE
         etas_norm=0.
      ENDIF

      uu(1,:,:)=eta
      uu(4,:,:)=j
      uu(5,:,:)=p/rho
      WHERE(p <=0. .OR. rho <= 0.)
         uu(5,:,:)=0.
      END WHERE

      SELECT CASE(eta_case)
      CASE("spitzer-chodura","spitzer","chodura")
c         CALL transport_seteta("spitzer-chodura",xyw(2,:,:),rho,p,j,
c     $        chod_const,etas_norm,etac_norm,v_chod_norm,r_eta,eta,
c     $        etavac,uu(1,:,:))
c         CALL transport_seteta("spitzer-chodura",xyw(2,:,:),rho,p,j,
c     $        chod_const,etas_norm,0._r8,v_chod_norm,r_eta,eta,etavac,
c     $        uu(2,:,:))
c         CALL transport_seteta("spitzer-chodura",xyw(2,:,:),rho,p,j,
c     $        chod_const,0._r8,etac_norm,v_chod_norm,r_eta,eta,etavac,
c     $        uu(3,:,:))
c      CASE("x-layer")
c         CALL transport_seteta(eta_case,xyw(2,:,:),rho,p,j,0._r8,0._r8,
c     $        0._r8,0._r8,r_eta,eta,etavac,uu(1,:,:))
      END SELECT

      uu(6,:,:)=uu(1,:,:)*j**2
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
      DO i=1,neta
         IF (MAXVAL(ABS(uu(i,:,:))) <= 1e-20)uu(i,:,:)=0
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postxmhd_geteta
c-----------------------------------------------------------------------
c     subprogram 3. postxmhd_ploteta.
c     generates contour plots of resistivity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postxmhd_ploteta(t,xyw,rho,p,j,header,footer,
     $     nt_next,ifile,stride,filename)

      LOGICAL, INTENT(IN) :: header,footer
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nt_next,ifile,stride
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,p,j
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw

      INTEGER, DIMENSION(3), SAVE :: xcont=0
      REAL(r8), DIMENSION(neta,SIZE(xyw,2),SIZE(xyw,3)) :: utemp
c-----------------------------------------------------------------------
c     if desired, compute and plot resistivity
c-----------------------------------------------------------------------
      IF(.NOT.flag_eta)RETURN

      CALL postxmhd_geteta(xyw,rho,p,j,utemp)

      CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $     nt_next,xcont,ifile,stride,filename,"eta")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postxmhd_ploteta
c-----------------------------------------------------------------------
c     subprogram 4. postxmhd_getkappa.
c     calculates quantities of interest for plotting thermal 
c     conductivity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postxmhd_getkappa(rho,p,bsq,uu)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,p,bsq
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: uu

      INTEGER :: i
      REAL(r8) :: tnorm,ke_norm,ki_norm,xe_norm,xi_norm
c-----------------------------------------------------------------------
c     zero uu.
c-----------------------------------------------------------------------
      IF(.NOT.flag_kappa)RETURN
      uu=0._r8
c-----------------------------------------------------------------------
c     compute thermal conductivity quantities:
c     Bsq, InvBsq       -->  uu(1:2,:,:)
c     T                 -->  uu(3,:,:)
c     kpare             -->  uu(4,:,:)
c     kpari             -->  uu(5,:,:)
c     kperpe            -->  uu(6,:,:)
c     kperpi            -->  uu(7,:,:)
c     kpar              -->  uu(8,:,:)
c     kperp             -->  uu(9,:,:)
c-----------------------------------------------------------------------
      uu(1,:,:) = bsq
      uu(2,:,:) = 1._r8/bsq
      WHERE(bsq == 0)
         uu(2,:,:)=HUGE(uu(2,:,:))
      END WHERE
      uu(3,:,:) = p/rho
      WHERE(p <= 0. .OR. rho <= 0.)
         uu(3,:,:) = 0.
      END WHERE
c-----------------------------------------------------------------------
c     calculate thermal conductivity.
c-----------------------------------------------------------------------
      SELECT CASE(kappa_case)
      CASE("braginskii")
         tnorm=L0*SQRT(mu0*n0*mi)/b0
         ke_norm = 3.56e21*(b0**2/(2*n0*mu0*qe))**2.5*tnorm/(2*n0*L0**2)
         xe_norm = 3.56e21*(b0**2/(2*n0*mu0*qe))**1.5*b0/n0
         ki_norm = 1.97e-7/SQRT(mi*mp)*(b0**2/(2*n0*mu0*qe))**2.5
     $        *tnorm/(2*n0*L0**2)
         xi_norm = 1.97e-7/SQRT(mi*mp)*(b0**2/(2*n0*mu0*qe))**1.5*b0/n0

         CALL transport_kbrag(rho,p,1._r8,bsq,ke_norm,0._r8,
     $        xe_norm,0._r8,kappa_min,kappa_max,uu(6,:,:),uu(4,:,:))
c-----------------------------------------------------------------------
c        kpare = (kpare - kperpe) + kperpe.
c-----------------------------------------------------------------------
         uu(4,:,:) = uu(4,:,:) + uu(6,:,:)

         CALL transport_kbrag(rho,p,0._r8,bsq,0._r8,ki_norm,
     $        0._r8,xi_norm,kappa_min,kappa_max,uu(7,:,:),uu(5,:,:))
c-----------------------------------------------------------------------
c        kpari = (kpari - kperpi) + kperpi
c-----------------------------------------------------------------------
         uu(5,:,:) = uu(5,:,:) + uu(7,:,:)
c-----------------------------------------------------------------------
c        kpar/perp = kpare+i or kperpe+i.
c-----------------------------------------------------------------------
         uu(8,:,:) = uu(4,:,:) + uu(5,:,:)
         uu(9,:,:) = uu(6,:,:) + uu(7,:,:)
      CASE("anisotropic")
         CALL transport_setkaniso(kappa_par,kappa_perp,bsq,
     $        uu(9,:,:),uu(8,:,:))
         uu(8,:,:) = uu(8,:,:) + uu(9,:,:)
         uu(9,:,:) = uu(9,:,:)
      CASE DEFAULT
         uu(8,:,:) = kappa_par
         uu(9,:,:) = kappa_par
      END SELECT
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
      DO i=1,nkappa
         IF (MAXVAL(ABS(uu(i,:,:))) <= 1e-20)uu(i,:,:)=0
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postxmhd_getkappa
c-----------------------------------------------------------------------
c     subprogram 5. postxmhd_plotkappa.
c     generates contour plots of thermal conductivity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postxmhd_plotkappa(t,xyw,rho,p,bsq,header,footer,
     $     nt_next,ifile,stride,filename)

      LOGICAL, INTENT(IN) :: header,footer
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nt_next,ifile,stride
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: rho,p,bsq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: xyw

      INTEGER, DIMENSION(3), SAVE :: xcont=0
      REAL(r8), DIMENSION(nkappa,SIZE(xyw,2),SIZE(xyw,3)) :: utemp
c-----------------------------------------------------------------------
c     if desired, calculate and plot thermal conductivity.
c-----------------------------------------------------------------------
      IF(.NOT. flag_kappa)RETURN
      CALL postxmhd_getkappa(rho,p,bsq,utemp)
      CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $     nt_next,xcont,ifile,stride,filename,"kappa")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postxmhd_plotkappa
c-----------------------------------------------------------------------
c     subprogram 6. postxmhd_screw.
c     plot various screw pinch quantities vs. time.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE postxmhd_screw(t,xyw,uw,pmin)

      REAL(r8), INTENT(IN) :: t,pmin
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw

      INTEGER :: ix,iy,nx,ny
      REAL(r8), DIMENSION(5) :: a
      REAL(r8), DIMENSION(2) :: value
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rbphi
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uw1
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      nx=SIZE(xyw,2)-1
      ny=SIZE(xyw,3)-1
      SELECT CASE(init_type)
      CASE("screw")
         ALLOCATE(uw1(SIZE(uw,1),0:nx,0:ny),rbphi(0:nx,0:ny))
      CASE DEFAULT
         RETURN
      END SELECT
      uw1 = uw
      value = zero
c-----------------------------------------------------------------------
c     eliminate momentum where p < 10*pmin.
c-----------------------------------------------------------------------
      WHERE(uw(8,:,:) < 10*pmin)
         uw1(4,:,:) = zero
      ENDWHERE
c-----------------------------------------------------------------------
c     find total axial momentum within given pressure contour.
c     SUM(mom*vol)=SUM(mom*dA*2pi*r)
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         DO iy=0,ny-1
            a(1)=(xyw(2,ix,iy)+xyw(2,ix,iy+1))*
     $           (xyw(1,ix,iy+1)-xyw(1,ix,iy))/2._r8
            a(2)=(xyw(2,ix,iy+1)+xyw(2,ix+1,iy+1))*
     $           (xyw(1,ix+1,iy+1)-xyw(1,ix,iy+1))/2._r8
            a(3)=(xyw(2,ix+1,iy+1)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix+1,iy)-xyw(1,ix+1,iy+1))/2._r8
            a(4)=(xyw(2,ix,iy)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix,iy)-xyw(1,ix+1,iy))/2._r8
            a(5)=ABS(SUM(a(1:4)))
            
            value(1)=value(1)+a(5)*twopi*SUM(xyw(2,ix:ix+1,iy:iy+1))*
     $           SUM(uw1(4,ix:ix+1,iy:iy+1))/16.
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     find the axial current in the region where p > 10*pmin.
c-----------------------------------------------------------------------
      rbphi = xyw(2,:,:)*uw(3,:,:)
      WHERE(uw(8,:,:) < 10*pmin)
         rbphi = zero
      ENDWHERE
      value(2)=MAXVAL(rbphi)*twopi

      DEALLOCATE(uw1,rbphi)
c-----------------------------------------------------------------------
c     open, write, and close screw.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      OPEN(UNIT=oned_unit,FILE="screw.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(oned_unit)REAL(t,4),REAL(value,4)
      CLOSE(UNIT=oned_unit)
c-----------------------------------------------------------------------
c     open, write, and close screw.asc file to store ascii data.
c-----------------------------------------------------------------------
      OPEN(UNIT=oned_unit,FILE="screw.asc",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="FORMATTED")
 10   FORMAT(9e15.5)
      WRITE(oned_unit,10)REAL(t,4),REAL(value,4)
      CLOSE(UNIT=oned_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE postxmhd_screw
      END MODULE postxmhd_mod
