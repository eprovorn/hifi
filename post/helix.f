c-----------------------------------------------------------------------
c     file helix.f.
c     post-processes output from sel code for helical symmetry.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. helix_mod.
c     1. helix_read.
c     2. helix_dUdt.
c     4. helix_contour.
c     5. helix_energy.
c     6. helix_maxUvsT.
c     7. helix_UxyT.
c     8. helix_Uprofile.
c     9. helix_q.
c     10. helix_Bspectr.
c     11. helix_swtenergy.
c-----------------------------------------------------------------------
c     subprogram 0. helix_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE helix_mod
      USE plotter_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: source
      CHARACTER(20), PRIVATE :: init_type
      REAL(r8), PRIVATE :: eta,mu,nu,di,r_s,Rinv,epsilon,x1,y1,x2,y2,
     $     t0_decay,t0_init,eta_0,mu_0,nu_0,beta0,beta_e,Tconst,
     $     q0_r,q0_delr,r_m
      REAL(r8), DIMENSION(0:100), PRIVATE :: x0,y0

      INTEGER, PRIVATE :: ncont
      INTEGER, DIMENSION(3), PRIVATE :: xcont=0
      REAL(r8), PRIVATE :: t_old
      REAL(r8), DIMENSION(1), PRIVATE :: u_old
      REAL(r8), DIMENSION(3), PRIVATE :: u2_old
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: uc_old

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. helix_read.
c     read necessary post-processing parameters from  post.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_read(indir)

      CHARACTER(*), INTENT(IN) :: indir

      CHARACTER(80) :: filename
      INTEGER :: itmax,itmax_incr,itmax_decr
      REAL(r8) :: dt_incr,dt_decr

      NAMELIST/helix_input/x1,y1,x2,y2,q0_r,q0_delr,r_m
      NAMELIST/kink_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,di,eta,mu,nu,eta_0,mu_0,nu_0,
     $     t0_init,t0_decay,beta0,beta_e,Tconst,r_s,Rinv,epsilon,
     $     source,init_type
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      CALL system("rm -f helix_cont.bin Bmodes.bin qprof.bin swten.bin")
      OPEN(UNIT=in_unit,FILE="post.in",STATUS="OLD")
      READ(in_unit,NML=helix_input)
      CLOSE(UNIT=in_unit)
      filename=TRIM(indir)//"/sel.in"
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
         OPEN(UNIT=in_unit,FILE=filename,STATUS="OLD")
         READ(in_unit,NML=kink_list)
         CLOSE(UNIT=in_unit)
         ncont=7
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="helix_cont.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
      CASE("XMHDhlx")
         OPEN(UNIT=in_unit,FILE=filename,STATUS="OLD")
         READ(in_unit,NML=kink_list)
         CLOSE(UNIT=in_unit)
         ncont=9
         IF(flag2d .AND. out_type/="hdf5")THEN
            OPEN(UNIT=Ucontour_unit,FILE="helix_cont.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,ncont
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_read
c-----------------------------------------------------------------------
c     subprogram 2. helix_dUdt.
c     generate du/dt(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_dUdt(t,xyw,uw,uxyw,first)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw
      LOGICAL, INTENT(IN) :: first

      INTEGER, DIMENSION(2) :: loc_index
      REAL(r8), DIMENSION(6) :: value
      REAL(r8), DIMENSION(3) :: u2
      REAL(r8), DIMENSION(0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1) :: gfac
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      u2=0
      loc_index=0
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
         gfac=1./(1. + (xyw(1,:,:)*Rinv)**2)
         CALL helix_energy(xyw,uw,uxyw,value)
         IF(value(1) /= 0)THEN
            u2(1)=LOG(ABS(value(1)))
         ELSE
            u2(1)=0
         ENDIF
         u2(2)=ABS(value(1))
         SELECT CASE(init_type)
         CASE("m1_park","m1")
            loc_index=MAXLOC(uw(5,:,:)
     $           - (2*xyw(1,:,:)**2-r_s**2)/gfac - 2*gfac)
         CASE("m1_axial","m1_axial_new")
            loc_index=MAXLOC(uw(5,:,:)
     $           - (2*xyw(1,:,:)**2-r_s**2)/gfac - 2
     $           + (Rinv*xyw(1,:,:))**2*(xyw(1,:,:)**2-r_s**2))
         END SELECT
         u2(3)=uw(1,loc_index(1),loc_index(2))
         CALL plotter_dUdt(t,t_old,u2,u2_old,first)
      CASE("XMHDhlx")
         gfac=1./(1. + (xyw(1,:,:)*Rinv)**2)
         CALL helix_energy(xyw,uw,uxyw,value)
         IF(value(1) /= 0)THEN
            u2(1)=LOG(ABS(value(1)))
         ELSE
            u2(1)=0
         ENDIF
         u2(2)=ABS(value(1))
         SELECT CASE(init_type)
         CASE("m1_axial")
            loc_index=MAXLOC(uw(4,:,:) - uw(7,:,:)*uw(1,:,:) 
     $           - di*SQRT(gfac)
     $           *((1.+ 2.*(Rinv*xyw(1,:,:))**2)*(xyw(1,:,:)**2-r_s**2)
     $           + xyw(1,:,:)**2/gfac - 2))
         END SELECT
         u2(3)=uw(5,loc_index(1),loc_index(2))
         CALL plotter_dUdt(t,t_old,u2,u2_old,first)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_dUdt
c-----------------------------------------------------------------------
c     subprogram 4. helix_contour.
c     generates contours of desired physical quantities
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_contour(t,xyw,uw,uxyw,header,
     $     footer,nt_next,ifile,stride,filename)

      LOGICAL, INTENT(IN) :: header,footer
      ChARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nt_next,ifile,stride
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw

      INTEGER :: ix
      REAL(r8) :: Eb,gfac_s
      REAL(r8), DIMENSION(0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: gfac,
     $     eta_local,Br,Btau,Be,Jr,Jtau,Je,Vr,Vtau,Ve,Jpar,Btot
      REAL(r8), DIMENSION(ncont,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) 
     $     :: utemp
c-----------------------------------------------------------------------
c     prepare auxiliary variables
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
         gfac=SQRT(1./(1. + (xyw(1,:,:)*Rinv)**2))
         gfac_s=1/(1+(r_s*Rinv)**2)
         SELECT CASE(init_type)
         CASE("m1_park")
            IF(.NOT. source)
     $           eta_local=eta*ABS((r_s**2/gfac_s - 2*gfac_s)
     $           /((2*xyw(1,:,:)**2-r_s**2)/gfac**2 - 2*gfac**2))
         CASE("m1_axial")
            Eb=eta*(r_s**2/gfac_s - 2)
            eta_local=ABS(Eb/((1+2*(xyw(1,:,:)*Rinv)**2)
     $           *(xyw(1,:,:)**2-r_s**2) + xyw(1,:,:)**2/gfac**2 - 2))
         CASE("m1_axial_new")
            Eb=((eta-eta_0)*EXP(-((t-t0_init)/t0_decay)**2) + eta_0)
     $           *(r_s**2/gfac_s - 2)
            eta_local=ABS(Eb/((1+2*(xyw(1,:,:)*Rinv)**2)
     $           *(xyw(1,:,:)**2-r_s**2) + xyw(1,:,:)**2/gfac**2 - 2))
         END SELECT
         Br=-uxyw(:,:,2,1)/xyw(1,:,:)
         Jr=uxyw(:,:,2,2)/xyw(1,:,:)
         Vr=-uxyw(:,:,2,6)/xyw(1,:,:)
         WHERE(xyw(1,:,:) == 0)
            Br(:,:)=0
            Jr(:,:)=0
            Vr(:,:)=0
         END WHERE
         Btau=gfac*uxyw(:,:,1,1)
         Jtau=-gfac*uxyw(:,:,1,2)
         Vtau=gfac*uxyw(:,:,1,6)
         Be=gfac*uw(2,:,:)
         Je=gfac*uw(5,:,:)
         Ve=gfac*uw(4,:,:)
c-----------------------------------------------------------------------
c     calculate |JxB|
c-----------------------------------------------------------------------
         utemp(1,:,:)=SQRT((Jtau*Be-Je*Btau)**2+(Je*Br-Jr*Be)**2
     $        +(Jr*Btau-Jtau*Br)**2)
c-----------------------------------------------------------------------
c     calculate Er
c-----------------------------------------------------------------------
         utemp(2,:,:)=eta_local*Jr + Ve*Btau - Vtau*Be
c-----------------------------------------------------------------------
c     calculate Etau
c-----------------------------------------------------------------------
         utemp(3,:,:)=eta_local*Jtau + Vr*Be - Ve*Br
c-----------------------------------------------------------------------
c     calculate Ee
c-----------------------------------------------------------------------
         utemp(4,:,:)=eta_local*Je + Vtau*Br - Vr*Btau
c-----------------------------------------------------------------------
c     calculate EdotB
c-----------------------------------------------------------------------
         utemp(5,:,:)=utemp(2,:,:)*Br + utemp(3,:,:)*Btau
     $        + utemp(4,:,:)*Be
c-----------------------------------------------------------------------
c     calculate Bz
c-----------------------------------------------------------------------
         utemp(6,:,:)=gfac*(Be + Rinv*xyw(1,:,:)*Btau)
c-----------------------------------------------------------------------
c     calculate B^2
c-----------------------------------------------------------------------
         utemp(7,:,:)=Br**2 + Btau**2 + Be**2
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
         DO ix=1,ncont
            IF (MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
         ENDDO
         CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $        nt_next,xcont,ifile,stride,filename,"helix_cont")
      CASE("XMHDhlx")
         gfac=SQRT(1./(1. + (xyw(1,:,:)*Rinv)**2))
         gfac_s=1/(1+(r_s*Rinv)**2)
         eta_local = di*uw(1,:,:)**1.5/(1.96*nu*uw(9,:,:)**1.5)
         Br=-uxyw(:,:,2,5)/xyw(1,:,:)
         Jr=uxyw(:,:,2,6)/xyw(1,:,:)
         Vr=uw(2,:,:)/uw(1,:,:)
         WHERE(xyw(1,:,:) == 0)
            Br(:,:)=0
            Jr(:,:)=0
         END WHERE
         Btau=gfac*uxyw(:,:,1,5)
         Jtau=-gfac*uxyw(:,:,1,6)
         Vtau=uw(3,:,:)/uw(1,:,:)
         Be=gfac*uw(6,:,:)
         Je=(uw(4,:,:)-uw(7,:,:)*uw(1,:,:))/di
         Ve=uw(4,:,:)/uw(1,:,:)
         Btot=SQRT(Br**2 + Btau**2 + Be**2)
         Jpar=(Jr*Br + Jtau*Btau + Je*Be)/Btot
c-----------------------------------------------------------------------
c     calculate |JxB|
c-----------------------------------------------------------------------
         utemp(1,:,:)=SQRT((Jtau*Be-Je*Btau)**2+(Je*Br-Jr*Be)**2
     $        +(Jr*Btau-Jtau*Br)**2)
c-----------------------------------------------------------------------
c     calculate Er
c-----------------------------------------------------------------------
         utemp(2,:,:)=eta_local*(1.96*Jr - .96*Jpar*Br/Btot) 
     $        + Ve*Btau - Vtau*Be
c-----------------------------------------------------------------------
c     calculate Etau
c-----------------------------------------------------------------------
         utemp(3,:,:)=eta_local*(1.96*Jtau - .96*Jpar*Btau/Btot) 
     $        + Vr*Be - Ve*Br
c-----------------------------------------------------------------------
c     calculate Ee
c-----------------------------------------------------------------------
         utemp(4,:,:)=eta_local*(1.96*Je - .96*Jpar*Be/Btot)
     $        + Vtau*Br - Vr*Btau
c-----------------------------------------------------------------------
c     calculate EdotB
c-----------------------------------------------------------------------
         utemp(5,:,:)=utemp(2,:,:)*Br + utemp(3,:,:)*Btau
     $        + utemp(4,:,:)*Be
c-----------------------------------------------------------------------
c     calculate div Vi_pol
c-----------------------------------------------------------------------
         utemp(6,:,:)=uw(2,:,:)/(uw(1,:,:)*xyw(1,:,:))
     $        + (uxyw(:,:,1,2) - uw(2,:,:)*uxyw(:,:,1,1)/uw(1,:,:))
     $        /uw(1,:,:)
     $        + (uxyw(:,:,2,3) - uw(3,:,:)*uxyw(:,:,2,1)/uw(1,:,:))
     $        /(uw(1,:,:)*xyw(1,:,:)*gfac)
c-----------------------------------------------------------------------
c     calculate curl Vi_pol
c-----------------------------------------------------------------------
         utemp(7,:,:)=uw(3,:,:)*gfac**2/(uw(1,:,:)*xyw(1,:,:))
     $        + (uxyw(:,:,1,3) - uw(3,:,:)*uxyw(:,:,1,1)/uw(1,:,:))
     $        /uw(1,:,:)
     $        - (uxyw(:,:,2,2) - uw(2,:,:)*uxyw(:,:,2,1)/uw(1,:,:))
     $        /(uw(1,:,:)*xyw(1,:,:)*gfac)
         WHERE(xyw(1,:,:) == 0)
            utemp(6,:,:)=0
            utemp(7,:,:)=0
         END WHERE
c-----------------------------------------------------------------------
c     calculate Ti
c-----------------------------------------------------------------------
         utemp(8,:,:)=uw(8,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     calculate Te
c-----------------------------------------------------------------------
         utemp(9,:,:)=uw(9,:,:)/uw(1,:,:)
c-----------------------------------------------------------------------
c     call plotter subroutine
c-----------------------------------------------------------------------
         DO ix=1,ncont
            IF (MAXVAL(ABS(utemp(ix,:,:))) <= 1e-12)utemp(ix,:,:)=0
         ENDDO
         CALL plotter_Ucontour(t,xyw,utemp,header,footer,
     $        nt_next,xcont,ifile,stride,filename,"helix_cont")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_contour
c-----------------------------------------------------------------------
c     subprogram 5. helix_energy.
c     generate a time plot for an integral of total energy
c     over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_energy(xyw,uw,uxyw,value)

      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw
      REAL(r8), DIMENSION(6), INTENT(OUT) :: value

      INTEGER :: ix,iy,nx,ny
      REAL(r8), DIMENSION(5) :: a
      REAL(r8), DIMENSION(0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1) :: 
     $     gfac,rad
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uw1,uw2,uw3
c-----------------------------------------------------------------------
c     prepare auxiliary parameters
c-----------------------------------------------------------------------
      value=0
      nx=SIZE(xyw,2)-1
      ny=SIZE(xyw,3)-1
c-----------------------------------------------------------------------
c     calculate energy integral
c-----------------------------------------------------------------------
      ALLOCATE(uw1(2,0:nx,0:ny),uw2(2,0:nx,0:ny),uw3(2,0:nx,0:ny))
      gfac=1/SQRT(1+(xyw(1,:,:)*Rinv)**2)
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
c-----------------------------------------------------------------------
c     kinetic and magnetic energy components.
c-----------------------------------------------------------------------
         uw1(1,:,:)=-uxyw(:,:,2,6)/xyw(1,:,:)
         uw1(2,:,:)=-uxyw(:,:,2,1)/xyw(1,:,:)
         WHERE(xyw(1,:,:)==0)
            uw1(1,:,:)=0
            uw1(2,:,:)=0
         END WHERE
         uw2(1,:,:)=uxyw(:,:,1,6)*gfac
         uw2(2,:,:)=uxyw(:,:,1,1)*gfac
         uw3(1,:,:)=uw(4,:,:)*gfac
         uw3(2,:,:)=uw(2,:,:)*gfac
      CASE("XMHDhlx")
c-----------------------------------------------------------------------
c     kinetic and magnetic energy components.
c-----------------------------------------------------------------------
         uw1(1,:,:)=uw(2,:,:)/SQRT(uw(1,:,:))
         uw1(2,:,:)=-uxyw(:,:,2,5)/xyw(1,:,:)
         WHERE(xyw(1,:,:)==0)
            uw1(1,:,:)=0
            uw1(2,:,:)=0
         END WHERE
         uw2(1,:,:)=uw(3,:,:)/SQRT(uw(1,:,:))
         uw2(2,:,:)=uxyw(:,:,1,5)*gfac
         uw3(1,:,:)=uw(4,:,:)/SQRT(uw(1,:,:))
         uw3(2,:,:)=uw(6,:,:)*gfac
      END SELECT
c-----------------------------------------------------------------------
c        calculate the integral of scalar energies over the domain.
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
            
            value(1)=value(1) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *(SUM(uw1(1,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw2(1,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw3(1,ix:ix+1,iy:iy+1))**2)
            value(2)=value(2) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *(SUM(uw1(2,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw2(2,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw3(2,ix:ix+1,iy:iy+1))**2)
            value(3)=value(3) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *(SUM(uw1(1,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw2(1,ix:ix+1,iy:iy+1))**2)
            value(4)=value(4) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *(SUM(uw1(2,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw2(2,ix:ix+1,iy:iy+1))**2)
            value(5)=value(5) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *SUM(uw3(1,ix:ix+1,iy:iy+1))**2
            value(6)=value(6) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *(SUM(uw3(2,ix:ix+1,iy:iy+1))**2)
         ENDDO
      ENDDO
      value=value/2
      DEALLOCATE(uw1,uw2,uw3)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_energy
c-----------------------------------------------------------------------
c     subprogram 6. helix_maxUvsT.
c     generate u(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_maxUvsT(t,xyw,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw

      INTEGER :: nqty
      REAL(r8), DIMENSION(0:SIZE(uw,2)-1,0:SIZE(uw,3)-1) :: gfac
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
      REAL(r8), DIMENSION(:), ALLOCATABLE :: value
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
         nqty=4
         gfac=SQRT(1./(1. + (xyw(1,:,:)*Rinv)**2))
         ALLOCATE(utemp(nqty,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),value(nqty))
         utemp(1,:,:)=-uxyw(:,:,2,6)/xyw(1,:,:)
         WHERE(xyw(1,:,:) == 0)
            utemp(1,:,:)=0
         END WHERE
         utemp(2,:,:)=gfac*uxyw(:,:,1,6)
         utemp(3,:,:)=gfac*uw(4,:,:)
         utemp(4,:,:)=SQRT(utemp(1,:,:)**2 + utemp(2,:,:)**2
     $        + utemp(3,:,:)**2)
         CALL plotter_maxUvsT("maxUvsT",t,utemp,value)
         DEALLOCATE(utemp,value)
      CASE("XMHDhlx")
         nqty=4
         gfac=SQRT(1./(1. + (xyw(1,:,:)*Rinv)**2))
         ALLOCATE(utemp(nqty,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),value(nqty))
         utemp(1,:,:)=uw(2,:,:)/uw(1,:,:)
         WHERE(xyw(1,:,:) == 0)
            utemp(1,:,:)=0
         END WHERE
         utemp(2,:,:)=uw(3,:,:)/uw(1,:,:)
         utemp(3,:,:)=uw(4,:,:)/uw(1,:,:)
         utemp(4,:,:)=SQRT(utemp(1,:,:)**2 + utemp(2,:,:)**2
     $        + utemp(3,:,:)**2)
         CALL plotter_maxUvsT("maxUvsT",t,utemp,value)
         DEALLOCATE(utemp,value)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_maxUvsT
c-----------------------------------------------------------------------
c     subprogram 7. helix_UxyT.
c     generate U vs. time plot at a point.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_UxyT(t,xyw,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw

      REAL(r8) :: value_q
      REAL(r8), DIMENSION(3) :: value_swten
      REAL(r8), DIMENSION(6) :: value_en
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      CALL helix_energy(xyw,uw,uxyw,value_en)
      CALL helix_swtenergy(r_m,xyw,uw,uxyw,value_swten)
      CALL helix_q(xyw,uw,uxyw,q0_r,q0_delr,value_q)

      OPEN(UNIT=UxyT_unit,FILE="UxyT.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(UxyT_unit)REAL(t,4),REAL(value_en(1),4),
     $     REAL(value_en(2),4),REAL(value_en(3),4),
     $     REAL(value_en(4),4),REAL(value_en(5),4),
     $     REAL(value_en(6),4),REAL(value_q,4)
      CLOSE(UNIT=UxyT_unit)

      OPEN(UNIT=swten_unit,FILE="swten.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(swten_unit)REAL(t,4),REAL(value_swten(1),4),
     $     REAL(value_swten(2),4),REAL(value_swten(3),4)
      CLOSE(UNIT=swten_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_UxyT
c-----------------------------------------------------------------------
c     subprogram 8. helix_Uprofile.
c     generate slices of U for each time-step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)

      INTEGER, INTENT(IN) :: nxs
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt,uxyw


      CHARACTER(16) :: filename
      INTEGER :: ipts
      REAL(r8) :: q0_r,q0_delr,value,y_plus
      REAL(r8), DIMENSION(0:nxs) :: q_prof
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: uprofile,uprof1,uprof2
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
      REAL(r8), DIMENSION(0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1) :: gfac,rad
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
         gfac=1/SQRT(1+(xyw(1,:,:)*Rinv)**2)
         ALLOCATE(utemp(10,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         utemp(1,:,:)=uw(1,:,:)
         utemp(2,:,:)=uw(2,:,:)
         utemp(3,:,:)=gfac*uxyw(:,:,1,1)
         utemp(4,:,:)=gfac**2*(uxyw(:,:,1,1)-Rinv*xyw(1,:,:)*uw(2,:,:))
         utemp(5,:,:)=gfac**2*(Rinv*xyw(1,:,:)*uxyw(:,:,1,1)+uw(2,:,:))
     $        *Rinv
         utemp(6,:,:)=uxyw(:,:,2,1)**2/xyw(1,:,:)**2 
     $        + gfac**2*(uxyw(:,:,1,1)**2+uw(2,:,:)**2)
         WHERE(xyw(1,:,:)==0)
            utemp(6,:,:)=gfac**2*(uxyw(:,:,1,1)**2+uw(2,:,:)**2)
         END WHERE
         utemp(7,:,:)=uw(5,:,:)
         utemp(8,:,:)=-uxyw(:,:,2,6)/xyw(1,:,:)
         WHERE(xyw(1,:,:)==0)
            utemp(8,:,:)=0
         END WHERE
         utemp(9,:,:)=uw(4,:,:)
         utemp(10,:,:)=uxyw(:,:,1,6)*gfac

         IF(x1 <= 1.e-2)THEN
            ALLOCATE(uprof1(11,0:nxs),uprof2(11,0:nxs))
            IF(y1 < pi)THEN
               y_plus=y1+pi
            ELSE
               y_plus=y1-pi
            ENDIF
            filename=" "
            CALL plotter_Uprofile(x1,y1,x2,y1,nxs,xyw,xyw_kt,utemp,
     $           uprof1,.FALSE.,filename)
            CALL plotter_Uprofile(x1,y_plus,x2,y_plus,nxs,xyw,xyw_kt,
     $           utemp,uprof2,.FALSE.,filename)
            OPEN(UNIT=Uprofile_unit,FILE="Uprofile.bin",
     $           STATUS="UNKNOWN",POSITION="APPEND",FORM="UNFORMATTED")
            DO ipts=nxs,0,-1
               WRITE(Uprofile_unit)
     $              REAL(-uprof2(1,ipts),4),REAL(uprof2(2:11,ipts),4)
            ENDDO
            OPEN(UNIT=qprof_unit,FILE="qprof.bin",
     $           STATUS="UNKNOWN",POSITION="APPEND",FORM="UNFORMATTED")
            DO ipts=0,nxs
               WRITE(Uprofile_unit)
     $              REAL(uprof1(1,ipts),4),REAL(uprof1(2:11,ipts),4)
               q0_r = MAX(1e-2_r8,x1 + (x2-x1)*ipts/REAL(nxs,8))
               q0_delr = MIN(0.25*q0_r,0.5*(x2-x1)/REAL(nxs,8))
               CALL helix_q(xyw,uw,uxyw,q0_r,q0_delr,value)
               IF(ABS(value) < 1e20)WRITE(qprof_unit)
     $              REAL(uprof1(1,ipts),4),REAL(value,4)
            ENDDO
            WRITE(Uprofile_unit)
            WRITE(qprof_unit)
            CLOSE(UNIT=qprof_unit)
            CLOSE(UNIT=Uprofile_unit)
            DEALLOCATE(uprof1,uprof2)
         ELSE
            ALLOCATE(uprofile(11,0:nxs))
            filename="Uprofile"
            CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $           uprofile,.TRUE.,filename)
c-----------------------------------------------------------------------
c     calculate and output axisymmetric q-profile.
c-----------------------------------------------------------------------
            OPEN(UNIT=qprof_unit,FILE="qprof.bin",
     $           STATUS="UNKNOWN",POSITION="APPEND",FORM="UNFORMATTED")
            DO ipts=0,nxs
               q0_r = MAX(1e-2_r8,x1 + (x2-x1)*ipts/REAL(nxs,8))
               q0_delr = MIN(0.25*q0_r,0.5*(x2-x1)/REAL(nxs,8))
               CALL helix_q(xyw,uw,uxyw,q0_r,q0_delr,value)
               IF(ABS(value) < 1e20)WRITE(qprof_unit)
     $              REAL(uprofile(1,ipts),4),REAL(value,4)
            ENDDO
            WRITE(qprof_unit)
            CLOSE(UNIT=qprof_unit)
            DEALLOCATE(uprofile)
         ENDIF
         DEALLOCATE(utemp)

      CASE("XMHDhlx")
         gfac=1./SQRT(1. + (xyw(1,:,:)*Rinv)**2)         
         ALLOCATE(utemp(15,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1),
     $        uprofile(16,0:nxs))
         utemp(1,:,:)=uw(1,:,:)
         utemp(2,:,:)=uw(2,:,:)
         utemp(3,:,:)=uw(3,:,:)
         utemp(4,:,:)=uw(4,:,:)
         utemp(5,:,:)=uw(5,:,:)
         utemp(6,:,:)=uw(6,:,:)
         utemp(7,:,:)=uw(7,:,:)
         utemp(8,:,:)=uw(8,:,:)
         utemp(9,:,:)=uw(9,:,:)
         utemp(10,:,:)=gfac*uxyw(:,:,1,5)
         utemp(11,:,:)=gfac**2*(Rinv*xyw(1,:,:)*uxyw(:,:,1,5)+uw(6,:,:))
         utemp(12,:,:)=(uw(4,:,:)-uw(1,:,:)*uw(7,:,:))/di
         utemp(13,:,:)=uw(2,:,:)/uw(1,:,:)
         utemp(14,:,:)=uw(8,:,:)/uw(1,:,:)
         utemp(15,:,:)=uw(9,:,:)/uw(1,:,:)
         filename="Uprofile"
         CALL plotter_Uprofile(x1,y1,x2,y2,nxs,xyw,xyw_kt,utemp,
     $        uprofile,.TRUE.,filename)
c-----------------------------------------------------------------------
c     calculate and output axisymmetric q-profile.
c-----------------------------------------------------------------------
         OPEN(UNIT=qprof_unit,FILE="qprof.bin",
     $        STATUS="UNKNOWN",POSITION="APPEND",FORM="UNFORMATTED")
         DO ipts=0,nxs
            q0_r = MAX(1e-2_r8,x1 + (x2-x1)*ipts/REAL(nxs,8))
            q0_delr = MIN(0.25*q0_r,0.5*(x2-x1)/REAL(nxs,8))
            CALL helix_q(xyw,uw,uxyw,q0_r,q0_delr,value)
            IF(ABS(value) < 1e20)
     $           WRITE(qprof_unit)REAL(uprofile(1,ipts),4),REAL(value,4)
         ENDDO
         WRITE(qprof_unit)
         CLOSE(UNIT=qprof_unit)
         DEALLOCATE(utemp,uprofile)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_Uprofile
c-----------------------------------------------------------------------
c     subprogram 9. helix_q.
c     generate a time plot for q value.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_q(xyw,uw,uxyw,r0,delr,value)

      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw
      REAL(r8), INTENT(IN) :: r0,delr
      REAL(r8), INTENT(OUT) :: value

      INTEGER :: ix,iy,nx,ny
      REAL(r8) :: r_c,Be0,phi,psi
      REAL(r8), DIMENSION(5) :: a
      REAL(r8), DIMENSION(0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1) :: gfac,rad
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uw_local
c-----------------------------------------------------------------------
c     prepare auxiliary parameters
c-----------------------------------------------------------------------
      value=zero
      phi=0
      psi=0
      Be0=1/Rinv
      nx=SIZE(xyw,2)-1
      ny=SIZE(xyw,3)-1
c-----------------------------------------------------------------------
c     calculate poloidal and toroidal fluxes
c-----------------------------------------------------------------------
      ALLOCATE(uw_local(2,0:nx,0:ny))
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
         IF(MAXVAL(ABS(uw(2,:,:))) > 1/(2*Rinv))Be0=0
         gfac=1/(1+(xyw(1,:,:)*Rinv)**2)
         uw_local(1,:,:)=Rinv*xyw(1,:,:)*gfac*
     $        (Rinv*xyw(1,:,:)*uxyw(:,:,1,1) + uw(2,:,:))
         uw_local(2,:,:)=gfac*(uxyw(:,:,1,1) 
     $        - Rinv*xyw(1,:,:)*uw(2,:,:))
      CASE("XMHDhlx")
         gfac=1/(1+(xyw(1,:,:)*Rinv)**2)
         uw_local(1,:,:)=Rinv*xyw(1,:,:)*gfac*
     $        (Rinv*xyw(1,:,:)*uxyw(:,:,1,5) + uw(6,:,:))
         uw_local(2,:,:)=gfac*(uxyw(:,:,1,5) 
     $        - Rinv*xyw(1,:,:)*uw(6,:,:))
      END SELECT
c-----------------------------------------------------------------------
c     calculate the integrals of fluxes .
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         DO iy=0,ny-1
            r_c=SUM(xyw(1,ix:ix+1,iy:iy+1))/4._r8
            IF(r_c < r0-delr .OR. r_c > r0+delr)CYCLE
            
            a(1)=(xyw(2,ix,iy)+xyw(2,ix,iy+1))*
     $           (xyw(1,ix,iy+1)-xyw(1,ix,iy))/2._r8
            a(2)=(xyw(2,ix,iy+1)+xyw(2,ix+1,iy+1))*
     $           (xyw(1,ix+1,iy+1)-xyw(1,ix,iy+1))/2._r8
            a(3)=(xyw(2,ix+1,iy+1)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix+1,iy)-xyw(1,ix+1,iy+1))/2._r8
            a(4)=(xyw(2,ix,iy)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix,iy)-xyw(1,ix+1,iy))/2._r8
            a(5)=ABS(SUM(a(1:4)))
            
            phi=phi + a(5)/4._r8
     $           *SUM(uw_local(1,ix:ix+1,iy:iy+1))
            psi=psi + a(5)/4._r8
     $           *SUM(uw_local(2,ix:ix+1,iy:iy+1))
         ENDDO
      ENDDO
      value=ABS(phi/psi)
      DEALLOCATE(uw_local)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_q
c-----------------------------------------------------------------------
c     subprogram 10. helix_Bspectr.
c     generate spectrum of B-field energy vs. time.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_Bspectr(t,xyw,xyw_kt,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt,uxyw


      CHARACTER(16) :: filename
      INTEGER :: ir,nr,ntau,nmode
      REAL(r8) :: rval,delr,deltau,eps,Be0
      REAL(r8), DIMENSION(1:9) :: Bmodes
      REAL(r8), DIMENSION(0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1) :: gfac,
     $     Br,Btau,Be
      REAL(r8), DIMENSION(6,0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1) :: utemp
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: utau_prof,ur_prof
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      filename=" "
      nr=SIZE(xyw,2)/4
      ntau=SIZE(xyw,3)/4
      delr=1._r8/REAL(nr,8)
      deltau=1._r8/REAL(ntau,8)
      Be0=1/Rinv
      eps=1e-5
      Bmodes=0

      gfac=1/SQRT(1+(xyw(1,:,:)*Rinv)**2)         
      ALLOCATE(utau_prof(7,0:ntau),ur_prof(6,1:nr))
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
         IF(MINVAL(ABS(uw(2,:,:))) < 1/(2*Rinv))Be0=0
         Br=-uxyw(:,:,2,1)/xyw(1,:,:)
         Btau=gfac*uxyw(:,:,1,1)
         Be=gfac*uw(2,:,:) - Be0
      CASE("XMHDhlx")
         IF(MINVAL(ABS(uw(6,:,:))) < 1/(2*Rinv))Be0=0
         Br=-uxyw(:,:,2,5)/xyw(1,:,:)
         Btau=gfac*uxyw(:,:,1,5)
         Be=gfac*uw(6,:,:) - Be0
      END SELECT
      WHERE(xyw(1,:,:) == 0)
         Br=0
      END WHERE
      DO nmode=1,9
         utemp = 0
         ur_prof = 0
         utau_prof = 0
         
         utemp(1,:,:)=Br*SIN(nmode*xyw(2,:,:))
         utemp(2,:,:)=Btau*SIN(nmode*xyw(2,:,:))
         utemp(3,:,:)=Be*SIN(nmode*xyw(2,:,:))
         utemp(4,:,:)=Br*COS(nmode*xyw(2,:,:))
         utemp(5,:,:)=Btau*COS(nmode*xyw(2,:,:))
         utemp(6,:,:)=Be*COS(nmode*xyw(2,:,:))
         DO ir=1,nr
            rval=(REAL(ir,8) - .5)*delr
            CALL plotter_Uprofile(rval,eps,rval,twopi-eps,ntau,
     $           xyw,xyw_kt,utemp,utau_prof,.FALSE.,filename)
            ur_prof(1,ir) = SUM(utau_prof(2,1:ntau))
            ur_prof(2,ir) = SUM(utau_prof(3,1:ntau))
            ur_prof(3,ir) = SUM(utau_prof(4,1:ntau))
            ur_prof(4,ir) = SUM(utau_prof(5,1:ntau))
            ur_prof(5,ir) = SUM(utau_prof(6,1:ntau))
            ur_prof(6,ir) = SUM(utau_prof(7,1:ntau))
            Bmodes(nmode) = Bmodes(nmode)
     $              + SUM(ur_prof(1:6,ir)**2)*rval
         ENDDO
      ENDDO
      Bmodes = Bmodes*twopi*delr*deltau**2
      DEALLOCATE(utau_prof,ur_prof)
      OPEN(UNIT=Bmodes_unit,FILE="Bmodes.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WHERE(Bmodes == 0)
         Bmodes=TINY(Bmodes)
      END WHERE
      WRITE(Bmodes_unit)REAL(t,4),REAL(LOG10(Bmodes),4)
      CLOSE(UNIT=Bmodes_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_Bspectr
c-----------------------------------------------------------------------
c     subprogram 11. helix_swtenergy.
c     generate a time plot for an integral of helical in- and 
c     out-of-plane magnetic energy inside some given radius.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_swtenergy(r_m,xyw,uw,uxyw,value)

      REAL(r8) :: r_m
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw
      REAL(r8), DIMENSION(3), INTENT(OUT) :: value

      INTEGER :: ix,iy,nx,ny
      REAL(r8), DIMENSION(5) :: a
      REAL(r8), DIMENSION(0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1) :: 
     $     gfac,rad
      REAL(r8), DIMENSION(1,0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1) :: 
     $     uw1,uw2,uw3
c-----------------------------------------------------------------------
c     prepare auxiliary parameters
c-----------------------------------------------------------------------
      value=0
      nx=SIZE(xyw,2)-1
      ny=SIZE(xyw,3)-1
c-----------------------------------------------------------------------
c     calculate energy integral
c-----------------------------------------------------------------------
      gfac=1/SQRT(1+(xyw(1,:,:)*Rinv)**2)
      SELECT CASE(job_type)
      CASE("fourfieldhlx")
c-----------------------------------------------------------------------
c     magnetic energy components.
c-----------------------------------------------------------------------
         uw1(1,:,:)=-uxyw(:,:,2,1)/xyw(1,:,:)
         WHERE(xyw(1,:,:)==0)
            uw1(1,:,:)=0
         END WHERE
         uw2(1,:,:)=uxyw(:,:,1,1)*gfac
         uw3(1,:,:)=(uw(2,:,:) - 1._r8/Rinv)*gfac
      CASE("XMHDhlx")
c-----------------------------------------------------------------------
c     magnetic energy components.
c-----------------------------------------------------------------------
         uw1(1,:,:)=-uxyw(:,:,2,5)/xyw(1,:,:)
         WHERE(xyw(1,:,:)==0)
            uw1(1,:,:)=0
         END WHERE
         uw2(1,:,:)=uxyw(:,:,1,5)*gfac
         uw3(1,:,:)=(uw(6,:,:) - 1._r8/Rinv)*gfac
      END SELECT
      WHERE(xyw(1,:,:) > r_m)
         uw1(1,:,:)=0
         uw2(1,:,:)=0
         uw3(1,:,:)=0
      END WHERE
c-----------------------------------------------------------------------
c        calculate the integral of scalar energies over the domain.
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
            
            value(1)=value(1) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *(SUM(uw1(1,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw2(1,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw3(1,ix:ix+1,iy:iy+1))**2)
            value(2)=value(2) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *(SUM(uw1(1,ix:ix+1,iy:iy+1))**2+
     $           SUM(uw2(1,ix:ix+1,iy:iy+1))**2)
            value(3)=value(3) + a(5)/4**3
     $           *SUM(xyw(1,ix:ix+1,iy:iy+1))
     $           *SUM(uw3(1,ix:ix+1,iy:iy+1))**2
         ENDDO
      ENDDO
      value=value/2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_swtenergy
      END MODULE helix_mod
