c-----------------------------------------------------------------------
c     file euler2D.f.
c     contains specifications for 2D Euler compressible flow.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module organization.
c-----------------------------------------------------------------------
c     0. euler2D_mod.
c     1. euler2D_openbc.
c     2. euler2D_nrbc.
c     3. euler2D_RDL.
c     4. euler2D_RDLt.
c     5. euler2D_bc_jac.
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
c     p. physics_main
c-----------------------------------------------------------------------
c     subprogram 0. euler2D_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE euler2D_mod
      USE local_mod
      IMPLICIT NONE

      CHARACTER(16) :: init_type="supersonic"
      LOGICAL :: ifbound_visc=.FALSE.
      REAL(r8), PARAMETER :: gamma=1.4_r8,eps=1.e-10_r8
c      REAL(r8), PARAMETER :: gamma=5._r8/3._r8,eps=1.e-10_r8
      REAL(r8) :: mach=1.75,mu=0.0,pin=0.7143,pout=0.1,p0=0.,
     $     lx=10.0,alpha=1.,rhoin=1.,vin=1.,rhoout=0.1,vout=1.,
     $     ddiff=0.,Ldiv=4.,Lconv=2.,ediff=0.,Adiv=0.,Aconv=0.,pmax=1.,
     $     rhomax=1.,expand=1.,gr_curve=0.,poutmax=0.

      REAL(r8) :: x1,x2,y1,y2,rs

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. physics_openbc.
c     computes normal flux at the boundary for an "open" boundary.
c     uses Roe's method per LeVeque Eq. 15.40.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler2D_openbc(lrtb,u,ux,uy,nhat,rot,c)

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy,nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: rot
      REAL(r8), DIMENSION(:,:,:), INTENT(INOUT) :: c

      INTEGER :: i,j
      REAL(r8), DIMENSION(2,2,SIZE(u,2),SIZE(u,3)) :: irot
      REAL(r8), DIMENSION(SIZE(u,2),SIZE(u,3),2) :: p
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: fcorr,
     $     fbound
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3),2) :: q,favg
      REAL(r8), DIMENSION(4,4,SIZE(u,2),SIZE(u,3)) :: RDL
c-----------------------------------------------------------------------
c     compute pressure and other required quantities.
c     momentum is rotated into boundary frame of reference with the x
c     direction normal to the surface.
c-----------------------------------------------------------------------
      p(:,:,1) = (gamma-one)*(u(4,:,:) - half*(u(2,:,:)**2
     $     + u(3,:,:)**2)/u(1,:,:))
      q(1,:,:,1) = u(1,:,:)
      q(2,:,:,1) = rot(1,1,:,:)*u(2,:,:) + rot(1,2,:,:)*u(3,:,:)
      q(3,:,:,1) = rot(2,1,:,:)*u(2,:,:) + rot(2,2,:,:)*u(3,:,:)
      q(4,:,:,1) = u(4,:,:)

      SELECT CASE(lrtb)
      CASE("right")
         q(1,:,:,2) = rhoout
         q(2,:,:,2) = vout*rhoout
         q(3,:,:,2) = zero
         p(:,:,2) = pout
         q(4,:,:,2) = p(:,:,2)/(gamma-one)
     $        + half*(q(2,:,:,2)**2 + q(3,:,:,2)**2)/q(1,:,:,2)
      CASE("left")
         q(1,:,:,2) = rhoin
         q(2,:,:,2) = vin*rhoin
         q(3,:,:,2) = zero
         p(:,:,2) = pin
         q(4,:,:,2) = p(:,:,2)/(gamma-one)
     $        + half*(q(2,:,:,2)**2 + q(3,:,:,2)**2)/q(1,:,:,2)
      END SELECT
c-----------------------------------------------------------------------
c     compute interior and exterior flux which will be averaged.
c-----------------------------------------------------------------------
      favg(1,:,:,:) = q(2,:,:,:)
      favg(2,:,:,:) = q(2,:,:,:)**2/q(1,:,:,:) + p
      favg(3,:,:,:) = q(2,:,:,:)*q(3,:,:,:)/q(1,:,:,:)
      favg(4,:,:,:) = q(2,:,:,:)/q(1,:,:,:)*(q(4,:,:,:) + p)
c-----------------------------------------------------------------------
c     compute the absolute value of the approximate flux jacobian.
c-----------------------------------------------------------------------
      CALL euler2D_RDL(q,"hat","abs",RDL)
c-----------------------------------------------------------------------
c        compute viscous correction term for Roe's method.
c-----------------------------------------------------------------------
      DO i=1,SIZE(u,2)
         DO j=1,SIZE(u,3)
            fcorr(:,i,j) = -half*MATMUL(RDL(:,:,i,j),
     $           (q(:,i,j,2) - q(:,i,j,1)))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute boundary flux.
c     fbound will now contain the normal fluxes for the variables in the
c     rotated c.s.
c-----------------------------------------------------------------------
      fbound = half*(favg(:,:,:,1) + favg(:,:,:,2)) + fcorr
c-----------------------------------------------------------------------
c     inverse rotate the normal momentum flux vector.
c-----------------------------------------------------------------------
      DO i=1,SIZE(u,2)
         DO j=1,SIZE(u,3)
            irot(:,:,i,j) = TRANSPOSE(rot(:,:,i,j))
         ENDDO
      ENDDO
      c(1,:,:) = fbound(1,:,:)
      c(2,:,:) = irot(1,1,:,:)*fbound(2,:,:) 
     $     + irot(1,2,:,:)*fbound(3,:,:)
      c(3,:,:) = irot(2,1,:,:)*fbound(2,:,:) 
     $     + irot(2,2,:,:)*fbound(3,:,:)
      c(4,:,:) = fbound(4,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler2D_openbc
c-----------------------------------------------------------------------
c     subprogram 2. euler2D_nrbc.
c     computes normal flux at the boundary for a non-reflecting boundary.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler2D_nrbc(u,ux,uy,nhat,rot_in,c)

      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy,nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: rot_in
      REAL(r8), DIMENSION(:,:,:), INTENT(INOUT) :: c

      INTEGER :: i,j
      REAL(r8), DIMENSION(2,SIZE(u,2),SIZE(u,3)) :: that
      REAL(r8), DIMENSION(4,4,SIZE(u,2),SIZE(u,3)) :: rot,irot
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: ddn_u,ddt_u
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3),2) :: q
      REAL(r8), DIMENSION(4,4,SIZE(u,2),SIZE(u,3)) :: RDL,RDL2,
     $     RDLt,RDLt2
c-----------------------------------------------------------------------
c     rotate into boundary frame of reference with the x
c     direction normal to the surface.
c-----------------------------------------------------------------------
      rot=0
      irot=0
      DO i=1,SIZE(u,2)
         DO j=1,SIZE(u,3)
            rot(1,1,i,j) = one
            rot(4,4,i,j) = one
            rot(2:3,2:3,i,j) = rot_in(:,:,i,j)
            irot(:,:,i,j) = TRANSPOSE(rot(:,:,i,j))
         ENDDO
      ENDDO
      q(1,:,:,1) = u(1,:,:)
      q(2,:,:,1) = rot_in(1,1,:,:)*u(2,:,:) + rot_in(1,2,:,:)*u(3,:,:)
      q(3,:,:,1) = rot_in(2,1,:,:)*u(2,:,:) + rot_in(2,2,:,:)*u(3,:,:)
      q(4,:,:,1) = u(4,:,:)
c-----------------------------------------------------------------------
c     non-reflecting boundary condition.
c     impose d/dt(u) + irot*Aplus*d/dn(un) + irot*B*d/dt(un) = 0,
c     where un is u rotated into the normal coordinate system and
c     irot is the inverse rotation matrix.
c     RDL = Aplus
c     RDLt = B
c     uses robin b.c. with c = -(RDL*d/dn(un) + RDLt*d/dt(ut)).
c-----------------------------------------------------------------------
      CALL euler2D_RDL(q,"q1","plus",RDL)
      CALL euler2D_RDLt(q(:,:,:,1),RDLt)
c-----------------------------------------------------------------------
c     c = -(irot*RDL*d/dn(un) + irot*RDLt*d/dt(un))
c       = -(irot*RDL*rot*d/dn(u) + irot*RDLt*rot*d/dt(u))
c     This assumes that d/dt(rot) = 0 which is only true for straight
c     boundaries.
c-----------------------------------------------------------------------
      DO i=1,SIZE(u,2)
         DO j=1,SIZE(u,3)
            RDL2(:,:,i,j) = MATMUL(MATMUL(irot(:,:,i,j),
     $           RDL(:,:,i,j)),rot(:,:,i,j))
            RDLt2(:,:,i,j) = MATMUL(MATMUL(irot(:,:,i,j),
     $           RDLt(:,:,i,j)),rot(:,:,i,j))
         ENDDO
      ENDDO

      ddn_u(1,:,:) = nhat(1,:,:)*ux(1,:,:) + nhat(2,:,:)*uy(1,:,:)
      ddn_u(2,:,:) = nhat(1,:,:)*ux(2,:,:) + nhat(2,:,:)*uy(2,:,:)
      ddn_u(3,:,:) = nhat(1,:,:)*ux(3,:,:) + nhat(2,:,:)*uy(3,:,:)
      ddn_u(4,:,:) = nhat(1,:,:)*ux(4,:,:) + nhat(2,:,:)*uy(4,:,:)

      that(1,:,:) = -nhat(2,:,:)
      that(2,:,:) = nhat(1,:,:)

      ddt_u(1,:,:) = that(1,:,:)*ux(1,:,:) + that(2,:,:)*uy(1,:,:)
      ddt_u(2,:,:) = that(1,:,:)*ux(2,:,:) + that(2,:,:)*uy(2,:,:)
      ddt_u(3,:,:) = that(1,:,:)*ux(3,:,:) + that(2,:,:)*uy(3,:,:)
      ddt_u(4,:,:) = that(1,:,:)*ux(4,:,:) + that(2,:,:)*uy(4,:,:)

      DO i=1,SIZE(u,2)
         DO j=1,SIZE(u,3)
            c(:,i,j) = -(MATMUL(RDL2(:,:,i,j),ddn_u(:,i,j))
     $           + MATMUL(RDLt2(:,:,i,j),ddt_u(:,i,j)))
         ENDDO
      ENDDO

      IF(mpi_rank==0.and.1==2)THEN
         write(*,*) "nhat: ",nhat(1,2,1),nhat(2,2,1)
         write(*,*) "rot 1,: ",rot(1,:,2,1)
         write(*,*) "rot 2,: ",rot(2,:,2,1)
         write(*,*) "rot 3,: ",rot(3,:,2,1)
         write(*,*) "rot 4,: ",rot(4,:,2,1)
         write(*,*) ".."

         write(*,*) "q1= ",q(1,2,1,:)
         write(*,*) "q2= ",q(2,2,1,:)
         write(*,*) "q3= ",q(3,2,1,:)
         write(*,*) "q4= ",q(4,2,1,:)
         write(*,*) ".."

         write(*,*) ".."
         write(*,*) ".."
      ENDIF

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler2D_nrbc
c-----------------------------------------------------------------------
c     subprogram 3. euler2D_RDL.
c     computes the absolute value of the approximate flux jacobian,
c     or the positive- or negative-going part of the flux jacobian.
c     R is a matrix of the +/- going right eigenvectors, D is a matrix
c     of the characteristics, and L is a matrix of +/- going left 
c     eigenvectors. q is the Roe-averaged vector of conserved variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler2D_RDL(q,type,direction,RDL)
      
      CHARACTER(*), INTENT(IN) :: type,direction
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: q
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: RDL

      INTEGER :: i,j
      REAL(r8), PARAMETER :: delta=1.e-6
      REAL(r8), DIMENSION(SIZE(q,2),SIZE(q,3)) :: mn,mt,rho,ke,
     $     p,tgp,a,denom,hi,hb,h,vn,vt
      REAL(r8), DIMENSION(4,4,SIZE(q,2),SIZE(q,3)) :: R,D,L
c-----------------------------------------------------------------------
c     initialize arrays.
c-----------------------------------------------------------------------
      D=0
      RDL=0

      SELECT CASE(type)
      CASE("hat")
         p = (gamma - 1)*(q(4,:,:,1) 
     $        - half*(q(2,:,:,1)**2 + q(3,:,:,1)**2)/q(1,:,:,1))
         hi = (q(4,:,:,1) + p)/q(1,:,:,1)
      
         p = (gamma - 1)*(q(4,:,:,2) 
     $        - half*(q(2,:,:,2)**2 + q(3,:,:,2)**2)/q(1,:,:,2))
         hb = (q(4,:,:,2) + p)/q(1,:,:,2)

         h = (SQRT(q(1,:,:,1))*hi + SQRT(q(1,:,:,2))*hb)
     $        /(SQRT(q(1,:,:,1)) + SQRT(q(1,:,:,2)))
         vn = (SQRT(q(1,:,:,1))*q(2,:,:,1)/q(1,:,:,1) 
     $        + SQRT(q(1,:,:,2))*q(2,:,:,2)/q(1,:,:,2))
     $        /(SQRT(q(1,:,:,1)) + SQRT(q(1,:,:,2)))
         vt = (SQRT(q(1,:,:,1))*q(3,:,:,1)/q(1,:,:,1) 
     $        + SQRT(q(1,:,:,2))*q(3,:,:,2)/q(1,:,:,2))
     $        /(SQRT(q(1,:,:,1)) + SQRT(q(1,:,:,2)))
         rho = SQRT(MAX(q(1,:,:,1),eps)*MAX(q(1,:,:,2),eps))
      CASE("q1")
         rho = q(1,:,:,1)
         p = (gamma - 1)*(q(4,:,:,1) 
     $        - half*(q(2,:,:,1)**2 + q(3,:,:,1)**2)/q(1,:,:,1))
         h = (q(4,:,:,1) + p)/q(1,:,:,1)
         vn = q(2,:,:,1)/rho
         vt = q(3,:,:,1)/rho
      END SELECT

      a = SQRT((gamma-1)*(h - half*(vn**2 + vt**2)))
      p = MAX((gamma-1)/gamma*(rho*h - half*rho*(vn**2 + vt**2)),eps)

      mn = rho*vn
      mt = rho*vt

      WHERE(ABS(mn) <= eps) mn = eps
      WHERE(ABS(mt) <= eps) mt = eps
      WHERE(mn == mt) mt = mt + eps

      ke = (mn**2 + mt**2)/(2*rho)
      tgp = two*gamma*p
      denom = (mn**2 - mt**2)*p*gamma
c-----------------------------------------------------------------------
c     compute D, abs(D), or D+
c     and implement Harten's entropy fix for abs(D).
c-----------------------------------------------------------------------
      D(1,1,:,:) = mn/rho
      D(2,2,:,:) = mn/rho
      D(3,3,:,:) = mn/rho - a
      D(4,4,:,:) = mn/rho + a
      SELECT CASE(direction)
      CASE("abs")
         D = ABS(D)
         DO i=1,4
            WHERE(D(i,i,:,:) < delta) 
     $           D(i,i,:,:) = (D(i,i,:,:)**2 + delta**2)/(2*delta)
         ENDDO
      CASE("plus")
         D = half*(D + ABS(D))
      CASE("minus")
         D = half*(D - ABS(D))
      CASE DEFAULT
         CALL program_stop("euler2D_RDL: bad case.")
      END SELECT
c-----------------------------------------------------------------------
c     compute right and left eigenvector matrices for the positive and 
c     negative waves.
c-----------------------------------------------------------------------
      R(1,1,:,:) = one
      R(2,1,:,:) = mn/rho
      R(3,1,:,:) = zero
      R(4,1,:,:) = (mn**2 - mt**2)/(2*rho**2)
      
      R(1,2,:,:) = one
      R(2,2,:,:) = mn/rho
      R(3,2,:,:) = (mt**2 - mn**2)/(2*mt*rho)
      R(4,2,:,:) = zero
      
      R(1,3,:,:) = one
      R(2,3,:,:) = mn/rho - a
      R(3,3,:,:) = mt/rho
      R(4,3,:,:) = h - a*mn/rho
      
      R(1,4,:,:) = one
      R(2,4,:,:) = mn/rho + a
      R(3,4,:,:) = mt/rho
      R(4,4,:,:) = h + a*mn/rho

      L(1,1,:,:) = one + two*ke**2*(1-gamma)*rho/denom
      L(2,1,:,:) = two*ke*mt**2*(gamma-1)/denom
      L(3,1,:,:) = (a*mn + ke*(gamma-1))/tgp
      L(4,1,:,:) = (-a*mn + ke*(gamma-1))/tgp
      
      L(1,2,:,:) = two*ke*mn*(gamma-1)*rho/denom
      L(2,2,:,:) = two*mn*mt**2*(1-gamma)/denom
      L(3,2,:,:) = (mn*(1-gamma) - a*rho)/tgp
      L(4,2,:,:) = (mn*(1-gamma) + a*rho)/tgp
      
      L(1,3,:,:) = two*mt*(ke*(gamma-1) + p*gamma)*rho/denom
      L(2,3,:,:) = -two*mt*(mt**2*(gamma-1) + p*gamma*rho)/denom
      L(3,3,:,:) = mt*(1-gamma)/tgp
      L(4,3,:,:) = L(3,3,:,:)
      
      L(1,4,:,:) = two*ke*(1-gamma)*rho**2/denom
      L(2,4,:,:) = two*mt**2*(gamma-1)*rho/denom
      L(3,4,:,:) = (gamma-1)*rho/tgp
      L(4,4,:,:) = L(3,4,:,:)

c-----------------------------------------------------------------------
c     compute R.abs(D).L.
c-----------------------------------------------------------------------
      DO i=1,SIZE(q,2)
         DO j=1,SIZE(q,3)
            RDL(:,:,i,j) = MATMUL(MATMUL(R(:,:,i,j),D(:,:,i,j)),
     $           L(:,:,i,j))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler2D_RDL
c-----------------------------------------------------------------------
c     subprogram 4. euler2D_RDLt.
c     computes the tangential flux jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler2D_RDLt(q,RDL)
      
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: q
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: RDL

      REAL(r8), DIMENSION(SIZE(q,2),SIZE(q,3)) :: e,rho,vn,vt,ke
c-----------------------------------------------------------------------
c     initialize arrays.
c-----------------------------------------------------------------------
      RDL=0

      rho = q(1,:,:)
      e = q(4,:,:)
      vn = q(2,:,:)/rho
      vt = q(3,:,:)/rho
      ke = half*rho*(vn**2 + vt**2)
c-----------------------------------------------------------------------
c     enter values for flux jacobian.
c-----------------------------------------------------------------------
      RDL(1,3,:,:) = one

      RDL(2,1,:,:) = -vn*vt
      RDL(2,2,:,:) = vt
      RDL(2,3,:,:) = vn
      
      RDL(3,1,:,:) = half*(vt**2*(gamma-3.) + vn**2*(gamma-1.))
      RDL(3,2,:,:) = vn*(1.-gamma)
      RDL(3,3,:,:) = vt*(3.-gamma)
      RDL(3,4,:,:) = gamma-1.
      
      RDL(4,1,:,:) = vt*(2.*ke*(gamma-1.) - gamma*e)/rho
      RDL(4,2,:,:) = vn*vt*(1.-gamma)
      RDL(4,3,:,:) = half*(1.-gamma)*(vn**2 + 3.*vt**2) + gamma*e/rho
      RDL(4,4,:,:) = gamma*vt
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler2D_RDLt
c-----------------------------------------------------------------------
c     subprogram 5. euler2D_bc_jac.
c     computes jacobian for a boundary condition.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE euler2D_bc_jac(lrtb,t,x,y,nhat,
     $     u,ux,uy,c_u,c_ux,c_uy)

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy,nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: c_u,c_ux,c_uy

      INTEGER :: i
      REAL(r8), PARAMETER :: du=1.e-6_r8
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: 
     $     u2,ux2,uy2,f,f2,nil
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

      nil=0
      CALL physics_edge_rhs(lrtb,t,x,y,nhat,u,ux,uy,nil,nil,nil,f)

      DO i=1,SIZE(u,1)
         u2(i,:,:) = u(i,:,:) + du
         CALL physics_edge_rhs(lrtb,t,x,y,nhat,u2,ux,uy,nil,nil,nil,f2)
         u2(i,:,:) = u(i,:,:)
         c_u(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

      DO i=1,SIZE(ux,1)
         ux2(i,:,:) = ux(i,:,:) + du
         CALL physics_edge_rhs(lrtb,t,x,y,nhat,u,ux2,uy,nil,nil,nil,f2)
         ux2(i,:,:) = ux(i,:,:)
         c_ux(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

      DO i=1,SIZE(uy,1)
         uy2(i,:,:) = uy(i,:,:) + du
         CALL physics_edge_rhs(lrtb,t,x,y,nhat,u,ux,uy2,nil,nil,nil,f2)
         uy2(i,:,:) = uy(i,:,:)
         c_uy(:,i,:,:) = (f2(:,:,:) - f(:,:,:))/du
      ENDDO

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE euler2D_bc_jac
      END MODULE euler2D_mod
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
      USE euler2D_mod
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
      NAMELIST/euler2D_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $  dtmax,tmax,nstep,init_type,mach,mu,pin,pout,lx,alpha,rhoin,vin,
     $  rhoout,vout,ddiff,Lconv,Ldiv,ediff,Aconv,Adiv,ifbound_visc,pmax,
     $  rhomax,expand,gr_curve,p0,poutmax

c-----------------------------------------------------------------------
c     definitions.
c-----------------------------------------------------------------------

c     mach - input mach number
c     mu - viscosity parameter
c     pin - pressure at the inflow
c     pout - pressure at the outflow
c     rhoout - density at the outflow
c     alpha - input variable ramp parameter
c     rhoin - density at inflow
c     vin - velocity at inflow
c     vout - velocity at ouflow

c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=euler2D_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      physics_type="euler2D"

      nqty=4
      nqty_schur=0
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
      USE euler2D_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     set PDE flags
c-----------------------------------------------------------------------
      ! none.
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(ifbound_visc,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(gr_curve,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mach,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(mu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pout,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(poutmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(p0,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(rhoout,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(vout,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(rhoin,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(vin,1,MPI_DOUBLE_PRECISION,0,comm,ierr) 
      CALL MPI_Bcast(lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ddiff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ediff,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Lconv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Ldiv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Adiv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Aconv,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(pmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(rhomax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(expand,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     set source/sink locations
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("source_sink")
         ! Sq_Pipe Source/Sink Location
ccc         x1=-0.75
ccc         x2=+0.75
ccc         y1=2.25
ccc         y2=0.25
ccc         rs=0.10

         ! S_Pipe Source/Sink Location
         x1=0.252
         x2=3.748
         y1=0.7755
         y2=1.122
         rs=0.1
      CASE("circle")
         x1=+0.75
         x2=-0.75
         y1=0.0
         y2=0.0
         rs=0.1
      END SELECT
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
      SUBROUTINE physics_init(xpi,ypi,ui)
      USE euler2D_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: ui
      
      REAL(r8), DIMENSION(SIZE(xpi,1),SIZE(xpi,2)) :: 
     $     area,area_s,dadx,p
      INTEGER :: i,j

c-----------------------------------------------------------------------
c     initial condition.
c-----------------------------------------------------------------------
      area = one + expand*half*(TANH(0.8*xpi - 3.) - TANH(-3.))

      SELECT CASE(init_type)
      CASE("supersonic")
         ui(1,:,:) = 1.0
         ui(2,:,:) = mach*SQRT(gamma*pin)
         ui(3,:,:) = 0.0
         ui(4,:,:) = pin/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("zero_normal")
         ui(1,:,:) = 1.0
         ui(2,:,:) = mach*SQRT(gamma*pin)
         ui(3,:,:) = mach*SQRT(gamma*pin)*ypi/area
     $        *.4*expand/COSH(3.-.8*xpi)**2
         ui(4,:,:) = pin/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("open_nozzle")
         ui(1,:,:) = rhoout
         ui(2,:,:) = mach
         ui(3,:,:) = 0.
         ui(4,:,:) = pout/(gamma-one)
      CASE("open_nozzle2")
         ui(1,:,:) = rhoout
         ui(2,:,:) = mach
         ui(3,:,:) = zero
         ui(4,:,:) = pout/(gamma-one)
      CASE("open_right")
         ui(1,:,:) = rhoin
         ui(2,:,:) = mach*SQRT(gamma*pin/rhoin)*rhoin
         ui(3,:,:) = mach*SQRT(gamma*pin/rhoin)*rhoin*ypi/area
     $        *.4*expand/COSH(3.-.8*xpi)**2
         ui(4,:,:) = pin/(gamma-one)
     $        + half*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("open_channel")
         ui(1,:,:) = .1_r8
         ui(2,:,:) = mach
         ui(3,:,:) = zero
         ui(4,:,:) = .01_r8/(gamma-one)
      CASE("phump")
         ui(1,:,:) = one
         p = one
         WHERE(xpi > lx/2.-.25 .and. xpi < lx/2+.25)
     $        p = p + mach*half*(1-COS(2*twopi*(xpi-(lx/2.-.25))))
     $        *half*(1-COS(twopi*ypi))
         ui(4,:,:) = p/(gamma-one)
      CASE("sharp_nozzle")
         DO i = 1,SIZE(xpi,1)
            DO j = 1,SIZE(xpi,2)
               IF (xpi(i,j) <= 5.7143) THEN
                  area_s(i,j) = 4e-20*xpi(i,j)**25 + 0.5256
                  dadx(i,j) = 1e-18*xpi(i,j)**24
               ELSE
                  area_s(i,j) = 0.861438080711430
                  dadx(i,j) = 0
               ENDIF
            ENDDO
         ENDDO
         ui(1,:,:) = 1.0
         ui(2,:,:) = mach*SQRT(gamma*pin)
         ui(3,:,:) = mach*SQRT(gamma*pin)*ypi/area_s*dadx
         ui(4,:,:) = pin/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("cosine_IC")
         ui(1,:,:) = 1.0
         ui(2,:,:) = 0.5*mach*SQRT(gamma*pin)*(COS(ypi/area*pi)+1)
         ui(3,:,:) = 0.0
         ui(4,:,:) = pin/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("parabolic_IC")
         ui(1,:,:) = 1.0
         ui(2,:,:) = mach*SQRT(gamma*pin)
         ui(3,:,:) = 0.0
         ui(4,:,:) = pin/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("pseudo1D_IC")
         ui(1,:,:) = 1.8396/((1.903 + 0.155*TANH(0.8*xpi - 3.67))*area)
         ui(2,:,:) = 1.8396/area
         ui(3,:,:) = 0.0
         ui(4,:,:) = 3.317/(1.428 + 0.429*TANH(0.8*xpi - 3.97))
      CASE("super_subsonic")
         ui(1,:,:) = 1.0
         ui(2,:,:) = mach*SQRT(gamma*pin)
         ui(3,:,:) = 0.0
         ui(4,:,:) = (pin + (pout-pin)*xpi/lx)/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("no_outflow")
         ui(1,:,:) = 1.0
         ui(2,:,:) = 0.0
         ui(3,:,:) = 0.0
         ui(4,:,:) = pin/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("source_sink")
         ui(1,:,:) = 1.0
         ui(2,:,:) = 0.0
         ui(3,:,:) = 0.0
         ui(4,:,:) = pin/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      CASE("circle")
         ui(1,:,:) = 1.0
         ui(2,:,:) = 0.0
         ui(3,:,:) = 0.0
         ui(4,:,:) = pin/(gamma-1)
     $        + 0.5*(ui(2,:,:)**2 + ui(3,:,:)**2)/ui(1,:,:)
      END SELECT

c-----------------------------------------------------------------------
c     add perturbation
c-----------------------------------------------------------------------
      ! none.

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
      USE euler2D_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
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
      CASE("supersonic")
         right%bc_type(1:4)="natural"

         top%bc_type(1)="natural"
         top%static(2:3)=.TRUE.
         top%bc_type(4)="zeroflux"
         bottom%bc_type(1)="natural"
         bottom%static(2:3)=.TRUE.
         bottom%bc_type(4)="zeroflux"

      CASE("zero_normal","cosine_IC","sharp_nozzle","parabolic_IC")

         right%bc_type(1)="natural"
         right%static(2:4)=.TRUE.

         top%bc_type(1)="natural"
         top%static(2:3)=.TRUE.
         top%bc_type(4)="zeroflux"
         bottom%bc_type(1)="natural"
         bottom%static(2:3)=.TRUE.
         bottom%bc_type(4)="zeroflux"

      CASE("open_right","open_nozzle","open_channel","open_nozzle2")
         IF(init_type /= "open_right")THEN
            left%bc_type(1:4)="normflux"
         ENDIF
         
         right%bc_type(1:4)="normflux"

         top%bc_type(1)="zeroflux"
         top%bc_type(4)="zeroflux"
         top%static(2:3)=.TRUE.

         bottom%bc_type(1)="zeroflux"
         bottom%bc_type(4)="zeroflux"
         bottom%static(2:3)=.TRUE.
      CASE("phump")
         top%bc_type="periodic"
         bottom%bc_type="periodic"

      CASE("pseudo1D_IC")

         right%bc_type(1)="natural"
         right%bc_type(2)="natural"
         right%bc_type(3:4)="natural"

         top%bc_type(1)="natural"
         top%static(2:3)=.TRUE.
         top%bc_type(4)="zeroflux"
         bottom%bc_type(1)="natural"
         bottom%static(2:3)=.TRUE.
         bottom%bc_type(4)="zeroflux"

      CASE("super_subsonic")
         left%bc_type(1:4)="robin"
         left%static(1:4)=.TRUE.

         right%bc_type(1)="natural"
         right%bc_type(2)="robin"
         right%static(2)=.TRUE.
         right%bc_type(3:4)="robin"
         right%static(3:4)=.TRUE.

         top%bc_type(1)="natural"
         top%bc_type(2)="robin"
         top%bc_type(3)="robin"
         top%bc_type(4)="zeroflux"
         top%static(2:3)=.TRUE.
         bottom%bc_type(1)="natural"
         bottom%bc_type(2)="robin"
         bottom%bc_type(3)="robin"
         bottom%bc_type(4)="zeroflux"
         bottom%static(2:3)=.TRUE.

      CASE("no_outflow")

         left%static(1)=.TRUE.
         left%static(3:4)=.TRUE.

         right%bc_type(1)="natural"
         right%static(2:3)=.TRUE.
         right%bc_type(4)="zeroflux"
         top%bc_type(1)="natural"
         top%static(2:3)=.TRUE.
         top%bc_type(4)="zeroflux"
         bottom%bc_type(1)="natural"
         bottom%static(2:3)=.TRUE.
         bottom%bc_type(4)="zeroflux"

      CASE("source_sink","circle")
         left%bc_type(1)="natural"
         left%static(2:3)=.TRUE.
         left%bc_type(4)="zeroflux"
         right%bc_type(1)="natural"
         right%static(2:3)=.TRUE.
         right%bc_type(4)="zeroflux"
         top%bc_type(1)="natural"
         top%static(2:3)=.TRUE.
         top%bc_type(4)="zeroflux"
         bottom%bc_type(1)="natural"
         bottom%static(2:3)=.TRUE.
         bottom%bc_type(4)="zeroflux"
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_boundary
c-----------------------------------------------------------------------
c     subprogram e. physics_edge_rhs.
c     computes rhs for edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_rhs(lrtb,t,x,y,nhat,u,ux,uy,uxx,uyy,uxy,c)
      USE euler2D_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c

      INTEGER :: i
      REAL(r8), PARAMETER :: beta=1.5
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: p,ramp,theta
      REAL(r8), DIMENSION(2,SIZE(x,1),SIZE(x,2)) :: vi,vi_x,vi_y
      REAL(r8), DIMENSION(2,2,SIZE(x,1),SIZE(x,2)) :: rot
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c  = 0
c-----------------------------------------------------------------------
c     boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("supersonic")
         SELECT CASE(lrtb)
         CASE("left")
         CASE("right")
            c(2,:,:) = 2*u(2,:,:)*ux(1,:,:)**2/u(1,:,:)**3
     $           - 2*ux(1,:,:)*ux(2,:,:)/u(1,:,:)**2
     $           - u(2,:,:)*uxx(1,:,:)/u(1,:,:)**2
     $           + uxx(2,:,:)/u(1,:,:)
            c(3,:,:) = ux(3,:,:)
            c(4,:,:) = ux(4,:,:)
         CASE("top","bottom")
            c(3,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
            c(2,:,:) = nhat(1,:,:)**2*ux(3,:,:) 
     $           - nhat(2,:,:)**2*uy(2,:,:)
     $           + nhat(1,:,:)*nhat(2,:,:)*(-ux(2,:,:) + uy(3,:,:))
         END SELECT

      CASE("zero_normal","cosine_IC","sharp_nozzle")
         SELECT CASE(lrtb)
         CASE("left")
         CASE("right")
            c(2,:,:) = 2*u(2,:,:)*ux(1,:,:)**2/u(1,:,:)**3
     $           - 2*ux(1,:,:)*ux(2,:,:)/u(1,:,:)**2
     $           - u(2,:,:)*uxx(1,:,:)/u(1,:,:)**2
     $           + uxx(2,:,:)/u(1,:,:)
            c(3,:,:) = ux(3,:,:)
            c(4,:,:) = ux(4,:,:)
         CASE("top","bottom")
            c(3,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
            c(2,:,:) = nhat(1,:,:)**2*ux(3,:,:)
     $           - nhat(2,:,:)**2*uy(2,:,:)
     $           + nhat(1,:,:)*nhat(2,:,:)*(-ux(2,:,:) + uy(3,:,:))
         END SELECT

      CASE("open_right","open_channel","open_nozzle","open_nozzle2")
         theta = ATAN2(nhat(2,:,:),nhat(1,:,:)) 
         rot(1,1,:,:) =  COS(theta)
         rot(1,2,:,:) =  SIN(theta)
         rot(2,1,:,:) = -SIN(theta)
         rot(2,2,:,:) =  COS(theta)
         WHERE(ABS(rot) < eps) rot = zero
         IF(init_type /= "open_right")THEN
            vin = zero
            IF(t < alpha)THEN
               pin = pout + (pmax - pout)
     $              *half*(1. - COS(pi*t/alpha))
               rhoin = rhoout + (rhomax - rhoout)
     $              *half*(1. - COS(pi*t/alpha))
            ELSEIF(t >= 2*alpha . AND. t < 3*alpha)THEN
               pout = p0 + (poutmax - p0)
     $              *half*(1. - COS(pi*t/alpha))
            ELSEIF(t >= 3*alpha)THEN
               pout = poutmax
            ENDIF
         ENDIF
         SELECT CASE(lrtb)
         CASE("left")
            IF(init_type /= "open_right")
     $           CALL euler2D_openbc(lrtb,u,ux,uy,nhat,rot,c)
         CASE("right")
            CALL euler2D_openbc(lrtb,u,ux,uy,nhat,rot,c)
         CASE("top","bottom")
            DO i=1,2
               vi(i,:,:)=u(i+1,:,:)/u(1,:,:)
               vi_x(i,:,:)=ux(i+1,:,:)/u(1,:,:)
     $              -u(i+1,:,:)*ux(1,:,:)/u(1,:,:)**2
               vi_y(i,:,:)=uy(i+1,:,:)/u(1,:,:)
     $              -u(i+1,:,:)*uy(1,:,:)/u(1,:,:)**2
            ENDDO
            c(2,:,:) = nhat(1,:,:)**2*vi_x(2,:,:)
     $           - nhat(2,:,:)**2*vi_y(1,:,:)
     $           + nhat(1,:,:)*nhat(2,:,:)*(-vi_x(1,:,:) + vi_y(2,:,:))
            c(3,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
         END SELECT
      CASE("phump")
         theta = ATAN2(nhat(2,:,:),nhat(1,:,:))
         rot(1,1,:,:) =  COS(theta)
         rot(1,2,:,:) =  SIN(theta)
         rot(2,1,:,:) = -SIN(theta)
         rot(2,2,:,:) =  COS(theta)
         SELECT CASE(lrtb)
         CASE("left","right")
            CALL euler2D_nrbc(u,ux,uy,nhat,rot,c)
         END SELECT

      CASE("parabolic_IC")
         SELECT CASE(lrtb)
         CASE("left")
         CASE("right")
            c(2,:,:) = 2*u(2,:,:)*ux(1,:,:)**2/u(1,:,:)**3
     $           - 2*ux(1,:,:)*ux(2,:,:)/u(1,:,:)**2
     $           - u(2,:,:)*uxx(1,:,:)/u(1,:,:)**2
     $           + uxx(2,:,:)/u(1,:,:)
            c(3,:,:) = ux(3,:,:)
            c(4,:,:) = ux(4,:,:)
         CASE("top","bottom")
            c(3,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
            c(2,:,:) = nhat(1,:,:)**2*ux(3,:,:)
     $           - nhat(2,:,:)**2*uy(2,:,:)
     $           + nhat(1,:,:)*nhat(2,:,:)*(-ux(2,:,:) + uy(3,:,:))
         END SELECT

      CASE("pseudo1D_IC")
         SELECT CASE(lrtb)
         CASE("left")
         CASE("right")
         CASE("top","bottom")
            c(3,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
            c(2,:,:) = nhat(1,:,:)**2*ux(3,:,:)
     $           - nhat(2,:,:)**2*uy(2,:,:)
     $           + nhat(1,:,:)*nhat(2,:,:)*(-ux(2,:,:) + uy(3,:,:))
         END SELECT

      CASE("super_subsonic")
         p = (gamma-1)*(u(4,:,:) - 0.5*(u(2,:,:)**2
     $        + u(3,:,:)**2)/u(1,:,:))
         SELECT CASE(lrtb)
         CASE("left")
            c(1,:,:) = u(1,:,:) - 1.0
            c(2,:,:) = u(2,:,:) - mach*SQRT(gamma*pin)
            c(3,:,:) = u(3,:,:)
            c(4,:,:) = u(4,:,:) - ( pin/(gamma-1)
     $           + 0.5*( mach**2*gamma*pin ) )
         CASE("right")
            c(2,:,:) = 2*u(2,:,:)*ux(1,:,:)**2/u(1,:,:)**3
     $           - 2*ux(1,:,:)*ux(2,:,:)/u(1,:,:)**2
     $           - u(2,:,:)*uxx(1,:,:)/u(1,:,:)**2
     $           + uxx(2,:,:)/u(1,:,:)
            c(4,:,:) = p - pout
         CASE("top","bottom")
            c(3,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
            c(2,:,:) = nhat(1,:,:)**2*ux(3,:,:)
     $           - nhat(2,:,:)**2*uy(2,:,:)
     $           + nhat(1,:,:)*nhat(2,:,:)*(-ux(2,:,:) + uy(3,:,:))
         END SELECT
      CASE("no_outflow")
         SELECT CASE(lrtb)
         CASE("left")
            IF (t <= alpha) THEN
               ramp=0.5*(1-COS(pi*t/alpha))*0.5*(COS(y/0.5256*pi)+1)
            ELSE
               ramp=0.5*(COS(y/0.5256*pi)+1)
            ENDIF
            c(1,:,:) = u(1,:,:) - (1.0*(1+(beta-1)*ramp))
            c(2,:,:) = u(2,:,:)**2*ux(1,:,:)/u(1,:,:)**2 
     $           - 2*u(2,:,:)*ux(2,:,:)/u(1,:,:)
     $           -(gamma-1)*ux(4,:,:)
     $           -(gamma-1)*0.5*u(2,:,:)**2*ux(1,:,:)/u(1,:,:)**2
     $           +(gamma-1)*u(2,:,:)*ux(2,:,:)/u(1,:,:)
            c(3,:,:) = u(3,:,:)
            c(4,:,:) = u(4,:,:) - ( pin*(1+(beta-1)*ramp)/(gamma-1)
     $           + 0.5*(u(2,:,:)**2 + u(3,:,:)**2)/u(1,:,:) )
         CASE("right")
            c(2,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
            c(3,:,:) = -nhat(1,:,:)*u(3,:,:) + nhat(2,:,:)*u(2,:,:)
         CASE("top","bottom")
            c(3,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
            c(2,:,:) = -nhat(1,:,:)*u(3,:,:) + nhat(2,:,:)*u(2,:,:)
         END SELECT

      CASE("source_sink","circle")
         c(2,:,:) = nhat(1,:,:)*u(2,:,:) + nhat(2,:,:)*u(3,:,:)
         c(3,:,:) = nhat(1,:,:)**2*ux(3,:,:) - nhat(2,:,:)**2*uy(2,:,:)
     $        +nhat(1,:,:)*nhat(2,:,:)*(-ux(2,:,:) + uy(3,:,:))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_edge_rhs
c-----------------------------------------------------------------------
c     subprogram f. physics_edge_drdu.
c     computes drdu for edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_drdu(lrtb,t,x,y,nhat,
     $     u,ux,uy,uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)
      USE euler2D_mod
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

      CALL euler2D_bc_jac(lrtb,t,x,y,nhat,u,ux,uy,c_u,c_ux,c_uy)
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
      USE euler2D_mod
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
      CASE("supersonic","zero_normal","cosine_IC","sharp_nozzle",
     $     "parabolic_IC","pseudo1D_IC","open_right")
         SELECT CASE(lrtb)
         CASE("left")
            DO iqty=1,4
               mass(iqty,iqty,:,:)=one
            ENDDO
         END SELECT
      CASE("no_outflow")
         SELECT CASE(lrtb)
         CASE("left")
            mass(2,2,:,:)=one
         END SELECT
      CASE("phump")
         SELECT CASE(lrtb)
         CASE("left","right")
            DO iqty=1,4
               mass(iqty,iqty,:,:)=one
            ENDDO
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
      USE euler2D_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8) :: ramp
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: p,sshape

c-----------------------------------------------------------------------
c     initialize
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
      p=(gamma-1)*(u(4,:,:) - 0.5*(u(2,:,:)**2 + u(3,:,:)**2)/u(1,:,:))
c-----------------------------------------------------------------------
c     flux and source terms
c-----------------------------------------------------------------------
      fx(1,:,:) = u(2,:,:) - ddiff*ux(1,:,:)
      fy(1,:,:) = u(3,:,:) - ddiff*uy(1,:,:)

      IF (init_type=="source_sink" .OR. init_type=="circle") THEN
         IF (t <= alpha) THEN
            ramp=0.5*(1-COS(2*pi*t/alpha))
         ELSE
            ramp=0.0
         ENDIF
         sshape = 0.10*ramp*(EXP(-((x-x1)**2 + (y-y1)**2)/rs**2)
     $        - EXP(-((x-x2)**2 + (y-y2)**2)/rs**2))
         s(1,:,:) = sshape
         s(2,:,:) = 0.0
         s(3,:,:) = 0.0
         s(4,:,:) = 1.5*p*sshape/u(1,:,:)
      ENDIF

      fx(2,:,:) = u(2,:,:)**2/u(1,:,:) + p 
     $     - mu*(ux(2,:,:)/u(1,:,:) - u(2,:,:)*ux(1,:,:)/u(1,:,:)**2)
      fy(2,:,:) = u(2,:,:)*u(3,:,:)/u(1,:,:)
     $     - mu*(uy(2,:,:)/u(1,:,:) - u(2,:,:)*uy(1,:,:)/u(1,:,:)**2)

      fx(3,:,:) = u(2,:,:)*u(3,:,:)/u(1,:,:)
     $     - mu*(ux(3,:,:)/u(1,:,:) - u(3,:,:)*ux(1,:,:)/u(1,:,:)**2)
      fy(3,:,:) = u(3,:,:)**2/u(1,:,:) + p 
     $     - mu*(uy(3,:,:)/u(1,:,:) - u(3,:,:)*uy(1,:,:)/u(1,:,:)**2)

      fx(4,:,:) = u(2,:,:)/u(1,:,:)*(u(4,:,:) + p)
     $     - u(2,:,:)/u(1,:,:)**3
     $     * mu*(ux(2,:,:)*u(1,:,:) - u(2,:,:)*ux(1,:,:))
     $     - u(3,:,:)/u(1,:,:)**3
     $     * mu*(ux(3,:,:)*u(1,:,:) - u(3,:,:)*ux(1,:,:))
     $     - ediff*ux(4,:,:)

      fy(4,:,:) = u(3,:,:)/u(1,:,:)*(u(4,:,:) + p)
     $     - u(2,:,:)/u(1,:,:)**3
     $     * mu*(uy(2,:,:)*u(1,:,:) - u(2,:,:)*uy(1,:,:))
     $     - u(3,:,:)/u(1,:,:)**3
     $     * mu*(uy(3,:,:)*u(1,:,:) - u(3,:,:)*uy(1,:,:))
     $     - ediff*uy(4,:,:)
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
      SUBROUTINE physics_drdu(t,x,y,u,ux,uy,fx_u,fx_ux,fx_uy,fy_u,
     $     fy_ux,fy_uy,s_u,s_ux,s_uy)
      USE euler2D_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8) :: ramp,length
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: p,sshape

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
      p=(gamma-1)*(u(4,:,:) - 0.5*(u(2,:,:)**2 + u(3,:,:)**2)/u(1,:,:))
c-----------------------------------------------------------------------
c     2D Euler flux Jacobian
c-----------------------------------------------------------------------

      fx_u(1,2,:,:) = 1.0
      fy_u(1,3,:,:) = 1.0
      fx_ux(1,1,:,:) = -ddiff
      fy_uy(1,1,:,:) = -ddiff

      IF (init_type=="source_sink" .OR. init_type=="circle") THEN
         IF (t <= alpha) THEN
            ramp=0.5*(1-COS(2*pi*t/alpha))
         ELSE
            ramp=0.0
         ENDIF
         sshape = 0.10*ramp*(EXP(-((x-x1)**2 + (y-y1)**2)/rs**2)
     $        - EXP(-((x-x2)**2 + (y-y2)**2)/rs**2))
         s_u(1,1,:,:) = 0.0
         s_u(2,2,:,:) = 0.0
         s_u(3,3,:,:) = 0.0
         s_u(4,1,:,:) = 1.5*sshape*(gamma-1)*(-u(4,:,:)
     $        + (u(2,:,:)**2+u(3,:,:)**2)/u(1,:,:))/u(1,:,:)**2
         s_u(4,2,:,:) = -1.5*sshape*(gamma-1)*u(2,:,:)/u(1,:,:)**2
         s_u(4,3,:,:) = -1.5*sshape*(gamma-1)*u(3,:,:)/u(1,:,:)**2
         s_u(4,4,:,:) = 1.5*sshape*(gamma-1)/u(1,:,:)
      END IF

      fx_u(2,1,:,:) = -u(2,:,:)**2/u(1,:,:)**2 + 0.5*(gamma-1)
     $     *(u(2,:,:)**2/u(1,:,:)**2 + u(3,:,:)**2/u(1,:,:)**2)
     $     + mu*(ux(2,:,:)/u(1,:,:)**2
     $     - 2*u(2,:,:)*ux(1,:,:)/u(1,:,:)**3)
      fx_u(2,2,:,:) = (3-gamma)*u(2,:,:)/u(1,:,:)
     $     + mu*ux(1,:,:)/u(1,:,:)**2
      fx_u(2,3,:,:) = -(gamma-1)*u(3,:,:)/u(1,:,:)
      fx_u(2,4,:,:) = (gamma-1)
      fx_ux(2,1,:,:) = mu*u(2,:,:)/u(1,:,:)**2
      fx_ux(2,2,:,:) = -mu/u(1,:,:)

      fy_u(2,1,:,:) = -u(2,:,:)*u(3,:,:)/u(1,:,:)**2
     $     + mu*(uy(2,:,:)/u(1,:,:)**2 
     $     - 2*u(2,:,:)*uy(1,:,:)/u(1,:,:)**3)
      fy_u(2,2,:,:) = u(3,:,:)/u(1,:,:) + mu*uy(1,:,:)/u(1,:,:)**2
      fy_u(2,3,:,:) = u(2,:,:)/u(1,:,:)
      fy_uy(2,1,:,:) = mu*u(2,:,:)/u(1,:,:)**2
      fy_uy(2,2,:,:) = -mu/u(1,:,:)

      fx_u(3,1,:,:) = -u(2,:,:)*u(3,:,:)/u(1,:,:)**2
     $     + mu*(ux(3,:,:)/u(1,:,:)**2 
     $     - 2*u(3,:,:)*ux(1,:,:)/u(1,:,:)**3)
      fx_u(3,2,:,:) = u(3,:,:)/u(1,:,:)
      fx_u(3,3,:,:) = u(2,:,:)/u(1,:,:) + mu*ux(1,:,:)/u(1,:,:)**2
      fx_ux(3,1,:,:) = mu*u(3,:,:)/u(1,:,:)**2
      fx_ux(3,3,:,:) = -mu/u(1,:,:)

      fy_u(3,1,:,:) = 
     $     -u(3,:,:)**2/u(1,:,:)**2 + 0.5*(gamma-1)
     $     *(u(2,:,:)**2/u(1,:,:)**2 + u(3,:,:)**2/u(1,:,:)**2)
     $     + mu*(uy(3,:,:)/u(1,:,:)**2 
     $     - 2*u(3,:,:)*uy(1,:,:)/u(1,:,:)**3)
      fy_u(3,2,:,:) = -(gamma-1)*u(2,:,:)/u(1,:,:)
      fy_u(3,3,:,:) = (3-gamma)*u(3,:,:)/u(1,:,:)
     $     + mu*uy(1,:,:)/u(1,:,:)**2
      fy_u(3,4,:,:) = (gamma-1)
      fy_uy(3,1,:,:) = mu*u(3,:,:)/u(1,:,:)**2
      fy_uy(3,3,:,:) = -mu/u(1,:,:)

      fx_u(4,1,:,:) = 
     $     -u(2,:,:)*u(4,:,:)/u(1,:,:)**2
     $     -(gamma-1)*u(2,:,:)*u(4,:,:)/u(1,:,:)**2
     $     +(gamma-1)*u(2,:,:)*(u(2,:,:)**2+u(3,:,:)**2)/u(1,:,:)**3
     $     +2*mu*u(2,:,:)*ux(2,:,:)/u(1,:,:)**3
     $     +2*mu*u(3,:,:)*ux(3,:,:)/u(1,:,:)**3
     $     -3*mu*u(2,:,:)**2*ux(1,:,:)/u(1,:,:)**4
     $     -3*mu*u(3,:,:)**2*ux(1,:,:)/u(1,:,:)**4

      fx_u(4,2,:,:) = 
     $     +u(4,:,:)/u(1,:,:)
     $     +(gamma-1)*u(4,:,:)/u(1,:,:)
     $     -(gamma-1)*3./2.*u(2,:,:)**2/u(1,:,:)**2
     $     -0.5*(gamma-1)*u(3,:,:)**2/u(1,:,:)**2
     $     -mu*ux(2,:,:)/u(1,:,:)**2
     $     +2*mu*u(2,:,:)*ux(1,:,:)/u(1,:,:)**3

      fx_u(4,3,:,:) =
     $     -(gamma-1)*u(2,:,:)*u(3,:,:)/u(1,:,:)**2
     $     -mu*ux(3,:,:)/u(1,:,:)**2
     $     +mu*2*ux(1,:,:)*u(3,:,:)/u(1,:,:)**3

      fx_u(4,4,:,:) = gamma*u(2,:,:)/u(1,:,:)
      fx_ux(4,1,:,:) = mu*(u(2,:,:)**2 + u(3,:,:)**2)/u(1,:,:)**3
      fx_ux(4,2,:,:) = -mu*u(2,:,:)/u(1,:,:)**2
      fx_ux(4,3,:,:) = -mu*u(3,:,:)/u(1,:,:)**2
      fx_ux(4,4,:,:) = -ediff

      fy_u(4,1,:,:) = 
     $     -u(3,:,:)*u(4,:,:)/u(1,:,:)**2
     $     -(gamma-1)*u(3,:,:)*u(4,:,:)/u(1,:,:)**2
     $     +(gamma-1)*u(3,:,:)*(u(2,:,:)**2+u(3,:,:)**2)/u(1,:,:)**3
     $     +2*mu*u(2,:,:)*uy(2,:,:)/u(1,:,:)**3
     $     +2*mu*u(3,:,:)*uy(3,:,:)/u(1,:,:)**3
     $     -3*mu*u(2,:,:)**2*uy(1,:,:)/u(1,:,:)**4
     $     -3*mu*u(3,:,:)**2*uy(1,:,:)/u(1,:,:)**4

      fy_u(4,2,:,:) =
     $     -(gamma-1)*u(3,:,:)*u(2,:,:)/u(1,:,:)**2
     $     -mu*uy(2,:,:)/u(1,:,:)**2
     $     +mu*2*uy(1,:,:)*u(2,:,:)/u(1,:,:)**3

      fy_u(4,3,:,:) =
     $     +u(4,:,:)/u(1,:,:)
     $     +(gamma-1)*u(4,:,:)/u(1,:,:)
     $     -(gamma-1)*3./2.*u(3,:,:)**2/u(1,:,:)**2
     $     -0.5*(gamma-1)*u(2,:,:)**2/u(1,:,:)**2
     $     -mu*uy(3,:,:)/u(1,:,:)**2
     $     +mu*2*u(3,:,:)*uy(1,:,:)/u(1,:,:)**3

      fy_u(4,4,:,:) = gamma*u(3,:,:)/u(1,:,:)
      fy_uy(4,1,:,:) = mu*(u(2,:,:)**2 + u(3,:,:)**2)/u(1,:,:)**3
      fy_uy(4,2,:,:) = -mu*u(2,:,:)/u(1,:,:)**2
      fy_uy(4,3,:,:) = -mu*u(3,:,:)/u(1,:,:)**2
      fy_uy(4,4,:,:) = -ediff
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_drdu
c-----------------------------------------------------------------------
c     subprogram l. physics_mass.
c     computes mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_mass(x,y,mass,mass_x,mass_y)
      USE euler2D_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     no modifications needed.
c-----------------------------------------------------------------------
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
      SUBROUTINE physics_grid(x,y,ksi,eta)
      USE euler2D_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,eta

      REAL(r8), PARAMETER :: c=1.45888 
      INTEGER :: ix,iy
      REAL(r8) :: a1,b1,f1,f2,d1,d2
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: xhat,Lhat
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("supersonic","zero_normal","cosine_IC","sharp_nozzle",
     $     "parabolic_IC","pseudo1D_IC","super_subsonic","no_outflow",
     $     "source_sink","open_right")
         ksi=lx*x
         eta=(2.*y-one)
     $           *(1. + expand*half*(TANH(0.8*ksi - 3.) - TANH(-3.)))
      CASE("phump")
         ksi=lx*x
         eta=y
      CASE("open_channel")
         ksi=lx*x
         eta=(y-0.5)
      CASE("open_nozzle2")
         lx = Ldiv + Lconv
         d1 = Ldiv/4.

         WHERE(x  > Lconv/lx)
            ksi = lx*x - Lconv
     $           + (one - ABS(two*y - one)**2)*d1*TAN(pi/8.)
     $           *(x*lx - Lconv)/Ldiv
         ELSEWHERE
c            xhat = x*lx/Lconv
c            xhat = (xhat**3 + gr_curve*xhat)/(one + gr_curve)
c            ksi = xhat*Lconv - Lconv
            ksi = lx*x - Lconv
         ENDWHERE
         xhat = -(ksi + SIN(Aconv))
         d1 = Lconv - SIN(Aconv)
         WHERE(ksi < -SIN(Aconv))
c            eta=(2.*y - 1.)*(2. - COS(Aconv)
c     $           - TAN(Aconv)*(ksi + SIN(Aconv)))
            eta = (2.*y - 1.)*(2. - COS(Aconv)
     $           + TAN(Aconv)*half*(xhat + d1/pi*(SIN(pi*xhat/d1))))
         ELSEWHERE(ksi < 0.)
            eta=(2.*y - 1.)*(2. - COS(ASIN(ABS(ksi))))
         ELSEWHERE
            eta=(2.*y - 1.)
     $           *(1. + expand*half*(TANH(0.8*ksi - 3.) - TANH(-3.)))
         ENDWHERE
      CASE("open_nozzle")
         lx = Ldiv + Lconv
         d1 = (2. - COS(Adiv))
     $        + TAN(Adiv)*(Ldiv - SIN(Adiv))
         WHERE(x  > Lconv/lx)
            ksi = lx*x - Lconv
     $           + (one - ABS(two*y - one)**2)*d1*TAN(pi/8.)
     $           *(x*lx - Lconv)/Ldiv
         ELSEWHERE
            ksi = lx*x - Lconv
         ENDWHERE
         WHERE(ksi < -SIN(Aconv))
            eta=(2.*y - 1.)*(2. - COS(Aconv)
     $           - TAN(Aconv)*(ksi + SIN(Aconv)))
         ELSEWHERE(ksi < SIN(Adiv))
            eta=(2.*y - 1.)*(2. - COS(ASIN(ABS(ksi))))
         ELSEWHERE
            eta=(2.*y - 1.)*(2. - COS(Adiv)
     $           + TAN(Adiv)*(ksi - SIN(Adiv)))
         ENDWHERE
      CASE("circle")
         DO ix=1,SIZE(x,1)
            DO iy=1,SIZE(x,2)

               a1=two*ABS(x(ix,iy)-half)
               b1=two*ABS(y(ix,iy)-half)

               d1=(SQRT((one+one/TAN(half*pi*a1))**2 
     $              - (COS(.25*pi*a1))**2)-SIN(.25*pi*a1))
               d2=(SQRT((one+one/TAN(half*pi*b1))**2 
     $              - (COS(.25*pi*b1))**2)-SIN(.25*pi*b1))

               f1=-(COS(half*pi*a1) 
     $              + SQRT((one+SIN(pi*a1))/(COS(.25*pi*a1))**2 
     $              - (SIN(half*pi*a1))**2))
               f2=-(COS(half*pi*b1) 
     $              + SQRT((one+SIN(pi*b1))/(COS(.25*pi*b1))**2 
     $              - (SIN(half*pi*b1))**2))

               IF(a1==one)THEN
                  ksi(ix,iy)=COS(.25*pi*b1)
                  eta(ix,iy)=SIN(.25*pi*b1)
               ELSEIF(b1==one)THEN
                  ksi(ix,iy)=SIN(.25*pi*a1)
                  eta(ix,iy)=COS(.25*pi*a1)
               ELSEIF(a1==zero .AND. b1==zero)THEN
                  ksi(ix,iy)=zero
                  eta(ix,iy)=zero
               ELSEIF(a1==zero)THEN
                  ksi(ix,iy)=zero
                  eta(ix,iy)=(one+one/TAN(half*pi*b1))-d2
               ELSEIF(b1==zero)THEN
                  ksi(ix,iy)=(one+one/TAN(half*pi*a1))-d1
                  eta(ix,iy)=zero
               ELSE
                  ksi(ix,iy)=(SQRT(4.*d1**2*(one-f1/d1**2-f2/d2**2)
     $                 - (f1-f2)**2/d2**2)-(two*d1+d1/d2**2*(f1-f2)))
     $                 /(two*(one+d1**2/d2**2))
                  eta(ix,iy)=ksi(ix,iy)*d1/d2 + half*(f1-f2)/d2
               ENDIF
               ksi(ix,iy)=ksi(ix,iy)*SIGN(one,x(ix,iy)-half)
               eta(ix,iy)=eta(ix,iy)*SIGN(one,y(ix,iy)-half)
            ENDDO
         ENDDO
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_grid
c-----------------------------------------------------------------------
c     subprogram n. physics_schur.
c     computes Schur complement.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_schur(t,hfac,x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
      USE euler2D_mod
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
c     deallocate variables
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_dealloc
      USE euler2D_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     terminate.
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
