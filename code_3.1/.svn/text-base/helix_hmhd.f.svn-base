c-----------------------------------------------------------------------
c     file helix_hmhd.f.
c     contains specifications for RMHD equations with helical symmetry.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. helix_hmhd_mod.
c     1. helix_hmhd_init.
c     2. helix_hmhd_boundary.
c     3. helix_hmhd_bound_rhs.
c     4. helix_hmhd_bound_drdu.
c     5. helix_hmhd_init_special.
c     6. helix_hmhd_rhs.
c     7. helix_hmhd_drdu.
c     8. helix_hmhd_mass.
c     9. helix_hmhd_equil.
c     10. helix_hmhd_grid.
c-----------------------------------------------------------------------
c     subprogram 0. helix_hmhd_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE helix_hmhd_mod
      USE local_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. helix_hmhd_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_init(ksi,eta,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: ksi,eta
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: ix,iy
      REAL(r8), DIMENSION(0:SIZE(ksi,1)-1,0:SIZE(ksi,2)-1) :: rad,
     $     cos_theta
      REAL(r8), DIMENSION(SIZE(ui,1),0:SIZE(ui,2)-1,0:SIZE(ui,3)-1) :: 
     $     ux,uy
c-----------------------------------------------------------------------
c     add perturbation to an equilibrium
c-----------------------------------------------------------------------
      rad=SQRT(ksi**2 + eta**2)
      DO ix=0,SIZE(rad,1)-1
         DO iy=0,SIZE(rad,2)-1
            IF(rad(ix,iy) /= 0)THEN
               cos_theta(ix,iy)=ksi(ix,iy)/rad(ix,iy)
            ELSE
               cos_theta(ix,iy)=1
            ENDIF
         ENDDO
      ENDDO
      CALL helix_hmhd_equil(ksi,eta,ui,ux,uy,.FALSE.)
      
      SELECT CASE(init_type)
      CASE("m1_park","m1_axial")
         ui(2,:,:)=ui(2,:,:) + epsilon*cos_theta
     $        *EXP(-((rad-r_s)/gr_curve)**2)
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_init
c-----------------------------------------------------------------------
c     subprogram 2. helix_hmhd_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_boundary(left,right,top,bottom,nqty)

      INTEGER, INTENT(IN) :: nqty
      TYPE(bc2_type) :: left,right,top,bottom

      INTEGER :: i,j
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      right%a=0
      right%b=0

      left%bc_type="natural"
      SELECT CASE(init_type)
      CASE("m1_park")
         right%bc_type(1)="robin"
         right%a(1,1)=1
         right%static(2:4)=.TRUE.
         right%bc_type(5)="natural"
         right%static(6)=.TRUE.
      CASE("m1_axial")
         right%bc_type(1:2)="robin"
         right%a(1,1)=1
         right%a(2,1)=Rinv
         right%static(3:4)=.TRUE.
         right%bc_type(5)="natural"
         right%static(6)=.TRUE.
      END SELECT
      top%bc_type = "periodic"
      bottom%bc_type = "periodic"
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_boundary
c-----------------------------------------------------------------------
c     subprogram 3. helix_hmhd_bound_rhs.
c     computes rhs for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_bound_rhs(edge_type,x,y,u,ux,uy,c)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: c

      REAL(r8) :: gfac,gfac_s,Eb,eta_local
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      c=0
      gfac=1/(1+Rinv**2)
      gfac_s=1/(1+(r_s*Rinv)**2)

      SELECT CASE(init_type)
      CASE("m1_park")
         Eb=eta*(r_s**2/gfac_s - 2*gfac_s)
         eta_local=ABS(Eb/((2 - r_s**2)/gfac - 2*gfac))
         IF(edge_type == "right")THEN
            IF(.NOT. source)c(1,:)=Eb
            IF(eta /= 0 .OR. di /= 0)THEN
               IF(source)THEN
                  c(2,:)=eta*(x*ux(2,:) + y*uy(2,:)) 
     $                 + di*u(2,:)*(x*uy(2,:) - y*ux(2,:))
               ELSE
                  c(2,:)=eta_local*(x*ux(2,:) + y*uy(2,:)) 
     $                 + di*u(2,:)*(x*uy(2,:) - y*ux(2,:))
               ENDIF
            ELSE
               c(2,:)=x*ux(2,:) + y*uy(2,:)
            ENDIF
            c(3,:)=u(3,:) + 2*Rinv*gfac*u(4,:) 
     $           - gfac*(x*ux(6,:) + y*uy(6,:))
            c(4,:)=x*ux(4,:) + y*uy(4,:) - gfac*Rinv**2*u(4,:)
            c(6,:)=u(6,:)
         ENDIF
      CASE("m1_axial")
         Eb=eta*(r_s**2/gfac_s - 2)
         eta_local=ABS(Eb/((1+2*Rinv**2)*(1-r_s**2) + 1/gfac - 2))
         IF(edge_type == "right")THEN
            c(1,:)=Eb/gfac + Rinv*(di*u(2,:)*(x*uy(2,:) - y*ux(2,:)) 
     $           + eta_local*(x*ux(2,:) + y*uy(2,:)))
            c(2,:)=-eta_local*(x*ux(2,:) + y*uy(2,:)) 
     $           - di*u(2,:)*(x*uy(2,:) - y*ux(2,:))
            c(3,:)=u(3,:) + 2*Rinv*gfac*u(4,:) 
     $           - gfac*(x*ux(6,:) + y*uy(6,:))
            c(4,:)=x*ux(4,:) + y*uy(4,:) - gfac*Rinv**2*u(4,:)
            c(6,:)=u(6,:)
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_bound_rhs
c-----------------------------------------------------------------------
c     subprogram 4. helix_hmhd_bound_drdu.
c     computes drdu for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_bound_drdu(edge_type,x,y,u,ux,uy,
     $     c_u,c_ux,c_uy)

      CHARACTER(*), INTENT(IN) :: edge_type
      REAL(r8), DIMENSION(:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy

      REAL(r8) :: gfac,gfac_s,Eb,eta_local
c-----------------------------------------------------------------------
c     zero out output.
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      gfac=1/(1+Rinv**2)
      gfac_s=1/(1+(r_s*Rinv)**2)

      SELECT CASE(init_type)
      CASE("m1_park")
         Eb=eta*(r_s**2/gfac_s - 2*gfac_s)
         eta_local=ABS(Eb/((2 - r_s**2)/gfac - 2*gfac))
         IF(edge_type == "right")THEN
            IF(eta /= 0 .OR. di /= 0)THEN
               c_u(2,2,:)=di*(x*uy(2,:) - y*ux(2,:))
               IF(source)THEN
                  c_ux(2,2,:)=eta*x - di*u(2,:)*y
                  c_uy(2,2,:)=eta*y + di*u(2,:)*x
               ELSE
                  c_ux(2,2,:)=eta_local*x - di*u(2,:)*y
                  c_uy(2,2,:)=eta_local*y + di*u(2,:)*x
               ENDIF
            ELSE
               c_ux(2,2,:)=x
               c_uy(2,2,:)=y
            ENDIF
            c_u(3,3,:)=1
            c_u(3,4,:)=2*Rinv*gfac
            c_ux(3,6,:)=-gfac*x
            c_uy(3,6,:)=-gfac*y
            c_u(4,4,:)=-gfac*Rinv**2
            c_ux(4,4,:)=x
            c_uy(4,4,:)=y
            c_u(6,6,:)=1
         ENDIF
      CASE("m1_axial")
         Eb=eta*(r_s**2/gfac_s - 2)
         eta_local=ABS(Eb/((1+2*Rinv**2)*(1-r_s**2) + 1/gfac - 2))
         IF(edge_type == "right")THEN
            c_u(1,2,:)=Rinv*di*(x*uy(2,:) - y*ux(2,:))
            c_ux(1,2,:)=Rinv*(eta_local*x - di*u(2,:)*y)
            c_uy(1,2,:)=Rinv*(eta_local*y + di*u(2,:)*x)
            c_u(2,2,:)=-di*(x*uy(2,:) - y*ux(2,:))
            c_ux(2,2,:)=-eta_local*x + di*u(2,:)*y
            c_uy(2,2,:)=-eta_local*y - di*u(2,:)*x
            c_u(3,3,:)=1
            c_u(3,4,:)=2*Rinv*gfac
            c_ux(3,6,:)=-gfac*x
            c_uy(3,6,:)=-gfac*y
            c_u(4,4,:)=-gfac*Rinv**2
            c_ux(4,4,:)=x
            c_uy(4,4,:)=y
            c_u(6,6,:)=1
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_bound_drdu
c-----------------------------------------------------------------------
c     subprogram 5. helix_hmhd_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_init_special(static,ground)

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground
c-----------------------------------------------------------------------
c     initialize run parameters
c-----------------------------------------------------------------------
      ground(6)=.TRUE.
      static(5:6)=.TRUE.
      IF(init_type == "m1_axial" .AND. eta == 0)
     $     CALL program_stop("Have to have non-zero resistivity!")
      IF(init_type == "m1_axial" .AND. source)
     $     CALL program_stop("Do not run m1_axial with source term!")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_init_special
c-----------------------------------------------------------------------
c     subprogram 6. helix_hmhd_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE helix_hmhd_rhs(x,y,u,ux,uy,fx,fy,s,first)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first

      REAL(r8) :: gfac_s,Eb
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac,xfac,yfac,pfac,
     $     Ve,Vex,Vey,eta_local
      REAL(r8), DIMENSION(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: u0,u0x,u0y,
     $     fx0,fy0,s0
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     calculate auxiliary variables.
c-----------------------------------------------------------------------
      gfac=1/(1+(x**2+y**2)*Rinv**2)
      xfac=1-Rinv**2*x**2*gfac
      yfac=1-Rinv**2*y**2*gfac
      pfac=x*y*Rinv**2*gfac
      Ve=u(4,:,:)-di*u(5,:,:)
      Vex=ux(4,:,:)-di*ux(5,:,:)
      Vey=uy(4,:,:)-di*uy(5,:,:)
c-----------------------------------------------------------------------
c     modify resistivity for boundary applied E-field(no current drive).
c-----------------------------------------------------------------------
      eta_local=eta
      gfac_s=1/(1+(r_s*Rinv)**2)
      SELECT CASE(init_type)
      CASE("m1_park")
         IF(.NOT. source)
     $        eta_local=eta*ABS((r_s**2/gfac_s - 2*gfac_s)
     $        /((2*(x**2+y**2)-r_s**2)/gfac - 2*gfac))
      CASE("m1_axial")
         Eb=eta*(r_s**2/gfac_s - 2)
         eta_local=ABS(Eb/((1+2*Rinv**2*(x**2+y**2))
     $        *(x**2+y**2-r_s**2) + (x**2+y**2)/gfac - 2))
      END SELECT
c-----------------------------------------------------------------------
c     Helical incomressible HMHD model in cartesian coordinates
c-----------------------------------------------------------------------
      fx(1,:,:)=-di*nu*(xfac*Vex - pfac*Vey)
      fy(1,:,:)=-di*nu*(yfac*Vey - pfac*Vex)

      s(1,:,:)=gfac*(ux(1,:,:)*(uy(6,:,:) + di*uy(2,:,:))
     $     - uy(1,:,:)*(ux(6,:,:) + di*ux(2,:,:))
     $     + eta_local*u(5,:,:) - 4*di*nu*Rinv**2*gfac**2*Ve)

      fx(2,:,:)=gfac*(-u(2,:,:)*(uy(6,:,:) + di*uy(2,:,:)) 
     $     + uy(1,:,:)*Ve - 2*di*nu*Rinv*(xfac*Vex - pfac*Vey)
     $     + 4*di*nu*Rinv**3*x*gfac**2*Ve)
     $     - eta_local*(xfac*ux(2,:,:) - pfac*uy(2,:,:))
  
      fy(2,:,:)=gfac*(u(2,:,:)*(ux(6,:,:) + di*ux(2,:,:)) 
     $     - ux(1,:,:)*Ve - 2*di*nu*Rinv*(yfac*Vey - pfac*Vex)
     $     + 4*di*nu*Rinv**3*y*gfac**2*Ve)
     $     - eta_local*(yfac*uy(2,:,:) - pfac*ux(2,:,:))

      fx(3,:,:)=gfac*(uy(1,:,:)*(u(5,:,:) + 2*Rinv*gfac*u(2,:,:))
     $     - uy(6,:,:)*(u(3,:,:) + 2*Rinv*gfac*u(4,:,:)) 
     $     + 2*nu*Rinv*(xfac*Vex - pfac*Vey)
     $     - 4*nu*Rinv**3*x*gfac**2*Ve)
     $     - mu*(xfac*ux(3,:,:) - pfac*uy(3,:,:)) 

      fy(3,:,:)=gfac*(ux(6,:,:)*(u(3,:,:) + 2*Rinv*gfac*u(4,:,:))
     $     - ux(1,:,:)*(u(5,:,:) + 2*Rinv*gfac*u(2,:,:))
     $     + 2*nu*Rinv*(yfac*Vey - pfac*Vex)
     $     - 4*nu*Rinv**3*y*gfac**2*Ve)
     $     - mu*(yfac*uy(3,:,:) - pfac*ux(3,:,:)) 
      
      s(3,:,:)=2*Rinv**2*gfac**2
     $     *(u(2,:,:)*(x*uy(2,:,:) - y*ux(2,:,:))
     $     - u(4,:,:)*(x*uy(4,:,:) - y*ux(4,:,:))
     $     - 4*Rinv*gfac
     $     *(u(2,:,:)*(x*uy(1,:,:) - y*ux(1,:,:))
     $     - u(4,:,:)*(x*uy(6,:,:) - y*ux(6,:,:))))
      
      fx(4,:,:)=-mu*(xfac*ux(4,:,:) - pfac*uy(4,:,:)) 
     $     - nu*(xfac*Vex - pfac*Vey)
      fy(4,:,:)=-mu*(yfac*uy(4,:,:) - pfac*ux(4,:,:)) 
     $     - nu*(yfac*Vey - pfac*Vex)

      s(4,:,:)=gfac*(ux(4,:,:)*uy(6,:,:) - uy(4,:,:)*ux(6,:,:)
     $     + ux(1,:,:)*uy(2,:,:) - uy(1,:,:)*ux(2,:,:)
     $     + 2*mu*Rinv*gfac*u(3,:,:)
     $     - 4*nu*Rinv**2*gfac**2*Ve)
      
      fx(5,:,:)=xfac*ux(1,:,:) - pfac*uy(1,:,:)
      fy(5,:,:)=yfac*uy(1,:,:) - pfac*ux(1,:,:)
      s(5,:,:)=gfac*(u(5,:,:) + 2*Rinv*gfac*u(2,:,:))
      
      fx(6,:,:)=xfac*ux(6,:,:) - pfac*uy(6,:,:)
      fy(6,:,:)=yfac*uy(6,:,:) - pfac*ux(6,:,:)
      s(6,:,:)=gfac*(u(3,:,:) + 2*Rinv*gfac*u(4,:,:))

      IF(source .AND. first)THEN
         SELECT CASE(init_type)
         CASE("m1_park")
            CALL helix_hmhd_equil(x,y,u0,u0x,u0y,.TRUE.)
            first=.FALSE.
            CALL helix_hmhd_rhs(x,y,u0,u0x,u0y,fx0,fy0,s0,first)
            fx=fx-fx0
            fy=fy-fy0
            s=s-s0
         CASE DEFAULT
            CALL program_stop("No sources for init_type = "//init_type)
         END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_rhs
c-----------------------------------------------------------------------
c     subprogram 7. helix_hmhd_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_drdu(x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy

      REAL(r8) :: gfac_s,Eb
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac,hfac,kfac,
     $     xfac,yfac,pfac,Ve,eta_local
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
c     calculate auxiliary variables.
c-----------------------------------------------------------------------
      gfac=1/(1+(x**2+y**2)*Rinv**2)
      hfac=2*Rinv*gfac
      kfac=4*Rinv**2*gfac**3
      xfac=1-Rinv**2*x**2*gfac
      yfac=1-Rinv**2*y**2*gfac
      pfac=x*y*Rinv**2*gfac
      Ve=u(4,:,:)-di*u(5,:,:)
c-----------------------------------------------------------------------
c     modify resistivity for boundary applied E-field(no current drive).
c-----------------------------------------------------------------------
      eta_local=eta
      gfac_s=1/(1+(r_s*Rinv)**2)
      SELECT CASE(init_type)
      CASE("m1_park")
         IF(.NOT. source)
     $        eta_local=eta*ABS((r_s**2/gfac_s - 2*gfac_s)
     $        /((2*(x**2+y**2)-r_s**2)/gfac - 2*gfac))
      CASE("m1_axial")
         Eb=eta*(r_s**2/gfac_s - 2)
         eta_local=ABS(Eb/((1+2*Rinv**2*(x**2+y**2))
     $        *(x**2+y**2-r_s**2) + (x**2+y**2)/gfac - 2))
      END SELECT
c-----------------------------------------------------------------------
c     Helical incomressible HMHD model in cartesian coordinates
c-----------------------------------------------------------------------
      fx_ux(1,4,:,:)=-di*nu*xfac
      fx_ux(1,5,:,:)=di**2*nu*xfac
      fx_uy(1,4,:,:)=di*nu*pfac
      fx_uy(1,5,:,:)=-di**2*nu*pfac
      fy_ux(1,4,:,:)=di*nu*pfac
      fy_ux(1,5,:,:)=-di**2*nu*pfac
      fy_uy(1,4,:,:)=-di*nu*yfac
      fy_uy(1,5,:,:)=di**2*nu*yfac

      s_u(1,4,:,:)=-di*nu*kfac
      s_u(1,5,:,:)=eta_local*gfac + di**2*nu*kfac
      s_ux(1,1,:,:)=gfac*(uy(6,:,:) + di*uy(2,:,:))
      s_ux(1,2,:,:)=-gfac*di*uy(1,:,:)
      s_ux(1,6,:,:)=-gfac*uy(1,:,:)
      s_uy(1,1,:,:)=-gfac*(ux(6,:,:) + di*ux(2,:,:))
      s_uy(1,2,:,:)=gfac*di*ux(1,:,:)
      s_uy(1,6,:,:)=gfac*ux(1,:,:)
      
      fx_u(2,2,:,:)=-gfac*(uy(6,:,:) + di*uy(2,:,:))
      fx_u(2,4,:,:)=gfac*uy(1,:,:) + di*nu*Rinv*x*kfac
      fx_u(2,5,:,:)=-di*(gfac*uy(1,:,:) + di*nu*Rinv*x*kfac)
      fx_ux(2,2,:,:)=-eta_local*xfac
      fx_ux(2,4,:,:)=-di*nu*xfac*hfac
      fx_ux(2,5,:,:)=di**2*nu*xfac*hfac
      fx_uy(2,1,:,:)=gfac*Ve
      fx_uy(2,2,:,:)=-di*gfac*u(2,:,:) + eta_local*pfac
      fx_uy(2,4,:,:)=di*nu*pfac*hfac
      fx_uy(2,5,:,:)=-di**2*nu*pfac*hfac
      fx_uy(2,6,:,:)=-gfac*u(2,:,:)

      fy_u(2,2,:,:)=gfac*(ux(6,:,:) + di*ux(2,:,:))
      fy_u(2,4,:,:)=-gfac*ux(1,:,:) + di*nu*Rinv*y*kfac
      fy_u(2,5,:,:)=di*(gfac*ux(1,:,:) - di*nu*Rinv*y*kfac)
      fy_ux(2,1,:,:)=-gfac*Ve
      fy_ux(2,2,:,:)=di*gfac*u(2,:,:) + eta_local*pfac
      fy_ux(2,4,:,:)=di*nu*pfac*hfac
      fy_ux(2,5,:,:)=-di**2*nu*pfac*hfac
      fy_ux(2,6,:,:)=gfac*u(2,:,:)
      fy_uy(2,2,:,:)=-eta_local*yfac
      fy_uy(2,4,:,:)=-di*nu*yfac*hfac
      fy_uy(2,5,:,:)=di**2*nu*yfac*hfac
      
      fx_u(3,2,:,:)=gfac*hfac*uy(1,:,:)
      fx_u(3,3,:,:)=-gfac*uy(6,:,:)
      fx_u(3,4,:,:)=-gfac*hfac*uy(6,:,:) - nu*Rinv*x*kfac
      fx_u(3,5,:,:)=gfac*uy(1,:,:) + di*nu*Rinv*x*kfac
      fx_ux(3,3,:,:)=-mu*xfac
      fx_ux(3,4,:,:)=nu*xfac*hfac
      fx_ux(3,5,:,:)=-di*nu*xfac*hfac
      fx_uy(3,1,:,:)=gfac*(u(5,:,:) + hfac*u(2,:,:))
      fx_uy(3,3,:,:)=mu*pfac
      fx_uy(3,4,:,:)=-nu*pfac*hfac
      fx_uy(3,5,:,:)=di*nu*pfac*hfac
      fx_uy(3,6,:,:)=-gfac*(u(3,:,:) + hfac*u(4,:,:))

      fy_u(3,2,:,:)=-gfac*hfac*ux(1,:,:)
      fy_u(3,3,:,:)=gfac*ux(6,:,:)
      fy_u(3,4,:,:)=gfac*hfac*ux(6,:,:) - nu*Rinv*y*kfac
      fy_u(3,5,:,:)=-gfac*ux(1,:,:) + di*nu*Rinv*y*kfac
      fy_ux(3,1,:,:)=-gfac*(u(5,:,:) + hfac*u(2,:,:))
      fy_ux(3,3,:,:)=mu*pfac
      fy_ux(3,4,:,:)=-nu*pfac*hfac
      fy_ux(3,5,:,:)=di*nu*pfac*hfac
      fy_ux(3,6,:,:)=gfac*(u(3,:,:) + hfac*u(4,:,:))
      fy_uy(3,3,:,:)=-mu*yfac
      fy_uy(3,4,:,:)=nu*yfac*hfac
      fy_uy(3,5,:,:)=-di*nu*yfac*hfac

      s_u(3,2,:,:)=hfac**2/2*(x*uy(2,:,:) - y*ux(2,:,:) 
     $     -2*hfac*(x*uy(1,:,:) - y*ux(1,:,:)))
      s_u(3,4,:,:)=-hfac**2/2*(x*uy(4,:,:) - y*ux(4,:,:) 
     $     -2*hfac*(x*uy(6,:,:) - y*ux(6,:,:)))
      s_ux(3,1,:,:)=hfac**3*u(2,:,:)*y
      s_ux(3,2,:,:)=-hfac**2/2*u(2,:,:)*y
      s_ux(3,4,:,:)=hfac**2/2*u(4,:,:)*y
      s_ux(3,6,:,:)=-hfac**3*u(4,:,:)*y
      s_uy(3,1,:,:)=-hfac**3*u(2,:,:)*x
      s_uy(3,2,:,:)=hfac**2/2*u(2,:,:)*x
      s_uy(3,4,:,:)=-hfac**2/2*u(4,:,:)*x
      s_uy(3,6,:,:)=hfac**3*u(4,:,:)*x

      fx_ux(4,4,:,:)=-(mu+nu)*xfac
      fx_ux(4,5,:,:)=di*nu*xfac
      fx_uy(4,4,:,:)=(mu+nu)*pfac
      fx_uy(4,5,:,:)=-di*nu*pfac

      fy_ux(4,4,:,:)=(mu+nu)*pfac
      fy_ux(4,5,:,:)=-di*nu*pfac
      fy_uy(4,4,:,:)=-(mu+nu)*yfac
      fy_uy(4,5,:,:)=di*nu*yfac

      s_u(4,3,:,:)=mu*gfac*hfac
      s_u(4,4,:,:)=-nu*kfac
      s_u(4,5,:,:)=di*nu*kfac
      s_ux(4,1,:,:)=gfac*uy(2,:,:)
      s_ux(4,2,:,:)=-gfac*uy(1,:,:)
      s_ux(4,4,:,:)=gfac*uy(6,:,:)
      s_ux(4,6,:,:)=-gfac*uy(4,:,:)
      s_uy(4,1,:,:)=-gfac*ux(2,:,:)
      s_uy(4,2,:,:)=gfac*ux(1,:,:)
      s_uy(4,4,:,:)=-gfac*ux(6,:,:)
      s_uy(4,6,:,:)=gfac*ux(4,:,:)
      
      fx_ux(5,1,:,:)=xfac
      fx_uy(5,1,:,:)=-pfac
      fy_ux(5,1,:,:)=-pfac
      fy_uy(5,1,:,:)=yfac
      s_u(5,2,:,:)=gfac*hfac
      s_u(5,5,:,:)=gfac

      fx_ux(6,6,:,:)=xfac
      fx_uy(6,6,:,:)=-pfac
      fy_ux(6,6,:,:)=-pfac
      fy_uy(6,6,:,:)=yfac
      s_u(6,3,:,:)=gfac
      s_u(6,4,:,:)=gfac*hfac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_drdu
c-----------------------------------------------------------------------
c     subprogram 8. helix_hmhd_mass.
c     computes mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_mass(x,y,mass)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass

      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac
c-----------------------------------------------------------------------
c     modify mass matrix
c-----------------------------------------------------------------------
      gfac=1/(1+(x**2+y**2)*Rinv**2)

      mass(1,1,:,:)=gfac
      mass(2,2,:,:)=gfac
      mass(2,1,:,:)=-2*Rinv*gfac**2
      mass(3,3,:,:)=gfac
      mass(3,4,:,:)=4*Rinv*gfac**2
      mass(4,4,:,:)=gfac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_mass
c-----------------------------------------------------------------------
c     subprogram 9. helix_hmhd_equil.
c     computes equilibrium for Helical RMHD model.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_equil(x,y,u,ux,uy,derivs)
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs

      INTEGER :: ix1,ix2
      REAL(r8) :: gfac_s,psi0
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: gfac,rad_sq
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      rad_sq=x**2 + y**2
      gfac=1/(1 + rad_sq*Rinv**2)
      SELECT CASE(init_type)
      CASE("m1_park")
         u(1,:,:)=rad_sq**3*Rinv**2/12 
     $        + rad_sq**2*(1 - r_s**2*Rinv**2)/8 - rad_sq*r_s**2/4
         u(2,:,:)=1/Rinv
         u(5,:,:)=(2*rad_sq - r_s**2)/gfac - 2*gfac
         IF(.NOT. derivs)RETURN
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
         ux(1,:,:)=x*rad_sq**2*Rinv**2/2 
     $        + x*rad_sq*(1 - r_s**2*Rinv**2)/2 - x*r_s**2/2
         ux(5,:,:)=2*x*((2*rad_sq - r_s**2)*Rinv**2 + 2/gfac 
     $        + 2*(gfac*Rinv)**2)
         uy(1,:,:)=y*rad_sq**2*Rinv**2/2 
     $        + y*rad_sq*(1 - r_s**2*Rinv**2)/2 - y*r_s**2/2
         uy(5,:,:)=2*y*((2*rad_sq - r_s**2)*Rinv**2 + 2/gfac 
     $        + 2*(gfac*Rinv)**2)
      CASE("m1_axial")
         IF(r_s <= 0)CALL program_stop("r_s has to be positive!")
         gfac_s=(1. + (Rinv*r_s)**2/3.)
         psi0=half/gfac_s
         u(1,:,:)=(rad_sq**2*Rinv**2/3. + rad_sq*(1. - (r_s*Rinv)**2)/2.
     $        - r_s**2)*rad_sq/4.
         u(2,:,:)=1/(gfac*Rinv) - Rinv*rad_sq/(2*gfac)*(rad_sq - r_s**2)
         u(5,:,:)=(2*rad_sq-r_s**2)/gfac 
     $        + Rinv**2*rad_sq*(rad_sq - r_s**2) - 2
         DO ix1 = 1,SIZE(x,1)
            DO ix2 = 1,SIZE(x,2)
               IF(rad_sq(ix1,ix2) <= r_s**2)THEN
                  u(1,ix1,ix2)=u(1,ix1,ix2) - psi0*gfac_s/(8*r_s**4)
     $                 *(rad_sq(ix1,ix2) - r_s**2)**4
                  u(5,ix1,ix2)=u(5,ix1,ix2) - 2*psi0*gfac_s/r_s**4
     $                 *(rad_sq(ix1,ix2) - r_s**2)**2
     $                 *(gfac(ix1,ix2)*(rad_sq(ix1,ix2) - r_s**2) 
     $                 + 3*rad_sq(ix1,ix2))
               ENDIF
            ENDDO
         ENDDO
      CASE DEFAULT
         CALL program_stop("No initial condition specified")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_equil
c-----------------------------------------------------------------------
c     subprogram 10. helix_hmhd_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE helix_hmhd_grid(x,y,ksi,eta,polar_axis)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta
      LOGICAL, INTENT(OUT) :: polar_axis

      REAL(r8) :: rs
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      rs=1/(1+((1-r_s)/r_s)**third)
      rs=rs+1.5*gr_curve*(r_s-rs)
      polar_axis=.TRUE.
      ksi=((x-rs)**3+5*gr_curve*x+rs**3)/(1+5*gr_curve+3*rs*(rs-1))
      eta=twopi*y
c-----------------------------------------------------------------------
c     convert to cartesian coordinates
c-----------------------------------------------------------------------
      ksi=ksi*COS(eta)
      eta=ksi*TAN(eta)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE helix_hmhd_grid
      END MODULE helix_hmhd_mod
