c-----------------------------------------------------------------------
c     file job2.f.
c     contains various 2D PDEs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. job2_mod.
c     1. job2_input.
c     2. job2_init.
c     3. job2_couplemat.
c     4. job2_edge_rhs.
c     5. job2_edge_drdu.
c     6. job2_edge_mass.
c     7. job2_rhs.
c     8. job2_drdu.
c     9. job2_mass.
c     10. job2_alloc.
c     11. job2_dealloc.
c     12. job2_bc_default.
c     13. job2_grid.
c     14. job2_schur.
c-----------------------------------------------------------------------
c     subprogram 0. job2_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE job2_mod
      USE beltrami_mod
      USE cubit_mod
      IMPLICIT NONE
 
      CHARACTER(8) :: grid_type="sel",grid_inv_type="jacobi"
      CHARACTER(16) :: job_type="physics"
      CHARACTER(80) :: cubit_file="."

      TYPE(block) :: blk
      INTEGER, DIMENSION(:), ALLOCATABLE :: grid_bpiv
      REAL(r8), DIMENSION(:), ALLOCATABLE :: inv_points
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: grid_bmat

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. job2_input.
c     set up input constants, executed only on mpi_rank=0.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $     nqty_schur,dt,dtmax,tmax,nstep,du_diagnose,physics_type,
     $     exit_flag)

      CHARACTER(*), INTENT(INOUT) :: physics_type
      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty,nqty_schur
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep,exit_flag
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,
     $        nqty,nqty_schur,dt,dtmax,tmax,nstep,du_diagnose,
     $        physics_type,exit_flag)
         USE local_mod
         IMPLICIT NONE
         CHARACTER(*), INTENT(INOUT) :: physics_type
         LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
         INTEGER, INTENT(OUT) :: nqty,nqty_schur
         INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep,exit_flag
         REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax
         END SUBROUTINE physics_input
      END INTERFACE
c-----------------------------------------------------------------------
c     call job-specific subroutines.
c-----------------------------------------------------------------------
      CALL physics_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $     nqty_schur,dt,dtmax,tmax,nstep,du_diagnose,physics_type,
     $     exit_flag)
c-----------------------------------------------------------------------
c     Read in CUBIT grid data.
c-----------------------------------------------------------------------
      IF(grid_type=="cubit")THEN
         blk%np=np
         CALL cubit_read_cdf(cubit_file,blk,exit_flag)
         IF (blk%nx /= nx)THEN
            WRITE(*,'(a)')
     $           'Switching nx and ny to accommodate CUBIT grid'
            nx = blk%nx
            ny = blk%ny
            nbx = mpi_size/nbx
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_input
c-----------------------------------------------------------------------
c     subprogram 2. job2_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_init(jtype,xperiodic,yperiodic,ixmin,iymin,
     $     ixmax,iymax,nx,ny,np,ground,couple_msize,couple_mlist,
     $     edge_order,edges,u)

      CHARACTER(*), INTENT(IN) :: jtype
      LOGICAL, INTENT(IN) :: xperiodic,yperiodic
      INTEGER, INTENT(IN) :: ixmin,iymin,ixmax,iymax,nx,ny,np
      LOGICAL, DIMENSION(:), INTENT(OUT) :: ground
      INTEGER, INTENT(OUT) :: couple_msize
      INTEGER, DIMENSION(:,:), POINTER :: couple_mlist
      INTEGER, DIMENSION(4), INTENT(OUT) :: edge_order
      TYPE(edge_type), DIMENSION(:) :: edges
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: u

      INTEGER :: iqty,jqty,nqty,jx,jy,info,ixx,iyy,ii
      LOGICAL, DIMENSION(SIZE(u,1)) :: static
      LOGICAL, DIMENSION(SIZE(u,1),SIZE(u,1)) :: couple_mmat
      REAL(r8), DIMENSION(0:np,0:np) :: xpi,ypi
      REAL(r8), DIMENSION(SIZE(u,1),0:np,0:np) :: ui
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_init_parameters(static,ground,adapt_qty)
         USE local_mod
         IMPLICIT NONE
         LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
         END SUBROUTINE physics_init_parameters
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_boundary(left,right,top,bottom,
     $        nqty,edge_order)
         USE local_mod
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: nqty
         INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
         TYPE(edge_type) :: left,right,top,bottom
         END SUBROUTINE physics_boundary
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_init(xpi,ypi,ui)
         USE local_mod
         IMPLICIT NONE
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: xpi,ypi
         REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: ui
         END SUBROUTINE physics_init
      END INTERFACE
c-----------------------------------------------------------------------
c     compute local sizes.
c-----------------------------------------------------------------------
      nqty=SIZE(u,1)
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      static=.FALSE.
      ground=.FALSE.
      couple_mmat=.TRUE.
c-----------------------------------------------------------------------
c     set default edge b.c.
c-----------------------------------------------------------------------
      edge_order=(/1,3,2,4/)

      DO ii=1,SIZE(edges)
         edges(ii)%bc_type="none"
         edges(ii)%static=.TRUE.
      ENDDO
      IF(xperiodic)THEN
         edges(1)%bc_type="periodic"
         edges(3)%bc_type="periodic"
      ENDIF
      IF(yperiodic)THEN
         edges(2)%bc_type="periodic"
         edges(4)%bc_type="periodic"
      ENDIF
c-----------------------------------------------------------------------
c     initialize parameters, boundary conditions and coupling matrices.
c-----------------------------------------------------------------------
      SELECT CASE(jtype)
      CASE("physics")
         ALLOCATE(adapt_qty(nqty))
         adapt_qty=.TRUE.
         CALL physics_init_parameters(static,ground,adapt_qty)
         CALL physics_boundary(edges(1),edges(3),edges(2),edges(4),
     $        nqty,edge_order)
         CALL job2_bc_default(edges,static)
         CALL job2_couplemat(edges,static,couple_mmat)
      CASE("beltrami")
         static=.TRUE.
         IF(bel_diagnose .AND. mpi_rank == 0)CALL beltrami_diagnose
         CALL beltrami_boundary(edges(1),edges(3),edges(2),edges(4))
         CALL job2_bc_default(edges,static)
         couple_mmat=.FALSE.
      END SELECT
c-----------------------------------------------------------------------
c     create mass matrix coupling list.
c-----------------------------------------------------------------------
      couple_msize = COUNT(couple_mmat)
      IF(couple_msize > 0)ALLOCATE(couple_mlist(2,couple_msize))

      ii=0
      DO iqty=1,nqty
         DO jqty=1,nqty
            IF(couple_mmat(jqty,iqty))THEN
               ii = ii+1
               couple_mlist(:,ii)=(/jqty,iqty/)
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     start loop over cell elements.
c-----------------------------------------------------------------------
      ui=0
      jx=0
      DO ixx=ixmin,ixmax-1
         jy=0
         DO iyy=iymin,iymax-1
            CALL job2_grid(jtype,nx,ny,np,ixx,iyy,xpi,ypi)
c-----------------------------------------------------------------------
c     initial conditions.
c-----------------------------------------------------------------------
            SELECT CASE(jtype)
            CASE("physics")
               CALL physics_init(xpi,ypi,ui)
            CASE("beltrami")
               CALL beltrami_init(xpi,ypi,ui)
            END SELECT
c-----------------------------------------------------------------------
c     solve for basis amplitudes.
c-----------------------------------------------------------------------
            DO iqty=1,nqty
               DO ii=1,2
                  CALL dgetrs('N',np+1,np+1,grid_bmat,np+1,grid_bpiv,
     $                 ui(iqty,:,:),np+1,info)
                  ui(iqty,:,:)=TRANSPOSE(ui(iqty,:,:))
               ENDDO
            ENDDO
            u(:,jx:jx+np,jy:jy+np)=ui
c-----------------------------------------------------------------------
c     finish loop over finite elements.
c-----------------------------------------------------------------------
            jy=jy+np+1
         ENDDO
         jx=jx+np+1
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_init
c-----------------------------------------------------------------------
c     subprogram 3. job2_couplemat.
c     determines the dependent variable coupling matrices
c     for the drdu & mass matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_couplemat(edges,static,couple_mmat)

      TYPE(edge_type), DIMENSION(:) :: edges
      LOGICAL, DIMENSION(:), INTENT(IN) :: static
      LOGICAL, DIMENSION(:,:), INTENT(OUT) :: couple_mmat

      INTEGER :: i,iqty
      REAL(r8) :: t_rand
      REAL(r8), DIMENSION(1) :: x_rand,y_rand
      REAL(r8), DIMENSION(ndim,1) :: nvec_rand
      REAL(r8), DIMENSION(1,SIZE(couple_mmat,1),SIZE(couple_mmat,2)) 
     $     :: mmat
      REAL(r8), DIMENSION(1,SIZE(couple_mmat,1),SIZE(couple_mmat,2),3)
     $     :: m_rand
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      couple_mmat=.TRUE.
      CALL RANDOM_NUMBER(t_rand)
      CALL RANDOM_NUMBER(x_rand)
      CALL RANDOM_NUMBER(y_rand)
      CALL RANDOM_NUMBER(nvec_rand)
c-----------------------------------------------------------------------
c     calculate the mass matrix coupling.
c-----------------------------------------------------------------------
      mmat = 0
      CALL job2_mass(x_rand,y_rand,m_rand(:,:,:,1),m_rand(:,:,:,2),
     $     m_rand(:,:,:,3))
      mmat = SUM(ABS(m_rand),4)
      DO i=1,SIZE(edges)
         CALL job2_edge_mass(edges(i),x_rand,y_rand,nvec_rand,
     $        m_rand(:,:,:,1),m_rand(:,:,:,2),m_rand(:,:,:,3))
         mmat = mmat + SUM(ABS(m_rand),4)
      ENDDO
      WHERE(mmat(1,:,:) == 0)
         couple_mmat = .FALSE.
      END WHERE
      
      DO iqty=1,SIZE(static)
         IF(static(iqty))couple_mmat(:,iqty) = .FALSE.
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_couplemat
c-----------------------------------------------------------------------
c     subprogram 4. job2_edge_rhs.
c     computes rhs for non-linear edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_edge_rhs(edge,t,x0,y0,nvec,u0,ux0,uy0,uxx0,uyy0,
     $     uxy0,c0)

      TYPE(edge_type) :: edge
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: nvec,u0,ux0,uy0,uxx0,uyy0,
     $     uxy0
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: c0

      LOGICAL :: nat_flag
      INTEGER :: iqty
      REAL(r8), DIMENSION(SIZE(x0),1) :: nnorm,fac,x,y
      REAL(r8), DIMENSION(ndim,SIZE(x0),1) :: nhat
      REAL(r8), DIMENSION(SIZE(x0),SIZE(u0,1)) :: fx0,fy0,s0
      REAL(r8), DIMENSION(SIZE(u0,1),SIZE(x0),1) :: c,u,ux,uy,uxx,uyy,
     $     uxy
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
c     initialize local variables
c-----------------------------------------------------------------------
      c = zero
      nnorm(:,1) = SQRT(SUM(nvec**2,1))
      DO iqty=1,ndim
         nhat(iqty,:,1) = nvec(iqty,:)/nnorm(:,1)
      ENDDO

      x(:,1) = x0
      y(:,1) = y0
      u(:,:,1) = u0
      ux(:,:,1) = ux0
      uy(:,:,1) = uy0
      uxx(:,:,1) = uxx0
      uyy(:,:,1) = uyy0
      uxy(:,:,1) = uxy0
c-----------------------------------------------------------------------
c     give values to c.
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("physics")
         CALL physics_edge_rhs(edge%edgename,t,x,y,nhat,u,ux,uy,uxx,uyy,
     $        uxy,c)
      CASE("beltrami")
         CALL beltrami_edge_rhs(edge%edgename,x,y,u,ux,uy,c)
      END SELECT
c-----------------------------------------------------------------------
c     reorder indices and re-scale "normflux" b.c. flux
c-----------------------------------------------------------------------
      DO iqty=1,SIZE(u0,1)
         IF(edge%bc_type(iqty)=="normflux")THEN
            fac = -nnorm
         ELSE
            fac = one
         ENDIF
         c0(:,iqty) = fac(:,1)*c(iqty,:,1)
      ENDDO
c-----------------------------------------------------------------------
c     if necessary, calculate "natural" b.c. fluxes
c-----------------------------------------------------------------------
      nat_flag=.FALSE.
      DO iqty=1,SIZE(u0,1)
         IF(edge%bc_type(iqty)=="natural")THEN
            nat_flag=.TRUE.
            EXIT
         ENDIF
      ENDDO
      IF(nat_flag)THEN
         CALL job2_rhs(t,x0,y0,u0,ux0,uy0,fx0,fy0,s0)
         DO iqty=1,SIZE(u0,1)
            IF(edge%bc_type(iqty)=="natural")THEN
               c0(:,iqty) = -(fx0(:,iqty)*nvec(1,:)
     $              + fy0(:,iqty)*nvec(2,:))
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_edge_rhs
c-----------------------------------------------------------------------
c     subprogram 5. job2_edge_drdu.
c     computes drdu for non-linear edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_edge_drdu(edge,t,x0,y0,nvec,u0,ux0,uy0,uxx0,uyy0,
     $     uxy0,c_u0,c_ux0,c_uy0,c_uxx0,c_uyy0,c_uxy0)

      TYPE(edge_type) :: edge
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: nvec,u0,ux0,uy0,uxx0,uyy0,
     $     uxy0
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c_u0,c_ux0,c_uy0,
     $     c_uxx0,c_uyy0,c_uxy0

      LOGICAL :: nat_flag
      INTEGER :: iqty,jqty
      REAL(r8), DIMENSION(SIZE(x0),1) :: nnorm,fac,x,y
      REAL(r8), DIMENSION(ndim,SIZE(x0),1) :: nhat
      REAL(r8), DIMENSION(SIZE(u0,1),SIZE(x0),1) :: u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(SIZE(x0),SIZE(u0,1),SIZE(u0,1)) :: 
     $     fx_u0,fx_ux0,fx_uy0,fy_u0,fy_ux0,fy_uy0,s_u0,s_ux0,s_uy0
      REAL(r8), DIMENSION(SIZE(u0,1),SIZE(u0,1),SIZE(x0),1) ::
     $     c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy
c-----------------------------------------------------------------------
c     interface block
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_edge_drdu(lrtb,t,x,y,nhat,
     $        u,ux,uy,uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)
         USE local_mod
         IMPLICIT NONE
         CHARACTER(*), INTENT(IN) :: lrtb
         REAL(r8), INTENT(IN) :: t
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
         REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,
     $        uxy
         REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $        c_uxx,c_uyy,c_uxy
         END SUBROUTINE physics_edge_drdu
      END INTERFACE
c-----------------------------------------------------------------------
c     initialize local variables
c-----------------------------------------------------------------------
      c_u = zero
      c_ux = zero
      c_uy = zero
      c_uxx = zero
      c_uyy = zero
      c_uxy = zero
      nnorm(:,1) = SQRT(SUM(nvec**2,1))
      DO iqty=1,ndim
         nhat(iqty,:,1) = nvec(iqty,:)/nnorm(:,1)
      ENDDO

      x(:,1) = x0
      y(:,1) = y0
      u(:,:,1) = u0
      ux(:,:,1) = ux0
      uy(:,:,1) = uy0
      uxx(:,:,1) = uxx0
      uyy(:,:,1) = uyy0
      uxy(:,:,1) = uxy0
c-----------------------------------------------------------------------
c     give values to c_u,c_ux,c_uy,c_uxx,c_uyy.
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("physics")
         CALL physics_edge_drdu(edge%edgename,t,x,y,nhat,u,ux,uy,
     $        uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)
      CASE("beltrami")
         CALL beltrami_edge_drdu(edge%edgename,x,y,u,ux,uy,
     $        c_u,c_ux,c_uy)
      END SELECT
c-----------------------------------------------------------------------
c     reorder indices and re-scale drdu for "normflux" b.c. flux
c-----------------------------------------------------------------------
      DO iqty=1,SIZE(u0,1)
         IF(edge%bc_type(iqty)=="normflux")THEN
            fac = -nnorm
         ELSE
            fac = one
         ENDIF
         
         DO jqty=1,SIZE(u0,1)
            c_u0(:,jqty,iqty) = fac(:,1)*c_u(iqty,jqty,:,1)
            c_ux0(:,jqty,iqty) = fac(:,1)*c_ux(iqty,jqty,:,1)
            c_uy0(:,jqty,iqty) = fac(:,1)*c_uy(iqty,jqty,:,1)
            c_uxx0(:,jqty,iqty) = fac(:,1)*c_uxx(iqty,jqty,:,1)
            c_uyy0(:,jqty,iqty) = fac(:,1)*c_uyy(iqty,jqty,:,1)
            c_uxy0(:,jqty,iqty) = fac(:,1)*c_uxy(iqty,jqty,:,1)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     if necessary, calculate "natural" b.c. fluxes
c-----------------------------------------------------------------------
      nat_flag=.FALSE.
      DO iqty=1,SIZE(u0,1)
         IF(edge%bc_type(iqty)=="natural")THEN
            nat_flag=.TRUE.
            EXIT
         ENDIF
      ENDDO
      IF(nat_flag)THEN
         CALL job2_drdu(t,x0,y0,u0,ux0,uy0,fx_u0,fx_ux0,fx_uy0,
     $        fy_u0,fy_ux0,fy_uy0,s_u0,s_ux0,s_uy0)
         DO iqty=1,SIZE(u0,1)
            IF(edge%bc_type(iqty)=="natural")THEN
               DO jqty=1,SIZE(u0,1)
                  c_u0(:,jqty,iqty) = -(fx_u0(:,jqty,iqty)*nvec(1,:)
     $                 + fy_u0(:,jqty,iqty)*nvec(2,:))
                  c_ux0(:,jqty,iqty) = -(fx_ux0(:,jqty,iqty)*nvec(1,:)
     $                 + fy_ux0(:,jqty,iqty)*nvec(2,:))
                  c_uy0(:,jqty,iqty) = -(fx_uy0(:,jqty,iqty)*nvec(1,:)
     $                 + fy_uy0(:,jqty,iqty)*nvec(2,:))
               ENDDO
               c_uxx0(:,:,iqty) = zero
               c_uyy0(:,:,iqty) = zero
               c_uxy0(:,:,iqty) = zero
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_edge_drdu
c-----------------------------------------------------------------------
c     subprogram 6. job2_edge_mass.
c     computes mass matrices for non-linear edge b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_edge_mass(edge,x0,y0,nvec,mass0,mass0_x,mass0_y)

      TYPE(edge_type) :: edge
      REAL(r8), DIMENSION(:), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: nvec
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: mass0,mass0_x,mass0_y

      INTEGER :: iqty,jqty,nqty
      REAL(r8), DIMENSION(SIZE(x0),1) :: x,y
      REAL(r8), DIMENSION(ndim,SIZE(x0),1) :: nhat
      REAL(r8), DIMENSION(SIZE(mass0,3),SIZE(mass0,3),SIZE(x0),1) ::
     $     mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_edge_mass(lrtb,x,y,nhat,mass,mass_x,mass_y)
         USE local_mod
         IMPLICIT NONE
         CHARACTER(*), INTENT(IN) :: lrtb
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
         REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
         REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y
         END SUBROUTINE physics_edge_mass
      END INTERFACE
c-----------------------------------------------------------------------
c     initialize local variables
c-----------------------------------------------------------------------
      nqty = SIZE(mass0,3)
      mass = zero
      mass_x = zero
      mass_y = zero

      DO iqty=1,ndim
         nhat(iqty,:,1) = nvec(iqty,:)/SQRT(SUM(nvec**2,1))
      ENDDO

      x(:,1) = x0
      y(:,1) = y0
c-----------------------------------------------------------------------
c     give values to mass, mass_x and mass_y.
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("physics")
         CALL physics_edge_mass(edge%edgename,x,y,nhat,
     $        mass,mass_x,mass_y)
      END SELECT
c-----------------------------------------------------------------------
c     reorder and reset mass matrices to zero for a static b.c.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         DO jqty=1,nqty
            mass0(:,jqty,iqty) = mass(iqty,jqty,:,1)
            mass0_x(:,jqty,iqty) = mass_x(iqty,jqty,:,1)
            mass0_y(:,jqty,iqty) = mass_y(iqty,jqty,:,1)
         ENDDO
         IF(edge%static(iqty))THEN
            mass0(:,:,iqty) = zero
            mass0_x(:,:,iqty) = zero
            mass0_y(:,:,iqty) = zero
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_edge_mass
c-----------------------------------------------------------------------
c     subprogram 7. job2_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_rhs(t,x0,y0,u0,ux0,uy0,fx0,fy0,s0)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u0,ux0,uy0
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: fx0,fy0,s0

      LOGICAL :: recurs
      REAL(r8), DIMENSION(SIZE(x0),1) :: x,y
      REAL(r8), DIMENSION(SIZE(u0,1),SIZE(x0),1) :: u,ux,uy,fx,fy,s
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         RECURSIVE SUBROUTINE physics_rhs(t,x,y,u,ux,uy,fx,fy,s,first)
         USE local_mod
         IMPLICIT NONE
         REAL(r8), INTENT(IN) :: t
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
         REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
         REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
         LOGICAL, INTENT(INOUT) :: first
         END SUBROUTINE physics_rhs
      END INTERFACE
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      x(:,1) = x0
      y(:,1) = y0
      u(:,:,1) = u0
      ux(:,:,1) = ux0
      uy(:,:,1) = uy0
      fx = zero
      fy = zero
      s = zero
      recurs = .TRUE.

      SELECT CASE(job_type)
      CASE("physics")
         CALL physics_rhs(t,x,y,u,ux,uy,fx,fy,s,recurs)
      CASE("beltrami")
         CALL beltrami_rhs(x,y,u,ux,uy,fx,fy,s)
      END SELECT
c-----------------------------------------------------------------------
c     reorder indices for optimization purposes
c-----------------------------------------------------------------------
      fx0 = TRANSPOSE(fx(:,:,1))
      fy0 = TRANSPOSE(fy(:,:,1))
      s0 = TRANSPOSE(s(:,:,1))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_rhs
c-----------------------------------------------------------------------
c     subprogram 8. job2_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_drdu(t,x0,y0,u0,ux0,uy0,
     $     fx_u0,fx_ux0,fx_uy0,fy_u0,fy_ux0,fy_uy0,s_u0,s_ux0,s_uy0)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u0,ux0,uy0
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) ::
     $     fx_u0,fx_ux0,fx_uy0,fy_u0,fy_ux0,fy_uy0,s_u0,s_ux0,s_uy0

      INTEGER :: iqty
      REAL(r8), DIMENSION(SIZE(x0),1) :: x,y
      REAL(r8), DIMENSION(SIZE(u0,1),SIZE(x0),1) :: u,ux,uy
      REAL(r8), DIMENSION(SIZE(u0,1),SIZE(u0,1),SIZE(x0),1) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_drdu(t,x,y,u,ux,uy,
     $        fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
         USE local_mod
         IMPLICIT NONE
         REAL(r8), INTENT(IN) :: t
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
         REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
         REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $        fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy
         END SUBROUTINE physics_drdu
      END INTERFACE
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      x(:,1) = x0
      y(:,1) = y0
      u(:,:,1) = u0
      ux(:,:,1) = ux0
      uy(:,:,1) = uy0
      fx_u = zero
      fx_ux = zero
      fx_uy = zero
      fy_u = zero
      fy_ux = zero
      fy_uy = zero
      s_u = zero
      s_ux = zero
      s_uy = zero

      SELECT CASE(job_type)
      CASE("physics")
         CALL physics_drdu(t,x,y,u,ux,uy,fx_u,fx_ux,fx_uy,
     $        fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
      CASE("beltrami")
         CALL beltrami_drdu(x,y,u,ux,uy,fx_u,fx_ux,fx_uy,fy_u,fy_ux,
     $        fy_uy,s_u,s_ux,s_uy)
      END SELECT
c-----------------------------------------------------------------------
c     reorder indices for optimization purposes
c-----------------------------------------------------------------------
      DO iqty=1,SIZE(u0,1)
         fx_u0(:,:,iqty) = TRANSPOSE(fx_u(iqty,:,:,1))
         fx_ux0(:,:,iqty) = TRANSPOSE(fx_ux(iqty,:,:,1))
         fx_uy0(:,:,iqty) = TRANSPOSE(fx_uy(iqty,:,:,1))
         fy_u0(:,:,iqty) = TRANSPOSE(fy_u(iqty,:,:,1))
         fy_ux0(:,:,iqty) = TRANSPOSE(fy_ux(iqty,:,:,1))
         fy_uy0(:,:,iqty) = TRANSPOSE(fy_uy(iqty,:,:,1))
         s_u0(:,:,iqty) = TRANSPOSE(s_u(iqty,:,:,1))
         s_ux0(:,:,iqty) = TRANSPOSE(s_ux(iqty,:,:,1))
         s_uy0(:,:,iqty) = TRANSPOSE(s_uy(iqty,:,:,1))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_drdu
c-----------------------------------------------------------------------
c     subprogram 9. job2_mass.
c     gives mass matrix couplings.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_mass(x0,y0,mass0,mass_x0,mass_y0)
      
      REAL(r8), DIMENSION(:), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: mass0,mass_x0,mass_y0
      
      INTEGER :: iqty
      REAL(r8), DIMENSION(SIZE(x0),1) :: x,y
      REAL(r8), DIMENSION(SIZE(mass0,3),SIZE(mass0,3),SIZE(x0),1) ::
     $     mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_mass(x,y,mass,mass_x,mass_y)
         USE local_mod
         IMPLICIT NONE
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
         REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,
     $        mass_y
         END SUBROUTINE physics_mass
      END INTERFACE
c-----------------------------------------------------------------------
c     zero out values.
c-----------------------------------------------------------------------
      x(:,1) = x0
      y(:,1) = y0
      mass = zero
      mass_x = zero
      mass_y = zero
c-----------------------------------------------------------------------
c     initialize mass matrix for the standard du/dt+del.F=S form
c     and correct for static equations
c-----------------------------------------------------------------------
      DO iqty=1,SIZE(mass0,3)
         mass(iqty,iqty,:,:) = one
      ENDDO
c-----------------------------------------------------------------------
c     special conditions.
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("physics")
         CALL physics_mass(x,y,mass,mass_x,mass_y)
      END SELECT
c-----------------------------------------------------------------------
c     reorder indices for optimization purposes
c-----------------------------------------------------------------------
      DO iqty=1,SIZE(mass0,3)
         mass0(:,:,iqty) =  TRANSPOSE(mass(iqty,:,:,1))
         mass_x0(:,:,iqty) = TRANSPOSE(mass_x(iqty,:,:,1))
         mass_y0(:,:,iqty) = TRANSPOSE(mass_y(iqty,:,:,1))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_mass
c-----------------------------------------------------------------------
c     subprogram 10. job2_alloc.
c     allocate job2 objects.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_alloc(nodal,np)

      LOGICAL, INTENT(IN) :: nodal
      INTEGER, INTENT(IN) :: np

      INTEGER :: n,info
      TYPE(jacobi_type) :: grid_jac
c-----------------------------------------------------------------------
c     compute and factor basis matrix.
c-----------------------------------------------------------------------
      ALLOCATE(grid_bmat(0:np,0:np),grid_bpiv(0:np),inv_points(0:np))

      CALL jacobi_alloc(grid_jac,np,nodal,.FALSE.,"gll")

      SELECT CASE(grid_inv_type)
      CASE("jacobi")
         inv_points=grid_jac%qzero
      CASE("uniform")
         inv_points=two*(/(n,n=0,np)/)/REAL(np,r8) - one
      END SELECT

      DO n=0,np
         CALL jacobi_basis(inv_points(n),grid_jac)
         grid_bmat(n,:) = grid_jac%pb
      ENDDO
      
      CALL dgetrf(np+1,np+1,grid_bmat,np+1,grid_bpiv,info)
      CALL jacobi_dealloc(grid_jac)
c-----------------------------------------------------------------------
c     initialize CUBIT objects.
c-----------------------------------------------------------------------
      IF(grid_type=="cubit")CALL cubit_init(blk)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_alloc
c-----------------------------------------------------------------------
c     subprogram 11. job2_dealloc.
c     deallocate objects allocated in the physical modules.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_dealloc(jtype)

      CHARACTER(*), INTENT(IN) :: jtype
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_dealloc
         IMPLICIT NONE
         END SUBROUTINE physics_dealloc
      END INTERFACE
c-----------------------------------------------------------------------
c     deallocate objects.
c-----------------------------------------------------------------------
      SELECT CASE(jtype)
      CASE("beltrami")
         CALL beltrami_dealloc
         RETURN
      CASE("physics")
         CALL physics_dealloc
         DEALLOCATE(adapt_qty)
      END SELECT

      DEALLOCATE(grid_bmat,grid_bpiv,inv_points)
      IF(grid_type=="cubit")CALL cubit_dealloc(blk)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_dealloc
c-----------------------------------------------------------------------
c     subprogram 12. job2_bc_default.
c     computes default boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_bc_default(edges,static)

      TYPE(edge_type), DIMENSION(:) :: edges
      LOGICAL, DIMENSION(:), INTENT(IN) :: static

      INTEGER :: i,iqty
c-----------------------------------------------------------------------
c     set overall b.c. types
c-----------------------------------------------------------------------
      DO i=1,SIZE(edges)
         IF(edges(i)%edgenum == 0)edges(i)%bc_type="none"
         DO iqty=1,SIZE(static)
            IF(edges(i)%bc_type(iqty) /= "robin"
     $           .AND. edges(i)%bc_type(iqty) /= "robin+")
     $           edges(i)%static(iqty)=static(iqty)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_bc_default
c-----------------------------------------------------------------------
c     subprogram 13. job2_grid.
c     computes initial 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_grid(jtype,nx,ny,np,ixx,iyy,ksi,eta)

      CHARACTER(*), INTENT(IN) :: jtype
      INTEGER, INTENT(IN) :: nx,ny,np,ixx,iyy
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta

      INTEGER :: ix
      REAL(r8), DIMENSION(0:np) :: xi
      REAL(r8), DIMENSION(0:np,0:np) :: x_temp,y_temp
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_grid(x,y,ksi,eta)
         USE local_mod
         IMPLICIT NONE
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
         REAL(r8), DIMENSION(:,:), INTENT(OUT) :: ksi,eta
         END SUBROUTINE physics_grid
      END INTERFACE
c-----------------------------------------------------------------------
c     calculate physical coordinates (ksi,eta) from logical (x,y)
c-----------------------------------------------------------------------
      xi = inv_points
      DO ix=0,np
         x_temp(ix,:) = (ixx + (one+xi(ix))*half)/REAL(nx,r8)
         y_temp(:,ix) = (iyy + (one+xi(ix))*half)/REAL(ny,r8)
      ENDDO
      ksi=x_temp
      eta=y_temp

      IF(jtype=="beltrami")RETURN

      IF(grid_type=="cubit")THEN
         CALL cubit_grid(ixx,iyy,ksi,eta,blk)
         RETURN
      ENDIF

      SELECT CASE(jtype)
      CASE("physics")
         CALL physics_grid(x_temp,y_temp,ksi,eta)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_grid
c-----------------------------------------------------------------------
c     subprogram 14. job2_schur.
c     computes Schur complement.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE job2_schur(t,hfac,x0,y0,u0,ux0,uy0,
     $     fx_u0,fx_ux0,fx_uy0,fy_u0,fy_ux0,fy_uy0,s_u0,s_ux0,s_uy0)

      REAL(r8), INTENT(IN) :: t,hfac
      REAL(r8), DIMENSION(:), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: u0,ux0,uy0
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) ::
     $     fx_u0,fx_ux0,fx_uy0,fy_u0,fy_ux0,fy_uy0,s_u0,s_ux0,s_uy0

      INTEGER :: iqty
      REAL(r8), DIMENSION(SIZE(x0),1) :: x,y
      REAL(r8), DIMENSION(SIZE(u0,1),SIZE(x0),1) :: u,ux,uy
      REAL(r8), DIMENSION(SIZE(fx_u0,3),SIZE(fx_u0,2),SIZE(x0),1) :: 
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE physics_schur(t,hfac,x,y,u,ux,uy,
     $        fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
         USE local_mod
         IMPLICIT NONE
         REAL(r8), INTENT(IN) :: t,hfac
         REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
         REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
         REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $        fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy
         END SUBROUTINE physics_schur
      END INTERFACE
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      x(:,1)=x0
      y(:,1)=y0
      u(:,:,1)=u0
      ux(:,:,1)=ux0
      uy(:,:,1)=uy0
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
c     get job-dependent schur complement.
c-----------------------------------------------------------------------
      SELECT CASE(job_type)
      CASE("physics")
         CALL physics_schur(t,hfac,x,y,u,ux,uy,
     $        fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
      END SELECT
c-----------------------------------------------------------------------
c     reorder indices for optimization purposes
c-----------------------------------------------------------------------
      DO iqty=1,SIZE(fx_u0,3)
         fx_u0(:,:,iqty) = TRANSPOSE(fx_u(iqty,:,:,1))
         fx_ux0(:,:,iqty) = TRANSPOSE(fx_ux(iqty,:,:,1))
         fx_uy0(:,:,iqty) = TRANSPOSE(fx_uy(iqty,:,:,1))
         fy_u0(:,:,iqty) = TRANSPOSE(fy_u(iqty,:,:,1))
         fy_ux0(:,:,iqty) = TRANSPOSE(fy_ux(iqty,:,:,1))
         fy_uy0(:,:,iqty) = TRANSPOSE(fy_uy(iqty,:,:,1))
         s_u0(:,:,iqty) = TRANSPOSE(s_u(iqty,:,:,1))
         s_ux0(:,:,iqty) = TRANSPOSE(s_ux(iqty,:,:,1))
         s_uy0(:,:,iqty) = TRANSPOSE(s_uy(iqty,:,:,1))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE job2_schur
      END MODULE job2_mod
