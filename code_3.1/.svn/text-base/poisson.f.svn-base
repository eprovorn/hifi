c-----------------------------------------------------------------------
c     file poisson.f.
c     contains specifications for Poisson's equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. poisson_mod.
c     1. poisson_input.
c     2. poisson_init.
c     3. poisson_init_special.
c     4. poisson_boundary.
c     5. poisson_rhs.
c     6. poisson_drdu.
c     7. poisson_dealloc.
c-----------------------------------------------------------------------
c     subprogram 0. poisson_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE poisson_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: poisson_static=.TRUE.
      CHARACTER(8), PRIVATE :: rho_type="random"
      CHARACTER(16), PRIVATE :: init_type=" "
      INTEGER, PRIVATE :: mmax=10,nmax=10,seed=0
      REAL(r8), PRIVATE :: d0=1,d1=0,d2=0,e0=0,e1=0,e2=0,
     $     f12=0,f21=0,f22=1,vx=0,vy=0

      REAL(r8), DIMENSION(:,:), POINTER, PRIVATE :: rho,u0
      REAL(r8), DIMENSION(2,2), PRIVATE :: dmat

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. poisson_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $     dt,dtmax,tmax,nstep,gr_curve,du_diagnose)

      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax,gr_curve

      INTEGER :: myios
      REAL(r8) :: tfac
      REAL(r8), DIMENSION(2) :: kvec
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/poisson_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,
     $     dtmax,tmax,nstep,gr_curve,du_diagnose,init_type,
     $     mmax,nmax,seed,rho_type,poisson_static,d0,d1,d2,e0,e1,e2,
     $     f12,f21,f22,vx,vy
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=poisson_list,IOSTAT=myios)

      nqty=1
c-----------------------------------------------------------------------
c     rescale dt-related input.
c-----------------------------------------------------------------------
      dmat(1,1)=d0
      dmat(2,2)=1/d0
      dmat(1,2)=d1+d2
      dmat(2,1)=d1-d2
      kvec=pi*(/mmax,nmax/)
      tfac=1/SUM(kvec*MATMUL(dmat,kvec))
      dt=dt*tfac
      dtmax=dtmax*tfac
      tmax=tmax*tfac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson_input
c-----------------------------------------------------------------------
c     subprogram 2. poisson_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson_init(xpi,ypi,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: m,n
      REAL(r8), DIMENSION(SIZE(xpi,1),SIZE(xpi,2),mmax) :: sinx
      REAL(r8), DIMENSION(SIZE(xpi,1),SIZE(xpi,2),nmax) :: siny
c-----------------------------------------------------------------------
c     compute sin factors.
c-----------------------------------------------------------------------
      DO m=1,mmax
         sinx(:,:,m)=SIN(m*pi*xpi)
      ENDDO
      DO n=1,nmax 
         siny(:,:,n)=SIN(n*pi*ypi)
      ENDDO
c-----------------------------------------------------------------------
c     compute initial conditions.
c-----------------------------------------------------------------------
      ui=0
      DO m=1,mmax
         DO n=1,nmax
            ui(1,:,:)=ui(1,:,:)+u0(m,n)*sinx(:,:,m)*siny(:,:,n)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson_init
c-----------------------------------------------------------------------
c     subprogram 3. poisson_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson_init_special(static)

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: m,n,i
      REAL(r8) :: kfac,ddet
      REAL(r8), DIMENSION(2) :: kvec
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/5x,"i",6x,"d1",9x,"d2"/)
 20   FORMAT(i6,1p,2e11.3)
 30   FORMAT(/5x,"m",5x,"n",5x,"rho",9x,"k1",9x,"k2",7x,"kfac",
     $     8x,"u0"/)
 40   FORMAT(2i6,1p,5e11.3)
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(rho_type,8,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(poisson_static,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast integer input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(mmax,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(nmax,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(seed,1,MPI_INTEGER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(d0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(d1,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(d2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(e0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(e1,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(e2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(f12,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(f21,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(f22,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(vx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(vy,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     set static.
c-----------------------------------------------------------------------
      static(1)=poisson_static
c-----------------------------------------------------------------------
c     define dmat.
c-----------------------------------------------------------------------
      dmat(1,1)=d0
      dmat(2,2)=1/d0
      dmat(1,2)=d1+d2
      dmat(2,1)=d1-d2
      ddet=dmat(1,1)*dmat(2,2)-dmat(1,2)*dmat(2,1)
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0)THEN
         WRITE(*,'(1p,3(a,e10.3))')" d0 = ",d0,", d1 = ",d1,", d2 = ",d2
         WRITE(*,'(1p,3(a,e10.3))')" vx = ",vx,", vy = ",vy,
     $        ", ddet = ",ddet
         WRITE(*,'(/5x,"i",2(4x,"dmati",i1,1x)/)')(i,i=1,2)
         WRITE(*,'(i6,1p,2e11.3)')(i,dmat(i,:),i=1,2)
      ENDIF
c-----------------------------------------------------------------------
c     allocate and fill rho.
c-----------------------------------------------------------------------
      ALLOCATE(rho(mmax,nmax))
      CALL put_seed(seed)
      SELECT CASE(rho_type)
      CASE("random")
         CALL RANDOM_NUMBER(rho)
         rho=2*rho-1
      CASE("one")
         rho=0
         rho(mmax,nmax)=1
      CASE DEFAULT
         CALL program_stop("poisson_init_special: cannot recognize "
     $        //"rho_type = "//TRIM(rho_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     allocate u0 and fill with solution.
c-----------------------------------------------------------------------
      ALLOCATE(u0(mmax,nmax))
      SELECT CASE(init_type)
      CASE("solve")
         IF(mpi_rank == 0 .AND. diagnose)THEN
            OPEN(UNIT=debug_unit,FILE="solve.out",STATUS="UNKNOWN")
            WRITE(debug_unit,'(4a,2(a,i1))')
     $           " rho_type = ",TRIM(rho_type),
     $           ", init_type = ",TRIM(init_type),
     $           ", mmax = ",mmax,", nmax = ",nmax
            WRITE(debug_unit,10)
            WRITE(debug_unit,20)(i,dmat(i,:),i=1,2)
            WRITE(debug_unit,10)
            WRITE(debug_unit,30)
         ENDIF
         DO m=1,mmax
            DO n=1,nmax
               kvec=pi*(/m,n/)
               kfac=SUM(kvec*MATMUL(dmat,kvec))
               u0(m,n)=rho(m,n)/kfac
               IF(mpi_rank == 0 .AND. diagnose)
     $              WRITE(debug_unit,40)m,n,rho(m,n),kvec,kfac,u0(m,n)
            ENDDO
         ENDDO
         IF(mpi_rank == 0 .AND. diagnose)THEN
            WRITE(debug_unit,30)
            CLOSE(UNIT=debug_unit)
         ENDIF
c-----------------------------------------------------------------------
c     compute other intializations.
c-----------------------------------------------------------------------
      CASE("random")
         CALL RANDOM_NUMBER(u0)
      CASE("zero")
         u0=0
      CASE DEFAULT
         CALL program_stop("poisson_special_init: cannot recognize "
     $        //"init_type = "//TRIM(init_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson_init_special
c-----------------------------------------------------------------------
c     subprogram 4. poisson_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson_boundary(left,right,top,bottom)

      TYPE(edge_type), INTENT(INOUT) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      left%a=1
      right%a=1
      top%a=1
      bottom%a=1
      left%b=0
      right%b=0
      top%b=0
      bottom%b=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson_boundary
c-----------------------------------------------------------------------
c     subprogram 5. poisson_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson_rhs(x,y,u,ux,uy,fx,fy,s)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s

      INTEGER :: m,n
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2),mmax) :: sinx
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2),nmax) :: siny
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     compute sin factors.
c-----------------------------------------------------------------------
      DO m=1,mmax
         sinx(:,:,m)=SIN(m*pi*x)
      ENDDO
      DO n=1,nmax 
         siny(:,:,n)=SIN(n*pi*y)
      ENDDO
c-----------------------------------------------------------------------
c     compute source terms.
c-----------------------------------------------------------------------
      s=0
      DO m=1,mmax
         DO n=1,nmax
            s(1,:,:)=s(1,:,:)+rho(m,n)*sinx(:,:,m)*siny(:,:,n)
         ENDDO
       ENDDO
c-----------------------------------------------------------------------
c     compute flux terms.
c-----------------------------------------------------------------------
      fx(1,:,:)=-dmat(1,1)*ux(1,:,:)-dmat(1,2)*uy(1,:,:)+vx*u(1,:,:)
      fy(1,:,:)=-dmat(2,1)*ux(1,:,:)-dmat(2,2)*uy(1,:,:)+vy*u(1,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson_rhs
c-----------------------------------------------------------------------
c     subprogram 6. poisson_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson_drdu(x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)

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
c     compute derivatives of fluxes.
c-----------------------------------------------------------------------
      fx_u(1,1,:,:)=vx
      fx_ux(1,1,:,:)=-dmat(1,1)
      fx_uy(1,1,:,:)=-dmat(1,2)
      fy_u(1,1,:,:)=vy
      fy_ux(1,1,:,:)=-dmat(2,1)
      fy_uy(1,1,:,:)=-dmat(2,2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson_drdu
c-----------------------------------------------------------------------
c     subprogram 7. poisson_dealloc.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson_dealloc
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      IF(ASSOCIATED(rho))DEALLOCATE(rho,u0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson_dealloc
      END MODULE poisson_mod
