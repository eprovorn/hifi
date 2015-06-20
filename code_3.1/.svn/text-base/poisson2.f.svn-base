c-----------------------------------------------------------------------
c     file poisson2.f.
c     contains specifications for Poisson's equation with nqty = 2.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. poisson2_mod.
c     1. poisson2_input.
c     2. poisson2_init.
c     3. poisson2_init_special.
c     4. poisson2_boundary.
c     5. poisson2_rhs.
c     6. poisson2_drdu.
c     7. poisson2_dealloc.
c     8. poisson2_solve.
c-----------------------------------------------------------------------
c     subprogram 0. poisson2_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE poisson2_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL, PRIVATE :: poisson_static=.TRUE.
      CHARACTER(8), PRIVATE :: rho_type="random"
      CHARACTER(16), PRIVATE :: init_type=" "
      INTEGER, PRIVATE :: mmax=10,nmax=10,seed=0
      REAL(r8), PRIVATE :: d0=1,d1=0,d2=0,e0=0,e1=0,e2=0,
     $     f12=0,f21=0,f22=1,vx=0,vy=0

      REAL(r8), DIMENSION(:,:,:), POINTER, PRIVATE :: rho,u0
      REAL(r8), DIMENSION(4,4), TARGET, PRIVATE :: dmat
      REAL(r8), DIMENSION(:,:), POINTER, PRIVATE :: d11,d12,d21,d22

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. poisson2_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson2_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,
     $     nqty,dt,dtmax,tmax,nstep,gr_curve,du_diagnose)

      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax,gr_curve

      INTEGER :: myios
      REAL(r8) :: tfac1,tfac2
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

      nqty=2
c-----------------------------------------------------------------------
c     rescale dt-related input.
c-----------------------------------------------------------------------
      d11(1,1)=d0
      d11(2,2)=1/d0
      d11(1,2)=d1+d2
      d11(2,1)=d1-d2
      d22=d11*f22
      kvec=pi*(/mmax,nmax/)
      tfac1=1/SUM(kvec*MATMUL(d11,kvec))
      tfac2=1/SUM(kvec*MATMUL(d22,kvec))
      tfac1=MIN(tfac1,tfac2)
      dt=dt*tfac1
      dtmax=dtmax*tfac1
      tmax=tmax*tfac1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson2_input
c-----------------------------------------------------------------------
c     subprogram 2. poisson2_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson2_init(xpi,ypi,ui)

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: m,n,iqty
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
      DO iqty=1,2
         DO m=1,mmax
            DO n=1,nmax
               ui(iqty,:,:)=ui(iqty,:,:)
     $              +u0(iqty,m,n)*sinx(:,:,m)*siny(:,:,n)
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson2_init
c-----------------------------------------------------------------------
c     subprogram 3. poisson2_init_special.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson2_init_special(static)

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static

      LOGICAL, PARAMETER :: diagnose=.TRUE.
      INTEGER :: m,n,i,info
      INTEGER, DIMENSION(2) :: ipiv
      REAL(r8) :: ddet
      REAL(r8), DIMENSION(4,4) :: dmat0
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
      static=poisson_static
c-----------------------------------------------------------------------
c     access individual quadrants.
c-----------------------------------------------------------------------
      d11 => dmat(1:2,1:2)
      d12 => dmat(1:2,3:4)
      d21 => dmat(3:4,1:2)
      d22 => dmat(3:4,3:4)
c-----------------------------------------------------------------------
c     define diagonal blocks.
c-----------------------------------------------------------------------
      d11(1,1)=d0
      d11(2,2)=1/d0
      d11(1,2)=d1+d2
      d11(2,1)=d1-d2
      d22=d11*f22
c-----------------------------------------------------------------------
c     define off-diagonal blocks.
c-----------------------------------------------------------------------
      d12(1,1)=e0
      d12(2,2)=1/e0
      d12(1,2)=e1+e2
      d12(2,1)=e1-e2
      d21=TRANSPOSE(d12)*f21
      d12=d12*f12
c-----------------------------------------------------------------------
c     compute ddet.
c-----------------------------------------------------------------------
      dmat0=dmat
      CALL dgetrf(4,4,dmat0,4,ipiv,info)
      ddet=PRODUCT((/(dmat0(i,i),i=1,4)/))
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0)THEN
         WRITE(*,'(1p,4(a,e10.3))')" d0 = ",d0,", d1 = ",d1,
     $        ", d2 = ",d2,", f12 = ",f12
         WRITE(*,'(1p,4(a,e10.3))')" e0 = ",e0,", e1 = ",e1,
     $        ", e2 = ",e2,", f21 = ",f21
         WRITE(*,'(1p,2(a,e10.3))')" f22 = ",f22,", ddet = ",ddet
         WRITE(*,'(/5x,"i",4(4x,"dmati",i1,1x)/)')(i,i=1,4)
         WRITE(*,'(i6,1p,4e11.3)')(i,dmat(i,:),i=1,4)
      ENDIF
c-----------------------------------------------------------------------
c     allocate and fill rho.
c-----------------------------------------------------------------------
      ALLOCATE(rho(2,mmax,nmax))
      CALL put_seed(seed)
      SELECT CASE(rho_type)
      CASE("random")
         CALL RANDOM_NUMBER(rho)
         rho=2*rho-1
      CASE("one")
         rho=0
         rho(:,mmax,nmax)=1
      CASE DEFAULT
         CALL program_stop("poisson2_init_special: cannot recognize "
     $        //"rho_type = "//TRIM(rho_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     allocate u0 and fill with random numbers.
c-----------------------------------------------------------------------
      ALLOCATE(u0(2,mmax,nmax))
      SELECT CASE(init_type)
      CASE("solve")
         DO m=1,mmax
            DO n=1,nmax
               u0(:,m,n)=poisson2_solve(m,n)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     compute other intializations.
c-----------------------------------------------------------------------
      CASE("random")
         CALL RANDOM_NUMBER(u0)
      CASE("zero")
         u0=0
      CASE DEFAULT
         CALL program_stop("poisson2_special_init: cannot recognize "
     $        //"init_type = "//TRIM(init_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson2_init_special
c-----------------------------------------------------------------------
c     subprogram 4. poisson2_boundary.
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson2_boundary(left,right,top,bottom)

      TYPE(edge_type), INTENT(INOUT) :: left,right,top,bottom
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      left%a=RESHAPE((/1,0,0,1/),(/2,2/))
      right%a=RESHAPE((/1,0,0,1/),(/2,2/))
      top%a=RESHAPE((/1,0,0,1/),(/2,2/))
      bottom%a=RESHAPE((/1,0,0,1/),(/2,2/))
      left%b=0
      right%b=0
      top%b=0
      bottom%b=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson2_boundary
c-----------------------------------------------------------------------
c     subprogram 5. poisson2_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson2_rhs(x,y,u,ux,uy,fx,fy,s)

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
            s(1,:,:)=s(1,:,:)+sinx(:,:,m)*siny(:,:,n)*rho(1,m,n)
            s(2,:,:)=s(2,:,:)+sinx(:,:,m)*siny(:,:,n)*rho(2,m,n)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute flux terms.
c-----------------------------------------------------------------------
      fx(1,:,:)=-dmat(1,1)*ux(1,:,:)-dmat(1,2)*uy(1,:,:)
     $     -dmat(1,3)*ux(2,:,:)-dmat(1,4)*uy(2,:,:)
      fy(1,:,:)=-dmat(2,1)*ux(1,:,:)-dmat(2,2)*uy(1,:,:)
     $     -dmat(2,3)*ux(2,:,:)-dmat(2,4)*uy(2,:,:)
      fx(2,:,:)=-dmat(3,1)*ux(1,:,:)-dmat(3,2)*uy(1,:,:)
     $     -dmat(3,3)*ux(2,:,:)-dmat(3,4)*uy(2,:,:)
      fy(2,:,:)=-dmat(4,1)*ux(1,:,:)-dmat(4,2)*uy(1,:,:)
     $     -dmat(4,3)*ux(2,:,:)-dmat(4,4)*uy(2,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson2_rhs
c-----------------------------------------------------------------------
c     subprogram 6. poisson2_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson2_drdu(x,y,u,ux,uy,
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
      fx_ux(1,1,:,:)=-dmat(1,1)
      fx_uy(1,1,:,:)=-dmat(1,2)
      fx_ux(1,2,:,:)=-dmat(1,3)
      fx_uy(1,2,:,:)=-dmat(1,4)
      fy_ux(1,1,:,:)=-dmat(2,1)
      fy_uy(1,1,:,:)=-dmat(2,2)
      fy_ux(1,2,:,:)=-dmat(2,3)
      fy_uy(1,2,:,:)=-dmat(2,4)
      fx_ux(2,1,:,:)=-dmat(3,1)
      fx_uy(2,1,:,:)=-dmat(3,2)
      fx_ux(2,2,:,:)=-dmat(3,3)
      fx_uy(2,2,:,:)=-dmat(3,4)
      fy_ux(2,1,:,:)=-dmat(4,1)
      fy_uy(2,1,:,:)=-dmat(4,2)
      fy_ux(2,2,:,:)=-dmat(4,3)
      fy_uy(2,2,:,:)=-dmat(4,4)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson2_drdu
c-----------------------------------------------------------------------
c     subprogram 7. poisson2_dealloc.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE poisson2_dealloc
c-----------------------------------------------------------------------
c     initialize grid coordinates
c-----------------------------------------------------------------------
      IF(ASSOCIATED(rho))DEALLOCATE(rho,u0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE poisson2_dealloc
c-----------------------------------------------------------------------
c     subprogram 8. poisson2_solve.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION poisson2_solve(m,n) RESULT(u)

      INTEGER, INTENT(IN) :: m,n
      REAL(r8), DIMENSION(2) :: u

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: info,i
      INTEGER, DIMENSION(2) :: ipiv
      REAL(r8), DIMENSION(4,2) :: amat
      REAL(r8), DIMENSION(2,2) :: bmat,bmat0
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/5x,"i",4(4x,"dmat",i1,2x)/)
 20   FORMAT(i6,1p,4e11.3)
 30   FORMAT(/5x,"i",2(4x,"amat",i1,2x)/)
 40   FORMAT(i6,1p,2e11.3)
 50   FORMAT(/5x,"i",2(4x,"bmat",i1,2x),5x,"rho"/)
 60   FORMAT(i6,1p,3e11.3)
 70   FORMAT(/5x,"i",2(4x,"bfac",i1,2x),6x,"u"/)
 80   FORMAT(i6,1p,3e11.3)
c-----------------------------------------------------------------------
c     solve.
c-----------------------------------------------------------------------
      amat=RESHAPE((/m,n,0,0,0,0,m,n/),(/4,2/))
      bmat=pi**2*MATMUL(TRANSPOSE(amat),MATMUL(dmat,amat))
      bmat0=bmat
      u=rho(:,m,n)
      CALL dgetrf(2,2,bmat,2,ipiv,info)
      CALL dgetrs('N',2,1,bmat,2,ipiv,u,2,info)
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0 .AND. diagnose)THEN
         OPEN(UNIT=debug_unit,FILE="solve.out",STATUS="UNKNOWN")
         WRITE(debug_unit,'(2(a,i2))')" m = ",m,", n = ",n
         WRITE(debug_unit,10)(i,i=1,4)
         WRITE(debug_unit,20)(i,dmat(i,:),i=1,4)
         WRITE(debug_unit,10)(i,i=1,4)
         WRITE(debug_unit,30)(i,i=1,2)
         WRITE(debug_unit,40)(i,amat(i,:),i=1,4)
         WRITE(debug_unit,30)(i,i=1,2)
         WRITE(debug_unit,50)(i,i=1,2)
         WRITE(debug_unit,60)(i,bmat0(i,:),rho(i,m,n),i=1,2)
         WRITE(debug_unit,50)(i,i=1,2)
         WRITE(debug_unit,70)(i,i=1,2)
         WRITE(debug_unit,80)(i,bmat(i,:),u(i),i=1,2)
         WRITE(debug_unit,70)(i,i=1,2)
         CLOSE(UNIT=debug_unit)
         CALL program_stop("Termination by poisson2_solve.")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION poisson2_solve
      END MODULE poisson2_mod
