c-----------------------------------------------------------------------
c     file local.F.
c     local defintions for most computers.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. local.
c     1. timer.
c     2. program_stop.
c     3. edge_alloc.
c     4. edge_dealloc.
c     5. gather.
c     6. log10_2d.
c     7. log10_3d.
c     8. put_seed.
c-----------------------------------------------------------------------
c     subprogram 0. local.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE local_mod
      USE io_mod
      IMPLICIT NONE

#include "finclude/petsc.h"
c-----------------------------------------------------------------------
c     numeric constants declaration
c-----------------------------------------------------------------------
      INTEGER, PARAMETER :: ndim=2,r8=SELECTED_REAL_KIND(13,307)
      REAL(r8), PARAMETER :: pi=3.1415926535897932385_r8,twopi=2*pi,
     $     pisq=pi*pi,rtod=180./pi,dtor=pi/180._r8,
     $     zero=0._r8,half=.5_r8,one=1._r8,two=2._r8,
     $     third=1./3._r8,neg_one=-1._r8,min_eps=1e-8
c-----------------------------------------------------------------------
c     variable declaration
c-----------------------------------------------------------------------
      MPI_Comm :: comm
      INTEGER :: mpi_rank=0,mpi_size=1,ierr=0

      TYPE :: edge_type
      INTEGER :: edgenum
      CHARACTER(6) :: edgename
      LOGICAL, DIMENSION(:), POINTER :: static
      CHARACTER(10), DIMENSION(:), POINTER :: bc_type
      END TYPE edge_type

      CHARACTER(10) :: date,time,zone
      CHARACTER(80) :: outdir=".",out_name="."
      LOGICAL :: monitor=.FALSE.
      INTEGER, DIMENSION(8) :: clock=0
      REAL(r8) :: wclock=0,wclock0=0

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. timer.
c     handles machine-dependent timing statistics.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      SUBROUTINE timer(out_mode,out_name,out_unit)
      
      CHARACTER(*), INTENT(IN) :: out_name
      INTEGER, INTENT(IN) :: out_mode,out_unit

      REAL(r8), SAVE :: seconds
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(a,1p,e10.3,a)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      CALL DATE_AND_TIME(date,time,zone,clock)

      IF(out_mode == 0)THEN
         seconds=(clock(5)*60+clock(6))*60+clock(7)+clock(8)/1000._r8
      ELSE
         seconds=(clock(5)*60+clock(6))*60+clock(7)+clock(8)/1000._r8
     $        -seconds
         WRITE(*,10)" Wallclock time = ",seconds," seconds"//CHAR(7)
         IF(monitor)THEN
            OPEN(UNIT=out_unit,FILE=TRIM(out_name),ACTION="WRITE",
     $           STATUS="OLD",POSITION="APPEND")
            WRITE(out_unit,10)" Wallclock time = ",seconds," seconds"
            CLOSE(UNIT=out_unit)
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     subprogram 2. program_stop.
c     terminates program with message, calls timer, closes output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE program_stop(message)

      CHARACTER(*), INTENT(IN) :: message
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0)THEN
         CALL timer(1,out_name,out_unit)
         WRITE(*,'(1x,a)')" PROGRAM_STOP => "//TRIM(message)
         IF(monitor)THEN         
            OPEN(UNIT=out_unit,FILE=TRIM(out_name),ACTION="WRITE",
     $           POSITION="APPEND")
            WRITE(out_unit,'(1x,a)')" PROGRAM_STOP => "//TRIM(message)
            CLOSE(UNIT=out_unit)
         ENDIF
      ENDIF
      CALL PetscFinalize(ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE program_stop
c-----------------------------------------------------------------------
c     subprogram 3. edge_alloc.
c     allocate edge_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE edge_alloc(nqty,n,edges)

      INTEGER, INTENT(IN) :: nqty,n
      TYPE(edge_type), DIMENSION(:), POINTER :: edges

      INTEGER :: i
c-----------------------------------------------------------------------
c     allocate fields of edges
c-----------------------------------------------------------------------
      ALLOCATE(edges(n))
      DO i=1,n
         edges(i)%edgenum=0
         SELECT CASE(i)
         CASE(1)
            edges(i)%edgename="left"
         CASE(2)
            edges(i)%edgename="top"
         CASE(3)
            edges(i)%edgename="right"
         CASE(4)
            edges(i)%edgename="bottom"
         END SELECT
         ALLOCATE(edges(i)%static(nqty),edges(i)%bc_type(nqty))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_alloc
c-----------------------------------------------------------------------
c     subprogram 4. edge_dealloc.
c     deallocate edge_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE edge_dealloc(edges)

      TYPE(edge_type), DIMENSION(:), POINTER :: edges

      INTEGER :: i
c-----------------------------------------------------------------------
c     deallocate fields of edge
c-----------------------------------------------------------------------
      DO i=1,SIZE(edges)
         DEALLOCATE(edges(i)%static,edges(i)%bc_type)
      ENDDO
      DEALLOCATE(edges)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_dealloc
c-----------------------------------------------------------------------
c     subprogram 5. gather.
c     gathers data from all processors.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE gather(scalar,vector)

      INTEGER, INTENT(IN) :: scalar
      INTEGER, DIMENSION(0:), INTENT(INOUT) :: vector

      INTEGER, DIMENSION(0:mpi_size-1) :: temp
c-----------------------------------------------------------------------
c     gather data from all processors.
c-----------------------------------------------------------------------
      temp=0
      vector=0
      temp(mpi_rank)=scalar
      CALL MPI_Allreduce(temp,vector,mpi_size,MPI_INTEGER,MPI_SUM,
     $     comm,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gather
c-----------------------------------------------------------------------
c     subprogram 6. log10_2d.
c     takes log10 of 2d array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION log10_2d(a) RESULT(b)

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: a
      REAL(r8), DIMENSION(SIZE(a,1),SIZE(a,2)) :: b
c-----------------------------------------------------------------------
c     compute log10(a).
c-----------------------------------------------------------------------
      b=a
      WHERE(b <= min_eps)
         b=min_eps
      ENDWHERE
      b=LOG10(b)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION log10_2d
c-----------------------------------------------------------------------
c     subprogram 7. log10_3d.
c     takes log10 of 3d array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION log10_3d(a) RESULT(b)

      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: a
      REAL(r8), DIMENSION(SIZE(a,1),SIZE(a,2),SIZE(a,3)) :: b
c-----------------------------------------------------------------------
c     compute log10(a).
c-----------------------------------------------------------------------
      b=a
      WHERE(b <= min_eps)
         b=min_eps
      ENDWHERE
      b=LOG10(b)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION log10_3d
c-----------------------------------------------------------------------
c     subprogram 8. put_seed.
c     deposits random number seed.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE put_seed(seed)

      INTEGER, INTENT(IN) :: seed

      INTEGER :: seed_size
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed_array
c-----------------------------------------------------------------------
c     deposit seed.
c-----------------------------------------------------------------------
      CALL RANDOM_SEED(SIZE=seed_size)
      ALLOCATE(seed_array(seed_size))
      CALL RANDOM_SEED(GET=seed_array)
      seed_array(1)=seed_array(1)*seed
      CALL RANDOM_SEED(PUT=seed_array)
      DEALLOCATE(seed_array)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE put_seed
      END MODULE local_mod
