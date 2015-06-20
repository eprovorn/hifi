c-----------------------------------------------------------------------
c     file cubit.f.
c     contains various cubit grid routines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. cubit_mod.
c     1. cubit_init.
c     2. cubit_read_inp.
c     3. cubit_read_cdf.
c     4. check.
c     5. cubit_write.
c     6. cubit_mapping.
c     7. cubit_grid.
c     8. cubit_alloc.
c     9. cubit_dealloc.
c-----------------------------------------------------------------------
c     subprogram 0. cubit_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE cubit_mod
      USE netcdf
      USE local_mod
      IMPLICIT NONE

      CHARACTER(8)   :: header1
      CHARACTER(100) :: header2
      CHARACTER(5)   :: header3
      CHARACTER(8)   :: header4
      CHARACTER(9)   :: header5
      CHARACTER(10)  :: header6

      TYPE :: block
      INTEGER :: blocks,nnodes,nels,cxns,nx,ny,np,align
      REAL(r8),DIMENSION(:,:), POINTER :: coords
      INTEGER, DIMENSION(:,:,:), POINTER :: connect,mapping
      END TYPE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. cubit_init.
c     initialize cubit data structures
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cubit_init(blk)
      
      TYPE(block) :: blk

      INTEGER :: bcast_size
c-----------------------------------------------------------------------
c     broadcast CUBIT grid sizes.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(blk%blocks,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(blk%cxns,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(blk%nnodes,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(blk%nels,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(blk%align,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(blk%np,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(blk%nx,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(blk%ny,1,MPI_INTEGER,0,comm,ierr)
c-----------------------------------------------------------------------
c     allocate CUBIT variables.
c-----------------------------------------------------------------------
      IF (mpi_rank /= 0)CALL cubit_alloc(blk)
c-----------------------------------------------------------------------
c     broadcast CUBIT grid data.
c-----------------------------------------------------------------------
      bcast_size = SIZE(blk%coords)
      CALL MPI_Bcast(blk%coords,bcast_size,MPI_DOUBLE_PRECISION,
     $     0,comm,ierr)
      bcast_size = SIZE(blk%connect)
      CALL MPI_Bcast(blk%connect,bcast_size,MPI_INTEGER,0,comm,ierr)
c-----------------------------------------------------------------------
c     generate mapping from 1D logical CUBIT data to 2D data
c-----------------------------------------------------------------------
      CALL cubit_mapping(blk)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cubit_init
c-----------------------------------------------------------------------
c     subprogram 2. cubit_read_inp.
c     read in a .inp (abaqus) file
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cubit_read_inp(cubit_file,blk)

      CHARACTER(*), INTENT(IN) :: cubit_file
      TYPE(block) :: blk

      INTEGER :: node,i,j,dummy
      REAL(r8) :: coord1, coord2
      INTEGER :: cubit_unit=28248
c-----------------------------------------------------------------------
c     specify parameters
c-----------------------------------------------------------------------
      blk%coords=0
      blk%connect=0
c-----------------------------------------------------------------------
c     open and read a ***.inp file (abaqus)
c-----------------------------------------------------------------------
      OPEN(UNIT=cubit_unit,FILE=TRIM(cubit_file),
     $     ACTION="READ",STATUS="OLD")
      REWIND(UNIT=cubit_unit)

      READ(cubit_unit,"(A)",IOSTAT=ierr) header1
      READ(cubit_unit,"(A)",IOSTAT=ierr) header2
      READ(cubit_unit,"(A)",IOSTAT=ierr) header3

      DO i = 1,blk%nnodes
          READ(cubit_unit,*,IOSTAT=ierr) node,coord1,coord2
          blk%coords(node,1) = coord1
          blk%coords(node,2) = coord2
      END DO
  
      DO j = 1,blk%blocks
      READ(cubit_unit,*,IOSTAT=ierr) header4, header5, header6
         DO i = 1,blk%nels
            IF (blk%cxns==4) THEN
              READ(cubit_unit,*,IOSTAT=ierr)
     $              dummy,blk%connect(j,1,i),blk%connect(j,2,i),
     $              blk%connect(j,3,i),blk%connect(j,4,i)
            ELSE
               CALL program_stop('These are not quadrilaterals')
            END IF

         END DO
      END DO

      CLOSE(cubit_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cubit_read_inp
c-----------------------------------------------------------------------
c     subprogram 3. cubit_read_cdf.
c     read in a .cdf (netcdf) file
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cubit_read_cdf(cubit_file,blk,exit_flag)

      CHARACTER(*), INTENT(IN) :: cubit_file
      INTEGER, INTENT(INOUT) :: exit_flag
      TYPE(block) :: blk

      INTEGER :: i,bcast_size
      INTEGER :: ncid,dimID,elemID,nodesID,npeID,nblkID
      INTEGER :: num_dim,num_elem,num_nodes,nodes_pe,num_blks
      INTEGER :: Cx1VarID,CoordVarID
c-----------------------------------------------------------------------
c     open a ***.cdf file (netcdf) of name 'cubit_file'
c-----------------------------------------------------------------------
      CALL check( NF90_OPEN(TRIM(cubit_file), NF90_NOWRITE, ncid) )
c-----------------------------------------------------------------------
c     get variable dimensions.
c-----------------------------------------------------------------------
      CALL check(NF90_INQ_DIMID(ncid,"num_dim",dimID))
      CALL check(NF90_INQ_DIMID(ncid,"num_elem",elemID))
      CALL check(NF90_INQ_DIMID(ncid,"num_nodes",nodesID))
      CALL check(NF90_INQ_DIMID(ncid,"num_nod_per_el1",npeID))
      CALL check(NF90_INQ_DIMID(ncid,"num_el_blk",nblkID))
      
      CALL check(NF90_INQUIRE_DIMENSION(ncid,dimID,len = num_dim))
      CALL check(NF90_INQUIRE_DIMENSION(ncid,elemID,len = num_elem))
      CALL check(NF90_INQUIRE_DIMENSION(ncid,nodesID,len=num_nodes))
      CALL check(NF90_INQUIRE_DIMENSION(ncid,npeID,len = nodes_pe))
      CALL check(NF90_INQUIRE_DIMENSION(ncid,nblkID,len = num_blks))
      
      IF (num_blks > 1)exit_flag=1
      IF (num_dim /= 2)exit_flag=2
c-----------------------------------------------------------------------
c     store grid data in 'blk' data structure.
c-----------------------------------------------------------------------
      blk%blocks = num_blks
      blk%cxns = nodes_pe
      blk%nnodes = num_nodes
      blk%nels = num_elem
c-----------------------------------------------------------------------
c     allocate variables for coordinate and grid connectivity
c-----------------------------------------------------------------------
      ALLOCATE(blk%connect(blk%blocks,blk%cxns,blk%nels),
     $     blk%coords(blk%nnodes,ndim))
      blk%coords=0
      blk%connect=0
c-----------------------------------------------------------------------
c     read in the coordinate and element numbering from netCDF file and
c     then save data to the 'blk' data structure.
c-----------------------------------------------------------------------
      CALL check( NF90_INQ_VARID(ncid,"connect1",Cx1VarID))
      CALL check( NF90_INQ_VARID(ncid,"coord",CoordVarID))         
      CALL check( NF90_GET_VAR(ncid,CoordVarID,blk%coords))
      CALL check( NF90_GET_VAR(ncid,Cx1VarID,blk%connect(1,:,:)))
c-----------------------------------------------------------------------
c     determine which way the elements are aligned.
c-----------------------------------------------------------------------
      IF(blk%connect(1,1,1)==1)THEN
         blk%align=1
      ELSE
         blk%align=2
      ENDIF
c-----------------------------------------------------------------------
c     determine which dimension the grid spans first.
c-----------------------------------------------------------------------
      i=1
      DO
         IF(COUNT((blk%connect(1,:,:)==blk%connect(1,2,i)))==1)EXIT
         i=i+1
      ENDDO
      IF (blk%align==1)THEN
         blk%nx=i/blk%np
         blk%ny=num_elem/(blk%np**ndim*blk%nx)
      ELSE
         blk%ny=i/blk%np
         blk%nx=num_elem/(blk%np**ndim*blk%ny)
      ENDIF
c-----------------------------------------------------------------------
c     close the netCDF file and print success statement to screen.
c-----------------------------------------------------------------------
      CALL check( NF90_CLOSE(ncid) )
      PRINT *,"Success Reading Grid File: ",TRIM(cubit_file),"! "
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cubit_read_cdf
c-----------------------------------------------------------------------
c     subprogram 4. check
c     checks the status of netCDF subroutines
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE check(status)

      INTEGER, INTENT(IN) :: status
    
      IF(status /= nf90_noerr) THEN
         PRINT *, trim(nf90_strerror(status))
         STOP "Stopped"
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE check
c-----------------------------------------------------------------------
c     subprogram 5. cubit_write.
c     write the .inp (abaqus) data to screen
c-----------------------------------------------------------------------
      SUBROUTINE cubit_write(blk)
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      TYPE(block) :: blk

      INTEGER :: i,j
c-----------------------------------------------------------------------
c     write (abaqus) data to screen
c-----------------------------------------------------------------------
      WRITE(*,"(A)") header1
      WRITE(*,"(A)") header2
      WRITE(*,"(A)") header3

      DO i = 1,blk%nnodes
         WRITE(*,"(I8,ES15.6,ES15.6)")
     $        i,blk%coords(i,1), blk%coords(i,2)
      END DO
  
      DO j = 1,blk%blocks
      WRITE(*,"(A,3x,A,3x,A)") header4, header5, header6
         DO i = 1,blk%nels
            IF (blk%cxns==4) THEN
              WRITE(*,"(I8,I8,I8,I8,I8)") 
     $              i,blk%connect(j,1,i),blk%connect(j,2,i),
     $              blk%connect(j,3,i),blk%connect(j,4,i)
            ELSE
               CALL program_stop('These are not quadrilaterals')
            END IF

         END DO
      END DO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cubit_write
c-----------------------------------------------------------------------
c     subprogram 6. cubit_mapping.
c     map the logically 1D grid data into a 2D array
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cubit_mapping(blk)

      TYPE(block) :: blk

      INTEGER :: i,j,nx,ny,np,elem_num
c-----------------------------------------------------------------------
c     create logical 1D to 2D mapping array
c-----------------------------------------------------------------------
      np=blk%np
      nx=blk%nx
      ny=blk%ny
      ALLOCATE(blk%mapping(nx*np+1,ny*np+1,2))
      blk%mapping=0
      elem_num=0
      SELECT CASE(blk%align)
      CASE(1)
         DO j = 1,ny*np+1
            DO i = 1,nx*np+1
               IF (i==(nx*np+1) .AND. j/=(ny*np+1))THEN
                  blk%mapping(i,j,1)=elem_num
                  blk%mapping(i,j,2)=2
               ELSEIF (j==(ny*np+1) .AND. i/=(nx*np+1))THEN
                  elem_num=elem_num+1
                  blk%mapping(i,j,1)=elem_num-nx*np
                  blk%mapping(i,j,2)=4
               ELSEIF (j==(ny*np+1) .AND. i==(nx*np+1))THEN
                  blk%mapping(i,j,1)=elem_num-nx*np
                  blk%mapping(i,j,2)=3
               ELSE
                  elem_num=elem_num+1
                  blk%mapping(i,j,1)=elem_num
                  blk%mapping(i,j,2)=1
               ENDIF
            ENDDO
         ENDDO
      CASE(2)
         DO i = 1,nx*np+1
            DO j = 1,ny*np+1
               IF (j==(ny*np+1) .AND. i/=(nx*np+1))THEN
                  blk%mapping(i,j,1)=elem_num
                  blk%mapping(i,j,2)=2
               ELSEIF (i==(nx*np+1) .AND. j/=(ny*np+1))THEN
                  elem_num=elem_num+1
                  blk%mapping(i,j,1)=elem_num-ny*np
                  blk%mapping(i,j,2)=4
               ELSEIF (i==(nx*np+1) .AND. j==(ny*np+1))THEN
                  blk%mapping(i,j,1)=elem_num-ny*np
                  blk%mapping(i,j,2)=1
               ELSE
                  elem_num=elem_num+1
                  blk%mapping(i,j,1)=elem_num
                  blk%mapping(i,j,2)=3
               ENDIF
            ENDDO
         ENDDO
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cubit_mapping
c-----------------------------------------------------------------------
c     subprogram 7. cubit_grid.
c     assign coordinates to SEL elements based on cubit grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cubit_grid(ixx,iyy,ksi,eta,blk)

      INTEGER, INTENT(IN) :: ixx,iyy
      TYPE(block) :: blk
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta

      INTEGER  :: block_num=1,node_num,elem_num,i,j,np
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      np=blk%np
c-----------------------------------------------------------------------
c     give physical coordinates to ksi and eta.
c-----------------------------------------------------------------------
      DO j=0,np
         DO i=0,np
            elem_num = blk%mapping(ixx*np+i+1,iyy*np+j+1,1)
            node_num = blk%mapping(ixx*np+i+1,iyy*np+j+1,2)

            ksi(i,j) = blk%coords(blk%connect(block_num,node_num,
     $           elem_num),1)
            eta(i,j) = blk%coords(blk%connect(block_num,node_num,
     $           elem_num),2)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cubit_grid
c-----------------------------------------------------------------------
c     subprogram 8. cubit_alloc
c     cubit data allocation
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cubit_alloc(blk)

      TYPE(block) :: blk
c-----------------------------------------------------------------------
c     allocate cubit variables
c-----------------------------------------------------------------------
      ALLOCATE(blk%connect(blk%blocks,blk%cxns,blk%nels))
      ALLOCATE(blk%coords(blk%nnodes,ndim))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cubit_alloc
c-----------------------------------------------------------------------
c     subprogram 9. cubit_dealloc
c     cubit data deallocation
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cubit_dealloc(blk)

      TYPE(block) :: blk
c-----------------------------------------------------------------------
c     deallocations
c-----------------------------------------------------------------------
      DEALLOCATE( blk%coords, blk%connect, blk%mapping )
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cubit_dealloc
c-----------------------------------------------------------------------
c     end cubit module
c-----------------------------------------------------------------------
      END MODULE cubit_mod
