c-----------------------------------------------------------------------
c     file p2_diagnose.F.
c     data I/O for the 2D code
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_diagnose_mod.
c     1. p2_diagnose_output.
c     2. p2_diagnose_output_p.
c     3. p2_diagnose_input.
c     4. p2_diagnose_input_p.
c     5. p2_diagnose_mpi_write.
c     6. p2_diagnose_mpi_read.
c     7. p2_diagnose_write.
c     8. p2_diagnose_read.
c     9. p2_diagnose_du_table.
c-----------------------------------------------------------------------
c     subprogram 0. p2_diagnose_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_diagnose_mod
      USE p2_sel_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_diagnose_output.
c     outputs calculated solution or grid after each time step to file 
c     with a local fortran input vector (uu).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_output(psv,outfile_type,data_flag,dmout,
     $     t,uu)

      TYPE(p2_sel_type) :: psv
      CHARACTER(*), INTENT(IN) :: outfile_type,data_flag
      INTEGER, INTENT(IN) :: dmout
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,psv%jxmin:,psv%jymin:), INTENT(IN) :: uu
      
      CHARACTER(160) :: filename
      INTEGER :: stepnum,istep
c-----------------------------------------------------------------------
c     output subsequent (adapted) grid to file.
c-----------------------------------------------------------------------
      stepnum = solver%istep
      SELECT CASE(data_flag)
      CASE("grid")
         istep = stepnum - solver%restart_step
         IF(istep > 0)
     $        stepnum=((istep-1)/dmout+1)*dmout + solver%restart_step
         SELECT CASE(outfile_type)
         CASE("hdf5")
            WRITE(filename,'(a,i5.5,a)')TRIM(outdir)//"/grid_",
     $           stepnum,".h5"
         CASE DEFAULT
            WRITE(filename,'(a,i5.5,a)')TRIM(outdir)//"/grid_",
     $           stepnum,".dat"
         END SELECT
      CASE("data")
         SELECT CASE(outfile_type)
         CASE("hdf5")
            WRITE(filename,'(a,i5.5,a)')TRIM(outdir)//"/sel_",
     $           stepnum,".h5"
         CASE DEFAULT
            WRITE(filename,'(a,i5.5,a)')TRIM(outdir)//"/sel_",
     $           stepnum,".dat"
         END SELECT
      END SELECT
      IF(parallel_write)THEN
         CALL p2_diagnose_mpi_write(psv,TRIM(filename),t,uu)
      ELSE
         CALL p2_diagnose_write(psv,TRIM(filename),outfile_type,t,uu)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_output
c-----------------------------------------------------------------------
c     subprogram 2. p2_diagnose_output_p.
c     outputs calculated solution or grid after each time step to file 
c     with a petsc input vector (petsc_uu).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_output_p(psv,outfile_type,data_flag,dmout,
     $     t,petsc_uu)

      TYPE(p2_sel_type) :: psv
      CHARACTER(*), INTENT(IN) :: outfile_type,data_flag
      INTEGER, INTENT(IN) :: dmout
      REAL(r8), INTENT(IN) :: t
      Vec :: petsc_uu

      REAL(r8), DIMENSION(psv%nqty,psv%jxmin:psv%jxmax,
     $     psv%jymin:psv%jymax) :: uu
c-----------------------------------------------------------------------
c     convert petsc_uu to local fortran uu and then write to file
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(psv%jscat,psv%jvec,petsc_uu,uu)
      CALL p2_diagnose_output(psv,outfile_type,data_flag,dmout,t,uu)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_output_p
c-----------------------------------------------------------------------
c     subprogram 3. p2_diagnose_input.
c     inputs data from restart files (grid or data files) 
c     with a local fortran input vector (uu).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_input(psv,restart_dir,outfile_type,
     $     data_flag,grid_step,t,uu)

      TYPE(p2_sel_type) :: psv
      CHARACTER(*), INTENT(IN) :: restart_dir,outfile_type,data_flag
      INTEGER, INTENT(IN) :: grid_step
      REAL(r8), INTENT(INOUT) :: t
      REAL(r8), DIMENSION(:,psv%jxmin:,psv%jymin:), INTENT(OUT) :: uu
   
      CHARACTER(160) :: filename
c-----------------------------------------------------------------------
c     input data to file.
c-----------------------------------------------------------------------
      SELECT CASE(data_flag)
      CASE("grid")
         SELECT CASE(outfile_type)
         CASE("hdf5")
            WRITE(filename,'(a,i5.5,a)')TRIM(restart_dir)//"/grid_",
     $           grid_step,".h5"
         CASE DEFAULT
            WRITE(filename,'(a,i5.5,a)')TRIM(restart_dir)//"/grid_",
     $           grid_step,".dat"
         END SELECT
      CASE("data")
         SELECT CASE(outfile_type)
         CASE("hdf5")
            WRITE(filename,'(a,i5.5,a)')TRIM(restart_dir)//"/sel_",
     $           solver%istep,".h5"
         CASE DEFAULT
            WRITE(filename,'(a,i5.5,a)')TRIM(restart_dir)//"/sel_",
     $           solver%istep,".dat"
         END SELECT
      END SELECT
      IF(parallel_read)THEN
         CALL p2_diagnose_mpi_read(psv,TRIM(filename),t,uu)
      ELSE
         CALL p2_diagnose_read(psv,TRIM(filename),outfile_type,t,uu)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_input
c-----------------------------------------------------------------------
c     subprogram 4. p2_diagnose_input_p.
c     inputs grid or data to file with a petsc vector (petsc_uu).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_input_p(psv,restart_dir,outfile_type,
     $     data_flag,grid_step,t,petsc_uu)

      TYPE(p2_sel_type) :: psv
      CHARACTER(*), INTENT(IN) :: restart_dir,outfile_type,data_flag
      INTEGER, INTENT(IN) :: grid_step
      REAL(r8), INTENT(INOUT) :: t
      Vec :: petsc_uu

      REAL(r8), DIMENSION(psv%nqty,psv%jxmin:psv%jxmax,
     $     psv%jymin:psv%jymax) :: uu
c-----------------------------------------------------------------------
c     convert petsc_uu to local fortran uu and then write to file
c-----------------------------------------------------------------------
      CALL p2_diagnose_input(psv,restart_dir,outfile_type,data_flag,
     $     grid_step,t,uu)
      CALL p2_sel_set(.TRUE.,uu,psv%jscat,psv%jvec,petsc_uu)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_input_p
c-----------------------------------------------------------------------
c     subprogram 5. p2_diagnose_mpi_write.
c     uses parallel HDF5 output to write data vector to file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_mpi_write(psv,filename,t,uu)

      TYPE(p2_sel_type) :: psv
      CHARACTER(*), INTENT(IN) :: filename
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,psv%jxmin:,psv%jymin:), INTENT(IN) :: uu

      INTEGER, PARAMETER :: u_rank=3,t_rank=1
      INTEGER :: info,h5_error,nqty
      INTEGER(HID_T) :: file_id,filespace_id,memspace_id,dset_id,
     $     attr_id,aspace_id,plist_id
      INTEGER(HSIZE_T), DIMENSION(1) :: t_dims
      INTEGER(HSIZE_T), DIMENSION(3) :: u_dims,ulocal_dims,h5offset
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      nqty=SIZE(uu,1)
      file_id=0
      filespace_id=0
      memspace_id=0
      dset_id=0
      attr_id=0
      aspace_id=0
      plist_id=0
      h5_error=0
      info=MPI_INFO_NULL
c-----------------------------------------------------------------------
c     write vector to file using parallel hdf5.
c-----------------------------------------------------------------------
      u_dims = (/nqty,psv%nxp+1,psv%nyp+1/)
      ulocal_dims = (/nqty,SIZE(uu,2),SIZE(uu,3)/)
      h5offset = (/0,psv%jxmin,psv%jymin/)

      CALL H5OPEN_F(h5_error)
c-----------------------------------------------------------------------
c     setup file access property list with parallel I/O access.
c-----------------------------------------------------------------------
      CALL H5PCREATE_F(H5P_FILE_ACCESS_F,plist_id,h5_error)
      CALL H5PSET_FAPL_MPIO_F(plist_id,comm,info,h5_error)
c-----------------------------------------------------------------------
c     create the file collectively
c-----------------------------------------------------------------------
      CALL H5FCREATE_F(filename,H5F_ACC_TRUNC_F,file_id,h5_error,
     $     access_prp = plist_id)
      CALL H5PCLOSE_F(plist_id,h5_error)
c-----------------------------------------------------------------------
c     create the dataspaces for the datasets
c-----------------------------------------------------------------------
      CALL H5SCREATE_SIMPLE_F(u_rank,u_dims,filespace_id,h5_error)
      CALL H5SCREATE_SIMPLE_F(u_rank,ulocal_dims,memspace_id,h5_error)
c-----------------------------------------------------------------------
c     create local datasets (for each processor)
c-----------------------------------------------------------------------
      CALL H5PCREATE_F(H5P_DATASET_CREATE_F,plist_id,h5_error)
      CALL H5DCREATE_F(file_id,"/AMP",H5T_NATIVE_DOUBLE,filespace_id,
     $     dset_id,h5_error,plist_id)
      CALL H5PCLOSE_F(plist_id,h5_error)
c-----------------------------------------------------------------------
c     select hyperslab in the file
c-----------------------------------------------------------------------
      CALL H5SSELECT_HYPERSLAB_F(filespace_id,H5S_SELECT_SET_F,
     $     h5offset,ulocal_dims,h5_error)
c-----------------------------------------------------------------------
c     create property list for collective dataset write
c-----------------------------------------------------------------------
      CALL H5PCREATE_F(H5P_DATASET_XFER_F,plist_id,h5_error)
      CALL H5PSET_DXPL_MPIO_F(plist_id,H5FD_MPIO_COLLECTIVE_F,h5_error)
c-----------------------------------------------------------------------
c     write the dataset collectively
c-----------------------------------------------------------------------
      CALL H5DWRITE_F(dset_id,H5T_NATIVE_DOUBLE,uu,u_dims,h5_error,
     $     file_space_id = filespace_id,mem_space_id = memspace_id,
     $     xfer_prp = plist_id)
c-----------------------------------------------------------------------
c     write "time", t as an attribute to the dataset.
c-----------------------------------------------------------------------
      t_dims = 1
      CALL H5SCREATE_SIMPLE_F(t_rank,t_dims,aspace_id,h5_error)
      CALL H5ACREATE_F(dset_id,"time",H5T_NATIVE_DOUBLE,aspace_id,
     $     attr_id,h5_error)
      CALL H5AWRITE_F(attr_id,H5T_NATIVE_DOUBLE,t,t_dims,h5_error)
      CALL H5ACLOSE_F(attr_id,h5_error)
      CALL H5SCLOSE_F(aspace_id,h5_error)
c-----------------------------------------------------------------------
c     close out hdf5 dataspaces,dataset,parallel property list, file, 
c     and fortran interface.
c-----------------------------------------------------------------------
      CALL H5DCLOSE_F(dset_id,h5_error)
      CALL H5SCLOSE_F(filespace_id,h5_error)
      CALL H5SCLOSE_F(memspace_id,h5_error)
      CALL H5PCLOSE_F(plist_id,h5_error)
      CALL H5FCLOSE_F(file_id,h5_error)
      CALL H5CLOSE_F(h5_error)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_mpi_write
c-----------------------------------------------------------------------
c     subprogram 6. p2_diagnose_mpi_read.
c     uses parallel HDF5 input to read data vector from file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_mpi_read(psv,filename,t,uu)

      TYPE(p2_sel_type) :: psv
      CHARACTER(*), INTENT(IN) :: filename
      REAL(r8), INTENT(OUT) :: t
      REAL(r8), DIMENSION(:,psv%jxmin:,psv%jymin:), INTENT(OUT) :: uu

      INTEGER, PARAMETER :: u_rank=3
      INTEGER :: info,h5_error,nqty
      INTEGER(HID_T) :: file_id,filespace_id,memspace_id,dset_id,
     $     dtype_id,attr_id,atype_id,plist_id
      INTEGER(HSIZE_T), DIMENSION(1) :: t_dims
      INTEGER(HSIZE_T), DIMENSION(3) :: u_dims,h5offset
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      nqty=SIZE(uu,1)
      file_id=0
      filespace_id=0
      memspace_id=0
      dset_id=0
      dtype_id=0
      attr_id=0
      atype_id=0
      plist_id=0
      h5_error=0
      info=MPI_INFO_NULL
      uu=0
c-----------------------------------------------------------------------
c     read vector from file using parallel hdf5.
c-----------------------------------------------------------------------
      u_dims = (/nqty,SIZE(uu,2),SIZE(uu,3)/)
      h5offset = (/0,psv%jxmin,psv%jymin/)

      CALL H5OPEN_F(h5_error)
c-----------------------------------------------------------------------
c     setup file access property list with parallel I/O access.
c-----------------------------------------------------------------------
      CALL H5PCREATE_F(H5P_FILE_ACCESS_F,plist_id,h5_error)
      CALL H5PSET_FAPL_MPIO_F(plist_id,comm,info,h5_error)    
c-----------------------------------------------------------------------
c     open the file collectively
c-----------------------------------------------------------------------
      CALL H5FOPEN_F(filename,H5F_ACC_RDONLY_F,file_id,h5_error,
     $     access_prp = plist_id)
      CALL H5PCLOSE_F(plist_id,h5_error)
c-----------------------------------------------------------------------
c     create the memory dataspace for the local dataset
c-----------------------------------------------------------------------
      CALL H5SCREATE_SIMPLE_F(u_rank,u_dims,memspace_id,h5_error)
c-----------------------------------------------------------------------
c     get local dataspace (for each processor)
c-----------------------------------------------------------------------
      CALL H5DOPEN_F(file_id,"/AMP",dset_id,h5_error)
      CALL H5DGET_TYPE_F(dset_id,dtype_id,h5_error)
      CALL H5DGET_SPACE_F(dset_id,filespace_id,h5_error)
c-----------------------------------------------------------------------
c     select hyperslab in the file
c-----------------------------------------------------------------------
      CALL H5SSELECT_HYPERSLAB_F(filespace_id,H5S_SELECT_SET_F,
     $     h5offset,u_dims,h5_error)
c-----------------------------------------------------------------------
c     create property list for collective dataset read
c-----------------------------------------------------------------------
      CALL H5PCREATE_F(H5P_DATASET_XFER_F,plist_id,h5_error)
      CALL H5PSET_DXPL_MPIO_F(plist_id,H5FD_MPIO_COLLECTIVE_F,h5_error)
c-----------------------------------------------------------------------
c     read the dataset collectively
c-----------------------------------------------------------------------      
      CALL H5DREAD_F(dset_id,dtype_id,uu,u_dims,h5_error,
     $     file_space_id = filespace_id,mem_space_id = memspace_id,
     $     xfer_prp = plist_id)
      CALL H5TCLOSE_F(dtype_id,h5_error)
c-----------------------------------------------------------------------
c     read "time", t as an attribute to the dataset.
c-----------------------------------------------------------------------
      t_dims = 1
      CALL H5AOPEN_NAME_F(dset_id,"time",attr_id,h5_error)
      CALL H5AGET_TYPE_F(attr_id,atype_id,h5_error)
      CALL H5AREAD_F(attr_id,atype_id,t,t_dims,h5_error)
      CALL H5TCLOSE_F(atype_id,h5_error)
      CALL H5ACLOSE_F(attr_id,h5_error)
c-----------------------------------------------------------------------
c     close out hdf5  dataset, and file.
c-----------------------------------------------------------------------
      CALL H5DCLOSE_F(dset_id,h5_error)
      CALL H5SCLOSE_F(filespace_id,h5_error)
      CALL H5SCLOSE_F(memspace_id,h5_error)
      CALL H5PCLOSE_F(plist_id,h5_error)
      CALL H5FCLOSE_F(file_id,h5_error)
      CALL H5CLOSE_F(h5_error)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_mpi_read
c-----------------------------------------------------------------------
c     subprogram 7. p2_diagnose_write.
c     writes the whole vector with one processor.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_write(psv,filename,outfile_type,t,uu)

      TYPE(p2_sel_type) :: psv
      CHARACTER(*), INTENT(IN) :: filename,outfile_type
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,psv%jxmin:,psv%jymin:), INTENT(IN) :: uu

      INTEGER, PARAMETER :: u_rank=3,t_rank=1
      INTEGER :: h5_error,nqty,np,jx,jy
      INTEGER(HID_T) :: file_id,filespace_id,dset_id,attr_id,aspace_id
      INTEGER(HSIZE_T), DIMENSION(1) :: t_dims
      INTEGER(HSIZE_T), DIMENSION(3) :: u_dims
      REAL(r8), DIMENSION(SIZE(uu,1),0:psv%np,0:psv%np) :: u0
      REAL(r8), DIMENSION(SIZE(uu,1),0:psv%nxp,0:psv%nyp) :: uglobal
c-----------------------------------------------------------------------
c     collect global vector.
c-----------------------------------------------------------------------
      nqty=SIZE(uu,1)
      np=psv%np
      uglobal=0
      DO jy=0,psv%nyp,np+1
         DO jx=0,psv%nxp,np+1
            u0 = 0
            IF(jx >= psv%jxmin .AND. jx < psv%jxmax 
     $           .AND. jy >= psv%jymin .AND. jy < psv%jymax)
     $           u0 = uu(:,jx:jx+np,jy:jy+np)
            CALL MPI_Allreduce(u0,uglobal(:,jx:jx+np,jy:jy+np),
     $           nqty*(np+1)**2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write vector to file.
c-----------------------------------------------------------------------
      IF(mpi_rank==0)THEN
         SELECT CASE(outfile_type)
c-----------------------------------------------------------------------
c     write vector to file using hdf5.
c-----------------------------------------------------------------------
         CASE("hdf5")
            u_dims = (/SIZE(uu,1),psv%nxp+1,psv%nyp+1/)
            CALL H5OPEN_F(h5_error)
            CALL H5FCREATE_F(filename,H5F_ACC_TRUNC_F,file_id,h5_error)
            CALL H5SCREATE_SIMPLE_F(u_rank,u_dims,filespace_id,h5_error)
            CALL H5DCREATE_F(file_id,"/AMP",H5T_NATIVE_DOUBLE,
     $           filespace_id,dset_id,h5_error)
            CALL H5DWRITE_F(dset_id,H5T_NATIVE_DOUBLE,uglobal,u_dims,
     $           h5_error)
c-----------------------------------------------------------------------
c     write "time", t as an attribute to the dataset.
c-----------------------------------------------------------------------
            t_dims = 1
            CALL H5SCREATE_SIMPLE_F(t_rank,t_dims,aspace_id,h5_error)
            CALL H5ACREATE_F(dset_id,"time",H5T_NATIVE_DOUBLE,aspace_id,
     $           attr_id,h5_error)
            CALL H5AWRITE_F(attr_id,H5T_NATIVE_DOUBLE,t,t_dims,h5_error)
            CALL H5ACLOSE_F(attr_id,h5_error)
            CALL H5SCLOSE_F(aspace_id,h5_error)
c-----------------------------------------------------------------------
c     close out hdf5 dataspaces, dataset, and file.
c-----------------------------------------------------------------------
            CALL H5DCLOSE_F(dset_id,h5_error)
            CALL H5SCLOSE_F(filespace_id,h5_error)
            CALL H5FCLOSE_F(file_id,h5_error)
            CALL H5CLOSE_F(h5_error)
c-----------------------------------------------------------------------
c     write vector to file using default binary format.
c-----------------------------------------------------------------------
         CASE DEFAULT
            OPEN(UNIT=dat_unit,FILE=TRIM(filename),ACTION="WRITE",
     $           STATUS="REPLACE",FORM="UNFORMATTED")
            WRITE(dat_unit)REAL(uglobal,8)
            WRITE(dat_unit)REAL(t,8)
            CLOSE(UNIT=dat_unit)
         END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_write
c-----------------------------------------------------------------------
c     subprogram 8. p2_diagnose_read.
c     reads the whole vector with one processor.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_read(psv,filename,outfile_type,t,uu_local)

      TYPE(p2_sel_type) :: psv
      CHARACTER(*), INTENT(IN) :: filename,outfile_type
      REAL(r8), INTENT(INOUT) :: t
      REAL(r8), DIMENSION(:,psv%jxmin:,psv%jymin:), INTENT(OUT) :: 
     $     uu_local

      INTEGER :: tflag,h5_error,nqty,np,jx,jy
      INTEGER(HID_T) :: file_id,dset_id,dtype_id,attr_id,atype_id
      INTEGER(HSIZE_T), DIMENSION(1) :: t_dims
      INTEGER(HSIZE_T), DIMENSION(3) :: u_dims
      REAL(r8) :: ttemp
      REAL(r8), DIMENSION(SIZE(uu_local,1),0:psv%nxp,0:psv%nyp) :: uu
c-----------------------------------------------------------------------
c     read vector from file.
c-----------------------------------------------------------------------
      nqty=SIZE(uu_local,1)
      np=psv%np
      uu=0
      IF(mpi_rank==0)THEN
         SELECT CASE(outfile_type)
c-----------------------------------------------------------------------
c     read vector to file using hdf5.
c-----------------------------------------------------------------------
         CASE("hdf5")
            CALL H5OPEN_F(h5_error)
            CALL H5FOPEN_F(filename,H5F_ACC_RDWR_F,file_id,h5_error)
            CALL H5DOPEN_F(file_id,"/AMP",dset_id,h5_error)
            u_dims = (/nqty,SIZE(uu,2),SIZE(uu,3)/)
            CALL H5DGET_TYPE_F(dset_id,dtype_id,h5_error)
            CALL H5DREAD_F(dset_id,dtype_id,uu,u_dims,h5_error)
            CALL H5TCLOSE_F(dtype_id,h5_error)
c-----------------------------------------------------------------------
c     read "time", t as an attribute to the dataset.
c-----------------------------------------------------------------------
            t_dims = 1
            CALL H5AOPEN_NAME_F(dset_id,"time",attr_id,h5_error)
            CALL H5AGET_TYPE_F(attr_id,atype_id,h5_error)
            CALL H5AREAD_F(attr_id,atype_id,t,t_dims,h5_error)
            CALL H5TCLOSE_F(atype_id,h5_error)
            CALL H5ACLOSE_F(attr_id,h5_error)
c-----------------------------------------------------------------------
c     close out hdf5  dataset, and file.
c-----------------------------------------------------------------------
            CALL H5DCLOSE_F(dset_id,h5_error)
            CALL H5FCLOSE_F(file_id,h5_error)
            CALL H5CLOSE_F(h5_error)
c-----------------------------------------------------------------------
c     read vector to file using default binary format.
c-----------------------------------------------------------------------
         CASE DEFAULT
            OPEN(UNIT=dat_unit,FILE=TRIM(filename),ACTION="READ",
     $           STATUS="OLD",FORM="UNFORMATTED")   
            READ(dat_unit)uu
            READ(dat_unit,IOSTAT=tflag)ttemp
            CLOSE(UNIT=dat_unit)
            IF(tflag==0)t=ttemp
         END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     broadcast solution vector to all processors.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(t,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      DO jy=0,psv%nyp,np+1
         DO jx=0,psv%nxp,np+1
            CALL MPI_Bcast(uu(:,jx:jx+np,jy:jy+np),nqty*(np+1)**2,
     $           MPI_DOUBLE_PRECISION,0,comm,ierr)
         ENDDO
      ENDDO
      uu_local = uu(:,psv%jxmin:psv%jxmax,psv%jymin:psv%jymax)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_read
c-----------------------------------------------------------------------
c     subprogram 9. p2_diagnose_du_table.
c     produce tabular display of du and related quantities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_diagnose_du_table(res,du)

      Vec, Intent(IN) :: res,du

      CHARACTER(80) :: filename,format2,format3
      CHARACTER(160) :: format1
      INTEGER :: ixp,iyp,iqty,nqty,np,ixpmin,iypmin,ixpmax,iypmax
      REAL(r8) :: errmax,errmin,rmax
      REAL(r8), DIMENSION(pt2%nqty) :: r_max,u_max,ju_max,err_max,
     $     r_min,u_min,ju_min,err_min
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: r,u,ju,err
      Vec :: jdu
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/3x,"iyp",3x,"ixp",',i1,
     $     '(7x,"du",i1,8x,"jdu",i1,7x,"res",i1,7x,"err",i1,2x)/)')
 20   FORMAT('(2i6,1p,',i2,'(2x,4e11.3))')
 30   FORMAT('(a6,6x,1p,',i2,'(2x,4e11.3))')
c-----------------------------------------------------------------------
c     compute local sizes.
c-----------------------------------------------------------------------
      nqty=pt2%nqty
      np=pt2%np
      ixpmin=pt2%jxmin
      iypmin=pt2%jymin
      ixpmax=pt2%jxmax
      iypmax=pt2%jymax
c-----------------------------------------------------------------------
c     create jduy and compute matrix-vector product.
c-----------------------------------------------------------------------
      CALL VecDuplicate(du,jdu,ierr)
      CALL MatMult(pt2%jac,du,jdu,ierr)
c-----------------------------------------------------------------------
c     get local arrays.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,res,r)
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,du,u)
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,jdu,ju)
c-----------------------------------------------------------------------
c     compute errors.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         rmax=MAXVAL(ABS(r(iqty,:,:)))
         IF(rmax == 0)THEN
            err(iqty,:,:)=0
         ELSE
            err(iqty,:,:)=ABS(ju(iqty,:,:)-r(iqty,:,:))/rmax
         ENDIF
      ENDDO
      errmax=MAXVAL(ABS(err))
      errmin=MINVAL(ABS(err))
c-----------------------------------------------------------------------
c     compute extrema.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         r_max(iqty)=MAXVAL(ABS(r(iqty,:,:)))
         r_min(iqty)=MINVAL(ABS(r(iqty,:,:)))
         u_max(iqty)=MAXVAL(ABS(u(iqty,:,:)))
         u_min(iqty)=MINVAL(ABS(u(iqty,:,:)))
         ju_max(iqty)=MAXVAL(ABS(ju(iqty,:,:)))
         ju_min(iqty)=MINVAL(ABS(ju(iqty,:,:)))
         err_max(iqty)=MAXVAL(ABS(err(iqty,:,:)))
         err_min(iqty)=MINVAL(ABS(err(iqty,:,:)))
      ENDDO
c-----------------------------------------------------------------------
c     create format statements.
c-----------------------------------------------------------------------
      WRITE(format1,10)nqty
      WRITE(format2,20)nqty
      WRITE(format3,30)nqty
c-----------------------------------------------------------------------
c     open output file.
c-----------------------------------------------------------------------
      WRITE(filename,'(a,2i1,a)')"du",mpi_rank,mpi_size,".out"
      OPEN(UNIT=debug_unit,FILE=TRIM(filename),STATUS="UNKNOWN")
c-----------------------------------------------------------------------
c     write main output.
c-----------------------------------------------------------------------
      DO iyp=iypmin,iypmax
         WRITE(debug_unit,format1)(iqty,iqty,iqty,iqty,iqty=1,nqty)
         DO ixp=ixpmin,ixpmax
            WRITE(debug_unit,format2)iyp,ixp,
     $           (u(iqty,ixp,iyp),ju(iqty,ixp,iyp),r(iqty,ixp,iyp),
     $           err(iqty,ixp,iyp),iqty=1,nqty)
         ENDDO
      ENDDO
      WRITE(debug_unit,format1)(iqty,iqty,iqty,iqty,iqty=1,nqty)
c-----------------------------------------------------------------------
c     write extrema.
c-----------------------------------------------------------------------
      WRITE(debug_unit,format3)"min",(u_min(iqty),ju_min(iqty),
     $     r_min(iqty),err_min(iqty),iqty=1,nqty)
      WRITE(debug_unit,format3)"max",(u_max(iqty),ju_max(iqty),
     $     r_max(iqty),err_max(iqty),iqty=1,nqty)
c-----------------------------------------------------------------------
c     write final output.
c-----------------------------------------------------------------------
      WRITE(debug_unit,'(/a,a4,1p,2(a,e10.3))')
     $     " solve_type = ",TRIM(pt2%solve_type),
     $     ", errmax = ",errmax,", errmin = ",errmin
c-----------------------------------------------------------------------
c     clean up.
c-----------------------------------------------------------------------
      CALL VecDestroy(jdu,ierr)
      IF(pt2%solve_type /= "full")CALL MatDestroy(pt2%jac,ierr)
c-----------------------------------------------------------------------
c     close output file and abort.
c-----------------------------------------------------------------------
      CLOSE(UNIT=debug_unit)
      CALL program_stop("p2_diagnose_du_table: "
     $     //"Termination after diagnosing du.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_diagnose_du_table
      END MODULE p2_diagnose_mod
