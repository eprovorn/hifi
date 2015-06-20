c-----------------------------------------------------------------------
c     file post.f.
c     reads and post-processes output from SEL code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. post_mod.
c     1. post_run.
c     2. post_ls.
c     4. post_write0.
c     5. post_write_amp.
c     6. post_griderr.
c     7. post_mass0.
c     9. post_read_header.
c     10. post_main.
c-----------------------------------------------------------------------
c     subprogram 0. post_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE post_mod
      USE post4field_mod
      USE post2fluid_mod
      USE postpn_mod
      USE postbreakout_mod
      USE postmast_mod
      USE helix_mod
      USE beltrami_mod
      USE posteuler_mod
      IMPLICIT NONE

      LOGICAL :: drawgrid,adapt_grid,contour=.TRUE.,cray_flag=.FALSE.
      INTEGER :: stride=1
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: mass0
      TYPE(jacobi_type) :: basis,quad

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. post_run.
c     main body of post-processing.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post_run

      LOGICAL :: xt_flag,yt_flag,init_flag,final_flag,first1d,
     $     first2d,last2d,polar_crd=.FALSE.,nodal=.FALSE.,
     $     header_dat=.FALSE.
      CHARACTER(160) :: filename
      CHARACTER(160), DIMENSION(:), POINTER :: file
      INTEGER :: nqty,nfile,mfile,firstfile,lastfile,nxw0,nyw0,nxw,nyw,
     $     nxp,nyp,ifile,gout,jfile,myios,nt_next,nxs,nys,ii
      INTEGER, DIMENSION(3) :: xcont
      REAL(r8) :: t_old,xs,ys,min_value=0.
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x,y,griderr
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: jac
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: u,uw,uw0,duw,xyw
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: uxyw,xyw_kt

      CHARACTER(160) :: h5_ufile
      INTEGER(HID_T) :: h5_uid,u_dset_id,attr_id,atype_id
      INTEGER(HSIZE_T), DIMENSION(1) :: t_dims
      INTEGER(HSIZE_T), DIMENSION(3) :: u_dims
      INTEGER :: h5_error

      NAMELIST/post_input/job_type,indir,postout,out_type,nxw0,nyw0,
     $     mfile,stride,contour,drawgrid,polar_crd,flag1d,flag2d,
     $     cray_flag,min_value
      NAMELIST/slice_input/xt_flag,yt_flag,init_flag,final_flag,
     $     nxs,nys,xs,ys
c-----------------------------------------------------------------------
c     declarations of unused input variables.
c-----------------------------------------------------------------------
      LOGICAL :: diagnose,restart_flag,maptest,monitor,fd_test,
     $     adapt_dt,read_grid,du_diagnose,parallel_write,parallel_read,
     $     always_pc_reset,dual_precon,restart_time,fd_bound_test
      CHARACTER(8) :: grid_inv_type,grid_type
      CHARACTER(160) :: outdir,restart_dir,quad_type,step_type,
     $     solve_type,cubit_file
      INTEGER :: restart_step,grid_step,ksp_restart,ksps_max,itmax,
     $     itmax_decr,itmax_incr,tflag
      REAL(r8) :: errtol,gridtol1,gridtol2,dtmin,dt_incr,dt_decr,theta,
     $     restart_t0
      
      NAMELIST/algorithm_input/maptest,monitor,nodal,quad_type,
     $     grid_inv_type,solve_type,fd_test,polar,adapt_grid,
     $     gridtol1,gridtol2,adapt_dt,itmax,itmax_incr,itmax_decr,
     $     dt_incr,dt_decr,step_type,theta,errtol,ksp_restart,ksps_max,
     $     always_pc_reset,parallel_write,parallel_read,du_diagnose,
     $     dual_precon,grid_type,outfile_type,fd_bound_test
      NAMELIST/universal_input/nx,ny,np,nq,nbx,xperiodic,yperiodic,
     $     dt,dtmax,dtmin,tmax,nstep,dmout,outdir,restart_flag,
     $     restart_dir,restart_step,read_grid,grid_step,restart_time,
     $     restart_t0,cubit_file
c-----------------------------------------------------------------------
c     read control files.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="post.in",STATUS="OLD")
      READ(in_unit,NML=post_input)
      READ(in_unit,NML=slice_input)
      CLOSE(UNIT=in_unit)

      filename=TRIM(indir)//"/sel.in"
      OPEN(UNIT=in_unit,FILE=filename,STATUS="OLD")
      READ(in_unit,NML=algorithm_input)
      READ(in_unit,NML=universal_input)
      CLOSE(UNIT=in_unit)

      WRITE(out_unit,'(2a)')"indir = ",TRIM(indir)
      IF(drawgrid .AND. adapt_grid)THEN
         filename=TRIM(indir)//"/grid"
         job_type="beltrami"
      ELSEIF(drawgrid)THEN
         filename=TRIM(indir)//"/header"
         job_type="beltrami"
      ELSE
         filename=TRIM(indir)//"/header"
      ENDIF
      WRITE(out_unit,'(2a)')"job_type = ",TRIM(job_type)
c-----------------------------------------------------------------------
c     check whether or not there is a grid to be read.
c-----------------------------------------------------------------------
      CALL system("rm -f ls.out")
      SELECT CASE(outfile_type)
      CASE("hdf5")
         CALL system("ls "//TRIM(indir)//"/grid_*.h5 | wc -l > ls.out")
      CASE DEFAULT
         CALL system("ls "//TRIM(indir)//"/grid_*.dat | wc -l > ls.out")
      END SELECT
      OPEN(UNIT=in_unit,FILE="ls.out",STATUS="OLD")
      READ(in_unit,*)nfile
      IF(nfile > 0)adapt_grid=.TRUE.
      CLOSE(UNIT=in_unit)
      CALL system("rm -f ls.out")
c-----------------------------------------------------------------------
c     job-specific read and diagnostics initialization.
c-----------------------------------------------------------------------
      CALL system("rm -f dUdt.bin Ucontour.bin VecSqInt.bin 
     $     maxUvsT.bin UxyT.bin Uprofile.bin contour.bin dcontour.bin
     $     amp.bin xt.bin yt.bin xy.bin layer.bin")
      SELECT CASE(job_type)
      CASE("euler1D","euler2D")
         CALL posteuler_read(indir,out_type)
      CASE("fourfield")
         CALL post4field_read
      CASE("fivefield","fourfieldplus","HallMHDiso","cmhd","HallMHD",
     $        "HallMHDicold","rmhd","emhd","epmhd","cmhdn")
         CALL post2fluid_read(indir)
      CASE("pn_ext")
         CALL postpn_read(indir)
      CASE("SolarHallMHD","SolarLnHMHD")
         CALL postbreakout_read(indir)
      CASE("MAST")
         CALL postmast_read(indir)
      CASE("fourfieldhlx","XMHDhlx")
         CALL helix_read(indir)
      CASE("beltrami")
         CALL beltrami_read
      END SELECT
c-----------------------------------------------------------------------
c     open time-dependent slice_files.
c-----------------------------------------------------------------------
      IF(xt_flag)OPEN(UNIT=xt_unit,FILE="xt.bin",STATUS="UNKNOWN",
     $     FORM="UNFORMATTED")
      IF(yt_flag)OPEN(UNIT=yt_unit,FILE="yt.bin",STATUS="UNKNOWN",
     $     FORM="UNFORMATTED")
c-----------------------------------------------------------------------
c     initialize variables and get listing of data files.
c-----------------------------------------------------------------------
      lastfile=-1
      jfile=0
      t_old=0
      CALL post_ls(indir,mfile,nfile,file)
      myios=0
      dmout=MAX(dmout,1)
      gout=1
      xcont=0
c-----------------------------------------------------------------------
c     read header file.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(filename)//".dat",EXIST=header_dat)
      IF(header_dat)THEN
         OPEN(UNIT=header_unit,FILE=TRIM(filename)//".dat",
     $        STATUS="UNKNOWN",FORM="UNFORMATTED")
         READ(header_unit)nqty,nx,ny,np,xperiodic,yperiodic,dt,t_old
      ELSE
         OPEN(UNIT=header_unit,FILE=TRIM(filename)//".txt",
     $        STATUS="UNKNOWN")
         READ(header_unit,'(i5,i5,i5,i5)')nqty,nx,ny,np
      ENDIF
      IF(drawgrid .AND. job_type=="beltrami")nqty=2
      ALLOCATE(x(0:nx),y(0:ny),mass0(0:np,0:np))
      IF(header_dat)THEN
         READ(header_unit)x,y
      ELSE
         x=(/(ii,ii=0,nx)/)/REAL(nx,r8)
         y=(/(ii,ii=0,ny)/)/REAL(ny,r8)
      ENDIF
c-----------------------------------------------------------------------
c     compute derived sizes.
c-----------------------------------------------------------------------
      nxp=nx*(np+1)-1
      nyp=ny*(np+1)-1
      nxw=nxw0*nx
      nyw=nyw0*ny
c-----------------------------------------------------------------------
c     allocate space for grid and physical variables.
c-----------------------------------------------------------------------
      ALLOCATE(u(nqty,0:nxp,0:nyp),uw(nqty,0:nxw,0:nyw),
     $     uxyw(0:nxw,0:nyw,2,nqty),duw(nqty,0:nxw,0:nyw),
     $     griderr(2*nqty),xyw(2,0:nxw,0:nyw),xyw_kt(0:nxw,0:nyw,2,2),
     $     jac(0:nxw,0:nyw))
c-----------------------------------------------------------------------
c     prepare for data output
c-----------------------------------------------------------------------
      SELECT CASE(out_type)
      CASE("hdf5")
         IF(contour .OR. flag2d)
     $        CALL system("mkdir -p "//TRIM(postout))
      CASE DEFAULT
         IF(contour)THEN
            OPEN(UNIT=Ucontour_unit,FILE="contour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,nqty
            CLOSE(UNIT=Ucontour_unit)
            OPEN(UNIT=Ucontour_unit,FILE="dcontour.bin",
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)1,0,nqty
            CLOSE(UNIT=Ucontour_unit)
         ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     allocate basis and calculate mass matrix.
c-----------------------------------------------------------------------
      CALL jacobi_alloc(basis,np,nodal,.FALSE.,"gll")
      CALL jacobi_alloc(quad,np,.TRUE.,.FALSE.,"gl0")
      CALL post_mass0
c-----------------------------------------------------------------------
c     start main loop
c-----------------------------------------------------------------------
      DO
c-----------------------------------------------------------------------
c     read in time-step/grid data from header.dat or header.txt file
c     and determine the number of output files with given grid
c-----------------------------------------------------------------------
         CALL post_read_header(header_dat,nfile,gout,jfile,myios,
     $        firstfile,lastfile,nt_next)
c-----------------------------------------------------------------------
c     cycle through datafiles
c-----------------------------------------------------------------------
         DO ifile=firstfile,lastfile
c-----------------------------------------------------------------------
c     set logical flags for producing post-processed output
c-----------------------------------------------------------------------
            IF(MOD(ifile,stride) /= 0)THEN
               IF(ifile == nfile)THEN
                  first1d=.FALSE.
                  first2d=(ifile == firstfile)
                  last2d=.TRUE.                  
               ELSEIF(ifile == firstfile 
     $                 .AND. (lastfile/stride-firstfile/stride) > 0)THEN
                  first1d=.FALSE.
                  first2d=.TRUE.
                  last2d=.FALSE.
               ELSE                  
                  CYCLE
               ENDIF
            ELSE
               first1d=(ifile == 0)
               first2d=(ifile == firstfile)
               last2d=(ifile == nfile) .AND. myios /= 0
            ENDIF
c-----------------------------------------------------------------------
c     read and close file.
c-----------------------------------------------------------------------
            WRITE(out_unit,'(a,i4,a,i4)')"ifile = ",ifile,"/",lastfile
            SELECT CASE(outfile_type)
            CASE("hdf5")
               h5_ufile = TRIM(file(ifile))
               CALL H5OPEN_F(h5_error)
               CALL H5FOPEN_F(h5_ufile,H5F_ACC_RDWR_F,h5_uid,h5_error)
               CALL H5DOPEN_F(h5_uid,"/AMP",u_dset_id,h5_error)
               u_dims = (/SIZE(u,1),SIZE(u,2),SIZE(u,3)/)
               CALL H5DREAD_F(u_dset_id,H5T_NATIVE_DOUBLE,u,u_dims,
     $              h5_error)
               t_dims = 1
               CALL H5AOPEN_NAME_F(u_dset_id,"time",attr_id,h5_error)
               CALL H5AGET_TYPE_F(attr_id,atype_id,h5_error)
               CALL H5AREAD_F(attr_id,atype_id,t_old,t_dims,h5_error)
               CALL H5TCLOSE_F(atype_id,h5_error)
               CALL H5ACLOSE_F(attr_id,h5_error)
               CALL H5DCLOSE_F(u_dset_id,h5_error)
               CALL H5FCLOSE_F(h5_uid,h5_error)
               CALL H5CLOSE_F(h5_error)
            CASE DEFAULT
               OPEN(UNIT=file_unit,FILE=TRIM(file(ifile)),
     $              STATUS="UNKNOWN",FORM="UNFORMATTED")
               READ(file_unit)u
               READ(file_unit,IOSTAT=tflag)t_old
               CLOSE(UNIT=file_unit)
            END SELECT
c-----------------------------------------------------------------------
c     interpolate.
c-----------------------------------------------------------------------
            CALL plotter_interp(nxw0,nyw0,nx,ny,basis,u,uw,uxyw)
c-----------------------------------------------------------------------
c     compute and output unstructured interpolated position array.
c-----------------------------------------------------------------------
            IF(first2d)CALL post_write0(nxp,nyp,x,y,file(ifile),
     $           xyw,xyw_kt,nxw0,nyw0,nxw,nyw,polar_crd)
c-----------------------------------------------------------------------
c     transform derivatives of u with respect to x and y.
c-----------------------------------------------------------------------
            jac=1.
            IF(adapt_grid .AND. .NOT. drawgrid)THEN
               CALL plotter_transf(xyw_kt,jac,uxyw)
            ENDIF
c-----------------------------------------------------------------------
c     draw slices.
c-----------------------------------------------------------------------
            IF(xt_flag)CALL slice_xt(basis,x,y,u,ifile,ys,nxs)
            IF(yt_flag)CALL slice_yt(basis,x,y,u,ifile,xs,nys)
            IF(init_flag .AND. ifile == firstfile)THEN
               CALL slice_xy(basis,x,y,u,nxs,nys,"initxy")
               CALL slice_yx(basis,x,y,u,nxs,nys,"inityx")
            ENDIF
            IF(final_flag .AND. ifile == lastfile)THEN
               CALL slice_xy(basis,x,y,u,nxs,nys,"finalxy")
               CALL slice_yx(basis,x,y,u,nxs,nys,"finalyx")
            ENDIF
c-----------------------------------------------------------------------
c     update initial u-vector if grid has changed.
c-----------------------------------------------------------------------
            IF(first1d)ALLOCATE(uw0(nqty,0:nxw,0:nyw))
            IF(first2d)uw0=uw
            duw=uw-uw0
c-----------------------------------------------------------------------
c     write diagnostic files for specific job_type
c-----------------------------------------------------------------------
            SELECT CASE(job_type)
            CASE("euler1D")
               CALL posteuler_Uprofile(nxs,xyw,xyw_kt,uw)
            CASE("euler2D")
               IF(flag1d)THEN
                  CALL posteuler_Uprofile(nxs,xyw,xyw_kt,uw)
                  CALL posteuler_energy(t_old,xyw,jac,uw)
               ENDIF
               IF(flag2d)THEN
                  CALL posteuler_contour(t_old,xyw,uw,first2d,
     $                 last2d,nt_next,ifile,stride,file(ifile))
               ENDIF
            CASE("fourfield")
               CALL post4field_UxyT(t_old,xyw,xyw_kt,uw)
               CALL post4field_VecSqInt(t_old,xyw,jac,uw,uxyw)
            CASE("fivefield","fourfieldplus","HallMHDiso","cmhd","cmhdn"
     $              ,"HallMHD","HallMHDicold","rmhd","emhd","epmhd")
               CALL post2fluid_main(t_old,jac,xyw,xyw_kt,uw,uxyw,
     $              flag1d,flag2d,first1d,first2d,last2d,nt_next,ifile,
     $              stride,file(ifile),nxs)
            CASE("pn_ext")
               CALL postpn_main(t_old,jac,xyw,xyw_kt,uw,uxyw,
     $              flag1d,flag2d,first1d,first2d,last2d,nt_next,ifile,
     $              stride,file(ifile),nxs)
            CASE("SolarHallMHD")
               IF(flag1d)THEN
                  CALL postbreakout_energy(t_old,xyw,jac,uw,uxyw)
               ENDIF
               IF(flag2d)THEN
                  CALL postbreakout_contour(t_old,xyw,uw,uxyw,
     $                 first2d,last2d,nt_next,ifile,stride,file(ifile))
               ENDIF
            CASE("SolarLnHMHD")
               IF(flag1d)THEN
                  CALL postbreakout_dUdt(t_old,xyw,uw,first1d,first2d,
     $                 last2d)
                  CALL postbreakout_maxUvsT(t_old,uw)
                  CALL postbreakout_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
                  CALL postbreakout_cstracking(t_old,xyw,uw,first1d)
                  CALL postbreakout_rTtracking(t_old,xyw,xyw_kt,uw)
               ENDIF
               IF(flag2d)THEN
                  CALL postbreakout_contour(t_old,xyw,uw,uxyw,
     $                 first2d,last2d,nt_next,ifile,stride,file(ifile))
               ENDIF
            CASE("fourfieldhlx","XMHDhlx")
               IF(flag1d)THEN
                  CALL helix_dUdt(t_old,xyw,uw,uxyw,first1d)
                  CALL helix_UxyT(t_old,xyw,uw,uxyw)
                  CALL helix_maxUvsT(t_old,xyw,uw,uxyw)
                  CALL helix_Uprofile(nxs,xyw,xyw_kt,uw,uxyw)
                  CALL helix_Bspectr(t_old,xyw,xyw_kt,uw,uxyw)
               ENDIF
               IF(flag2d)THEN
                  CALL helix_contour(t_old,xyw,uw,uxyw,first2d,
     $                 last2d,nt_next,ifile,stride,file(ifile))
               ENDIF
            CASE("MAST")
               IF(flag1d)THEN
                  CALL postmast_energy(t_old,xyw,jac,uw,uxyw)
               ENDIF
               IF(flag2d)THEN
                  CALL postmast_contour(t_old,xyw,uw,uxyw,
     $                 first2d,last2d,nt_next,ifile,stride,file(ifile))
               ENDIF
            CASE("beltrami")
               IF(ifile > nfile - stride .AND. myios /= 0)THEN
                  CALL beltrami_grid(basis,x,y,u,polar_crd)
               ENDIF
            END SELECT
c-----------------------------------------------------------------------
c     write out dependent variable data 
c-----------------------------------------------------------------------
            IF(contour)THEN
               WHERE(ABS(uw) < min_value)
                  uw=0
               END WHERE
               SELECT CASE(out_type)
               CASE("hdf5")
                  CALL plotter_Ucontour(t_old,xyw,uw,first2d,last2d,
     $                 nt_next,xcont,ifile,stride,file(ifile),"post")
               CASE DEFAULT
                  CALL plotter_Ucontour(t_old,xyw,uw,first2d,last2d,
     $                 nt_next,xcont,ifile,stride,file(ifile),
     $                 "contour")
                  CALL plotter_Ucontour(t_old,xyw,duw,first2d,last2d,
     $                 nt_next,xcont,ifile,stride,file(ifile),
     $                 "dcontour")
               END SELECT
            ENDIF
            IF(ifile > 0 .AND. MOD(ifile,stride) == 0)THEN
               CALL post_griderr(nx,ny,np,nqty,u,griderr)
               CALL post_write_amp(t_old,griderr,duw)
            ENDIF
c-----------------------------------------------------------------------
c     end ifile loop 
c-----------------------------------------------------------------------
         ENDDO
         IF(myios /= 0)EXIT
      ENDDO
c-----------------------------------------------------------------------
c     clean up.
c-----------------------------------------------------------------------
      DEALLOCATE(x,y,u,uw,uxyw,duw,xyw,xyw_kt,jac,griderr,file,uw0,
     $     mass0)
      CALL jacobi_dealloc(basis)
      CALL jacobi_dealloc(quad)
      CLOSE(UNIT=header_unit)
      IF(xt_flag)CLOSE(UNIT=xt_unit)
      IF(yt_flag)CLOSE(UNIT=yt_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post_run
c-----------------------------------------------------------------------
c     subprogram 2. post_ls.
c     get directory listing.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post_ls(indir,mfile,nfile,file)

      CHARACTER(*), INTENT(IN) :: indir
      INTEGER, INTENT(IN) :: mfile
      INTEGER, INTENT(OUT) :: nfile
      CHARACTER(*), DIMENSION(:), POINTER, INTENT(OUT) :: file

      CHARACTER(11) :: label
      INTEGER :: ifile
c-----------------------------------------------------------------------
c     list data files.
c-----------------------------------------------------------------------
      IF(drawgrid)THEN
         label="/grid_*."
      ELSE
         label="/sel_*."
      ENDIF

      IF(outfile_type=="hdf5")THEN
         label=TRIM(label)//"h5"
      ELSE
         label=TRIM(label)//"dat"
      ENDIF

      CALL system("rm -f ls.out")
      IF(cray_flag)THEN
         CALL system("ls "//TRIM(indir)//TRIM(label)
     $        //" | ~/wc -l > ls.out")
         CALL system("ls "//TRIM(indir)//TRIM(label)
     $        //" | ~/sort >> ls.out")
      ELSE
         CALL system("ls "//TRIM(indir)//TRIM(label)
     $        //" | wc -l > ls.out")
         CALL system("ls "//TRIM(indir)//TRIM(label)
     $        //" >> ls.out")
      ENDIF
c-----------------------------------------------------------------------
c     read number of files and allocate space.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="ls.out",STATUS="OLD")
      READ(in_unit,*)nfile
      nfile=nfile-1
      nfile=MIN(mfile,nfile)
      ALLOCATE(file(0:nfile))
c-----------------------------------------------------------------------
c     read filenames.
c-----------------------------------------------------------------------
      DO ifile=0,nfile
         READ(in_unit,'(a)')file(ifile)
         file(ifile)=TRIM(file(ifile))
      ENDDO
c-----------------------------------------------------------------------
c     close and delete file.
c-----------------------------------------------------------------------
      CLOSE(UNIT=in_unit)
      CALL system("rm -f ls.out")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post_ls
c-----------------------------------------------------------------------
c     subprogram 4. post_write0.
c     write coordinate and header information
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post_write0(nxp,nyp,xx,yy,filename,xyw,xyw_kt,
     $     nxw0,nyw0,nxw,nyw,polar_crd)

      CHARACTER(*), INTENT(IN) :: filename
      LOGICAL, INTENT(IN) :: polar_crd
      INTEGER, INTENT(IN) :: nxp,nyp,nxw0,nyw0,nxw,nyw
      REAL(r8), DIMENSION(0:), INTENT(IN) :: xx,yy
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: xyw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(OUT) :: xyw_kt

      CHARACTER(160) :: gridname
      INTEGER :: ix,iqty,msize,nx,ny,h5_error
      INTEGER, DIMENSION(3) :: temp
      INTEGER(HSIZE_T), DIMENSION(3) :: xy_rdims
      INTEGER(HID_T) :: h5_xyid,xy_dset_id,h5_gid
      REAL(r8) :: xi,dx,dy
      REAL(r8), DIMENSION(SIZE(xyw,1),SIZE(xyw,2),SIZE(xyw,3)) :: 
     $     xyw_loc
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xw,yw
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uu
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      nx=SIZE(xx)-1
      ny=SIZE(yy)-1
      xyw_kt=0
c-----------------------------------------------------------------------
c     compute unstructured interpolated position array.
c-----------------------------------------------------------------------
      IF(adapt_grid .AND. .NOT. drawgrid)THEN
         ALLOCATE(uu(2,0:nxp,0:nyp))
         msize=LEN_TRIM(filename)
         SELECT CASE(outfile_type)
         CASE("hdf5")
            gridname = TRIM(filename(1:msize-13))//"/grid"//
     $           TRIM(filename(msize-8:msize))
            CALL H5OPEN_F(h5_error)
            CALL H5FOPEN_F(gridname,H5F_ACC_RDWR_F,h5_xyid,h5_error)
            CALL H5DOPEN_F(h5_xyid,"/AMP",xy_dset_id,h5_error)
            xy_rdims = (/SIZE(uu,1),SIZE(uu,2),SIZE(uu,3)/)
            CALL H5DREAD_F(xy_dset_id,H5T_NATIVE_DOUBLE,uu,xy_rdims,
     $           h5_error)
            CALL H5DCLOSE_F(xy_dset_id,h5_error)
            CALL H5FCLOSE_F(h5_xyid,h5_error)
            CALL H5CLOSE_F(h5_error)
         CASE DEFAULT
            gridname=TRIM(filename(1:msize-13))//"grid"//
     $           TRIM(filename(msize-9:msize))
            OPEN(UNIT=file_unit,FILE=TRIM(gridname),
     $           STATUS="UNKNOWN",FORM="UNFORMATTED")   
            READ(file_unit)uu
            CLOSE(UNIT=file_unit)
         END SELECT
         CALL plotter_interp(nxw0,nyw0,nx,ny,basis,uu,xyw,xyw_kt)
         DEALLOCATE(uu)
      ELSE
c-----------------------------------------------------------------------
c     compute structured interpolated position array.
c-----------------------------------------------------------------------
         dx=1._r8/REAL(nx,r8)
         dy=1._r8/REAL(ny,r8)
         ALLOCATE(xw(0:nxw),yw(0:nyw))
         xi=-1._r8
         DO ix=0,nxw0
            xw(ix:(nx-1)*nxw0+ix:nxw0) = xx(0:nx-1) + 0.5*dx*(1.+xi)
            xi = xi + 2./REAL(nxw0,r8)
         ENDDO
         xi=-1._r8
         DO ix=0,nyw0
            yw(ix:(ny-1)*nyw0+ix:nyw0) = yy(0:ny-1) + 0.5*dy*(1.+xi)
            xi = xi + 2./REAL(nyw0,r8)
         ENDDO

         DO ix=0,nxw
            xyw(1,ix,:)=xw(ix)
            xyw(2,ix,:)=yw
         ENDDO
         DEALLOCATE(xw,yw)
      ENDIF
c-----------------------------------------------------------------------
c     write HDF5 grid data (x and y)
c     if necessary, transform the coordinate system.
c-----------------------------------------------------------------------
      IF((contour .OR. flag2d) .AND. out_type == "hdf5")THEN
         IF(polar_crd)THEN
            xyw_loc(1,:,:)=xyw(1,:,:)*COS(xyw(2,:,:))
            xyw_loc(2,:,:)=xyw(1,:,:)*SIN(xyw(2,:,:))
         ELSE
            xyw_loc=xyw
         ENDIF
         CALL plotter_Ucontour(0._r8,xyw,xyw_loc,.FALSE.,.FALSE.,
     $        0,temp,0,1,filename,"grid")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post_write0
c-----------------------------------------------------------------------
c     subprogram 5. post_write_amp
c     write solution amplitude data for XDRAW.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post_write_amp(t,griderr,duw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:), INTENT(IN) :: griderr
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: duw

      INTEGER :: iqty,nqty
      REAL(r8), DIMENSION(SIZE(duw,1)) :: amp
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      nqty=SIZE(duw,1)
c-----------------------------------------------------------------------
c     write amp  output.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         amp(iqty)=MAXVAL(ABS(duw(iqty,:,:)))
      ENDDO
      OPEN(UNIT=bin_unit,FILE="amp.bin",STATUS="UNKNOWN",
     $     POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(bin_unit)REAL(t,4),REAL(amp,4),REAL(LOG10(griderr),4)
      CLOSE(UNIT=bin_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post_write_amp
c-----------------------------------------------------------------------
c     subprogram 6. post_griderr
c     computes grid error.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post_griderr(nx,ny,np,nqty,uu,griderr)

      INTEGER, INTENT(IN) :: nx,ny,np,nqty
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(r8), DIMENSION(:), INTENT(OUT) :: griderr

      INTEGER :: ix,jx,iy,jy,iqty
      REAL(r8) :: maxnorm
      REAL(r8), DIMENSION(0:np,0:np) :: u,delu
      REAL(r8), DIMENSION(nqty,nx,ny) :: norm,errx,erry
c-----------------------------------------------------------------------
c     compute norm and error.
c-----------------------------------------------------------------------
      jy=0
      DO iy=1,ny
         jx=0
         DO ix=1,nx
            DO iqty=1,nqty
               delu=0
               u=uu(iqty,jx:jx+np,jy:jy+np)
               norm(iqty,ix,iy)=SUM(u*MATMUL(MATMUL(mass0,u),mass0))
               delu=0
               delu(np-1,:)=u(np-1,:)
               errx(iqty,ix,iy)=SUM(delu*MATMUL(MATMUL(mass0,delu),
     $              mass0))
               delu=0
               delu(:,np-1)=u(:,np-1)
               erry(iqty,ix,iy)=SUM(delu*MATMUL(MATMUL(mass0,delu),
     $              mass0))
            ENDDO
            jx=jx+np+1
         ENDDO
         jy=jy+np+1
      ENDDO

      DO iqty=1,nqty
         maxnorm=MAXVAL(norm(iqty,:,:))
         WHERE(norm(iqty,:,:) <= maxnorm*min_eps**.25)
            norm(iqty,:,:)=1./min_eps
         ENDWHERE
      ENDDO
      WHERE(norm <= min_eps)
         norm=1./min_eps
      ENDWHERE

      griderr=0
      DO iqty=1,nqty
         griderr(iqty)=
     $        MAXVAL(SQRT(ABS(errx(iqty,:,:)/norm(iqty,:,:))))
         griderr(iqty+nqty)=
     $        MAXVAL(SQRT(ABS(erry(iqty,:,:)/norm(iqty,:,:))))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post_griderr
c-----------------------------------------------------------------------
c     subprogram 7. post_mass0.
c     computes elementary mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post_mass0

      INTEGER :: ix,m,n
      REAL(r8) :: x,wt
c-----------------------------------------------------------------------
c     compute block mass matrix, upper half.
c-----------------------------------------------------------------------
      mass0=0
      DO ix=1,quad%np
         x=quad%pzero(ix)
         wt=quad%weight(ix)
         CALL jacobi_basis(x,basis)
         DO m=0,basis%np
            DO n=m,basis%np
               mass0(m,n)=mass0(m,n) + wt*basis%pb(m)*basis%pb(n)
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     symmetrize.
c-----------------------------------------------------------------------
      DO m=0,basis%np
         DO n=0,m-1
            mass0(m,n)=mass0(n,m)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post_mass0
c-----------------------------------------------------------------------
c     subprogram 9. post_read_header.
c     reads in time-step/grid data from header.dat or header.txt file
c     and determines the number of output files with given grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post_read_header(header_dat,nfile,gout,jfile,myios,
     $     firstfile,lastfile,nt_next)

      LOGICAL, INTENT(IN) :: header_dat
      INTEGER, INTENT(IN) :: nfile
      INTEGER, INTENT(INOUT) :: gout,jfile,myios,lastfile
      INTEGER, INTENT(OUT) :: firstfile,nt_next

      INTEGER :: iout,regrid_step
      REAL(r8) :: dt_temp
c-----------------------------------------------------------------------
c     read timestep vector
c-----------------------------------------------------------------------
      dt_temp=0.
      firstfile=lastfile+1
      IF(header_dat)THEN
         DO
            DO iout=gout,dmout
               READ(header_unit,IOSTAT=myios)dt_temp
               IF(dt_temp == -1.0 .OR. myios /= 0)EXIT
            ENDDO
            IF(dt_temp == -1.0)THEN
               gout=iout
               EXIT
            ENDIF
            IF(myios /= 0 .OR. jfile == nfile)EXIT
            gout=1
            jfile=jfile+1
         ENDDO
      ELSE
         jfile=nfile+1
      ENDIF
c-----------------------------------------------------------------------
c     determine how many output files there are with given grid
c-----------------------------------------------------------------------
      IF(jfile == nfile)THEN
         lastfile=nfile
         myios=-1
      ELSE
         IF(header_dat)THEN
            READ(header_unit,IOSTAT=myios)regrid_step
         ELSE
            READ(header_unit,'(i8)',IOSTAT=myios)regrid_step
         ENDIF
         IF(myios /= 0)THEN
            lastfile=nfile
         ELSE
            IF(MOD(regrid_step,dmout) == 0)THEN
               regrid_step=regrid_step/dmout
            ELSE
               regrid_step=regrid_step/dmout+1
            ENDIF
            lastfile=MIN(regrid_step-1,jfile)
         ENDIF
      ENDIF
      nt_next=lastfile/stride-firstfile/stride
      IF(MOD(firstfile,stride)==0)nt_next=nt_next+1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post_read_header
      END MODULE post_mod
c-----------------------------------------------------------------------
c     subprogram 10. post_main.
c     trivial main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM post_main
      USE post_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     do it.
c-----------------------------------------------------------------------
      OPEN(UNIT=out_unit,FILE="post.out",STATUS="UNKNOWN")
      CALL timer(0,out_unit)
      CALL post_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination.")
      END PROGRAM post_main
      
