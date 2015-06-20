c-----------------------------------------------------------------------
c     file plotter.f.
c     subroutines for various ways to plot post-processed data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  plotter_mod.
c     1.  plotter_UxyT.
c     2.  plotter_VecSqInt.
c     2b.  plotter_integral.
c     3.  plotter_Ucontour.
c     4.  plotter_Uprofile.
c     5.  plotter_dUdt.
c     6.  plotter_maxUvsT.
c     8.  plotter_width.
c     10. plotter_write_xdmf.
c     11. plotter_interp.
c     12. plotter_transf.
c     13. plotter_newton_search.
c     14. plotter_newton_residual.
c     15. plotter_evaluate.
c-----------------------------------------------------------------------
c     subprogram 0. plotter_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE plotter_mod
      USE slice_mod
      USE HDF5
      IMPLICIT NONE

      CHARACTER(160), PRIVATE :: h5_gfile="."

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. plotter_UxyT.
c     generate u(x,y) vs. time plot for a given physical quantity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_UxyT(t,x0,y0,xyw,xyw_kt,uw,draw,value,err)

      REAL(r8), INTENT(IN) :: t,x0,y0
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt
      LOGICAL, INTENT(IN) :: draw
      REAL(r8), DIMENSION(:), INTENT(OUT) :: value
      LOGICAL, INTENT(OUT) :: err

      LOGICAL :: check
      INTEGER :: ix,iy
      REAL(r8), DIMENSION(2) :: xy0,kt0
c-----------------------------------------------------------------------
c     find the location of (x0,y0) point in xyw matrix.
c-----------------------------------------------------------------------
      ix=0
      iy=0
      xy0=(/x0,y0/)
      kt0=.5
      value=0.
      CALL plotter_newton_search(xy0,xyw,xyw_kt,kt0,check)

      IF(.NOT. check)THEN
         err=.TRUE.
         RETURN
      ELSE
         err=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c     find the value of u at (x0,y0).
c-----------------------------------------------------------------------
      ix = INT(kt0(1)*(SIZE(xyw,2)-1))
      iy = INT(kt0(2)*(SIZE(xyw,3)-1))

      value = plotter_evaluate(x0,y0,xyw(:,ix:ix+1,iy:iy+1),
     $     uw(:,ix:ix+1,iy:iy+1))
c-----------------------------------------------------------------------
c     open, write, and close UxyT.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      IF(draw)THEN
         OPEN(UNIT=UxyT_unit,FILE="UxyT.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WHERE(value == 0)
            value=TINY(value)
         END WHERE
         WRITE(UxyT_unit)REAL(t,4),REAL(value,4)
         CLOSE(UNIT=UxyT_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_UxyT
c-----------------------------------------------------------------------
c     subprogram 2. plotter_VecSqInt.
c     generate a time plot for an integral of a square of a vector
c     quantity over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_VecSqInt(t,xyw,uw1,uw2,uw3,draw,value)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw1,uw2,uw3
      LOGICAL, INTENT(IN) :: draw
      REAL(r8), DIMENSION(:), INTENT(OUT) :: value

      INTEGER :: ix,iy,iqty
      REAL(r8), DIMENSION(5) :: a
c-----------------------------------------------------------------------
c     calculate the integral over the domain.
c-----------------------------------------------------------------------
      value=0
      DO ix=0,SIZE(xyw,2)-2
         DO iy=0,SIZE(xyw,3)-2
            a(1)=0.5*(xyw(2,ix,iy)+xyw(2,ix,iy+1))*
     $           (xyw(1,ix,iy+1)-xyw(1,ix,iy))
            a(2)=0.5*(xyw(2,ix,iy+1)+xyw(2,ix+1,iy+1))*
     $           (xyw(1,ix+1,iy+1)-xyw(1,ix,iy+1))
            a(3)=0.5*(xyw(2,ix+1,iy+1)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix+1,iy)-xyw(1,ix+1,iy+1))
            a(4)=0.5*(xyw(2,ix,iy)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix,iy)-xyw(1,ix+1,iy))
            a(5)=ABS(SUM(a(1:4)))
            DO iqty=1,SIZE(uw1,1)
               value(iqty)=value(iqty) + 0.25*a(5)*
     $              SUM(uw1(iqty,ix:ix+1,iy:iy+1)**2 +
     $              uw2(iqty,ix:ix+1,iy:iy+1)**2 +
     $              uw3(iqty,ix:ix+1,iy:iy+1)**2)
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     open, write, and close VecSqInt.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      IF(draw)THEN
         OPEN(UNIT=VecSqInt_unit,FILE="VecSqInt.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(VecSqInt_unit)REAL(t,4),REAL(value,4)
         CLOSE(UNIT=VecSqInt_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_VecSqInt
c-----------------------------------------------------------------------
c     subprogram 2b. plotter_integral.
c     generate a time plot for an integral of a quantity 
c     over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_integral(t,xyw,uw1,draw,value)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw1
      LOGICAL, INTENT(IN) :: draw
      REAL(r8), DIMENSION(:), INTENT(OUT) :: value

      INTEGER :: ix,iy,iqty
      REAL(r8), DIMENSION(5) :: a
c-----------------------------------------------------------------------
c     calculate the integral over the domain.
c-----------------------------------------------------------------------
      value=0
      DO ix=0,SIZE(xyw,2)-2
         DO iy=0,SIZE(xyw,3)-2
            a(1)=0.5*(xyw(2,ix,iy)+xyw(2,ix,iy+1))*
     $           (xyw(1,ix,iy+1)-xyw(1,ix,iy))
            a(2)=0.5*(xyw(2,ix,iy+1)+xyw(2,ix+1,iy+1))*
     $           (xyw(1,ix+1,iy+1)-xyw(1,ix,iy+1))
            a(3)=0.5*(xyw(2,ix+1,iy+1)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix+1,iy)-xyw(1,ix+1,iy+1))
            a(4)=0.5*(xyw(2,ix,iy)+xyw(2,ix+1,iy))*
     $           (xyw(1,ix,iy)-xyw(1,ix+1,iy))
            a(5)=ABS(SUM(a(1:4)))
            DO iqty=1,SIZE(uw1,1)
               value(iqty)=value(iqty) + 0.25*a(5)*
     $              SUM(uw1(iqty,ix:ix+1,iy:iy+1))
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     open, write, and close VecSqInt.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      IF(draw)THEN
         OPEN(UNIT=VecSqInt_unit,FILE="integral.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(VecSqInt_unit)REAL(t,4),REAL(value,4)
         CLOSE(UNIT=VecSqInt_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_integral
c-----------------------------------------------------------------------
c     subprogram 3. plotter_Ucontour.
c     generates file(s) with 2D data for a set of physical variables
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_Ucontour(t,xyw,uw,first,last,
     $     nt_next,xcont_prev,ifile,stride,fname_in,fname_out)

      LOGICAL, INTENT(IN) :: first,last
      CHARACTER(*), INTENT(IN) :: fname_in,fname_out
      INTEGER, INTENT(IN) :: nt_next,ifile,stride
      INTEGER, DIMENSION(3), INTENT(INOUT) :: xcont_prev
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw

      CHARACTER(4) :: fname
      CHARACTER(10) :: var_name
      CHARACTER(160) :: filename
      INTEGER, PARAMETER :: uw_rank=2
      INTEGER :: i,iqty,nqty,nxw,nyw,msize,h5_error
      INTEGER(HID_T) :: h5_uid,dset_id,uw_dspace
      INTEGER(HSIZE_T), DIMENSION(2) :: uw_dims
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uu
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      nxw=SIZE(xyw,2)-1
      nyw=SIZE(xyw,3)-1
c-----------------------------------------------------------------------
c     construct the filename.
c-----------------------------------------------------------------------
      SELECT CASE(out_type)
      CASE("hdf5")
         msize=LEN_TRIM(fname_in)
         fname = "____"
         fname(1:MIN(4,LEN_TRIM(fname_out))) 
     $        = fname_out(1:MIN(4,LEN_TRIM(fname_out)))
         SELECT CASE(outfile_type)
         CASE("hdf5")
            WRITE(filename,'(A)')TRIM(postout)//'/'//
     $           TRIM(fname)//TRIM(fname_in(msize-8:msize))
         CASE DEFAULT
            WRITE(filename,'(A)')TRIM(postout)//'/'//
     $           TRIM(fname)//TRIM(fname_in(msize-9:msize-4))//'.h5'
         END SELECT
      CASE DEFAULT
         WRITE(filename,'(A)')TRIM(fname_out)//'.bin'
      END SELECT         
c-----------------------------------------------------------------------
c     when needed, write initial header for XDRAW.
c-----------------------------------------------------------------------
      IF(first .AND. out_type /= "hdf5")THEN
         OPEN(UNIT=Ucontour_unit,FILE=TRIM(filename),STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(Ucontour_unit)zero,nxw,nyw,nt_next,xcont_prev
c-----------------------------------------------------------------------
c     if necessary, modify coordinate system for the grid. 
c-----------------------------------------------------------------------
         IF(polar)THEN
            ALLOCATE(uu(2,0:SIZE(xyw,2)-1,0:SIZE(xyw,3)-1))
            uu(1,:,:)=xyw(1,:,:)*COS(xyw(2,:,:))
            uu(2,:,:)=xyw(1,:,:)*SIN(xyw(2,:,:))
            WRITE(Ucontour_unit)(REAL(uu(iqty,:,:),4),iqty=1,2)
            DEALLOCATE(uu)
         ELSE
            WRITE(Ucontour_unit)(REAL(xyw(iqty,:,:),4),iqty=1,2)
         ENDIF
         xcont_prev(1)=nxw
         xcont_prev(2)=nyw
         xcont_prev(3)=nt_next
         CLOSE(UNIT=Ucontour_unit)
      ENDIF
c-----------------------------------------------------------------------
c     write the data.
c-----------------------------------------------------------------------
      IF(MOD(ifile,stride) == 0)THEN
         nqty=SIZE(uw,1)

         SELECT CASE(out_type)
         CASE("hdf5")
            uw_dims=(/nxw+1,nyw+1/)
            CALL H5OPEN_F(h5_error)
            CALL H5FCREATE_F(filename,H5F_ACC_TRUNC_F,h5_uid,h5_error)
            CALL H5SCREATE_SIMPLE_F(uw_rank,uw_dims,uw_dspace,h5_error)
            DO iqty = 1,nqty
               WRITE(var_name,'(A,i2.2)')'/U',iqty
               CALL H5DCREATE_F(h5_uid,var_name,H5T_NATIVE_DOUBLE,
     $              uw_dspace,dset_id,h5_error)
               CALL H5DWRITE_F(dset_id,H5T_NATIVE_DOUBLE,
     $              uw(iqty,:,:),uw_dims,h5_error)
               CALL H5DCLOSE_F(dset_id,h5_error)
            ENDDO                  
            CALL H5SCLOSE_F(uw_dspace,h5_error)
            CALL H5FCLOSE_F(h5_uid,h5_error)
            CALL H5CLOSE_F(h5_error)
            CALL plotter_write_xdmf(nqty,nxw+1,nyw+1,t,filename,fname)
         CASE DEFAULT
            OPEN(UNIT=Ucontour_unit,FILE=TRIM(filename),
     $           STATUS="UNKNOWN",POSITION="APPEND",FORM="UNFORMATTED")
            WRITE(Ucontour_unit)one
            DO iqty=1,nqty
               WRITE(Ucontour_unit)REAL(uw(iqty,:,:),4)
            ENDDO
            CLOSE(UNIT=Ucontour_unit)
         END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     when needed, write closing header for XDRAW.
c-----------------------------------------------------------------------
      IF(last .AND. out_type /= "hdf5")THEN
         OPEN(UNIT=Ucontour_unit,FILE=TRIM(filename),STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WRITE(Ucontour_unit)zero,zero,zero,zero,xcont_prev
         CLOSE(UNIT=Ucontour_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_Ucontour
c-----------------------------------------------------------------------
c     subprogram 4. plotter_Uprofile.
c     generate u(x,y) vs. (x',y') for a given physical quantity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_Uprofile(x0,y0,x1,y1,nxs,xyw,xyw_kt,uw,
     $     uprofile,wflag,filename)

      INTEGER, INTENT(IN) :: nxs
      REAL(r8), INTENT(IN) :: x0,y0,x1,y1
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt
      REAL(r8), DIMENSION(:,0:), INTENT(OUT) :: uprofile
      LOGICAL, INTENT(IN) :: wflag
      CHARACTER(*), INTENT(IN) :: filename

      INTEGER :: ipts
      LOGICAL :: err
      REAL(r8) :: x,y,dist,dl,dx,dy,t
      REAL(r8), DIMENSION(SIZE(uw,1)) :: value,value_old
c-----------------------------------------------------------------------
c     open Uprofile.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      IF(wflag)OPEN(UNIT=Uprofile_unit,FILE=TRIM(filename)//".bin",
     $     STATUS="UNKNOWN",POSITION="APPEND",FORM="UNFORMATTED")
c-----------------------------------------------------------------------
c     prepare for constructing the profile
c-----------------------------------------------------------------------
      dl=SQRT((x1-x0)**2+(y1-y0)**2)/nxs
      dx=(x1-x0)/nxs
      dy=(y1-y0)/nxs
      dist=0
      x=x0
      y=y0
      t=0
      value=0
      value_old=0
      DO ipts=0,nxs
         CALL plotter_UxyT(t,x,y,xyw,xyw_kt,uw,.FALSE.,value,err)
         IF(err)value=value_old
c-----------------------------------------------------------------------
c        write to Uprofile.bin
c-----------------------------------------------------------------------
         IF(wflag)WRITE(Uprofile_unit)REAL(dist,4),REAL(value,4)
         uprofile(1,ipts)=dist
         uprofile(2:SIZE(value)+1,ipts)=value
c-----------------------------------------------------------------------
c        set up for the next point along the profile
c-----------------------------------------------------------------------
         dist=dist+dl
         x=x+dx
         y=y+dy
         value_old=value
      ENDDO
      IF(wflag)WRITE(Uprofile_unit)
c-----------------------------------------------------------------------
c     close Uprofile.bin file.
c-----------------------------------------------------------------------
      IF(wflag)CLOSE(UNIT=Uprofile_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_Uprofile
c-----------------------------------------------------------------------
c     subprogram 5. plotter_dUdt.
c     generate du/dt vs. time plot for a given physical quantity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_dUdt(t,t_old,uw,uw_old,first)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), INTENT(INOUT) :: t_old
      REAL(r8), DIMENSION(:), INTENT(IN) :: uw
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: uw_old
      LOGICAL, INTENT(IN) :: first

      REAL(r8), DIMENSION(SIZE(uw)) :: value
c-----------------------------------------------------------------------
c     calculate value of the time derivative
c-----------------------------------------------------------------------
      IF(.NOT. first)THEN
         value=(uw-uw_old)/(t-t_old)
c-----------------------------------------------------------------------
c     open, write, and close UxyT.bin file to store data for xdraw.
c-----------------------------------------------------------------------
         OPEN(UNIT=dUdt_unit,FILE="dUdt.bin",STATUS="UNKNOWN",
     $        POSITION="APPEND",FORM="UNFORMATTED")
         WHERE(value == 0)
            value=TINY(value)
         END WHERE
         WRITE(dUdt_unit)REAL((t+t_old)/2.,4),REAL(value,4)
         CLOSE(UNIT=dUdt_unit)
      ENDIF
      uw_old=uw
      t_old=t
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_dUdt
c-----------------------------------------------------------------------
c     subprogram 6. plotter_maxUvsT.
c     generate MAX(u(x,y)) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_maxUvsT(filename,t,uw,value)

      CHARACTER(*), INTENT(IN) :: filename
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uw
      REAL(r8), DIMENSION(:), INTENT(OUT) :: value

      INTEGER :: i
c-----------------------------------------------------------------------
c     find the maximum of u.
c-----------------------------------------------------------------------
      value=0
      DO i=1,SIZE(uw,1)
         value(i)=MAXVAL(ABS(uw(i,:,:)))
      ENDDO
c-----------------------------------------------------------------------
c     open, write, and close UxyT.bin file to store data for xdraw.
c-----------------------------------------------------------------------
      OPEN(UNIT=maxUvsT_unit,FILE=TRIM(filename)//".bin",
     $     STATUS="UNKNOWN",POSITION="APPEND",FORM="UNFORMATTED")
      WRITE(maxUvsT_unit)REAL(t,4),REAL(value,4)
      CLOSE(UNIT=maxUvsT_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_maxUvsT
c-----------------------------------------------------------------------
c     subprogram 8. plotter_width.
c     finds half-width of a peak from (x0,y0) to (x1,y1).
c------------------------------------------------------------------ -----
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_width(x0,y0,x1,y1,nxs,xyw,xyw_kt,uw,wtype,
     $     width)

      INTEGER, INTENT(IN) :: nxs
      REAL(r8), INTENT(IN) :: x0,y0,x1,y1
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uw,xyw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt
      CHARACTER(*), INTENT(IN) :: wtype
      REAL(r8), INTENT(OUT) :: width

      INTEGER :: ipts
      LOGICAL :: err
      REAL(r8) :: x,y,dist,dl,dx,dy,t,val_width
      REAL(r8), DIMENSION(1,0:nxs) :: value
c-----------------------------------------------------------------------
c     prepare for constructing the profile
c-----------------------------------------------------------------------
      dl=SQRT((x1-x0)**2+(y1-y0)**2)/nxs
      dx=(x1-x0)/nxs
      dy=(y1-y0)/nxs
      dist=0
      x=x0
      y=y0
      t=0
      value=0
      DO ipts=0,nxs
         CALL plotter_UxyT(t,x,y,xyw,xyw_kt,uw,.FALSE.,value(:,ipts),
     $        err)
         IF(err)value(:,ipts)=value(:,MAX(ipts-1,0))
c-----------------------------------------------------------------------
c        set up for the next point along the profile
c-----------------------------------------------------------------------
         dist=dist+dl
         x=x+dx
         y=y+dy
      ENDDO
c-----------------------------------------------------------------------
c     find half-width
c-----------------------------------------------------------------------
      value=ABS(value)
      width=0

      SELECT CASE(wtype)
      CASE("half_max")
         val_width=0.5*MAXVAL(value)
         DO ipts=1,nxs
            IF((value(1,ipts) <= val_width) .AND.
     $           (value(1,ipts-1) > val_width))
     $           width=dl*(ipts - 1 + (value(1,ipts-1) - val_width)/
     $           (value(1,ipts-1) - value(1,ipts)))
         ENDDO
      CASE("local_max")
         val_width=0
         DO ipts=1,nxs
            IF((value(1,ipts) <= val_width) .AND.
     $           (value(1,ipts-1) <= val_width)
     $           .AND. val_width > MAXVAL(value)/10.)THEN
               width=dl*ipts
               RETURN
            ENDIF
            val_width=MAX(val_width,value(1,ipts))
         ENDDO
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_width
c-----------------------------------------------------------------------
c     subprogram 10. plotter_write_xdmf
c     write XDMF (.xmf) file that corresponds to the HDF5 data files.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_write_xdmf(nqty,nxw,nyw,t,h5_ufile,post_name)

      CHARACTER(4), INTENT(IN) :: post_name
      CHARACTER(*), INTENT(IN) :: h5_ufile
      INTEGER, INTENT(IN) :: nqty,nxw,nyw
      REAL(r8), INTENT(IN) :: t

      CHARACTER(10) :: var_name
      CHARACTER(160) :: xdmf_file,h5_gfile_loc,h5_ufile_loc
      INTEGER :: iqty,fsize
c-----------------------------------------------------------------------
c     set full path/name for the grid file.
c-----------------------------------------------------------------------
      IF(post_name == "grid")THEN
         h5_gfile=h5_ufile
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     get file names.
c-----------------------------------------------------------------------
      fsize = LEN_TRIM(h5_ufile)
      h5_ufile_loc = h5_ufile(fsize-12:fsize)
      WRITE(xdmf_file,'(A)')TRIM(h5_ufile(1:fsize-3))//".xmf"
      fsize = LEN_TRIM(h5_gfile)
      h5_gfile_loc = h5_gfile(fsize-12:fsize)
c-----------------------------------------------------------------------
c     write XDMF (.xmf) file for VisIt visualization
c-----------------------------------------------------------------------
      OPEN(UNIT=xdmf_unit,FILE=TRIM(xdmf_file),STATUS="REPLACE",
     $     ACTION="WRITE",POSITION="APPEND")
      WRITE(xdmf_unit,'(A)')'<?xml version="1.0" ?>'
      WRITE(xdmf_unit,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      WRITE(xdmf_unit,'(A)')'<Xdmf Version="2.0">'
      WRITE(xdmf_unit,'(A)')' <Domain>'
      WRITE(xdmf_unit,'(A)')'   <Grid Name="mesh2D" GridType="Uniform">'
      WRITE(xdmf_unit,'(A,e12.5,A)')'     <Time Value="',t,'" />'
      WRITE(xdmf_unit,'(A,i7,1x,i7,A)')'     <Topology 
     $     TopologyType="2DSMesh" NumberOfElements="',nyw,nxw,'"/>'
      WRITE(xdmf_unit,'(A)')'     <Geometry GeometryType="X_Y">'
      WRITE(xdmf_unit,'(A,i7,1x,i7,A)')'       <DataItem 
     $     Dimensions="',nyw,nxw,'" NumberType="Float" Precision="4" 
     $     Format="HDF">'
      WRITE(xdmf_unit,'(A)')'         '//TRIM(h5_gfile_loc)//':/U01'
      WRITE(xdmf_unit,'(A)')'       </DataItem>'
      WRITE(xdmf_unit,'(A,i7,1x,i7,A)')'       <DataItem Dimensions="',
     $     nyw,nxw,'" NumberType="Float" Precision="4" Format="HDF">'
      WRITE(xdmf_unit,'(A)')'         '//TRIM(h5_gfile_loc)//':/U02'
      WRITE(xdmf_unit,'(A)')'       </DataItem>'
      WRITE(xdmf_unit,'(A)')'     </Geometry>'
      DO iqty = 1,nqty
         WRITE(xdmf_unit,'(A,i2.2,A)')'     <Attribute Name="U',iqty,
     $        '" AttributeType="Scalar" Center="Node">'
         WRITE(xdmf_unit,'(A,i7,1x,i7,A)')'       <DataItem 
     $        Dimensions="',nyw,nxw,'" NumberType="Float" Precision="4" 
     $        Format="HDF">'
         WRITE(var_name,'(A,i2.2)')'/U',iqty
         WRITE(xdmf_unit,'(A)')'          '//TRIM(h5_ufile_loc)//':'
     $        //TRIM(var_name)
         WRITE(xdmf_unit,'(A)')'       </DataItem>'
         WRITE(xdmf_unit,'(A)')'     </Attribute>'
      ENDDO
      WRITE(xdmf_unit,'(A)')'   </Grid>'
      WRITE(xdmf_unit,'(A)')' </Domain>'
      WRITE(xdmf_unit,'(A)')'</Xdmf>'
      CLOSE(UNIT=xdmf_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_write_xdmf
c-----------------------------------------------------------------------
c     subprogram 11. plotter_interp.
c     interpolate solution vector in the spectral element representation
c     onto a grid uniformly spaced in the logical space.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_interp(nxw0,nyw0,nx,ny,basis,uu,uw,uxyw)

      INTEGER, INTENT(IN) :: nxw0,nyw0,nx,ny
      TYPE(jacobi_type) :: basis
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(OUT) :: uxyw

      INTEGER :: ix,iy,ixx,iyy,jx,jy,np,imx,imy,iqty
      REAL(r8) :: xi,yi,dxi,dyi,dxfac,dyfac
      REAL(r8), DIMENSION(0:basis%np) :: pbx,pby,qbx,qby
      REAL(r8), DIMENSION(SIZE(uw,1),0:nxw0,0:nyw0) :: luw
      REAL(r8), DIMENSION(SIZE(uw,1),2,0:nxw0,0:nyw0) :: luxyw
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      np=basis%np
c-----------------------------------------------------------------------
c     interpolate solution.
c-----------------------------------------------------------------------
      uw=0
      uxyw=0

      dxfac=two*nx
      dyfac=two*ny
      dxi=two/REAL(nxw0,r8)
      dyi=two/REAL(nyw0,r8)

      jy=0
      DO iyy=0,ny-1
         jx=0
         DO ixx=0,nx-1

            luw=0
            luxyw=0

            yi=-1._r8
            DO iy=0,nyw0
c-----------------------------------------------------------------------
c     compute y positions, basis functions.
c-----------------------------------------------------------------------
               CALL jacobi_basis(yi,basis)
               pby=basis%pb
               qby=basis%qb*dyfac

               xi=-1._r8
               DO ix=0,nxw0
c-----------------------------------------------------------------------
c     compute x positions, basis functions.
c-----------------------------------------------------------------------
                  CALL jacobi_basis(xi,basis)
                  pbx=basis%pb
                  qbx=basis%qb*dxfac
c-----------------------------------------------------------------------
c     interpolate solutions and derivatives.
c-----------------------------------------------------------------------
                  DO imy=0,np
                     DO imx=0,np
                        luw(:,ix,iy) = luw(:,ix,iy)
     $                       + pbx(imx)*pby(imy)*uu(:,imx+jx,imy+jy)
                        luxyw(:,1,ix,iy) = luxyw(:,1,ix,iy)
     $                       + qbx(imx)*pby(imy)*uu(:,imx+jx,imy+jy)
                        luxyw(:,2,ix,iy) = luxyw(:,2,ix,iy)
     $                       + pbx(imx)*qby(imy)*uu(:,imx+jx,imy+jy)
                     ENDDO
                  ENDDO
                  xi=xi+dxi
               ENDDO
               yi=yi+dyi
            ENDDO

            uw(:,ixx*nxw0:(1+ixx)*nxw0,iyy*nyw0:(1+iyy)*nyw0) = luw
            DO iqty=1,SIZE(uw,1)
               uxyw(ixx*nxw0:(1+ixx)*nxw0,iyy*nyw0:(1+iyy)*nyw0,1,iqty) 
     $              = luxyw(iqty,1,:,:)
               uxyw(ixx*nxw0:(1+ixx)*nxw0,iyy*nyw0:(1+iyy)*nyw0,2,iqty) 
     $              = luxyw(iqty,2,:,:)
            ENDDO
c-----------------------------------------------------------------------
c     finish outer loops.
c-----------------------------------------------------------------------
            jx=jx+np+1
         ENDDO
         jy=jy+np+1
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_interp
c-----------------------------------------------------------------------
c     subprogram 12. plotter_transf.
c     transform from a gradient w.r.t logical coordinates
c     to a gradint w.r.t physical coordinates.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_transf(xyw_kt,jac,uxyw)

      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: jac
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(INOUT) :: uxyw

      INTEGER :: iqty
      REAL(r8), DIMENSION(0:SIZE(uxyw,1)-1,0:SIZE(uxyw,2)-1) :: u_temp
c-----------------------------------------------------------------------
c     transform.
c-----------------------------------------------------------------------
      jac = xyw_kt(:,:,1,1)*xyw_kt(:,:,2,2)
     $     - xyw_kt(:,:,1,2)*xyw_kt(:,:,2,1)

      WHERE(ABS(jac) <= 1e-20)
         jac=SIGN(1._r8,jac)
      END WHERE

      DO iqty=1,SIZE(uxyw,4)
         u_temp=uxyw(:,:,1,iqty)

         uxyw(:,:,1,iqty)=(u_temp*xyw_kt(:,:,2,2) 
     $        - uxyw(:,:,2,iqty)*xyw_kt(:,:,1,2))/jac

         uxyw(:,:,2,iqty)=(-u_temp*xyw_kt(:,:,2,1)
     $        + uxyw(:,:,2,iqty)*xyw_kt(:,:,1,1))/jac
      ENDDO 
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_transf
c-----------------------------------------------------------------------
c     subprogram 13. plotter_newton_search.
c     finds the logical coordinates kt0 of a point with coordinates 
c     xy0 in the physical space via Newton iteration.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_newton_search(xy0,xyw,xyw_kt,kt0,check)

      REAL(r8), DIMENSION(2), INTENT(IN) :: xy0
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt
      REAL(r8), DIMENSION(2), INTENT(INOUT) :: kt0
      LOGICAL, INTENT(OUT) :: check

      INTEGER, PARAMETER :: max_it=100
      INTEGER :: ii,jj
      INTEGER, DIMENSION(2) :: nkt,ixy0
      REAL(r8), PARAMETER :: errtol = 1.e-4
      REAL(r8) :: norm
      REAL(r8), DIMENSION(2) :: res0
      REAL(r8), DIMENSION(2,2) :: jac
      REAL(r8), DIMENSION(1,2,2) :: temp_array
      REAL(r8), DIMENSION(2,2,2) :: kt_array
      REAL(r8), DIMENSION(2,2,2,2) :: identity,transf_array
c-----------------------------------------------------------------------
c     set initial values.
c-----------------------------------------------------------------------
      check=.FALSE.
      identity(:,:,1,1) = 1.
      identity(:,:,1,2) = 0.
      identity(:,:,2,1) = 0.
      identity(:,:,2,2) = 1.
      nkt=(/SIZE(xyw,2)-1,SIZE(xyw,3)-1/)
c-----------------------------------------------------------------------
c     calculate a norm for the Newton convergence criteria.
c-----------------------------------------------------------------------
      ixy0 = INT(kt0*nkt)
      kt_array = xyw(:,ixy0(1):ixy0(1)+1,ixy0(2):ixy0(2)+1)
      CALL plotter_newton_residual(nkt,ixy0,kt0,xy0,kt_array,res0)
      norm = MAX(SQRT(SUM(res0**2)),min_eps)
c-----------------------------------------------------------------------
c     check for some errors.
c-----------------------------------------------------------------------
      IF(MAXVAL(xyw(1,:,:)) < xy0(1) .OR. MAXVAL(xyw(2,:,:)) < xy0(2) 
     $     .OR. MINVAL(xyw(1,:,:)) > xy0(1)
     $     .OR. MINVAL(xyw(2,:,:)) > xy0(2))THEN
         WRITE(*,'(a,1p,2e12.3,a)')
     $        "Location (x,y)=(",xy0,") is outside of the domain."
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     iterate until convergence or maximum number of iterations.
c-----------------------------------------------------------------------
      DO jj=1,max_it
         transf_array=identity
         CALL plotter_transf
     $        (xyw_kt(ixy0(1):ixy0(1)+1,ixy0(2):ixy0(2)+1,:,:),jac,
     $        transf_array)
         DO ii=1,2
            temp_array(1,:,:) = transf_array(:,:,1,ii) 
            jac(1:1,ii) = plotter_evaluate(kt0(1),kt0(2),kt_array,
     $           temp_array)
            temp_array(1,:,:) = transf_array(:,:,2,ii) 
            jac(2:2,ii) = plotter_evaluate(kt0(1),kt0(2),kt_array,
     $           temp_array)
         ENDDO

         kt0 = kt0 - MATMUL(TRANSPOSE(jac),res0)
         kt0 = MIN(MAX(kt0,0._r8),1._r8-min_eps)

         ixy0 = INT(kt0*nkt)
         kt_array = xyw(:,ixy0(1):ixy0(1)+1,ixy0(2):ixy0(2)+1)
         CALL plotter_newton_residual(nkt,ixy0,kt0,xy0,kt_array,res0)

         IF(SQRT(SUM(res0**2)) < errtol*norm)THEN
            check=.TRUE.
            RETURN
         ENDIF
      ENDDO
      WRITE(*,'(a,1p,2e12.3,a)')
     $     "Could not locate (x,y)=(",xy0,") in the domain."
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_newton_search
c-----------------------------------------------------------------------
c     subprogram 14. plotter_newton_residual.
c     evaluates the Newton residual for subroutine plotter_newton_search
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plotter_newton_residual(nkt,ixy0,kt0,xy0,xyw,res0)

      INTEGER, DIMENSION(2), INTENT(IN) :: nkt,ixy0
      REAL(r8), DIMENSION(2), INTENT(INOUT) :: kt0
      REAL(r8), DIMENSION(2), INTENT(IN) :: xy0
      REAL(r8), DIMENSION(2,2,2), INTENT(INOUT) :: xyw
      REAL(r8), DIMENSION(2), INTENT(OUT) :: res0

      REAL(r8), DIMENSION(2,2,2) :: kt_array
c-----------------------------------------------------------------------
c     set initial values.
c-----------------------------------------------------------------------
      WHERE(kt0 < 0.)
         kt0 = zero
      END WHERE
      WHERE(kt0 >= 1.)
         kt0 = one-min_eps
      END WHERE

      kt_array(:,1,1)=ixy0/REAL(nkt,r8)
      kt_array(:,2,1)=(ixy0 + (/1,0/))/REAL(nkt,r8)
      kt_array(:,1,2)=(ixy0 + (/0,1/))/REAL(nkt,r8)
      kt_array(:,2,2)=(ixy0 + 1)/REAL(nkt,r8)
      
      res0 = plotter_evaluate(kt0(1),kt0(2),kt_array,xyw) - xy0
      xyw = kt_array
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plotter_newton_residual
c-----------------------------------------------------------------------
c     function 15. plotter_evaluate.
c     evaluate quantity u at some point (x0,y0) within a quadrilateral
c     based on its values at the corners of the quadrilateral.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION plotter_evaluate(x0,y0,xyw,uw) RESULT(u0)

      REAL(r8), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(2,2,2), INTENT(IN) :: xyw
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: uw
      REAL(r8), DIMENSION(SIZE(uw,1)) :: u0

      INTEGER :: count,info
      INTEGER, DIMENSION(4) :: ipiv
      REAL(r8), DIMENSION(4,SIZE(uw,1)) :: u
      REAL(r8), DIMENSION(4,4) :: bmat
c-----------------------------------------------------------------------
c     find the value of uw at (x0,y0).
c-----------------------------------------------------------------------
      bmat(:,1)=RESHAPE(xyw(1,:,:)*xyw(2,:,:),(/4/))
      bmat(:,2)=RESHAPE(xyw(1,:,:),(/4/))
      bmat(:,3)=RESHAPE(xyw(2,:,:),(/4/))
      bmat(:,4)=one
      DO count=1,SIZE(u,2)
         u(:,count)=RESHAPE(uw(count,:,:),(/4/))
      ENDDO
      CALL dgetrf(4,4,bmat,4,ipiv,info)
      CALL dgetrs('N',4,SIZE(u,2),bmat,4,ipiv,u,4,info)
      u0 = x0*y0*u(1,:) + x0*u(2,:) + y0*u(3,:) + u(4,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION plotter_evaluate
      END MODULE plotter_mod
