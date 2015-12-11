c-----------------------------------------------------------------------
c     file p2_ct.F.
c     coordinate transformation between logical cartesian and 
c     arbitrary physical coordinate systems.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_ct_mod.
c     1. p2_ct_alloc.
c     2. p2_ct_dealloc.
c     3. p2_ct_transf.
c     4. p2_ct_calc.
c     5. p2_ct_default.
c-----------------------------------------------------------------------
c     subprogram 0. p2_ct_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_ct_mod
      USE p2_diagnose_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      TYPE :: coord_type
      REAL(r8), DIMENSION(:), POINTER :: x_ksi,y_ksi,x_eta,y_eta,jac_xy
      REAL(r8), DIMENSION(:,:), POINTER :: ksita
      REAL(r8), DIMENSION(:,:,:), POINTER :: crd_mat
      Vec :: crd_vec
      END TYPE coord_type

      INTEGER :: grid_step
      LOGICAL :: read_grid=.FALSE.
      TYPE(coord_type) :: coord
      
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_ct_alloc.
c     allocate objects for coordinate transformations
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_ct_alloc(crd,restart_dir,outfile_type)

      TYPE(coord_type) :: crd
      CHARACTER(*), INTENT(IN) :: restart_dir,outfile_type

      INTEGER :: nx,ny,np,jx,jy,ix,iy,info,iqty,j
      REAL(r8) :: dummy
      REAL(r8), DIMENSION(ndim,0:pt2%np,0:pt2%np) :: ui
c-----------------------------------------------------------------------
c     allocate and initialize coordinate mappings
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      ALLOCATE(crd%crd_mat(ndim,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax))
      IF(read_grid)THEN
         CALL p2_diagnose_input(pt2,restart_dir,outfile_type,"grid",
     $        grid_step,dummy,crd%crd_mat)
      ELSE
         jx=pt2%jxmin
c-----------------------------------------------------------------------
c     begin loop over finite elements.
c-----------------------------------------------------------------------
         DO ix=pt2%ixmin,pt2%ixmax-1
            jy=pt2%jymin
            DO iy=pt2%iymin,pt2%iymax-1
               CALL job2_grid("physics",nx,ny,np,ix,iy,
     $              ui(1,:,:),ui(2,:,:))
c-----------------------------------------------------------------------
c     solve for basis amplitudes.
c-----------------------------------------------------------------------
               DO iqty=1,ndim
                  DO j=1,2
                     CALL dgetrs('N',np+1,np+1,grid_bmat,np+1,grid_bpiv,
     $                    ui(iqty,:,:),np+1,info)
                     ui(iqty,:,:)=TRANSPOSE(ui(iqty,:,:))
                  ENDDO
               ENDDO
               crd%crd_mat(:,jx:jx+np,jy:jy+np)=ui
c-----------------------------------------------------------------------
c     finish loop over finite elements.
c-----------------------------------------------------------------------
               jy=jy+np+1
            ENDDO
            jx=jx+np+1
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     allocate transformation objects
c-----------------------------------------------------------------------
      ALLOCATE(crd%ksita(ndim,pt2%qsize),crd%jac_xy(pt2%qsize),
     $     crd%x_ksi(pt2%qsize),crd%y_ksi(pt2%qsize),
     $     crd%x_eta(pt2%qsize),crd%y_eta(pt2%qsize))
c-----------------------------------------------------------------------
c     store initial transformation in PETSc framework.
c-----------------------------------------------------------------------
      IF(adapt_grid)THEN
         CALL VecDuplicate(grid%global,crd%crd_vec,ierr)
         CALL p2_sel_set(.TRUE.,crd%crd_mat,grid%jscat,grid%jvec,
     $        crd%crd_vec)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_ct_alloc
c-----------------------------------------------------------------------
c     subprogram 2. p2_ct_dealloc.
c     deallocates objects for coordinate transformations
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_ct_dealloc(crd)

      TYPE(coord_type) :: crd
c-----------------------------------------------------------------------
c     deallocate objects
c-----------------------------------------------------------------------
      IF(adapt_grid)CALL VecDestroy(crd%crd_vec,ierr)
      IF(ASSOCIATED(crd%jac_xy))DEALLOCATE(crd%ksita,crd%x_ksi,
     $     crd%x_eta,crd%y_ksi,crd%y_eta,crd%crd_mat,crd%jac_xy)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_ct_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. p2_ct_transf.
c     get transformation between current and pre-fixed
c     coordinate systems
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_ct_transf(ii,crd,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $     jac_xy)

      INTEGER, INTENT(IN) :: ii
      TYPE(coord_type) :: crd
      REAL(r8), DIMENSION((pt2%nq+1)**2), INTENT(OUT) :: 
     $     ksi,eta,x_ksi,y_ksi,x_eta,y_eta,jac_xy

      INTEGER :: jj,nqq
c-----------------------------------------------------------------------
c     get coord. system transformation over the given cell (ix,iy)
c-----------------------------------------------------------------------
      nqq=(pt2%nq+1)**2
      jj=ii*nqq

      ksi=crd%ksita(1,jj+1:jj+nqq)
      eta=crd%ksita(2,jj+1:jj+nqq)
      x_ksi=crd%x_ksi(jj+1:jj+nqq)
      y_ksi=crd%y_ksi(jj+1:jj+nqq)
      x_eta=crd%x_eta(jj+1:jj+nqq)
      y_eta=crd%y_eta(jj+1:jj+nqq)
      jac_xy=crd%jac_xy(jj+1:jj+nqq)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_ct_transf
c-----------------------------------------------------------------------
c     subprogram 4. p2_ct_calc.
c     calculate coordinate transformations matrices
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_ct_calc(crd)

      TYPE(coord_type) :: crd

      INTEGER :: ix,iy,jj,np,nqq
      REAL(r8), DIMENSION(0:pt2%nq) :: quad,x,y
      REAL(r8), DIMENSION((pt2%nq+1)**2) :: jac_xy
      REAL(r8), DIMENSION(ndim,(pt2%nq+1)**2) :: kt,kt_x,kt_y
c-----------------------------------------------------------------------
c     zero out transformation matrices
c-----------------------------------------------------------------------
      np=pt2%np
      nqq=(pt2%nq+1)**2
      crd%jac_xy=0.
      crd%ksita=0.
      crd%x_ksi=0.
      crd%y_ksi=0.
      crd%x_eta=0.
      crd%y_eta=0.
c-----------------------------------------------------------------------
c     get local grid.
c-----------------------------------------------------------------------
      IF(adapt_grid)CALL p2_sel_get_local(grid%jscat,grid%jvec,
     $     crd%crd_vec,crd%crd_mat)
c-----------------------------------------------------------------------
c     select quadrature points.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%quad%quadr)
      CASE("gl0")
         quad=pt2%quad%pzero
      CASE("gll")
         quad=pt2%quad%qzero
      END SELECT
      jj = 0
c-----------------------------------------------------------------------
c     loop over cells.
c-----------------------------------------------------------------------
      DO iy=pt2%iymin,pt2%iymax-1
         DO ix=pt2%ixmin,pt2%ixmax-1
c-----------------------------------------------------------------------
c     interpolate physical positions and coordinates.
c-----------------------------------------------------------------------
            CALL p2_sel_interp_all(ix,iy,quad,crd%crd_mat
     $           (:,ix*(np+1):ix*(np+1)+np,iy*(np+1):iy*(np+1)+np),
     $           x,y,kt,kt_x,kt_y)
c-----------------------------------------------------------------------
c     calculate jacobian of the coordinate system
c-----------------------------------------------------------------------
            jac_xy = kt_x(1,:)*kt_y(2,:) - kt_x(2,:)*kt_y(1,:)
c-----------------------------------------------------------------------
c     calculate derivatives of current w.r.t. fixed coord. system
c-----------------------------------------------------------------------
            crd%ksita(:,jj+1:jj+nqq) = kt
            crd%x_ksi(jj+1:jj+nqq) = kt_y(2,:)/jac_xy
            crd%y_ksi(jj+1:jj+nqq) = -kt_x(2,:)/jac_xy
            crd%x_eta(jj+1:jj+nqq) = -kt_y(1,:)/jac_xy
            crd%y_eta(jj+1:jj+nqq) = kt_x(1,:)/jac_xy
            crd%jac_xy(jj+1:jj+nqq) = jac_xy
            jj = jj + nqq
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_ct_calc
c-----------------------------------------------------------------------
c     subprogram 5. p2_ct_default.
c     default values for coordinate transformations
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_ct_default(x,y,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $     jac_xy)
      
      REAL(r8), DIMENSION(0:pt2%nq), INTENT(IN) :: x,y
      REAL(r8), DIMENSION((pt2%nq+1)**2), INTENT(OUT) :: 
     $     x_ksi,y_ksi,x_eta,y_eta,ksi,eta,jac_xy

      INTEGER :: i,nq
c-----------------------------------------------------------------------
c     give default values to the transformations
c-----------------------------------------------------------------------
      nq=pt2%nq

      x_ksi=one
      x_eta=zero
      y_ksi=zero
      y_eta=one
      jac_xy=one

      DO i=0,(nq+1)**2-1
         ksi(i+1)=x(MOD(i,nq+1))
         eta(i+1)=y(i/(nq+1))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_ct_default
      END MODULE p2_ct_mod

