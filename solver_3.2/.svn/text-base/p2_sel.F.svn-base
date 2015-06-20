c-----------------------------------------------------------------------
c     file p2_sel.F.
c     spatial discretization of 2D PDEs with spectral element method.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_sel_mod.
c     1. p2_sel_alloc.
c     2. p2_sel_dealloc.
c     3. p2_sel_global_alloc.
c     4. p2_sel_mat_alloc.
c     5. p2_sel_proc_limits.
c     6. p2_sel_map_init.
c     7. p2_sel_scatter.
c     8. p2_sel_get_local.
c     9. p2_sel_get_global.
c     10. p2_sel_set
c     11. p2_sel_interp.
c     12. p2_sel_interp_all.
c     13. p2_sel_elementmap.
c-----------------------------------------------------------------------
c     subprogram 0. p2_sel_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_sel_mod
      USE job2_mod
      USE debug_mod
      USE HDF5
      IMPLICIT NONE

#include "finclude/petscdef.h"

      TYPE :: p2_sel_type
      LOGICAL :: xperiodic,yperiodic,curve_grid,nodal
      LOGICAL, DIMENSION(:), POINTER :: ground,polar_axis,stat_qty
      CHARACTER(8) :: step_type,solve_type
      INTEGER :: nx,ny,nxp,nyp,np,nq,nqty,ixmin,ixmax,iymin,iymax,
     $     jxmin,jxmax,jymin,jymax,nstep,jsize,istep,ksps,jacs,newts,
     $     itmax,itmax_decr,itmax_incr,lsize,elsize,qsize,restart_step,
     $     nqty_schur,couple_msize
      REAL(r8) :: griderr,t,h,h_old,theta,dtmax,dtmin,bdfc,
     $     dt_decr,dt_incr,tmax
      LOGICAL, DIMENSION(:,:,:), POINTER :: theta_fac
      INTEGER, DIMENSION(:), POINTER :: bpiv,edge_order
      INTEGER, DIMENSION(:,:), POINTER :: couple_mlist
      INTEGER, DIMENSION(:,:,:), POINTER :: idmap,elmap
      REAL(r8), DIMENSION(:), POINTER :: weight
      REAL(r8), DIMENSION(:,:), POINTER :: bmat,mass0,tinv,
     $     ppbq,pqbq,qpbq
      REAL(r8), DIMENSION(:,:,:), POINTER :: bq_amp
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: mass
      TYPE(jacobi_type) :: basis,quad
      TYPE(edge_type), DIMENSION(:), POINTER :: edges
      Vec :: global,jvec
      VecScatter :: rscat,jscat
      KSP :: ksp
      Mat :: jac
      END TYPE p2_sel_type

      LOGICAL :: parallel_write=.FALSE.,parallel_read=.FALSE.,
     $     fd_test=.FALSE.,fd_bound_test=.FALSE.,adapt_grid=.FALSE.,
     $     grid_solver=.FALSE.,polar=.FALSE.

      TYPE(p2_sel_type), POINTER :: pt2
      TYPE(p2_sel_type), TARGET :: solver,grid

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_sel_alloc.
c     allocates p2_sel_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_alloc(psv,nx,ny,nbx,np,nq,nqty,nqty_schur,nodal,
     $     xperiodic,yperiodic,quad_type,step_type,solve_type,
     $     restart_step)

      TYPE(p2_sel_type) :: psv
      INTEGER, INTENT(IN) :: nx,ny,nbx,np,nq,nqty,nqty_schur,
     $     restart_step
      LOGICAL, INTENT(IN) :: nodal,xperiodic,yperiodic
      CHARACTER(*), INTENT(IN) :: quad_type,step_type,solve_type

      INTEGER :: info,n,imm,imx,imy,np1,nq1
      INTEGER, DIMENSION(0:np) :: tpiv 
      REAL(r8), DIMENSION(0:np,0:np) :: tmat 
      TYPE(jacobi_type) :: modal
c-----------------------------------------------------------------------
c     transfer to members.
c-----------------------------------------------------------------------
      np1=np+1
      nq1=nq+1
      psv%ksps=0
      psv%jacs=0
      psv%newts=0
      psv%nstep=0
      psv%istep=restart_step
      psv%restart_step=restart_step
      psv%nx=nx
      psv%ny=ny
      psv%nxp=nx*np1-1
      psv%nyp=ny*np1-1
      psv%np=np
      psv%nq=nq
      psv%nqty=nqty
      psv%nqty_schur=nqty_schur
      psv%itmax=100
      psv%itmax_incr=3
      psv%itmax_decr=5
      psv%dt_decr=.8
      psv%dt_incr=1.2
      psv%bdfc=0
      psv%theta=.5
      psv%griderr=0
      psv%t=0
      psv%h=0
      psv%h_old=0
      psv%dtmax=1.
      psv%dtmin=0
      psv%tmax=0
c-----------------------------------------------------------------------
c     set logical variables.
c-----------------------------------------------------------------------
      psv%xperiodic=xperiodic
      psv%yperiodic=yperiodic
      psv%curve_grid=.FALSE.
      psv%step_type=step_type
      psv%solve_type=solve_type
      psv%nodal=nodal
c-----------------------------------------------------------------------
c     allocate local arrays.
c-----------------------------------------------------------------------
      ALLOCATE(psv%mass0(0:np,0:np),psv%bq_amp(0:nq,0:np,3),
     $     psv%ppbq(nq1**2,np1**2),psv%pqbq(nq1**2,np1**2),
     $     psv%qpbq(nq1**2,np1**2),psv%weight(nq1**2),
     $     psv%bmat(0:np,0:np),psv%bpiv(0:np),psv%stat_qty(nqty),
     $     psv%ground(nqty),psv%polar_axis(nqty),
     $     psv%edge_order(4))
      psv%stat_qty=.TRUE.
      psv%ground=.FALSE.
      psv%polar_axis=.FALSE.
      psv%edge_order=(/1,3,2,4/)
c-----------------------------------------------------------------------
c     allocate jacobi_types.
c-----------------------------------------------------------------------
      CALL jacobi_alloc(psv%basis,np,nodal,.FALSE.,"gll")
      CALL jacobi_alloc(psv%quad,nq,.TRUE.,.FALSE.,quad_type)
c-----------------------------------------------------------------------
c     compute and factor basis matrix.
c-----------------------------------------------------------------------
      DO n=0,np
         CALL jacobi_basis(psv%basis%qzero(n),psv%basis)
         psv%bmat(n,:)=psv%basis%pb
      ENDDO
      CALL dgetrf(np1,np1,psv%bmat,np1,psv%bpiv,info)
c-----------------------------------------------------------------------
c     compute modal mass matrix.
c-----------------------------------------------------------------------
      CALL jacobi_alloc(modal,np,.FALSE.,.FALSE.,"gll")
      CALL jacobi_product(modal,psv%quad,psv%mass0)
c-----------------------------------------------------------------------
c     compute transformation matrix from nodal to modal basis.
c-----------------------------------------------------------------------
      IF(nodal)THEN
         ALLOCATE(psv%tinv(0:np,0:np))
         DO n=0,np
            CALL jacobi_basis(psv%basis%qzero(n),modal)
            tmat(n,:)=modal%pb
         ENDDO
         psv%tinv=RESHAPE
     $        ((/(one,(zero,n=0,np),imm=0,np-1),one/),(/np1,np1/))
         CALL dgetrf(np1,np1,tmat,np1,tpiv,info)
         CALL dgetrs('N',np1,np1,tmat,np1,tpiv,psv%tinv,np1,info)
      ENDIF
      CALL jacobi_dealloc(modal)
c-----------------------------------------------------------------------
c     compute matrix of values of basis functions and their 
c     derivatives at quadrature points.
c-----------------------------------------------------------------------
      DO n=0,nq
         IF(psv%quad%quadr=="gl0")THEN
            CALL jacobi_basis(psv%quad%pzero(n),psv%basis)
         ELSEIF(psv%quad%quadr=="gll")THEN
            CALL jacobi_basis(psv%quad%qzero(n),psv%basis)
         ENDIF
         psv%bq_amp(n,:,1)=psv%basis%pb
         psv%bq_amp(n,:,2)=psv%basis%qb
         psv%bq_amp(n,:,3)=psv%basis%qb2
      ENDDO
      DO imy=0,np
         DO imx=0,np
            imm = imx + imy*np1 + 1
            psv%ppbq(:,imm) =  RESHAPE(MATMUL
     $           (RESHAPE(psv%bq_amp(:,imx,1),(/nq1,1/)),
     $           RESHAPE(psv%bq_amp(:,imy,1),(/1,nq1/))),(/nq1**2/))
            psv%pqbq(:,imm) =  RESHAPE(MATMUL
     $           (RESHAPE(psv%bq_amp(:,imx,1),(/nq1,1/)),
     $           RESHAPE(psv%bq_amp(:,imy,2),(/1,nq1/))),(/nq1**2/))
            psv%qpbq(:,imm) =  RESHAPE(MATMUL
     $           (RESHAPE(psv%bq_amp(:,imx,2),(/nq1,1/)),
     $           RESHAPE(psv%bq_amp(:,imy,1),(/1,nq1/))),(/nq1**2/))
         ENDDO
      ENDDO
      psv%qpbq=psv%qpbq*two*nx
      psv%pqbq=psv%pqbq*two*ny

      psv%weight = RESHAPE(MATMUL(RESHAPE(psv%quad%weight,(/nq1,1/)),
     $     RESHAPE(psv%quad%weight,(/1,nq1/))),(/nq1**2/))
     $     /(4._r8*nx*ny)
c-----------------------------------------------------------------------
c     create global vectors.
c-----------------------------------------------------------------------
      CALL p2_sel_global_alloc(psv,nbx)
c-----------------------------------------------------------------------
c     compute processor boundaries and set up global mapping and scatter
c-----------------------------------------------------------------------
      CALL p2_sel_proc_limits(psv)
      ALLOCATE(psv%idmap(nqty,psv%jxmin-1:psv%jxmax,
     $     psv%jymin-1:psv%jymax))
      CALL p2_sel_map_init(psv)

      psv%elsize=nqty*np1**2
      ALLOCATE(psv%elmap(nqty,0:np,0:np))
      CALL p2_sel_elementmap(nqty,np,psv%elmap)

      CALL VecCreateSeq(PETSC_COMM_SELF,psv%jsize,psv%jvec,ierr)
      CALL p2_sel_scatter(psv,.FALSE.,psv%jscat)
c-----------------------------------------------------------------------
c     allocate arrays for boundary conditions.
c-----------------------------------------------------------------------
      CALL edge_alloc(nqty,4,psv%edges)
      IF(psv%ixmin == 0 .AND. .NOT. xperiodic)psv%edges(1)%edgenum=1
      IF(psv%iymin == 0 .AND. .NOT. yperiodic)psv%edges(4)%edgenum=4
      IF(psv%ixmax == nx .AND. .NOT. xperiodic)psv%edges(3)%edgenum=3
      IF(psv%iymax == ny .AND. .NOT. yperiodic)psv%edges(2)%edgenum=2
c-----------------------------------------------------------------------
c     create theta_fac vector for theta time-stepping.
c-----------------------------------------------------------------------
      ALLOCATE(psv%theta_fac(nqty,psv%jxmin:psv%jxmax,
     $     psv%jymin:psv%jymax))
      psv%theta_fac=.FALSE.
c-----------------------------------------------------------------------
c     create full linear solver and PETSc matrix.
c-----------------------------------------------------------------------
      IF(solve_type == "full")THEN
         CALL KSPCreate(comm,psv%ksp,ierr)
         CALL KSPSetComputeSingularValues(psv%ksp,PETSC_TRUE,ierr)
         CALL KSPSetFromOptions(psv%ksp,ierr)
         CALL p2_sel_mat_alloc(psv,psv%jac)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_alloc
c-----------------------------------------------------------------------
c     subprogram 2. p2_sel_dealloc.
c     deallocates p2_sel_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_dealloc(psv)

      TYPE(p2_sel_type) :: psv
c-----------------------------------------------------------------------
c     deallocate arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(psv%mass0,psv%bq_amp,psv%bmat,psv%bpiv,
     $     psv%ppbq,psv%qpbq,psv%pqbq,psv%weight,psv%idmap,
     $     psv%elmap,psv%stat_qty,psv%ground,psv%polar_axis,
     $     psv%edge_order,psv%theta_fac)
      IF(ASSOCIATED(psv%mass))DEALLOCATE(psv%mass)
      CALL edge_dealloc(psv%edges)
      IF(ASSOCIATED(psv%couple_mlist))DEALLOCATE(psv%couple_mlist)
      IF(psv%nodal)DEALLOCATE(psv%tinv)
      CALL jacobi_dealloc(psv%basis)
      CALL jacobi_dealloc(psv%quad)
c-----------------------------------------------------------------------
c     destroy petsc objects.
c-----------------------------------------------------------------------
      CALL VecDestroy(psv%global,ierr)
      CALL VecDestroy(psv%jvec,ierr)
      CALL VecScatterDestroy(psv%jscat,ierr)
      CALL VecScatterDestroy(psv%rscat,ierr)
      IF(psv%solve_type == "full")THEN
         CALL MatDestroy(psv%jac,ierr)
         CALL KSPDestroy(psv%ksp,ierr)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. p2_sel_global_alloc.
c     allocates distributed array and associated vectors.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_global_alloc(psv,nbx)

      TYPE(p2_sel_type) :: psv
      INTEGER, INTENT(IN) :: nbx

      CHARACTER(100) :: message
      INTEGER :: nby,ip,l0,l1,last,nl,ib,lsize
      INTEGER, DIMENSION(:), ALLOCATABLE :: lx,ly
c-----------------------------------------------------------------------
c     choose number of blocks.
c-----------------------------------------------------------------------
      IF(MOD(mpi_size,nbx) /= 0)THEN
         WRITE(message,'(a,i6,a,i4)')"mpi_size = ",mpi_size,
     $        " is not a multiple of nbx = ",nbx
         CALL program_stop(message)
      ENDIF
      nby=mpi_size/nbx
      IF(psv%nx < nbx .OR. psv%ny < nby)THEN
         WRITE(message,'(a,i6,3(a,i4))')
     $        "Too many processors, mpi_size=",mpi_size,
     $        " for grid size nx,ny=",psv%nx,
     $        ",",psv%ny," given nbx=",nbx
         CALL program_stop(message)
      ENDIF
c-----------------------------------------------------------------------
c     define x sizes.
c-----------------------------------------------------------------------
      ALLOCATE(lx(nbx))
      l0=psv%np*(psv%nx/nbx)
      nl=MOD(psv%nx,nbx)
      l1=l0+psv%np*MIN(nl,1)
      IF(psv%xperiodic)THEN
         last=0
      ELSE
         last=1
      ENDIF
      lx=(/(l1,ip=1,nl),(l0,ip=1,nbx-nl-1),l0+last/)
c-----------------------------------------------------------------------
c     define y sizes.
c-----------------------------------------------------------------------
      ALLOCATE(ly(nby))
      l0=psv%np*(psv%ny/nby)
      nl=MOD(psv%ny,nby)
      l1=l0+psv%np*MIN(nl,1)
      IF(psv%yperiodic)THEN
         last=0
      ELSE
         last=1
      ENDIF
      ly=(/(l1,ip=1,nl),(l0,ip=1,nby-nl-1),l0+last/)
c-----------------------------------------------------------------------
c     create global vector and define limits.
c-----------------------------------------------------------------------
      lsize=ly(mpi_rank/nbx+1)*lx(MOD(mpi_rank,nbx)+1)*psv%nqty
      CALL VecCreateMPI(comm,lsize,PETSC_DETERMINE,psv%global,ierr)
      psv%lsize=lsize
c-----------------------------------------------------------------------
c     calculate min/max cell indices.
c-----------------------------------------------------------------------
      ib = mpi_rank/nbx
      psv%iymin=SUM(ly(1:ib))/psv%np
      psv%iymax=MIN(SUM(ly(1:ib+1))/psv%np,psv%ny)
      ib = mpi_rank - nbx*ib
      psv%ixmin=SUM(lx(1:ib))/psv%np
      psv%ixmax=MIN(SUM(lx(1:ib+1))/psv%np,psv%nx)
c-----------------------------------------------------------------------
c     deallocate objects.
c-----------------------------------------------------------------------
      DEALLOCATE(lx,ly)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_global_alloc
c-----------------------------------------------------------------------
c     subprogram 4. p2_sel_mat_alloc.
c     allocates matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_mat_alloc(psv,new_mat)

      TYPE(p2_sel_type) :: psv
      Mat :: new_mat

      INTEGER :: lsize,nysize,np,nqty,gsize
      INTEGER, DIMENSION(psv%lsize) :: d_z,o_z
c-----------------------------------------------------------------------
c     set local constants.
c-----------------------------------------------------------------------
      np=psv%np
      nqty=psv%nqty
      lsize=psv%lsize
      CALL VecGetSize(psv%global,gsize,ierr)
c-----------------------------------------------------------------------
c     allocate matrix.
c-----------------------------------------------------------------------
      d_z=(2*np+1)**2
      o_z=(2*np+1)**2
      IF(polar .AND. mpi_rank == 0)THEN
         nysize=psv%iymax-psv%iymin
         IF(psv%ixmax > 1)THEN
            d_z(1:nqty)=np**2*nysize+1
            o_z(1:nqty)=np**2*(psv%ny-nysize)
         ELSE
            d_z(1:nqty)=np*(np-1)*nysize+1
            o_z(1:nqty)=np**2*psv%ny-(d_z(1:nqty)-1)
         ENDIF
      ENDIF
      d_z=d_z*nqty
      o_z=o_z*nqty
      d_z=MIN(d_z,lsize)
      o_z=MIN(o_z,gsize)
c-----------------------------------------------------------------------
c     create sparse matrix.
c-----------------------------------------------------------------------
      CALL MatCreate(comm,new_mat,ierr)
      CALL MatSetSizes(new_mat,lsize,lsize,PETSC_DETERMINE,
     $     PETSC_DETERMINE,ierr)
      CALL MatSetFromOptions(new_mat,ierr)
      CALL MatMPIAIJSetPreallocation(new_mat,PETSC_NULL_INTEGER,d_z,
     $     PETSC_NULL_INTEGER,o_z,ierr)
      CALL MatSetOption(new_mat,MAT_NEW_NONZERO_ALLOCATION_ERR,
     $     PETSC_FALSE,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_mat_alloc
c-----------------------------------------------------------------------
c     subprogram 5. p2_sel_proc_limits.
c     computes processor boundaries.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_proc_limits(psv)

      TYPE(p2_sel_type) :: psv

      INTEGER :: np
c-----------------------------------------------------------------------
c     define scalar boundaries.
c-----------------------------------------------------------------------
      np=psv%np
c-----------------------------------------------------------------------
c     define y limits.
c-----------------------------------------------------------------------
      psv%jymin=psv%iymin*(np+1)
      psv%jymax=psv%iymax*(np+1)-1
c-----------------------------------------------------------------------
c     define x limits.
c-----------------------------------------------------------------------
      psv%jxmin=psv%ixmin*(np+1)
      psv%jxmax=psv%ixmax*(np+1)-1
c-----------------------------------------------------------------------
c     compute local size of Fortran vectors.
c-----------------------------------------------------------------------
      psv%jsize=(psv%jxmax-psv%jxmin+1)*(psv%jymax-psv%jymin+1)*psv%nqty
c-----------------------------------------------------------------------
c     compute local number of quadrature points.
c-----------------------------------------------------------------------
      psv%qsize=(psv%ixmax-psv%ixmin)*(psv%iymax-psv%iymin)
     $     *(psv%nq+1)**2
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_proc_limits
c-----------------------------------------------------------------------
c     subprogram 6. p2_sel_map_init.
c     create and store global 3d to 1d mapping.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_map_init(psv)
      
      TYPE(p2_sel_type) :: psv

      INTEGER :: i,iqty,ix,iy,ixp,iyp,jxp,jyp,jxpmax,jypmax, 
     $     nx,ny,np,nqty,irank,iix,iiy
      INTEGER, DIMENSION(0:mpi_size-1) :: ixmin_v,ixmax_v,iymin_v,
     $     iymax_v
      INTEGER, DIMENSION(psv%nqty,
     $     MAX(psv%ixmin*psv%np-1,0):psv%ixmax*psv%np,
     $     MAX(psv%iymin*psv%np-1,0):psv%iymax*psv%np) :: jdmap
c-----------------------------------------------------------------------
c     initialize variables  
c-----------------------------------------------------------------------
      np=psv%np
      nx=psv%nx
      ny=psv%ny
      nqty=psv%nqty
      CALL gather(psv%ixmin,ixmin_v)
      CALL gather(psv%ixmax,ixmax_v)
      CALL gather(psv%iymin,iymin_v)
      CALL gather(psv%iymax,iymax_v)
c-----------------------------------------------------------------------
c     create 3d to 1d mapping of continuous full vector
c-----------------------------------------------------------------------
      i = 0
      jdmap=0
      DO irank=0,mpi_size-1
         DO iy = iymin_v(irank)+1,iymax_v(irank)
            jypmax = np - 1
            IF (iy == ny .AND. (.NOT. psv%yperiodic))jypmax = np
            DO ix = ixmin_v(irank)+1,ixmax_v(irank)
               jxpmax = np - 1
               IF (ix == nx .AND. (.NOT. psv%xperiodic))jxpmax = np

               DO jyp = 0,jypmax
                  iyp = (iy-1)*np+jyp
                  iiy = -1
                  IF(iyp == 0)THEN
                     iiy=MAX(psv%iymin*np-1,0)
                  ELSEIF(iyp >= psv%iymin*np
     $                    .AND. iyp <= psv%iymax*np)THEN
                     iiy=iyp
                  ENDIF

                  DO jxp = 0,jxpmax
                     ixp = (ix-1)*np+jxp
                     iix = -1
                     IF(ixp == 0)THEN
                        iix=MAX(psv%ixmin*np-1,0)
                     ELSEIF(ixp >= psv%ixmin*np
     $                       .AND. ixp <= psv%ixmax*np)THEN
                        iix=ixp
                     ENDIF

                     DO iqty = 1,nqty
                        IF(iix >= 0 .AND. iiy >= 0)
     $                       jdmap(iqty,iix,iiy) = i
                        i = i + 1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     transfer to the discontinuous representation
c-----------------------------------------------------------------------
      psv%idmap=0
      DO iyp = psv%jymin-1,psv%jymax
         jyp = iyp/(np+1)*np + MOD(iyp,np+1)
         IF((psv%yperiodic .AND. iyp == psv%nyp)
     $        .OR. (iyp < psv%jymin))jyp=MAX(psv%iymin*np-1,0)
         DO ixp =  psv%jxmin-1,psv%jxmax
            jxp = ixp/(np+1)*np + MOD(ixp,np+1)
            IF((psv%xperiodic .AND. ixp == psv%nxp)
     $           .OR. (ixp < psv%jxmin))jxp=MAX(psv%ixmin*np-1,0)
            psv%idmap(:,ixp,iyp) = jdmap(:,jxp,jyp)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_map_init
c-----------------------------------------------------------------------
c     subprogram 7. p2_sel_scatter.
c     creates scatter context for full continuous vectors.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_scatter(psv,rhs_flag,vscatter)

      TYPE(p2_sel_type) :: psv
      LOGICAL, INTENT(IN) :: rhs_flag
      VecScatter :: vscatter

      INTEGER :: ixp,iyp,jxp,jyp,iqty,counter
      INTEGER, DIMENSION(0:psv%jsize-1) :: is0
      IS :: ris0,ris1
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      counter=0
      is0=-1
c-----------------------------------------------------------------------
c     compute index set.
c-----------------------------------------------------------------------
      DO iyp=psv%jymin,psv%jymax
         DO ixp=psv%jxmin,psv%jxmax
            DO iqty=1,psv%nqty

               jxp=ixp
               jyp=iyp
               IF(rhs_flag)THEN
                  IF(ixp == psv%nxp .AND. 
     $                 psv%edges(3)%bc_type(iqty)=="periodic")
     $                 jxp=psv%jxmin-1
                  IF((psv%polar_axis(iqty) .AND. ixp == 0)
     $                 .OR. (iyp == psv%nyp .AND. 
     $                 psv%edges(2)%bc_type(iqty)=="periodic"))
     $                 jyp=psv%jymin-1
               ENDIF
               is0(counter)=psv%idmap(iqty,jxp,jyp)
               counter=counter+1

            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     create global scatter context over processor.
c-----------------------------------------------------------------------
      CALL ISCreateGeneral(PETSC_COMM_SELF,psv%jsize,is0,
     $     PETSC_COPY_VALUES,ris0,ierr)
      CALL ISCreateStride(PETSC_COMM_SELF,psv%jsize,0,1,ris1,ierr)
      CALL VecScatterCreate(psv%global,ris0,psv%jvec,ris1,vscatter,ierr)
      CALL ISDestroy(ris0,ierr)
      CALL ISDestroy(ris1,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_scatter
c-----------------------------------------------------------------------
c     subprogram 8. p2_sel_get_local.
c     gets a local copy of the vector.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_get_local(jscat,jvec,global,uu)

      VecScatter :: jscat
      Vec :: jvec,global
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: uu

      INTEGER :: i,jsize
      INTEGER, DIMENSION(0:SIZE(uu)-1) :: idx
c-----------------------------------------------------------------------
c     scatter.
c-----------------------------------------------------------------------
      jsize=SIZE(uu)
      CALL VecScatterBegin(jscat,global,jvec,
     $     INSERT_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterEnd(jscat,global,jvec,
     $     INSERT_VALUES,SCATTER_FORWARD,ierr)
c-----------------------------------------------------------------------
c     transfer to Fortran array.
c-----------------------------------------------------------------------
      uu=0
      idx=(/(i,i=0,jsize-1)/)
      CALL VecGetValues(jvec,jsize,idx,uu,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_get_local
c-----------------------------------------------------------------------
c     subprogram 9. p2_sel_get_global.
c     get a global copy of a PETSc vector into Fortran array.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_get_global(psv,global,uglobal)

      TYPE(p2_sel_type) :: psv
      Vec :: global
      REAL(r8), DIMENSION(psv%nqty,0:psv%nxp,0:psv%nyp), INTENT(OUT) 
     $     :: uglobal

      INTEGER :: np,jx,jy
      REAL(r8), DIMENSION(psv%nqty,0:psv%np,0:psv%np) :: ulocal
      REAL(r8), DIMENSION(psv%nqty,psv%jxmin:psv%jxmax,
     $     psv%jymin:psv%jymax) :: u0
c-----------------------------------------------------------------------
c     get local solution.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(psv%jscat,psv%jvec,global,u0)
c-----------------------------------------------------------------------
c     create global version of solution vector.
c-----------------------------------------------------------------------
      uglobal=0
      np=psv%np
      DO jy=0,psv%nyp,np+1
         DO jx=0,psv%nxp,np+1
            ulocal=0
            IF(jx >= psv%jxmin .AND. jx < psv%jxmax 
     $           .AND. jy >= psv%jymin .AND. jy < psv%jymax)
     $           ulocal = u0(:,jx:jx+np,jy:jy+np)
            CALL MPI_Allreduce(ulocal,uglobal(:,jx:jx+np,jy:jy+np),
     $           psv%nqty*(np+1)**2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,
     $           ierr)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_get_global
c-----------------------------------------------------------------------
c     subprogram 10. p2_sel_set.
c     store data from Fortran array into PETSc vector.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_set(insert,uu,jscat,jvec,global)

      LOGICAL, INTENT(IN) :: insert
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: uu
      VecScatter :: jscat
      Vec :: jvec,global

      InsertMode :: flag
      INTEGER :: i,jsize
      INTEGER, DIMENSION(0:SIZE(uu)-1) :: idx
c-----------------------------------------------------------------------
c     insert 3D array into global vector.
c-----------------------------------------------------------------------
      jsize=SIZE(uu)
      idx=(/(i,i=0,jsize-1)/)
      CALL VecSet(jvec,zero,ierr)
      CALL VecSetValues(jvec,jsize,idx,uu,INSERT_VALUES,ierr)
      CALL VecAssemblyBegin(jvec,ierr)
      CALL VecAssemblyEnd(jvec,ierr)
c-----------------------------------------------------------------------
c     scatter.
c-----------------------------------------------------------------------
      IF(insert)THEN
         flag=INSERT_VALUES
      ELSE
         flag=ADD_VALUES
      ENDIF
      CALL VecSet(global,zero,ierr)
      CALL VecScatterBegin(jscat,jvec,global,flag,SCATTER_REVERSE,ierr)
      CALL VecScatterEnd(jscat,jvec,global,flag,SCATTER_REVERSE,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_set
c-----------------------------------------------------------------------
c     subprogram 11. p2_sel_interp.
c     interpolate solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_interp(x,y,uu,u)

      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(r8), DIMENSION(:), INTENT(OUT) :: u

      INTEGER :: jx,jy,inx,iny,np
      REAL(r8), DIMENSION(0:pt2%np) :: px,py,qx,qy
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      np=pt2%np
c-----------------------------------------------------------------------
c     find values of jacobi polynomials at interpolatory points.
c-----------------------------------------------------------------------
      CALL jacobi_interp(pt2%nx,pt2%ny,x,y,jx,jy,pt2%basis,
     $     px,py,qx,qy)
      jx=(jx-1)*(np+1)
      jy=(jy-1)*(np+1)
c-----------------------------------------------------------------------
c     interpolate solutions.
c-----------------------------------------------------------------------
      u=0
      DO iny=0,np
         DO inx=0,np
            u=u+uu(:,jx+inx,jy+iny)*px(inx)*py(iny)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_interp
c-----------------------------------------------------------------------
c     subprogram 12. p2_sel_interp_all.
c     interpolate solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_interp_all(ix,iy,quad,uu,x,y,u,ux,uy)

      INTEGER, INTENT(IN) :: ix,iy
      REAL(r8), DIMENSION(0:), INTENT(IN) :: quad
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(r8), DIMENSION(0:), INTENT(OUT) :: x,y
      REAL(r8), DIMENSION(:,0:), INTENT(OUT) :: u,ux,uy

      INTEGER :: imx,jx,imy,jy,np,nq,jjy
      REAL(r8) :: dxf,dyf
      REAL(r8), DIMENSION(0:pt2%nq,0:pt2%np,2) :: bq_amp
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      np=pt2%np
      nq=pt2%nq
      bq_amp=pt2%bq_amp(:,:,1:2)
c-----------------------------------------------------------------------
c     compute x,y positions and increments.
c-----------------------------------------------------------------------
      dxf = two*pt2%nx
      dyf = two*pt2%ny
      x = (REAL(ix,r8) + 0.5*(one+quad))/REAL(pt2%nx,r8)
      y = (REAL(iy,r8) + 0.5*(one+quad))/REAL(pt2%ny,r8)
c-----------------------------------------------------------------------
c     interpolate solutions and derivatives.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
      DO imy=0,np
         DO imx=0,np
            DO jy=0,nq
               jjy = jy*(nq+1)
               DO jx=0,nq
                  u(:,jx+jjy)=u(:,jx+jjy) + uu(:,imx,imy)
     $                 *bq_amp(jx,imx,1)*bq_amp(jy,imy,1)
                  ux(:,jx+jjy)=ux(:,jx+jjy) + uu(:,imx,imy)
     $                 *bq_amp(jx,imx,2)*bq_amp(jy,imy,1)*dxf
                  uy(:,jx+jjy)=uy(:,jx+jjy) + uu(:,imx,imy)
     $                 *bq_amp(jx,imx,1)*bq_amp(jy,imy,2)*dyf
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_interp_all
c-----------------------------------------------------------------------
c     subprogram 13. p2_sel_elementmap
c     create and store 1d to 3d mapping of element degrees of freedom  
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_sel_elementmap(nqty,np,elmap)
      
      INTEGER, INTENT(IN) :: nqty,np
      INTEGER, DIMENSION(nqty,0:np,0:np), INTENT(OUT) :: elmap

      INTEGER :: iqty,ix,iy,iel
c-----------------------------------------------------------------------
c     initialize variables  
c-----------------------------------------------------------------------
      iel=0
      elmap=0
c-----------------------------------------------------------------------
c     1d to 3d mapping of element layout
c-----------------------------------------------------------------------
      DO iy = 0,np
         DO ix = 0,np
            DO iqty = 1,nqty
               iel = iel+1
               elmap(iqty,ix,iy)=iel
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_sel_elementmap
      END MODULE p2_sel_mod
