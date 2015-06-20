c-----------------------------------------------------------------------
c     file p2_rj.F.
c     evaluation of 2D rhs, jacobian ans mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_rj_mod.
c     1. p2_rj_rhs.
c     2. p2_rj_jac.
c     3. p2_rj_j2petsc.
c     4. p2_rj_mass.
c     5. p2_rj_stop.
c-----------------------------------------------------------------------
c     subprogram 0. p2_rj_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_rj_mod
      USE p2_schur_mod
      USE p2_interior_mod
      USE p2_edge_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_rj_rhs.
c     computes right-hand side of differential equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_rj_rhs(t_var,crd,u,rr)
      
      REAL(r8), INTENT(IN) :: t_var
      TYPE(coord_type) :: crd
      Vec :: u
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(OUT) :: rr

      INTEGER :: iqty,nqty,nx,ny,np,ixmin,ixmax,iymin,iymax,i
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: uu
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      nqty=pt2%nqty
      ixmin=pt2%ixmin
      ixmax=pt2%ixmax
      iymin=pt2%iymin
      iymax=pt2%iymax
c-----------------------------------------------------------------------
c     get local u.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,u,uu)
c-----------------------------------------------------------------------
c     integrate rhs over block interior.
c-----------------------------------------------------------------------
      rr=0

      CALL p2_interior_rhs(t_var,crd,uu,rr)
c-----------------------------------------------------------------------
c     integrate b.c. over block edges.
c-----------------------------------------------------------------------
      DO i=1,4
         SELECT CASE(pt2%edges(pt2%edge_order(i))%edgenum)
         CASE(1)
            CALL p2_edge_rhs_lr(t_var,pt2%edges(1),
     $           uu(:,0:np,:),rr(:,0,:),crd%crd_mat(:,0:np,:))
         CASE(2)
            CALL p2_edge_rhs_tb(t_var,pt2%edges(2),
     $           uu(:,:,pt2%nyp-np:pt2%nyp),rr(:,:,pt2%nyp),
     $           crd%crd_mat(:,:,pt2%nyp-np:pt2%nyp))
         CASE(3)
            CALL p2_edge_rhs_lr(t_var,pt2%edges(3),
     $           uu(:,pt2%nxp-np:pt2%nxp,:),rr(:,pt2%nxp,:),
     $           crd%crd_mat(:,pt2%nxp-np:pt2%nxp,:))
         CASE(4)
            CALL p2_edge_rhs_tb(t_var,pt2%edges(4),
     $           uu(:,:,0:np),rr(:,:,0),crd%crd_mat(:,:,0:np))
         END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     treat grounded potentials.
c-----------------------------------------------------------------------
      IF(pt2%xperiodic .AND. pt2%yperiodic)THEN
         DO iqty=1,nqty
            IF(.NOT. pt2%ground(iqty))CYCLE
            IF(ixmin==0 .AND. iymin==0)THEN
               rr(iqty,0,0)=0
            ELSEIF(ixmax==nx .AND. iymin==0)THEN
               rr(iqty,pt2%nxp,0)=0
            ELSEIF(ixmin==0 .AND. iymax==ny)THEN
               rr(iqty,0,pt2%nyp)=0
            ELSEIF(ixmax==nx .AND. iymax==ny)THEN
               rr(iqty,pt2%nxp,pt2%nyp)=0
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_rj_rhs
c-----------------------------------------------------------------------
c     subprogram 2. p2_rj_jac.
c     computes derivatives of rhs with respect to dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_rj_jac(t_var,crd,u,drdu,mmat)

      REAL(r8), INTENT(IN) :: t_var
      TYPE(coord_type) :: crd
      Vec :: u
      REAL(r8), DIMENSION(pt2%elsize,pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(OUT) :: drdu
      REAL(r8), DIMENSION((pt2%np+1)**2,pt2%couple_msize,
     $     pt2%jxmin:pt2%jxmax,pt2%jymin:pt2%jymax), INTENT(IN), 
     $     OPTIONAL :: mmat

      INTEGER :: iqty,jqty,nqty,nx,ny,np,ixmin,ixmax,iymin,iymax,i,ii
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: uu
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      nx=pt2%nx
      ny=pt2%ny
      np=pt2%np
      nqty=pt2%nqty
      ixmin=pt2%ixmin
      ixmax=pt2%ixmax
      iymin=pt2%iymin
      iymax=pt2%iymax
c-----------------------------------------------------------------------
c     get local u.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,u,uu)
c-----------------------------------------------------------------------
c     integrate drdu over block interior.
c-----------------------------------------------------------------------
      CALL p2_interior_jac(t_var,crd,uu,drdu)
c-----------------------------------------------------------------------
c     apply r=0 boundary correction.
c-----------------------------------------------------------------------
      IF(ixmin == 0 .AND. polar)THEN
         DO iqty=1,nqty
            IF(.NOT. pt2%polar_axis(iqty))CYCLE
            DO i=1,np-1
               ii=pt2%elmap(iqty,0,i)
               drdu(ii,:,0:np,:)=0
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     integrate b.c. over block edges.
c-----------------------------------------------------------------------
      DO i=1,4
         SELECT CASE(pt2%edges(pt2%edge_order(i))%edgenum)
         CASE(1)
            CALL p2_edge_drdu_lr(t_var,pt2%edges(1),uu(:,0:np,:),
     $           drdu(:,:,0,:),crd%crd_mat(:,0:np,:))
         CASE(2)
            CALL p2_edge_drdu_tb(t_var,pt2%edges(2),
     $           uu(:,:,pt2%nyp-np:pt2%nyp),drdu(:,:,:,pt2%nyp),
     $           crd%crd_mat(:,:,pt2%nyp-np:pt2%nyp))
         CASE(3)
            CALL p2_edge_drdu_lr(t_var,pt2%edges(3),
     $           uu(:,pt2%nxp-np:pt2%nxp,:),drdu(:,:,pt2%nxp,:),
     $           crd%crd_mat(:,pt2%nxp-np:pt2%nxp,:))
         CASE(4)
            CALL p2_edge_drdu_tb(t_var,pt2%edges(4),uu(:,:,0:np),
     $           drdu(:,:,:,0),crd%crd_mat(:,:,0:np))
         END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     combine with mass matrix.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%step_type)
      CASE("theta","theta2")
         drdu=-drdu*(pt2%h*pt2%theta)
      CASE("bdf2")
         drdu=-drdu*pt2%bdfc
      CASE DEFAULT
         CALL p2_rj_stop("Cannot recognize step_type = "
     $        //TRIM(pt2%step_type)//".")
      END SELECT
      IF(PRESENT(mmat))THEN
         DO i=1,pt2%couple_msize
            jqty=pt2%couple_mlist(1,i)
            iqty=pt2%couple_mlist(2,i)
            drdu(jqty:pt2%elsize:nqty,iqty,:,:) 
     $           = drdu(jqty:pt2%elsize:nqty,iqty,:,:) + mmat(:,i,:,:)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     treat grounded potentials.
c-----------------------------------------------------------------------
      IF(pt2%xperiodic .AND. pt2%yperiodic)THEN
         DO iqty=1,nqty
            IF(.NOT. pt2%ground(iqty))CYCLE
            IF(ixmin==0 .AND. iymin==0)THEN
               drdu(:,iqty,0,0)=0
               drdu(iqty,iqty,0,0)=1
            ELSEIF(ixmax==nx .AND. iymin==0)THEN
               drdu(:,iqty,pt2%nxp,0)=0
            ELSEIF(ixmin==0 .AND. iymax==ny)THEN
               drdu(:,iqty,0,pt2%nyp)=0
            ELSEIF(ixmax==nx .AND. iymax==ny)THEN
               drdu(:,iqty,pt2%nxp,pt2%nyp)=0
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_rj_jac
c-----------------------------------------------------------------------
c     subprogram 3. p2_rj_j2petsc.
c     transfers jacobian to petsc storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_rj_j2petsc(drdu,jac)

      REAL(r8), DIMENSION(pt2%elsize,pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(IN) :: drdu
      Mat :: jac

      INTEGER :: n,nn,ixp,iyp,jxp,jyp,iqty,jqty,np,nqty,kxp,kyp
      INTEGER, DIMENSION(1) :: idxm
      INTEGER, DIMENSION(2) :: idxn_p
      REAL(r8), DIMENSION(1,2), PARAMETER :: 
     $     values_p=RESHAPE((/one,-one/),(/1,2/))
      INTEGER, DIMENSION(pt2%elsize) :: idxn
      REAL(r8), DIMENSION(1,pt2%elsize) :: values
c-----------------------------------------------------------------------
c     set local constants.
c-----------------------------------------------------------------------
      np=pt2%np
      nqty=pt2%nqty
c-----------------------------------------------------------------------
c     start loop over matrix rows.
c-----------------------------------------------------------------------
      DO iyp=pt2%jymin,pt2%jymax
         DO ixp=pt2%jxmin,pt2%jxmax
            DO iqty=1,nqty
c-----------------------------------------------------------------------
c     find matrix row index for jacobian data.
c-----------------------------------------------------------------------
               idxm(1) = pt2%idmap(iqty,ixp,iyp)
c-----------------------------------------------------------------------
c     special call for polar/periodic values.
c-----------------------------------------------------------------------
               jxp=ixp
               jyp=iyp
               IF(pt2%polar_axis(iqty) .AND. ixp == 0 .AND.
     $              MOD(MOD(iyp,np+1),np) == 0)jyp=pt2%jymin-1
               IF(ixp == pt2%nxp .AND.
     $              pt2%edges(3)%bc_type(iqty)=="periodic")
     $              jxp=pt2%jxmin-1
               IF(iyp == pt2%nyp .AND.
     $              pt2%edges(2)%bc_type(iqty)=="periodic")
     $              jyp=pt2%jymin-1
               IF(ixp /= jxp .OR. iyp /= jyp)THEN
                  idxn_p(1)=idxm(1)
                  idxn_p(2)=pt2%idmap(iqty,jxp,jyp)
                  CALL MatSetValues(jac,1,idxm,2,idxn_p,values_p,
     $                 ADD_VALUES,ierr)
                  idxm(1)=idxn_p(2)
               ENDIF

               n=0
               idxn=0
               values=0
               kxp=ixp-MOD(ixp,np+1)
               kyp=iyp-MOD(iyp,np+1)
c-----------------------------------------------------------------------
c     gather indices and values.
c-----------------------------------------------------------------------
               DO jyp=kyp,kyp+np
                  DO jxp=kxp,kxp+np
                     DO jqty=1,nqty
                        nn=pt2%elmap(jqty,jxp-kxp,jyp-kyp)
                        IF(drdu(nn,iqty,ixp,iyp) == 0)CYCLE
                        n=n+1
                        idxn(n) = pt2%idmap(jqty,jxp,jyp)
                        values(1,n) = drdu(nn,iqty,ixp,iyp)
                     ENDDO
                  ENDDO
               ENDDO
c-----------------------------------------------------------------------
c     transfer to petsc storage.
c-----------------------------------------------------------------------
               IF(n > 0)THEN
                  CALL MatSetValues(jac,1,idxm,n,idxn(1:n),
     $                 values(:,1:n),ADD_VALUES,ierr)
               ENDIF
c-----------------------------------------------------------------------
c     finish loop over matrix rows.
c-----------------------------------------------------------------------
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     assemble jacobian matrix.
c-----------------------------------------------------------------------
      CALL MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_rj_j2petsc
c-----------------------------------------------------------------------
c     subprogram 4. p2_rj_mass.
c     computes mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_rj_mass(crd,stat_qty,theta_fac,mmat)

      TYPE(coord_type) :: crd
      LOGICAL, DIMENSION(pt2%nqty), INTENT(OUT) :: stat_qty
      LOGICAL, DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(OUT) :: theta_fac
      REAL(r8), DIMENSION((pt2%np+1)**2,pt2%couple_msize,
     $     pt2%jxmin:pt2%jxmax,pt2%jymin:pt2%jymax), INTENT(OUT) :: mmat
      
      LOGICAL, DIMENSION(pt2%nqty) :: lstat_qty
      INTEGER :: iqty,nqty,np,i,ii,jj,ix,iy
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      IF(pt2%couple_msize <= 0)
     $     CALL program_stop("p2_rj_mass is called with no mass matrix")
      nqty=pt2%nqty
      np=pt2%np
c-----------------------------------------------------------------------
c     integrate mass matrix over block interior.
c-----------------------------------------------------------------------
      CALL p2_interior_mass(crd,mmat)
c-----------------------------------------------------------------------
c     apply r=0 boundary correction.
c-----------------------------------------------------------------------
      IF(pt2%ixmin == 0 .AND. polar)THEN
         DO ii=1,pt2%couple_msize
            iqty=pt2%couple_mlist(1,ii)
            IF(.NOT. pt2%polar_axis(iqty))CYCLE 
            DO i=1,np-1
               jj=pt2%elmap(iqty,0,i)/nqty+1
               mmat(jj,ii,0:np,:)=0
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     integrate b.c. over block edges.
c-----------------------------------------------------------------------
      DO i=1,4
         SELECT CASE(pt2%edges(pt2%edge_order(i))%edgenum)
         CASE(1)
            CALL p2_edge_mass_lr(pt2%edges(1),mmat(:,:,0,:),
     $           crd%crd_mat(:,0:np,:))
         CASE(2)
            CALL p2_edge_mass_tb(pt2%edges(2),mmat(:,:,:,pt2%nyp),
     $           crd%crd_mat(:,:,pt2%nyp-np:pt2%nyp))
         CASE(3)
            CALL p2_edge_mass_lr(pt2%edges(3),mmat(:,:,pt2%nxp,:),
     $           crd%crd_mat(:,pt2%nxp-np:pt2%nxp,:))
         CASE(4)
            CALL p2_edge_mass_tb(pt2%edges(4),mmat(:,:,:,0),
     $           crd%crd_mat(:,:,0:np))
         END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     set theta_fac parameter for the "theta" time-step 
c-----------------------------------------------------------------------
      lstat_qty=.TRUE.
      theta_fac=.FALSE.
      DO iy=pt2%jymin,pt2%jymax
         DO ix=pt2%jxmin,pt2%jxmax
            DO ii=1,pt2%couple_msize
               iqty = pt2%couple_mlist(2,ii)
               lstat_qty(iqty) = .FALSE.
               IF(SUM(ABS(mmat(:,ii,ix,iy))) /= 0.)
     $              theta_fac(iqty,ix,iy) = .TRUE.
            ENDDO
         ENDDO
      ENDDO
      CALL MPI_Allreduce(lstat_qty,stat_qty,nqty,MPI_LOGICAL,MPI_LAND,
     $     comm,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_rj_mass
c-----------------------------------------------------------------------
c     subprogram 5. p2_rj_stop.
c     stops execution and deallocates all objects.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_rj_stop(message)

      CHARACTER(*), INTENT(IN) :: message
c-----------------------------------------------------------------------
c     free work space and finalize PETSc.
c-----------------------------------------------------------------------
      CALL p2_sel_dealloc(solver)
      SELECT CASE(solver%solve_type)
      CASE("condense")
         CALL p2_condense_dealloc(ct_solver)
      CASE("schur")
         CALL p2_schur_dealloc(scv)
      END SELECT
      IF(.NOT. (fd_test .OR. fd_bound_test))THEN
         IF(adapt_grid)THEN
            CALL p2_sel_dealloc(grid)
            CALL job2_dealloc("beltrami")
            SELECT CASE(grid%solve_type)
            CASE("condense")
               CALL p2_condense_dealloc(ct_grid)
            END SELECT
         ENDIF
         CALL p2_ct_dealloc(coord)
      ENDIF
      CALL job2_dealloc("physics")
      CALL program_stop(message)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_rj_stop
      END MODULE p2_rj_mod
