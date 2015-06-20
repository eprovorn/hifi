c-----------------------------------------------------------------------
c     file p2_condense.F.
c     performs static condensation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_condense_mod.
c     1. p2_condense_rhs.
c     2. p2_condense_jac.
c     3. p2_condense_j2petsc.
c     4. p2_condense_solve.
c     5. p2_condense_alloc.
c     6. p2_condense_dealloc.
c     7. p2_condense_mat_alloc.
c     8. p2_condense_scatter.
c     9. p2_condense_map_init.
c     10. p2_condense_map_cell.
c-----------------------------------------------------------------------
c     subprogram 0. p2_condense_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_condense_mod
      USE p2_sel_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      TYPE :: condense_type
      INTEGER, DIMENSION(:), POINTER :: ipiv
      REAL(r8), DIMENSION(:), POINTER :: r2
      REAL(r8), DIMENSION(:,:), POINTER :: jac12,jac22,jac21
      END TYPE condense_type

      TYPE :: condense_obj_type
      LOGICAL :: xperiodic,yperiodic
      LOGICAL, DIMENSION(:), POINTER :: polar_axis
      INTEGER :: s1,s2,lsize,elsize,nx,ny,nxp,nyp,nqty,np,ixmin,ixmax,
     $     iymin,iymax,jxmin,jxmax,jymin,jymax,jsize,its
      INTEGER, DIMENSION(:), POINTER :: id2iqty_i,id2ixp_i,id2iyp_i,
     $     id2iqty_g,id2ixp_g,id2iyp_g
      INTEGER, DIMENSION(:,:,:), POINTER :: idmap,elmap
      REAL(r8) :: emax,emin,condno
      TYPE(condense_type), DIMENSION(:,:), POINTER :: ct
      TYPE(edge_type), DIMENSION(:), POINTER :: edges
      VecScatter :: cjscat,crscat
      Vec :: cdu,cres,jvec
      Mat :: jac
      KSP :: ksp
      END TYPE condense_obj_type

      TYPE(condense_obj_type), POINTER :: ct_obj
      TYPE(condense_obj_type), TARGET :: ct_solver,ct_grid

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_condense_rhs.
c     applies static condensation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_rhs(ctv,rhs,cres)

      TYPE(condense_obj_type) :: ctv
      REAL(r8), DIMENSION(ctv%nqty,ctv%jxmin:ctv%jxmax,
     $     ctv%jymin:ctv%jymax), INTENT(INOUT) :: rhs
      Vec :: cres

      INTEGER :: nqty,np1,ix,iy,ixp,iyp,ixp0,iyp0,iqty,i,info,s1,s2
      INTEGER, DIMENSION(ctv%s2) :: ipiv
      REAL(r8), DIMENSION(ctv%s1) :: r1
      REAL(r8), DIMENSION(ctv%s2) :: r2
      REAL(r8), DIMENSION(ctv%s2,ctv%s2) :: jac22
c-----------------------------------------------------------------------
c     set local constants.
c-----------------------------------------------------------------------
      np1=ctv%np+1
      nqty=ctv%nqty
      s1=ctv%s1
      s2=ctv%s2
c-----------------------------------------------------------------------
c     start loop over elements.
c-----------------------------------------------------------------------
      iyp0=ctv%jymin
      DO iy=ctv%iymin+1,ctv%iymax
         ixp0=ctv%jxmin
         DO ix=ctv%ixmin+1,ctv%ixmax
            r2=0
c-----------------------------------------------------------------------
c     assemble interior rhs.
c-----------------------------------------------------------------------
            DO i=1,s2
               iqty=ctv%id2iqty_i(i)
               ixp=ctv%id2ixp_i(i)+ixp0
               iyp=ctv%id2iyp_i(i)+iyp0

               r2(i)=rhs(iqty,ixp,iyp)
            ENDDO
c-----------------------------------------------------------------------
c     compute composite rhs.
c-----------------------------------------------------------------------
            jac22 = ctv%ct(ix,iy)%jac22
            ipiv = ctv%ct(ix,iy)%ipiv

            CALL dgetrs('N',s2,1,jac22,s2,ipiv,r2,s2,info)
            r1=MATMUL(ctv%ct(ix,iy)%jac12,r2)
            ctv%ct(ix,iy)%r2=r2
c-----------------------------------------------------------------------
c     assemble interior rhs.
c-----------------------------------------------------------------------
            DO i=1,s1
               iqty=ctv%id2iqty_g(i)
               ixp=ctv%id2ixp_g(i)+ixp0
               iyp=ctv%id2iyp_g(i)+iyp0

               rhs(iqty,ixp,iyp) = rhs(iqty,ixp,iyp) - r1(i)
            ENDDO
c-----------------------------------------------------------------------
c     finish loop over elements.
c-----------------------------------------------------------------------
            ixp0=ixp0+np1
         ENDDO
         iyp0=iyp0+np1
      ENDDO
c-----------------------------------------------------------------------
c     transfer condensed rhs to petsc storage.
c-----------------------------------------------------------------------
      CALL p2_sel_set(.FALSE.,rhs,ctv%crscat,ctv%jvec,cres)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_rhs
c-----------------------------------------------------------------------
c     subprogram 2. p2_condense_jac.
c     applies static condensation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_jac(ctv,drdu,jac)

      TYPE(condense_obj_type) :: ctv
      REAL(r8), DIMENSION(ctv%elsize,ctv%nqty,ctv%jxmin:ctv%jxmax,
     $     ctv%jymin:ctv%jymax), INTENT(IN) :: drdu
      Mat :: jac

      INTEGER :: np,ix,iy,iqty,jqty,ixp,jxp,iyp,jyp,ixp0,iyp0,
     $     s1,s2,i,j,info,mm,nn
      INTEGER, DIMENSION(ctv%s2) :: ipiv
      REAL(r8), DIMENSION(ctv%s1,ctv%s1) :: jac11
      REAL(r8), DIMENSION(ctv%s1,ctv%s2) :: jac12
      REAL(r8), DIMENSION(ctv%s2,ctv%s1) :: jac21
      REAL(r8), DIMENSION(ctv%s2,ctv%s2) :: jac22
c-----------------------------------------------------------------------
c     set local constants.
c-----------------------------------------------------------------------
      np=ctv%np
      s1=ctv%s1
      s2=ctv%s2
c-----------------------------------------------------------------------
c     start loop over elements.
c-----------------------------------------------------------------------
      iyp0=ctv%jymin
      DO iy=ctv%iymin+1,ctv%iymax
         ixp0=ctv%jxmin
         DO ix=ctv%ixmin+1,ctv%ixmax
            jac12=0
            jac21=0
            jac22=0
            ipiv=0
            DO i=1,s2
               iqty=ctv%id2iqty_i(i)
               ixp=ctv%id2ixp_i(i)
               iyp=ctv%id2iyp_i(i)

               mm = ctv%elmap(iqty,ixp,iyp)
c-----------------------------------------------------------------------
c     assemble interior matrix.
c-----------------------------------------------------------------------               
               DO j=1,s2
                  nn = ctv%elmap(ctv%id2iqty_i(j),ctv%id2ixp_i(j),
     $                 ctv%id2iyp_i(j))
                  jac22(i,j)=drdu(nn,iqty,ixp+ixp0,iyp+iyp0)
               ENDDO
c-----------------------------------------------------------------------
c     assemble coupling to the edges.
c-----------------------------------------------------------------------
               DO j=1,s1
                  jqty=ctv%id2iqty_g(j)
                  jxp=ctv%id2ixp_g(j)
                  jyp=ctv%id2iyp_g(j)
                  nn = ctv%elmap(jqty,jxp,jyp)
                  
                  jac21(i,j)=drdu(nn,iqty,ixp+ixp0,iyp+iyp0)
                  jac12(j,i)=drdu(mm,jqty,jxp+ixp0,jyp+iyp0)
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     compute composite matrices.
c-----------------------------------------------------------------------
            CALL dgetrf(s2,s2,jac22,s2,ipiv,info)
            CALL dgetrs('N',s2,s1,jac22,s2,ipiv,jac21,s2,info)

            jac11=MATMUL(jac12,jac21)
            ctv%ct(ix,iy)%jac12 = jac12
            ctv%ct(ix,iy)%jac21 = jac21
            ctv%ct(ix,iy)%jac22 = jac22
            ctv%ct(ix,iy)%ipiv = ipiv
c-----------------------------------------------------------------------
c     subtract composite matrices from the edges.
c-----------------------------------------------------------------------
            DO i=1,s1
               iqty=ctv%id2iqty_g(i)
               ixp=ctv%id2ixp_g(i)
               iyp=ctv%id2iyp_g(i)
                                 
               DO j=1,s1
                  nn = ctv%elmap(ctv%id2iqty_g(j),ctv%id2ixp_g(j),
     $                 ctv%id2iyp_g(j))
                  jac11(i,j)
     $                 = drdu(nn,iqty,ixp+ixp0,iyp+iyp0) - jac11(i,j)
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     transfer condensed jacobian to petsc storage.
c-----------------------------------------------------------------------
            CALL p2_condense_j2petsc(ctv,ixp0,iyp0,jac11,jac)
c-----------------------------------------------------------------------
c     finish loop over elements.
c-----------------------------------------------------------------------
            ixp0=ixp0+np+1
         ENDDO
         iyp0=iyp0+np+1
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
      END SUBROUTINE p2_condense_jac
c-----------------------------------------------------------------------
c     subprogram 3. p2_condense_j2petsc.
c     transfers jacobian to petsc storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_j2petsc(ctv,ixp0,iyp0,jac11,jac)

      TYPE(condense_obj_type) :: ctv
      INTEGER, INTENT(IN) :: ixp0,iyp0
      REAL(r8), DIMENSION(ctv%s1,ctv%s1), INTENT(IN) :: jac11
      Mat :: jac

      INTEGER :: n,ixp,iyp,jxp,jyp,iqty,jqty,np,i,j,s1
      INTEGER, DIMENSION(1) :: idxm
      INTEGER, DIMENSION(2) :: idxn_p
      REAL(r8), DIMENSION(1,2), PARAMETER ::
     $     values_p=RESHAPE((/one,-one/),(/1,2/))
      INTEGER, DIMENSION(ctv%s1) :: idxn
      REAL(r8), DIMENSION(1,ctv%s1) :: values
c-----------------------------------------------------------------------
c     set local constants.
c-----------------------------------------------------------------------
      np=ctv%np
      s1=ctv%s1
c-----------------------------------------------------------------------
c     start loop over matrix rows.
c-----------------------------------------------------------------------
      DO i=1,s1
         iqty=ctv%id2iqty_g(i)
         ixp=ctv%id2ixp_g(i)+ixp0
         iyp=ctv%id2iyp_g(i)+iyp0
c-----------------------------------------------------------------------
c     find matrix row index for jacobian data.
c-----------------------------------------------------------------------
         idxm(1)=ctv%idmap(iqty,ixp,iyp)
c-----------------------------------------------------------------------
c     special call for polar/periodic values.
c-----------------------------------------------------------------------
         jxp=ixp
         jyp=iyp
         IF(ixp == ctv%nxp .AND.
     $        ctv%edges(3)%bc_type(iqty)=="periodic")
     $        jxp=ctv%jxmin-1
         IF((ctv%polar_axis(iqty) .AND. ixp == 0 .AND. 
     $        MOD(MOD(iyp,np+1),np) == 0) .OR. 
     $        (iyp == ctv%nyp .AND.
     $        ctv%edges(2)%bc_type(iqty)=="periodic"))jyp=ctv%jymin-1
         IF(ixp /= jxp .OR. iyp /= jyp)THEN
            idxn_p(1)=idxm(1)
            idxn_p(2)=ctv%idmap(iqty,jxp,jyp)
            CALL MatSetValues(jac,1,idxm,2,idxn_p,values_p,
     $           ADD_VALUES,ierr)
            idxm(1)=idxn_p(2)
         ENDIF
         
         n=0
         idxn=0
         values=0
c-----------------------------------------------------------------------
c     accumulate values over matrix columns.
c-----------------------------------------------------------------------
         DO j=1,s1
            jqty=ctv%id2iqty_g(j)
            jxp=ctv%id2ixp_g(j)+ixp0
            jyp=ctv%id2iyp_g(j)+iyp0
            
            IF(jac11(i,j) == 0)CYCLE
            n=n+1
            idxn(n)=ctv%idmap(jqty,jxp,jyp)
            values(1,n)=jac11(i,j)
         ENDDO
c-----------------------------------------------------------------------
c     store values.
c-----------------------------------------------------------------------
         IF(n > 0)CALL MatSetValues(jac,1,idxm,n,idxn(1:n),
     $        values(:,1:n),ADD_VALUES,ierr)
c-----------------------------------------------------------------------
c     finish loop over matrix rows.
c-----------------------------------------------------------------------
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_j2petsc
c-----------------------------------------------------------------------
c     subprogram 4. p2_condense_solve.
c     reconsitutes full solution from condensed solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_solve(ctv,cdu,duu)

      TYPE(condense_obj_type) :: ctv
      Vec :: cdu
      REAL(r8), DIMENSION(ctv%nqty,ctv%jxmin:ctv%jxmax,
     $     ctv%jymin:ctv%jymax), INTENT(OUT) :: duu
     
      INTEGER :: np,nqty,ix,iy,iqty,ixp0,iyp0,i,ixp,iyp,s1,s2
      REAL(r8), DIMENSION(ctv%s1) :: du1
      REAL(r8), DIMENSION(ctv%s2) :: du2
c-----------------------------------------------------------------------
c     set local constants.
c-----------------------------------------------------------------------
      np=ctv%np
      nqty=ctv%nqty
      s1=ctv%s1
      s2=ctv%s2
c-----------------------------------------------------------------------
c     scatter from condensed PETSc to Fortran full vector.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(ctv%cjscat,ctv%jvec,cdu,duu)
c-----------------------------------------------------------------------
c     start loop over elements to compute interior solutions.
c-----------------------------------------------------------------------
      iyp0=ctv%jymin
      DO iy=ctv%iymin+1,ctv%iymax
         ixp0=ctv%jxmin
         DO ix=ctv%ixmin+1,ctv%ixmax
            du1=0
            DO i=1,s1
               iqty=ctv%id2iqty_g(i)
               ixp=ctv%id2ixp_g(i)
               iyp=ctv%id2iyp_g(i)

               du1(i)=duu(iqty,ixp+ixp0,iyp+iyp0)
            ENDDO
c-----------------------------------------------------------------------
c     compute interior solutions and transfer to global storage.
c-----------------------------------------------------------------------
            du2=ctv%ct(ix,iy)%r2-MATMUL(ctv%ct(ix,iy)%jac21,du1)
            DO i=1,s2
               iqty=ctv%id2iqty_i(i)
               ixp=ctv%id2ixp_i(i)
               iyp=ctv%id2iyp_i(i)

               duu(iqty,ixp+ixp0,iyp+iyp0)=du2(i)
            ENDDO
c-----------------------------------------------------------------------
c     finish loop over elements.
c-----------------------------------------------------------------------
            ixp0=ixp0+np+1
         ENDDO
         iyp0=iyp0+np+1
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_solve
c-----------------------------------------------------------------------
c     subprogram 5. p2_condense_alloc.
c     allocates arrays for static condensation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_alloc(psv,ctv,rnge)

      TYPE(p2_sel_type) :: psv
      TYPE(condense_obj_type) :: ctv
      INTEGER, DIMENSION(:), INTENT(IN) :: rnge

      INTEGER :: np,nqty,ix,iy,nx,ny,lsize,xsize,ysize,iedge,nedge
c-----------------------------------------------------------------------
c     set local sizes.
c-----------------------------------------------------------------------
      np=psv%np
      nx=psv%nx
      ny=psv%ny
      nqty=SIZE(rnge)
c-----------------------------------------------------------------------
c     copy values from psv to ctv.
c-----------------------------------------------------------------------
      ctv%nqty=nqty
      ctv%nx=psv%nx
      ctv%ny=psv%ny
      ctv%np=psv%np
      ctv%nxp=psv%nxp
      ctv%nyp=psv%nyp
      ctv%ixmax=psv%ixmax
      ctv%ixmin=psv%ixmin
      ctv%iymax=psv%iymax
      ctv%iymin=psv%iymin
      ctv%jxmax=psv%jxmax
      ctv%jxmin=psv%jxmin
      ctv%jymax=psv%jymax
      ctv%jymin=psv%jymin
      ctv%xperiodic=psv%xperiodic
      ctv%yperiodic=psv%yperiodic
      ctv%jsize=(ctv%jxmax-ctv%jxmin+1)*(ctv%jymax-ctv%jymin+1)*nqty
c-----------------------------------------------------------------------
c     allocate and copy polar_axis.
c-----------------------------------------------------------------------
      ALLOCATE(ctv%polar_axis(nqty))
      ctv%polar_axis=psv%polar_axis(rnge)
c-----------------------------------------------------------------------
c     duplicate edges.
c-----------------------------------------------------------------------
      nedge=SIZE(psv%edges)
      CALL edge_alloc(nqty,nedge,ctv%edges)
      DO iedge=1,nedge
         ctv%edges(iedge)%edgenum=psv%edges(iedge)%edgenum
         ctv%edges(iedge)%edgename=psv%edges(iedge)%edgename
         ctv%edges(iedge)%static=psv%edges(iedge)%static(rnge)
         ctv%edges(iedge)%bc_type=psv%edges(iedge)%bc_type(rnge)
      ENDDO
c-----------------------------------------------------------------------
c     allocate condense_type local arrays.
c-----------------------------------------------------------------------
      ALLOCATE(ctv%ct(ctv%ixmin+1:ctv%ixmax,ctv%iymin+1:ctv%iymax))
      ctv%s1=4*np*nqty
      ctv%s2=(np-1)**2*nqty
      DO iy=ctv%iymin+1,ctv%iymax
         DO ix=ctv%ixmin+1,ctv%ixmax
            ALLOCATE(ctv%ct(ix,iy)%r2(ctv%s2))
            ALLOCATE(ctv%ct(ix,iy)%jac22(ctv%s2,ctv%s2))
            ALLOCATE(ctv%ct(ix,iy)%jac12(ctv%s1,ctv%s2))
            ALLOCATE(ctv%ct(ix,iy)%jac21(ctv%s2,ctv%s1))
            ALLOCATE(ctv%ct(ix,iy)%ipiv(ctv%s2))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     define sizes for global arrays.
c-----------------------------------------------------------------------
      xsize=(ctv%ixmax-ctv%ixmin)*np
      ysize=(ctv%iymax-ctv%iymin)*np
      IF(ctv%ixmax == nx .AND. .NOT. ctv%xperiodic)xsize=xsize+1
      IF(ctv%iymax == ny .AND. .NOT. ctv%yperiodic)ysize=ysize+1
      lsize=nqty*(xsize*ysize - (ctv%ixmax-ctv%ixmin)
     $     *(ctv%iymax-ctv%iymin)*(np-1)**2)
c-----------------------------------------------------------------------
c     create global condensed vectors.
c-----------------------------------------------------------------------
      CALL VecCreateMPI(comm,lsize,PETSC_DETERMINE,ctv%cres,ierr)
      ctv%lsize=lsize
      CALL VecDuplicate(ctv%cres,ctv%cdu,ierr)
c-----------------------------------------------------------------------
c     create global skeleton and cell-wise mappings
c-----------------------------------------------------------------------
      ALLOCATE(ctv%idmap(nqty,ctv%jxmin-1:ctv%jxmax,
     $     ctv%jymin-1:ctv%jymax))
      CALL p2_condense_map_init(ctv)

      ALLOCATE(ctv%id2iqty_i(ctv%s2),ctv%id2ixp_i(ctv%s2), 
     $     ctv%id2iyp_i(ctv%s2))
      ALLOCATE(ctv%id2iqty_g(ctv%s1),ctv%id2ixp_g(ctv%s1), 
     $     ctv%id2iyp_g(ctv%s1))
      CALL p2_condense_map_cell(ctv)

      ctv%elsize=nqty*(np+1)**2
      ALLOCATE(ctv%elmap(nqty,0:np,0:np))
      CALL p2_sel_elementmap(nqty,np,ctv%elmap)
c-----------------------------------------------------------------------
c     create scatterers and sparse matrix.
c-----------------------------------------------------------------------
      CALL VecCreateSeq(PETSC_COMM_SELF,ctv%jsize,ctv%jvec,ierr)
      CALL p2_condense_scatter(ctv,.FALSE.,ctv%cjscat)
      CALL p2_condense_scatter(ctv,.TRUE.,ctv%crscat)
      CALL p2_condense_mat_alloc(ctv)
c-----------------------------------------------------------------------
c     create linear solver.
c-----------------------------------------------------------------------
      CALL KSPCreate(comm,ctv%ksp,ierr)
      CALL KSPSetComputeSingularValues(ctv%ksp,PETSC_TRUE,ierr)
      CALL KSPSetFromOptions(ctv%ksp,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_alloc
c-----------------------------------------------------------------------
c     subprogram 6. p2_condense_dealloc.
c     deallocates arrays for static condensation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_dealloc(ctv)

      TYPE(condense_obj_type) :: ctv

      INTEGER :: ix,iy
c-----------------------------------------------------------------------
c     destroy petsc objects.
c-----------------------------------------------------------------------
      CALL VecDestroy(ctv%cres,ierr)
      CALL VecDestroy(ctv%cdu,ierr)
      CALL VecDestroy(ctv%jvec,ierr)
      CALL VecScatterDestroy(ctv%cjscat,ierr)
      CALL VecScatterDestroy(ctv%crscat,ierr)
      CALL MatDestroy(ctv%jac,ierr)
      CALL KSPDestroy(ctv%ksp,ierr)
c-----------------------------------------------------------------------
c     deallocate global objects.
c-----------------------------------------------------------------------
      CALL edge_dealloc(ctv%edges)
      DEALLOCATE(ctv%polar_axis)
      DEALLOCATE(ctv%idmap,ctv%elmap)
      DEALLOCATE(ctv%id2iqty_i,ctv%id2ixp_i,ctv%id2iyp_i)
      DEALLOCATE(ctv%id2iqty_g,ctv%id2ixp_g,ctv%id2iyp_g)
c-----------------------------------------------------------------------
c     deallocate cell-wise objects.
c-----------------------------------------------------------------------
      DO iy=ctv%iymin+1,ctv%iymax
         DO ix=ctv%ixmin+1,ctv%ixmax
            DEALLOCATE(ctv%ct(ix,iy)%r2)
            DEALLOCATE(ctv%ct(ix,iy)%jac22)
            DEALLOCATE(ctv%ct(ix,iy)%jac12)
            DEALLOCATE(ctv%ct(ix,iy)%jac21)
            DEALLOCATE(ctv%ct(ix,iy)%ipiv)
         ENDDO
      ENDDO
      DEALLOCATE(ctv%ct)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_dealloc
c-----------------------------------------------------------------------
c     subprogram 7. p2_condense_mat_alloc.
c     allocates matrix for static condensation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_mat_alloc(ctv)

      TYPE(condense_obj_type) :: ctv

      INTEGER :: np,nqty,lsize,nysize,gsize,i,ix,iy,ixp,iyp,jxp,jyp,
     $     jxpmax,jypmax
      INTEGER, DIMENSION(ctv%lsize) :: d_z,o_z
c-----------------------------------------------------------------------
c     set local sizes.
c-----------------------------------------------------------------------
      np=ctv%np
      nqty=ctv%nqty
      lsize=ctv%lsize
      nysize=ctv%iymax-ctv%iymin
      CALL VecGetSize(ctv%cres,gsize,ierr)
c-----------------------------------------------------------------------
c     prepare for allocation.
c-----------------------------------------------------------------------
      i = 0
      DO iy = ctv%iymin+1,ctv%iymax
         jypmax = np - 1
         IF(iy == ctv%ny .AND. (.NOT. ctv%yperiodic))jypmax = np
         DO ix = ctv%ixmin+1,ctv%ixmax
            jxpmax = np - 1
            IF(ix == ctv%nx .AND. (.NOT. ctv%xperiodic))jxpmax = np
            DO jyp = 0,jypmax
               iyp = (iy-1)*np+jyp
               DO jxp = 0,jxpmax
                  ixp = (ix-1)*np+jxp
                  IF(polar .AND. ixp == 0)THEN
                     IF(iyp == 0)THEN
                        d_z(i+1:i+nqty) = (2*np-1)*nysize+1
                        o_z(i+1:i+nqty) = (2*np-1)*ctv%ny 
     $                       - (np-1)*nysize
                     ELSE
                        d_z(i+1:i+nqty) = 2
                        o_z(i+1:i+nqty) = 1
                     ENDIF
                  ELSEIF(MOD(ixp,np)==0 .AND. MOD(iyp,np)==0)THEN
                     d_z(i+1:i+nqty) = 12*np - 3
                     o_z(i+1:i+nqty) = 10*np - 2
                  ELSEIF(MOD(ixp,np)==0 .OR. MOD(iyp,np)==0)THEN
                     d_z(i+1:i+nqty) = 7*np - 1
                     o_z(i+1:i+nqty) = 5*np
                  ELSE
                     CYCLE   
                  ENDIF
                  
                  i = i + nqty
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      d_z=d_z*nqty
      o_z=o_z*nqty
      d_z=MIN(d_z,lsize)
      o_z=MIN(o_z,gsize)
c-----------------------------------------------------------------------
c     create sparse matrix.
c-----------------------------------------------------------------------
      CALL MatCreate(comm,ctv%jac,ierr)
      CALL MatSetSizes(ctv%jac,lsize,lsize,PETSC_DETERMINE,
     $     PETSC_DETERMINE,ierr)
      CALL MatSetFromOptions(ctv%jac,ierr)
      CALL MatMPIAIJSetPreallocation(ctv%jac,PETSC_NULL_INTEGER,d_z,
     $     PETSC_NULL_INTEGER,o_z,ierr)
      CALL MatSetOption(ctv%jac,MAT_NEW_NONZERO_ALLOCATION_ERR,
     $     PETSC_FALSE,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_mat_alloc
c-----------------------------------------------------------------------
c     subprogram 8. p2_condense_scatter.
c     creates scatter context between full and condensed vectors.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_scatter(ctv,rhs_flag,vscatter)

      TYPE(condense_obj_type) :: ctv
      LOGICAL, INTENT(IN) :: rhs_flag
      VecScatter :: vscatter

      INTEGER :: ixp,iyp,jxp,jyp,iqty,counter,jcount
      INTEGER, DIMENSION(0:ctv%jsize-1) :: is0,is1
      IS :: ris0,ris1
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      counter=0
      jcount=-1
      is0=-1
      is1=-1
c-----------------------------------------------------------------------
c     compute index set.
c-----------------------------------------------------------------------
      DO iyp=ctv%jymin,ctv%jymax
         DO ixp=ctv%jxmin,ctv%jxmax
            DO iqty=1,ctv%nqty

               jcount=jcount+1
               IF(ctv%idmap(iqty,ixp,iyp) < 0)CYCLE

               jxp=ixp
               jyp=iyp
               IF(rhs_flag)THEN
                  IF(ixp == ctv%nxp .AND.
     $                 ctv%edges(3)%bc_type(iqty)=="periodic")
     $                 jxp=ctv%jxmin-1
                  IF((ctv%polar_axis(iqty) .AND. ixp == 0) .OR.
     $                 (iyp == ctv%nyp .AND.
     $                 ctv%edges(2)%bc_type(iqty)=="periodic"))
     $                 jyp=ctv%jymin-1
               ENDIF

               is0(counter) = ctv%idmap(iqty,jxp,jyp)
               is1(counter) = jcount
               counter = counter+1

            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     create global scatter context over processor.
c-----------------------------------------------------------------------
      CALL ISCreateGeneral(PETSC_COMM_SELF,counter,is0(0:counter-1),
     $     PETSC_COPY_VALUES,ris0,ierr)
      CALL ISCreateGeneral(PETSC_COMM_SELF,counter,is1(0:counter-1),
     $     PETSC_COPY_VALUES,ris1,ierr)
      CALL VecScatterCreate(ctv%cres,ris0,ctv%jvec,ris1,vscatter,
     $     ierr)
      CALL ISDestroy(ris0,ierr)
      CALL ISDestroy(ris1,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_scatter
c-----------------------------------------------------------------------
c     subprogram 9. p2_condense_map_init
c     create and store condensed 3d to 1d mapping.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_map_init(ctv)
      
      TYPE(condense_obj_type) :: ctv

      INTEGER :: i,iqty,ix,iy,ixp,iyp,jxp,jyp,jxpmax,jypmax, 
     $     nx,ny,np,nqty,irank,iix,iiy
      INTEGER, DIMENSION(0:mpi_size-1) :: ixmin_v,ixmax_v,iymin_v,
     $     iymax_v
      INTEGER, DIMENSION(ctv%nqty,
     $     MAX(ctv%ixmin*ctv%np-1,0):ctv%ixmax*ctv%np,
     $     MAX(ctv%iymin*ctv%np-1,0):ctv%iymax*ctv%np) :: jdmap
c-----------------------------------------------------------------------
c     initialize variables  
c-----------------------------------------------------------------------
      np=ctv%np
      nx=ctv%nx
      ny=ctv%ny
      nqty=ctv%nqty
      CALL gather(ctv%ixmin,ixmin_v)
      CALL gather(ctv%ixmax,ixmax_v)
      CALL gather(ctv%iymin,iymin_v)
      CALL gather(ctv%iymax,iymax_v)
c-----------------------------------------------------------------------
c     create 3d to 1d mapping of continuous skeleton  
c-----------------------------------------------------------------------
      i = 0
      jdmap=0
      DO irank=0,mpi_size-1
         DO iy = iymin_v(irank)+1,iymax_v(irank)
            jypmax = np - 1
            IF (iy == ny .AND. (.NOT. ctv%yperiodic))jypmax = np
            DO ix = ixmin_v(irank)+1,ixmax_v(irank)
               jxpmax = np - 1
               IF (ix == nx .AND. (.NOT. ctv%xperiodic))jxpmax = np
               DO jyp = 0,jypmax
                  iyp = (iy-1)*np+jyp
                  iiy = -1
                  IF(iyp == 0)THEN
                     iiy=MAX(ctv%iymin*np-1,0)
                  ELSEIF(iyp >= ctv%iymin*np
     $                    .AND. iyp <= ctv%iymax*np)THEN
                     iiy=iyp
                  ENDIF

                  DO jxp = 0,jxpmax
                     IF (MOD(jxp,np) > 0 .AND. MOD(jyp,np) > 0)CYCLE
                     ixp = (ix-1)*np+jxp
                     iix = -1
                     IF(ixp == 0)THEN
                        iix=MAX(ctv%ixmin*np-1,0)
                     ELSEIF(ixp >= ctv%ixmin*np
     $                       .AND. ixp <= ctv%ixmax*np)THEN
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
      ctv%idmap=-1

      DO iyp = ctv%jymin-1,ctv%jymax
         jyp = iyp/(np+1)*np + MOD(iyp,np+1)
         IF((ctv%yperiodic .AND. iyp == ctv%nyp)
     $        .OR. (iyp < ctv%jymin))jyp=MAX(ctv%iymin*np-1,0)
         DO ixp = ctv%jxmin-1,ctv%jxmax
            jxp = ixp/(np+1)*np + MOD(ixp,np+1)
            IF((ctv%xperiodic .AND. ixp == ctv%nxp)
     $           .OR. (ixp < ctv%jxmin))jxp=MAX(ctv%ixmin*np-1,0)

            IF (MOD(MOD(ixp,np+1),np) > 0 
     $           .AND. MOD(MOD(iyp,np+1),np) > 0)CYCLE
            ctv%idmap(:,ixp,iyp) = jdmap(:,jxp,jyp)

         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_map_init
c-----------------------------------------------------------------------
c     subprogram 10. p2_condense_map_cell.
c     create and store 3d to 1d mapping of single cell indices 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_condense_map_cell(ctv)

      TYPE(condense_obj_type) :: ctv

      INTEGER :: i1,i2,ixp,iyp,iqty,nqty,np
c-----------------------------------------------------------------------
c     create and store 3d to 1d mapping of single cell indices 
c-----------------------------------------------------------------------
      np=ctv%np
      nqty=ctv%nqty

      i1 = 0
      i2 = 0
      DO iyp = 0,np
         DO ixp = 0,np
c-----------------------------------------------------------------------
c     interior points
c-----------------------------------------------------------------------
            IF(MOD(ixp,np) > 0 .AND. MOD(iyp,np) > 0)THEN
               DO iqty = 1,nqty
                  i2 = i2 + 1
                  ctv%id2iqty_i(i2) = iqty
                  ctv%id2ixp_i(i2) = ixp
                  ctv%id2iyp_i(i2) = iyp
               ENDDO
c-----------------------------------------------------------------------
c     interface points
c-----------------------------------------------------------------------
            ELSE
               DO iqty = 1,nqty
                  i1 = i1 + 1
                  ctv%id2iqty_g(i1) = iqty
                  ctv%id2ixp_g(i1) = ixp
                  ctv%id2iyp_g(i1) = iyp
               ENDDO
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_condense_map_cell
      END MODULE p2_condense_mod
