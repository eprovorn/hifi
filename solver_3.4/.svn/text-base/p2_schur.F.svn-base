c-----------------------------------------------------------------------
c     file p2_schur.F.
c     physics-based preconditioning.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_schur_mod.
c     1. p2_schur_alloc.
c     2. p2_schur_dealloc.
c     3. p2_schur_mat.
c     4. p2_schur_schur.
c     5. p2_schur_matmul.
c     6. p2_schur_solve.
c     7. p2_schur_diagnose.
c     8. p2_schur_assemble.
c     9. p2_schur_diagnose_alloc.
c     10. p2_schur_diagnose_dealloc.
c     11. p2_schur_map_init.
c-----------------------------------------------------------------------
c     subprogram 0. p2_schur_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_schur_mod
      USE p2_ct_mod
      USE p2_condense_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      TYPE :: schur_diagnose_type
      INTEGER :: size
      INTEGER, DIMENSION(:), POINTER :: index
      Vec :: local,global
      VecScatter :: scat
      END TYPE schur_diagnose_type

      TYPE :: schur_type
      LOGICAL :: xperiodic,yperiodic
      INTEGER :: nqty,nqty1,nqty2,elsize1,elsize2,np,nq,
     $     ixmin,ixmax,iymin,iymax,jxmin,jxmax,jymin,jymax,nqty_schur,
     $     lsize,nx,ny,nxp,nyp
      INTEGER, DIMENSION(:), POINTER :: rn1,rn2
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: mat12,mat21,mat22,
     $     mat110,mat220
      REAL(r8), DIMENSION(:,:,:), POINTER :: r1,r2,s1,s2,u1,u2,t1,t2
      TYPE(condense_obj_type) :: m11c,schurc
      TYPE(schur_diagnose_type) :: sdt1,sdt2
      END TYPE schur_type

      INTEGER, PRIVATE :: schur_solve_00,schur_solve_s1,
     $     schur_solve_u2,schur_solve_u1,schur_solve_uu
      CHARACTER(8) :: schur_solve_type="condense",
     $     mass_solve_type="condense"
      LOGICAL :: schur_diagnose=.FALSE.
      TYPE(schur_type) :: scv

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_schur_alloc.
c     alocates storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_alloc(psv,scv)

      TYPE(p2_sel_type) :: psv
      TYPE(schur_type), INTENT(OUT) :: scv

      INTEGER :: iqty
      INTEGER, DIMENSION(:), ALLOCATABLE :: rnge
c$$$c-----------------------------------------------------------------------
c$$$c     register events, schur_solve.
c$$$c-----------------------------------------------------------------------
c$$$      CALL PetscLogEventRegister(schur_solve_00,"schur_solve_00",0,ierr)
c$$$      CALL PetscLogEventRegister(schur_solve_s1,"schur_solve_s1",0,ierr)
c$$$      CALL PetscLogEventRegister(schur_solve_u2,"schur_solve_u2",0,ierr)
c$$$      CALL PetscLogEventRegister(schur_solve_u1,"schur_solve_u1",0,ierr)
c$$$      CALL PetscLogEventRegister(schur_solve_uu,"schur_solve_uu",0,ierr)
c-----------------------------------------------------------------------
c     abort or set sizes.
c----------------------------------------------------------------------- 
      IF(psv%nqty_schur == 0)THEN
         CALL program_stop("This job type does not have a Schur"
     $        //" complement defined.")
      ELSE
         scv%nqty=psv%nqty
         scv%nqty1=psv%nqty-psv%nqty_schur
         scv%nqty2=psv%nqty_schur
      ENDIF
c-----------------------------------------------------------------------
c     set other scalar integers.
c-----------------------------------------------------------------------
      scv%np=psv%np
      scv%nq=psv%nq
      scv%nx=psv%nx
      scv%ny=psv%ny
      scv%nxp=psv%nxp
      scv%nyp=psv%nyp
      scv%ixmin=psv%ixmin
      scv%ixmax=psv%ixmax
      scv%iymin=psv%iymin
      scv%iymax=psv%iymax
      scv%jxmin=psv%jxmin
      scv%jxmax=psv%jxmax
      scv%jymin=psv%jymin
      scv%jymax=psv%jymax
      scv%lsize=psv%lsize
      scv%elsize1=scv%nqty1*(scv%np+1)**2
      scv%elsize2=scv%nqty2*(scv%np+1)**2
      scv%xperiodic=psv%xperiodic
      scv%yperiodic=psv%yperiodic
c-----------------------------------------------------------------------
c     allocate native vectors.
c-----------------------------------------------------------------------
      ALLOCATE(
     $     scv%r1(scv%nqty1,scv%jxmin:scv%jxmax,scv%jymin:scv%jymax),
     $     scv%s1(scv%nqty1,scv%jxmin:scv%jxmax,scv%jymin:scv%jymax),
     $     scv%u1(scv%nqty1,scv%jxmin:scv%jxmax,scv%jymin:scv%jymax),
     $     scv%t1(scv%nqty1,scv%jxmin:scv%jxmax,scv%jymin:scv%jymax),
     $     scv%t2(scv%nqty1,scv%jxmin:scv%jxmax,scv%jymin:scv%jymax),
     $     scv%r2(scv%nqty2,scv%jxmin:scv%jxmax,scv%jymin:scv%jymax),
     $     scv%s2(scv%nqty2,scv%jxmin:scv%jxmax,scv%jymin:scv%jymax),
     $     scv%u2(scv%nqty2,scv%jxmin:scv%jxmax,scv%jymin:scv%jymax))
c-----------------------------------------------------------------------
c     allocate native matrices.
c-----------------------------------------------------------------------
      ALLOCATE(scv%mat12(scv%elsize2,scv%nqty1,scv%jxmin:scv%jxmax,
     $     scv%jymin:scv%jymax))
      ALLOCATE(scv%mat21(scv%elsize1,scv%nqty2,scv%jxmin:scv%jxmax,
     $     scv%jymin:scv%jymax))
      ALLOCATE(scv%mat22(scv%elsize2,scv%nqty2,scv%jxmin:scv%jxmax,
     $     scv%jymin:scv%jymax))
c-----------------------------------------------------------------------
c     allocate diagonal block matrices for diagnostics.
c-----------------------------------------------------------------------
      IF(schur_diagnose)THEN
         ALLOCATE(scv%mat110(scv%elsize1,scv%nqty1,
     $        scv%jxmin:scv%jxmax,scv%jymin:scv%jymax))
         ALLOCATE(scv%mat220(scv%elsize2,scv%nqty2,
     $        scv%jxmin:scv%jxmax,scv%jymin:scv%jymax))
      ENDIF
c-----------------------------------------------------------------------
c     allocate mass objects.
c-----------------------------------------------------------------------
      SELECT CASE(mass_solve_type)
      CASE("condense")
         ALLOCATE(rnge(scv%nqty1))
         rnge=(/(iqty,iqty=1,scv%nqty1)/)
         CALL p2_condense_alloc(psv,scv%m11c,rnge)
         DEALLOCATE(rnge)
      CASE DEFAULT
         CALL program_stop
     $        ("p2_schur_alloc: cannot recognize mass_solve_type = "
     $        //TRIM(mass_solve_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     allocate schur objects.
c-----------------------------------------------------------------------
      SELECT CASE(schur_solve_type)
      CASE("condense")
         ALLOCATE(rnge(scv%nqty2))
         rnge=(/(iqty,iqty=scv%nqty1+1,scv%nqty)/)
         CALL p2_condense_alloc(psv,scv%schurc,rnge)
         DEALLOCATE(rnge)
      CASE DEFAULT
         CALL program_stop
     $        ("p2_schur_alloc: cannot recognize schur_solve_type = "
     $        //TRIM(schur_solve_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     allocate and define ranges.
c-----------------------------------------------------------------------
      ALLOCATE(scv%rn1(scv%nqty1),scv%rn2(scv%nqty2))
      scv%rn1=(/(iqty,iqty=1,scv%nqty1)/)
      scv%rn2=(/(iqty,iqty=scv%nqty1+1,scv%nqty)/)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_alloc
c-----------------------------------------------------------------------
c     subprogram 2. p2_schur_dealloc.
c     dealocates storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_dealloc(scv)

      TYPE(schur_type), INTENT(INOUT) :: scv

      INTEGER :: iqty
c-----------------------------------------------------------------------
c     deallocate mass objects.
c-----------------------------------------------------------------------
      SELECT CASE(mass_solve_type)
      CASE("condense")
         CALL p2_condense_dealloc(scv%m11c)
      END SELECT
c-----------------------------------------------------------------------
c     deallocate schur objects.
c-----------------------------------------------------------------------
      SELECT CASE(schur_solve_type)
      CASE("condense")
         CALL p2_condense_dealloc(scv%schurc)
      END SELECT
c-----------------------------------------------------------------------
c     deallocate native arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(scv%r1,scv%s1,scv%u1,scv%r2,scv%s2,scv%u2,scv%t1,
     $     scv%t2,scv%mat12,scv%mat21,scv%mat22,scv%rn1,scv%rn2)
c-----------------------------------------------------------------------
c     deallocate diagonal block matrices for diagnostics.
c-----------------------------------------------------------------------
      IF(schur_diagnose) DEALLOCATE(scv%mat110,scv%mat220)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. p2_schur_mat.
c     assembles matrix blocks.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_mat(psv,scv,crd,matrix,u)

      TYPE(p2_sel_type), INTENT(IN) :: psv
      TYPE(schur_type), INTENT(INOUT) :: scv
      TYPE(coord_type) :: crd
      REAL(r8), DIMENSION(:,:,:,:), INTENT(IN) :: matrix
      Vec, INTENT(IN) :: u

      INTEGER :: iqty
      REAL(r8), DIMENSION(scv%nqty,scv%jxmin:scv%jxmax,
     $     scv%jymin:scv%jymax) :: uu
      REAL(r8), DIMENSION(scv%elsize1,scv%nqty1,
     $     scv%jxmin:scv%jxmax,scv%jymin:scv%jymax) :: mat11
c-----------------------------------------------------------------------
c     construct approximate schur complement.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(psv%jscat,psv%jvec,u,uu)
      CALL p2_schur_schur(psv,scv,crd,uu)
c-----------------------------------------------------------------------
c     set up matrix blocks.
c-----------------------------------------------------------------------
      DO iqty=1,SIZE(scv%rn1)
         mat11(iqty:scv%elsize1:scv%nqty1,:,:,:)
     $        = matrix(scv%rn1(iqty):psv%elsize:scv%nqty,scv%rn1,:,:)
         scv%mat21(iqty:scv%elsize1:scv%nqty1,:,:,:)
     $        = matrix(scv%rn1(iqty):psv%elsize:scv%nqty,scv%rn2,:,:)
      ENDDO
      DO iqty=1,SIZE(scv%rn2)
         scv%mat12(iqty:scv%elsize2:scv%nqty2,:,:,:)
     $        = matrix(scv%rn2(iqty):psv%elsize:scv%nqty,scv%rn1,:,:)
         scv%mat22(iqty:scv%elsize2:scv%nqty2,:,:,:)
     $        = matrix(scv%rn2(iqty):psv%elsize:scv%nqty,scv%rn2,:,:)
     $        - scv%mat22(iqty:scv%elsize2:scv%nqty2,:,:,:)
      ENDDO
c-----------------------------------------------------------------------
c     copy diagonal block matrices for diagnostics.
c-----------------------------------------------------------------------
      IF(schur_diagnose)THEN
         scv%mat110=mat11
         scv%mat220=scv%mat22
      ENDIF
c-----------------------------------------------------------------------
c     set up mass objects.
c-----------------------------------------------------------------------
      SELECT CASE(mass_solve_type)
      CASE("condense")
         CALL p2_condense_jac(scv%m11c,mat11,scv%m11c%jac)
         CALL KSPSetOperators(scv%m11c%ksp,scv%m11c%jac,scv%m11c%jac,
     $        DIFFERENT_NONZERO_PATTERN,ierr)
      END SELECT
c-----------------------------------------------------------------------
c     set up schur objects.
c-----------------------------------------------------------------------
      SELECT CASE(schur_solve_type)
      CASE("condense")
         CALL p2_condense_jac(scv%schurc,scv%mat22,scv%schurc%jac)
         CALL KSPSetOperators(scv%schurc%ksp,scv%schurc%jac,
     $        scv%schurc%jac,DIFFERENT_NONZERO_PATTERN,ierr)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_mat
c-----------------------------------------------------------------------
c     subprogram 4. p2_schur_schur.
c     constructs Schur complement matrix block. 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_schur(psv,scv,crd,uu)

      TYPE(p2_sel_type), INTENT(IN) :: psv
      TYPE(schur_type), INTENT(INOUT) :: scv
      TYPE(coord_type) :: crd
      REAL(r8), DIMENSION(psv%nqty,psv%jxmin:psv%jxmax,
     $     psv%jymin:psv%jymax), INTENT(IN) :: uu

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      CHARACTER(80) :: filename,format1,format2
      INTEGER :: imm,mm,inn,nn,iqty,jqty,nqty,ix,iy,np,ii
      REAL(r8), DIMENSION(0:psv%nq) :: quad,x,y
      REAL(r8), DIMENSION((psv%nq+1)**2) :: ksi,eta,x_ksi,y_ksi,
     $     x_eta,y_eta,jac_xy
      REAL(r8), DIMENSION(psv%nqty,(psv%nq+1)**2) :: u,ux,uy,uk,ue
      REAL(r8), DIMENSION(scv%elsize2,scv%elsize2) :: mat_1d
      REAL(r8), DIMENSION((psv%nq+1)**2,scv%nqty2,scv%nqty2) ::
     $     fk_u,fk_uk,fk_ue,fe_u,fe_uk,fe_ue,s_u,s_uk,s_ue
      REAL(r8), DIMENSION((psv%nq+1)**2,(psv%np+1)**2) :: 
     $     kap1,kap2,ppbq,pqbq,qpbq
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT('(/3x,"inx",3x,"iny",4x,"jx",4x,"jy",2x,"iqty",'
     $     i1,'(4x,"jqty",i1,2x)/)')
 20   FORMAT('(4i6,1p,',i1,'e11.3)')
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      np=psv%np
      nqty=scv%nqty2
      ppbq=psv%ppbq
      pqbq=psv%pqbq
      qpbq=psv%qpbq
c-----------------------------------------------------------------------
c     select quadrature points.
c-----------------------------------------------------------------------
      SELECT CASE(psv%quad%quadr)
      CASE("gl0")
         quad=psv%quad%pzero
      CASE("gll")
         quad=psv%quad%qzero
      END SELECT
c-----------------------------------------------------------------------
c     start loop over cells.
c-----------------------------------------------------------------------
      ii=0
      scv%mat22=0
      DO iy=psv%iymin,psv%iymax-1
         DO ix=psv%ixmin,psv%ixmax-1
c-----------------------------------------------------------------------
c     interpolate positions.
c-----------------------------------------------------------------------
            CALL p2_sel_interp_all(ix,iy,quad,
     $           uu(:,ix*(np+1):ix*(np+1)+np,iy*(np+1):iy*(np+1)+np),
     $           x,y,u,ux,uy)
            IF(.NOT. psv%curve_grid)THEN
               CALL p2_ct_default(x,y,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $              jac_xy)
            ELSE
               CALL p2_ct_transf(ii,crd,ksi,eta,x_ksi,x_eta,y_ksi,y_eta,
     $              jac_xy)
            ENDIF
            DO iqty=1,scv%nqty
               uk(iqty,:)=ux(iqty,:)*x_ksi+uy(iqty,:)*y_ksi
               ue(iqty,:)=ux(iqty,:)*x_eta+uy(iqty,:)*y_eta
            ENDDO
c-----------------------------------------------------------------------
c     compute fluxes and sources.
c-----------------------------------------------------------------------
            CALL job2_schur(psv%t,psv%h*psv%theta,ksi,eta,u,uk,ue,
     $           fk_u,fk_uk,fk_ue,fe_u,fe_uk,fe_ue,s_u,s_uk,s_ue)
            jac_xy = jac_xy*psv%weight
c-----------------------------------------------------------------------
c     calculate kappa transformation matrix.
c-----------------------------------------------------------------------
            kap1=0
            kap2=0
            DO imm=1,(np+1)**2
               kap1(:,imm) = x_ksi*qpbq(:,imm) + y_ksi*pqbq(:,imm)
               kap2(:,imm) = x_eta*qpbq(:,imm) + y_eta*pqbq(:,imm)
            ENDDO
c-----------------------------------------------------------------------
c     loop over trial polynomials.
c-----------------------------------------------------------------------
            DO imm=0,scv%elsize2-1
               iqty=MOD(imm,nqty)+1
               mm=imm/nqty+1
c-----------------------------------------------------------------------
c     loop over basis polynomials.
c-----------------------------------------------------------------------
               DO inn=0,scv%elsize2-1
                  jqty=MOD(inn,nqty)+1
                  nn=inn/nqty+1
c-----------------------------------------------------------------------
c     compute derivatives of rhs with respect to dependent variables.
c-----------------------------------------------------------------------
                  mat_1d(inn+1,imm+1) = SUM(jac_xy
     $                 *(ppbq(:,mm)*(ppbq(:,nn)*s_u(:,jqty,iqty)
     $                 + kap1(:,nn)*s_uk(:,jqty,iqty)
     $                 + kap2(:,nn)*s_ue(:,jqty,iqty))
     $                 + kap1(:,mm)*(ppbq(:,nn)*fk_u(:,jqty,iqty)
     $                 + kap1(:,nn)*fk_uk(:,jqty,iqty)
     $                 + kap2(:,nn)*fk_ue(:,jqty,iqty))
     $                 + kap2(:,mm)*(ppbq(:,nn)*fe_u(:,jqty,iqty)
     $                 + kap1(:,nn)*fe_uk(:,jqty,iqty)
     $                 + kap2(:,nn)*fe_ue(:,jqty,iqty))))
               ENDDO
            ENDDO
            scv%mat22(:,:,ix*(np+1):ix*(np+1)+np,iy*(np+1):iy*(np+1)+np)
     $           = RESHAPE(mat_1d,(/scv%elsize2,nqty,np+1,np+1/))
            ii = ii + 1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         WRITE(format1,10)nqty
         WRITE(format2,20)nqty
         WRITE(filename,'(a,2i1,a)')"smatrix",mpi_rank,mpi_size,".out"
         OPEN(UNIT=debug_unit,FILE=TRIM(filename),STATUS="UNKNOWN")
         DO iy=scv%jxmin,scv%jxmax
            DO ix=scv%jymin,scv%jymax
               DO inn=0,scv%elsize2-1,nqty
                  nn=inn/nqty+1
                  WRITE(debug_unit,format1)(jqty,jqty=1,nqty)
                  DO iqty=1,nqty
                     WRITE(debug_unit,format2)inn,ix,iy,iqty,
     $                    scv%mat22(inn+1:inn+nqty,iqty,ix,iy)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         WRITE(debug_unit,format1)(jqty,jqty=1,nqty)
         CLOSE(UNIT=debug_unit)
         CALL program_stop("Termination by p2_schur_schur.")
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_schur
c-----------------------------------------------------------------------
c     subprogram 5. p2_schur_matmul.
c     multiplies matrix times vector.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_matmul(scv,matrix,vector,matvec)

      TYPE(schur_type), INTENT(IN) :: scv
      REAL(r8), DIMENSION(:,:,0:,0:), INTENT(IN) :: matrix
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: vector
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: matvec

      INTEGER :: iqty,nqty,np,ix,iy,jx,jy,jx0,jy0
c-----------------------------------------------------------------------
c     set local sizes.
c-----------------------------------------------------------------------
      np=scv%np
      nqty=SIZE(matvec,1)
c-----------------------------------------------------------------------
c     multiply matrix times vector.
c-----------------------------------------------------------------------
      matvec=0
      DO iy=0,scv%iymax-scv%iymin-1
         jy0=iy*(np+1)
         DO ix=0,scv%ixmax-scv%ixmin-1
            jx0=ix*(np+1)
            DO jy=0,np
               DO jx=0,np
                  DO iqty=1,nqty
                     matvec(iqty,jx0+jx,jy0+jy)
     $                    =SUM(matrix(:,iqty,jx0+jx,jy0+jy)
     $                    *RESHAPE(vector(:,jx0:jx0+np,jy0:jy0+np),
     $                    (/SIZE(matrix,1)/)))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_matmul
c-----------------------------------------------------------------------
c     subprogram 6. p2_schur_solve.
c     solves the preconditioning equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_solve(scv,res,uu,emax,emin,its)

      TYPE(schur_type), INTENT(INOUT) :: scv
      REAL(r8), DIMENSION(scv%nqty,scv%jxmin:scv%jxmax,
     $     scv%jymin:scv%jymax), INTENT(IN) :: res
      REAL(r8), DIMENSION(scv%nqty,scv%jxmin:scv%jxmax,
     $     scv%jymin:scv%jymax), INTENT(OUT) :: uu
      REAL(r8), INTENT(OUT) :: emax,emin
      INTEGER, INTENT(OUT) :: its

      INTEGER :: iqty
c-----------------------------------------------------------------------
c     begin schur_solve event.
c-----------------------------------------------------------------------
c$$$      CALL PetscLogEventBegin(schur_solve_00,0,0,0,0,ierr)
c-----------------------------------------------------------------------
c     get r1 and r2.
c-----------------------------------------------------------------------
c$$$      CALL PetscLogEventBegin(schur_solve_s1,0,0,0,0,ierr)
      scv%r1=res(scv%rn1,:,:)
      scv%r2=res(scv%rn2,:,:)
c-----------------------------------------------------------------------
c     solve for s1.
c-----------------------------------------------------------------------
      SELECT CASE(mass_solve_type)
      CASE("condense")
         CALL p2_condense_rhs(scv%m11c,scv%r1,scv%m11c%cres)
         CALL KSPSolve(scv%m11c%ksp,scv%m11c%cres,scv%m11c%cdu,ierr)
         CALL KSPComputeExtremeSingularValues(scv%m11c%ksp,
     $        scv%m11c%emax,scv%m11c%emin,ierr)
         CALL KSPGetIterationNumber(scv%m11c%ksp,scv%m11c%its,ierr)
         scv%m11c%condno=scv%m11c%emax/scv%m11c%emin
         CALL p2_condense_solve(scv%m11c,scv%m11c%cdu,scv%s1)
      END SELECT
c$$$      CALL PetscLogEventEnd(schur_solve_s1,0,0,0,0,ierr)
c-----------------------------------------------------------------------
c     compute s2 and solve for u2.
c-----------------------------------------------------------------------
c$$$      CALL PetscLogEventBegin(schur_solve_u2,0,0,0,0,ierr)
      CALL p2_schur_matmul(scv,scv%mat21,scv%s1,scv%s2)
      scv%s2=scv%r2-scv%s2
      SELECT CASE(schur_solve_type)
      CASE("condense")
         CALL p2_condense_rhs(scv%schurc,scv%s2,scv%schurc%cres)
         CALL KSPSolve(scv%schurc%ksp,scv%schurc%cres,
     $        scv%schurc%cdu,ierr)
         CALL KSPComputeExtremeSingularValues(scv%schurc%ksp,
     $        scv%schurc%emax,scv%schurc%emin,ierr)
         scv%schurc%condno=scv%schurc%emax/scv%schurc%emin
         CALL KSPGetIterationNumber(scv%schurc%ksp,scv%schurc%its,ierr)
         CALL p2_condense_solve(scv%schurc,scv%schurc%cdu,scv%u2)
         emax=scv%schurc%emax
         emin=scv%schurc%emin
         its=scv%schurc%its
      END SELECT
c$$$      CALL PetscLogEventEnd(schur_solve_u2,0,0,0,0,ierr)
c-----------------------------------------------------------------------
c     compute t1 and solve for u1.
c-----------------------------------------------------------------------
c$$$      CALL PetscLogEventBegin(schur_solve_u1,0,0,0,0,ierr)
      CALL p2_schur_matmul(scv,scv%mat12,scv%u2,scv%t1)
      SELECT CASE(mass_solve_type)
      CASE("condense")
         CALL p2_condense_rhs(scv%m11c,scv%t1,scv%m11c%cres)
         CALL KSPSolve(scv%m11c%ksp,scv%m11c%cres,scv%m11c%cdu,ierr)
         CALL KSPComputeExtremeSingularValues(scv%m11c%ksp,
     $        scv%m11c%emax,scv%m11c%emin,ierr)
         CALL KSPGetIterationNumber(scv%m11c%ksp,scv%m11c%its,ierr)
         scv%m11c%condno=scv%m11c%emax/scv%m11c%emin
         CALL p2_condense_solve(scv%m11c,scv%m11c%cdu,scv%u1)
      END SELECT
c$$$      CALL PetscLogEventEnd(schur_solve_u1,0,0,0,0,ierr)
c-----------------------------------------------------------------------
c     assemble full solution uu.
c-----------------------------------------------------------------------
c$$$      CALL PetscLogEventBegin(schur_solve_uu,0,0,0,0,ierr)
      uu(scv%rn1,:,:)=scv%s1-scv%u1
      uu(scv%rn2,:,:)=scv%u2
c$$$      CALL PetscLogEventEnd(schur_solve_uu,0,0,0,0,ierr)
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(schur_diagnose)CALL p2_schur_diagnose(scv,res,uu)
c-----------------------------------------------------------------------
c     end schur_solve event.
c-----------------------------------------------------------------------
c$$$      CALL PetscLogEventEnd(schur_solve_00,0,0,0,0,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_solve
c-----------------------------------------------------------------------
c     subprogram 7. p2_schur_diagnose.
c     solves the preconditioning equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_diagnose(scv,res,uu)

      TYPE(schur_type), INTENT(INOUT) :: scv
      REAL(r8), DIMENSION(scv%nqty,scv%jxmin:scv%jxmax,
     $     scv%jymin:scv%jymax), INTENT(IN) :: res,uu

      CHARACTER(80) :: filename
      CHARACTER(32) :: name
      INTEGER :: iqty
      REAL(r8), DIMENSION(SIZE(scv%s1,1),SIZE(scv%s1,2),SIZE(scv%s1,3))
     $     :: ls1,lu1
      REAL(r8), DIMENSION(SIZE(scv%u2,1),SIZE(scv%u2,2),SIZE(scv%u2,3))
     $     :: su2,su2_err
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/2x,"solve",3x,"its",4x,"condno"/)
 20   FORMAT(1x,a6,i6,1p,e11.3)
c-----------------------------------------------------------------------
c     allocate diagnostic objects.
c-----------------------------------------------------------------------
      CALL p2_schur_diagnose_alloc(scv,scv%sdt1,scv%nqty1)
      CALL p2_schur_diagnose_alloc(scv,scv%sdt2,scv%nqty2)
c-----------------------------------------------------------------------
c     open output file.
c-----------------------------------------------------------------------
      WRITE(filename,'(a,2i1,a)')"schur",mpi_rank,mpi_size,".out"
      OPEN(UNIT=debug_unit,FILE=TRIM(filename),STATUS="UNKNOWN")
c-----------------------------------------------------------------------
c     diagnose unassembled res and r1.
c-----------------------------------------------------------------------
      CALL debug2(res,"res",debug_unit)
      CALL debug2(scv%r1,"ur1",debug_unit)
      CALL debug2(scv%r2,"ur2",debug_unit)
c-----------------------------------------------------------------------
c     compute and assemble ls1.
c-----------------------------------------------------------------------
      CALL p2_schur_assemble(scv%sdt1,scv%r1)
      CALL p2_schur_matmul(scv,scv%mat110,scv%s1,ls1)
      CALL p2_schur_assemble(scv%sdt1,ls1)
c-----------------------------------------------------------------------
c     diagnose r1, s1, and ls1.
c-----------------------------------------------------------------------
      CALL debug2(scv%r1,"ar1",debug_unit)
      CALL debug2(scv%s1,"s1",debug_unit)
      CALL debug2(ls1,"ls1",debug_unit)
c-----------------------------------------------------------------------
c     compute and assemble su2.
c-----------------------------------------------------------------------
      CALL p2_schur_assemble(scv%sdt2,scv%s2)
      CALL p2_schur_matmul(scv,scv%mat220,scv%u2,su2)
      CALL p2_schur_assemble(scv%sdt2,su2)
c-----------------------------------------------------------------------
c     compute su2_err.
c-----------------------------------------------------------------------
      WHERE(scv%s2 /= 0)
         su2_err=su2/scv%s2-1
      ELSEWHERE
         su2_err=0
      END WHERE
      WRITE(name,'(a,1p,e9.3)')"su2_err, max = ",MAXVAL(ABS(su2_err))
c-----------------------------------------------------------------------
c     diagnose s2, su2, and su2_err.
c-----------------------------------------------------------------------
      CALL debug2(scv%s2,"s2",debug_unit)
      CALL debug2(su2,"su2",debug_unit)
      CALL debug2(su2_err,name,debug_unit)
c-----------------------------------------------------------------------
c     compute and assemble lu1.
c-----------------------------------------------------------------------
      CALL p2_schur_assemble(scv%sdt1,scv%t1)
      CALL p2_schur_matmul(scv,scv%mat110,scv%u1,lu1)
      CALL p2_schur_assemble(scv%sdt1,lu1)
c-----------------------------------------------------------------------
c     diagnose t1, u2, and su2.
c-----------------------------------------------------------------------
      CALL debug2(scv%t1,"t1",debug_unit)
      CALL debug2(lu1,"lu1",debug_unit)
c-----------------------------------------------------------------------
c     diagnose uu.
c-----------------------------------------------------------------------
      scv%u1=scv%s1-scv%u1
      CALL debug2(scv%u1,"u1",debug_unit)
      CALL debug2(scv%u2,"u2",debug_unit)
      CALL debug2(uu,"uu",debug_unit)
c-----------------------------------------------------------------------
c     write iterations and condition numbers, mass matrix.
c-----------------------------------------------------------------------
      WRITE(debug_unit,10)
      SELECT CASE(mass_solve_type)
      CASE("condense")
         WRITE(debug_unit,20)"m11c",scv%m11c%its,scv%m11c%condno
      END SELECT
c-----------------------------------------------------------------------
c     write iterations and condition numbers, schur complement.
c-----------------------------------------------------------------------
      SELECT CASE(schur_solve_type)
      CASE("condense")
         WRITE(debug_unit,20)" schur",scv%schurc%its,scv%schurc%condno
      END SELECT
c-----------------------------------------------------------------------
c     deallocate diagnostic objects.
c-----------------------------------------------------------------------
      CALL p2_schur_diagnose_dealloc(scv%sdt1)
      CALL p2_schur_diagnose_dealloc(scv%sdt2)
c-----------------------------------------------------------------------
c     close output file and terminate.
c-----------------------------------------------------------------------
      CLOSE(UNIT=debug_unit)
      CALL program_stop("Termination after p2_schur_diagnose.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_diagnose
c-----------------------------------------------------------------------
c     subprogram 8. p2_schur_assemble.
c     assembles a vector.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_assemble(sdt,uu)

      TYPE(schur_diagnose_type), INTENT(INOUT) :: sdt
      REAL(r8), DIMENSION(:,:,:), INTENT(INOUT) :: uu
c-----------------------------------------------------------------------
c     transfer to flat petsc local vector.
c-----------------------------------------------------------------------
      CALL VecSetValues(sdt%local,sdt%size,sdt%index,uu,
     $     INSERT_VALUES,ierr)
      CALL VecAssemblyBegin(sdt%local,ierr)
      CALL VecAssemblyEnd(sdt%local,ierr)
c-----------------------------------------------------------------------
c     scatter from local to global, adding cell boundary values.
c-----------------------------------------------------------------------
      CALL VecSet(sdt%global,zero,ierr)
      CALL VecScatterBegin(sdt%scat,sdt%local,sdt%global,
     $     ADD_VALUES,SCATTER_FORWARD,ierr)
      CALL VecScatterEnd(sdt%scat,sdt%local,sdt%global,
     $     ADD_VALUES,SCATTER_FORWARD,ierr)
c-----------------------------------------------------------------------
c     scatter assembled values back from global to local.
c-----------------------------------------------------------------------
      CALL VecScatterBegin(sdt%scat,sdt%global,sdt%local,
     $     INSERT_VALUES,SCATTER_REVERSE,ierr)
      CALL VecScatterEnd(sdt%scat,sdt%global,sdt%local,
     $     INSERT_VALUES,SCATTER_REVERSE,ierr)
c-----------------------------------------------------------------------
c     transfer values back to Fortran array.
c-----------------------------------------------------------------------
      CALL VecGetValues(sdt%local,sdt%size,sdt%index,uu,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_assemble
c-----------------------------------------------------------------------
c     subprogram 9. p2_schur_diagnose_alloc.
c     allocates schur-diagnose types.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_diagnose_alloc(scv,sdt,nqty)

      TYPE(schur_type), INTENT(IN) :: scv
      TYPE(schur_diagnose_type), INTENT(OUT) :: sdt
      INTEGER, INTENT(IN) :: nqty

      INTEGER :: i,ix,ixp,ixp0,iy,iyp,iyp0,iqty,counter,np,s1,lsize
      INTEGER, DIMENSION(:), ALLOCATABLE :: is0
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: idmap
      IS :: ris0,ris1
c-----------------------------------------------------------------------
c     create native and petsc arrays.
c-----------------------------------------------------------------------
      lsize=(scv%jxmax-scv%jxmin+1)*(scv%jymax-scv%jymin+1)*nqty
      sdt%size=lsize
      CALL VecCreateMPI(comm,lsize,PETSC_DETERMINE,sdt%global,ierr)
      CALL VecCreateSeq(PETSC_COMM_SELF,lsize,sdt%local,ierr)
      ALLOCATE(sdt%index(lsize))
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      ALLOCATE(idmap(nqty,0:scv%nxp,0:scv%nyp),is0(0:lsize-1))
      sdt%index=(/(i,i=0,lsize-1)/)
      CALL p2_schur_map_init(scv,nqty,idmap)
      counter=0
      is0=-1
c-----------------------------------------------------------------------
c     compute index set.
c-----------------------------------------------------------------------
      DO iyp=scv%jymin,scv%jymax
         DO ixp=scv%jxmin,scv%jxmax
            DO iqty=1,nqty
               is0(counter)=idmap(iqty,ixp,iyp)
               counter = counter+1
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     create global scatter context over processor.
c-----------------------------------------------------------------------
      CALL ISCreateGeneral(PETSC_COMM_SELF,lsize,is0,
     $     PETSC_COPY_VALUES,ris0,ierr)
      CALL ISCreateStride(PETSC_COMM_SELF,lsize,0,1,ris1,ierr)
      CALL VecScatterCreate(sdt%local,ris1,sdt%global,ris0,sdt%scat,
     $     ierr)
      CALL ISDestroy(ris0,ierr)
      CALL ISDestroy(ris1,ierr)
      DEALLOCATE(is0,idmap)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_diagnose_alloc
c-----------------------------------------------------------------------
c     subprogram 10. p2_schur_diagnose_dealloc.
c     deallocates schur-diagnose types.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_diagnose_dealloc(sdt)

      TYPE(schur_diagnose_type), INTENT(INOUT) :: sdt
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DEALLOCATE(sdt%index)
      CALL VecDestroy(sdt%local,ierr)
      CALL VecDestroy(sdt%global,ierr)
      CALL VecScatterDestroy(sdt%scat,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_diagnose_dealloc
c-----------------------------------------------------------------------
c     subprogram 11. p2_schur_map_init
c     create and store global 3d to 1d mapping.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_schur_map_init(scv,nqty,idmap)
      
      TYPE(schur_type) :: scv
      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(:,0:,0:), INTENT(OUT) :: idmap

      INTEGER :: i,iqty,ix,iy,ixp,iyp,jxp,jyp,jxpmax,jypmax, nx,ny,np,
     $     irank
      INTEGER, DIMENSION(0:mpi_size-1) ::
     $     ixmin_v,ixmax_v,iymin_v,iymax_v
      INTEGER, DIMENSION(nqty,0:scv%nx*scv%np,0:scv%ny*scv%np) :: jdmap
c-----------------------------------------------------------------------
c     initialize local variables  
c-----------------------------------------------------------------------
      np=scv%np
      nx=scv%nx
      ny=scv%ny
c-----------------------------------------------------------------------
c     gather indices over processors.
c-----------------------------------------------------------------------
      CALL gather(scv%ixmin,ixmin_v)
      CALL gather(scv%ixmax,ixmax_v)
      CALL gather(scv%iymin,iymin_v)
      CALL gather(scv%iymax,iymax_v)
c-----------------------------------------------------------------------
c     create 3d to 1d mapping of continuous full vector
c-----------------------------------------------------------------------
      i = 0
      jdmap=0
      DO irank=0,mpi_size-1
         DO iy = iymin_v(irank)+1,iymax_v(irank)
            jypmax = np - 1
            IF (iy == ny .AND. (.NOT. scv%yperiodic))jypmax = np
            DO ix = ixmin_v(irank)+1,ixmax_v(irank)
               jxpmax = np - 1
               IF (ix == nx .AND. (.NOT. scv%xperiodic))jxpmax = np
               DO jyp = 0,jypmax
                  iyp = (iy-1)*np+jyp
                  DO jxp = 0,jxpmax
                     ixp = (ix-1)*np+jxp
                     DO iqty = 1,nqty
                        jdmap(iqty,ixp,iyp) = i
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
      idmap=0
      DO iyp = 0,scv%nyp
         jyp = iyp/(np+1)*np + MOD(iyp,np+1)
         IF(scv%yperiodic .AND. iyp == scv%nyp)jyp=0
         DO ixp = 0,scv%nxp
            jxp = ixp/(np+1)*np + MOD(ixp,np+1)
            IF(scv%xperiodic .AND. ixp == scv%nxp)jxp=0
            idmap(:,ixp,iyp) = jdmap(:,jxp,jyp)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_schur_map_init
      END MODULE p2_schur_mod
