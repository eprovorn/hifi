c-----------------------------------------------------------------------
c     file p2_snes.F.
c     uses petsc snes nonlinear solvers for time step advance.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_snes_mod.
c     1. p2_snes_init.
c     2. p2_snes_stop.
c     3. p2_snes_step.
c     4. p2_snes_function.
c     5. p2_snes_norm.
c     6. p2_snes_error.
c     7. p2_snes_jacobian.
c     8. p2_snes_precon.
c     9. p2_snes_reduce_dt.
c     10. p2_snes_diagnose1.
c     11. p2_snes_diagnose2.
c-----------------------------------------------------------------------
c     subprogram 0. p2_snes_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_snes_mod
      USE p2_rj_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      TYPE :: snes_type
      INTEGER :: jacs,newjacs
      REAL(r8) :: dt_dfac
      REAL(r8), DIMENSION(:), POINTER :: errnorm
      REAL(r8), DIMENSION(:,:,:), POINTER :: rhs,res
      SNES :: snes
      KSP :: ksp
      PC :: pc
      Vec :: u,du,f,ug0,ug00,ug000
      Mat :: fullmat
      END TYPE snes_type


      CHARACTER(16), PARAMETER, PRIVATE :: norm_type="max"
      INTEGER :: ksp_restart=30,ksps_max = 5000
      REAL(r8) :: errtol=1e-4,condno=0
      LOGICAL :: adapt_dt=.FALSE.,always_pc_reset=.TRUE.,
     $     du_diagnose=.FALSE.,newjac=.FALSE.,pc_reset=.FALSE.
      TYPE(snes_type), POINTER :: stv
      TYPE(snes_type), TARGET :: st_solver,st_grid

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_snes_init.
c     initializes nonlinear solver.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_init(psv,stv)

      TYPE(p2_sel_type) :: psv
      TYPE(snes_type) :: stv

      INTEGER :: ctx
c-----------------------------------------------------------------------
c     create vectors.
c-----------------------------------------------------------------------
      CALL VecDuplicate(psv%global,stv%u,ierr)
      CALL VecDuplicate(psv%global,stv%du,ierr)
      CALL VecDuplicate(psv%global,stv%f,ierr)
      CALL VecDuplicate(psv%global,stv%ug0,ierr)
      IF(psv%step_type=="bdf2")THEN
         CALL VecDuplicate(psv%global,stv%ug00,ierr)
         CALL VecDuplicate(psv%global,stv%ug000,ierr)
      ENDIF
c-----------------------------------------------------------------------
c     create snes objects.
c-----------------------------------------------------------------------
      CALL SNESCreate(comm,stv%snes,ierr)
      CALL SNESSetType(stv%snes,SNESNEWTONLS,ierr)
      CALL SNESSetFunction(stv%snes,stv%f,p2_snes_function,ctx,ierr)
      CALL MatCreateSNESMF(stv%snes,stv%fullmat,ierr)
      CALL SNESSetJacobian(stv%snes,stv%fullmat,stv%fullmat,
     $     p2_snes_jacobian,ctx,ierr)
      CALL SNESSetTolerances(stv%snes,min_eps*errtol,errtol,zero,
     $     PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr)
      CALL SNESSetConvergenceTest(stv%snes,p2_snes_error,ctx,
     $     PETSC_NULL_FUNCTION,ierr)
      CALL SNESSetFromOptions(stv%snes,ierr)
c-----------------------------------------------------------------------
c     set up preconditioner.
c-----------------------------------------------------------------------
      CALL SNESGetKSP(stv%snes,stv%ksp,ierr)
      CALL KSPSetType(stv%ksp,KSPPREONLY,ierr)
      CALL KSPGetPC(stv%ksp,stv%pc,ierr)
      CALL PCSetType(stv%pc,PCSHELL,ierr)
      CALL PCShellSetContext(stv%pc,ctx,ierr)
      CALL PCShellSetApply(stv%pc,p2_snes_precon,ierr)
c-----------------------------------------------------------------------
c     allocate native arrays.
c-----------------------------------------------------------------------
      ALLOCATE(
     $     stv%rhs(psv%nqty,psv%jxmin:psv%jxmax,psv%jymin:psv%jymax),
     $     stv%res(psv%nqty,psv%jxmin:psv%jxmax,psv%jymin:psv%jymax),
     $     stv%errnorm(psv%nqty))
      stv%rhs=0
      stv%res=0
      stv%errnorm=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_init
c-----------------------------------------------------------------------
c     subprogram 2. p2_snes_stop.
c     stops execution and deallocates all objects.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_stop(message)

      CHARACTER(*), INTENT(IN) :: message
c-----------------------------------------------------------------------
c     deallocate all snes objects.
c-----------------------------------------------------------------------
      CALL VecDestroy(st_solver%u,ierr)
      CALL VecDestroy(st_solver%du,ierr)
      CALL VecDestroy(st_solver%f,ierr)
      CALL VecDestroy(st_solver%ug0,ierr)
      IF(solver%step_type=="bdf2")THEN
         CALL VecDestroy(st_solver%ug00,ierr)
         CALL VecDestroy(st_solver%ug000,ierr)
      ENDIF
      CALL MatDestroy(st_solver%fullmat,ierr)
      CALL SNESDestroy(st_solver%snes,ierr)
      DEALLOCATE(st_solver%rhs,st_solver%res,st_solver%errnorm)
      IF(.NOT. (fd_test .OR. fd_bound_test) .AND. adapt_grid)THEN
         CALL VecDestroy(st_grid%u,ierr)
         CALL VecDestroy(st_grid%du,ierr)
         CALL VecDestroy(st_grid%f,ierr)
         CALL VecDestroy(st_grid%ug0,ierr)
         CALL MatDestroy(st_grid%fullmat,ierr)
         CALL SNESDestroy(st_grid%snes,ierr)
         DEALLOCATE(st_grid%rhs,st_grid%res,st_grid%errnorm)
      ENDIF
c-----------------------------------------------------------------------
c     free the rest of the work space.
c-----------------------------------------------------------------------
      CALL p2_rj_stop(message)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_stop
c-----------------------------------------------------------------------
c     subprogram 3. p2_snes_step.
c     takes a time step with nonlinear solver.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_step(crd,stv)

      TYPE(coord_type) :: crd
      TYPE(snes_type) :: stv

      CHARACTER(80) :: message
      INTEGER :: itry,ntry,newts,dnewts
      REAL(r8) :: hsum,denom,bdfa,bdfb,error
      SNESConvergedReason :: reason
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      CALL VecCopy(pt2%global,stv%u,ierr)
      ntry=pt2%itmax/pt2%itmax_decr
      stv%dt_dfac=pt2%dt_decr
      stv%jacs=0
      stv%newjacs=0
c-----------------------------------------------------------------------
c     initialize step, theta scheme.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%step_type)
      CASE("theta","theta2")
         CALL p2_rj_rhs(pt2%t,crd,stv%u,stv%rhs)
         CALL VecCopy(stv%u,stv%ug0,ierr)
         CALL p2_snes_norm(stv%rhs,stv%errnorm)
c-----------------------------------------------------------------------
c     initialize step, bdf2.
c-----------------------------------------------------------------------
      CASE("bdf2")
         hsum=pt2%h_old+pt2%h
         denom=pt2%h_old+2*pt2%h
         bdfa=hsum**2/(pt2%h_old*denom)
         bdfb=pt2%h**2/(pt2%h_old*denom)
         pt2%bdfc=pt2%h*hsum/denom
         CALL VecCopy(stv%ug00,stv%ug000,ierr)
         CALL VecCopy(stv%u,stv%ug00,ierr)
         CALL VecCopy(stv%ug00,stv%ug0,ierr)
         CALL VecScale(stv%ug0,bdfa,ierr)
         CALL VecAXPY(stv%ug0,-bdfb,stv%ug000,ierr)
         CALL p2_snes_norm(stv%rhs,stv%errnorm)
c-----------------------------------------------------------------------
c     initialize step, abort.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL p2_snes_stop("p2_snes_step: cannot recognize step_type = 
     $        "//TRIM(pt2%step_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     nonlinear solve.
c-----------------------------------------------------------------------
      pt2%t=pt2%t+pt2%h
      newts=0
      DO itry=0,ntry
         CALL SNESSolve(stv%snes,PETSC_NULL_OBJECT,stv%u,ierr)
c-----------------------------------------------------------------------
c     abort on failure.
c-----------------------------------------------------------------------
         IF(ierr /= 0)THEN
            CALL SNESGetConvergedReason(stv%snes,reason,ierr)
            WRITE(message,'(a,i2)')
     $           "SNES solver did not converge. Reason=",reason
            CALL p2_snes_stop(message)
         ENDIF
c-----------------------------------------------------------------------
         CALL SNESGetIterationNumber(stv%snes,dnewts,ierr)
         CALL SNESGetConvergedReason(stv%snes,reason,ierr)
         newts=newts+dnewts
         IF(reason > 0 .AND. reason /= 7)EXIT
         IF(stv%jacs >= pt2%itmax_decr)THEN
            newjac=.TRUE.
            IF(pt2%couple_msize > 0 .AND. adapt_dt .AND. 
     $           stv%newjacs > 0)CALL p2_snes_reduce_dt
            stv%jacs=0
         ENDIF
      ENDDO
      CALL VecCopy(stv%u,pt2%global,ierr)
      pt2%newts=pt2%newts+newts
c-----------------------------------------------------------------------
c     if desired, increase or decrease time-step size.
c-----------------------------------------------------------------------
      pt2%h_old=pt2%h
      IF(adapt_dt)THEN
         IF(pt2%ksps/pt2%newts >= ksps_max)pt2%h=pt2%h*pt2%dt_decr
         IF(newts < pt2%itmax_incr .AND. itry == 0
     $     .AND. pt2%ksps/pt2%newts < ksps_max/2)
     $        pt2%h=MIN(pt2%h*pt2%dt_incr,pt2%dtmax)
         IF(pt2%h < pt2%dtmin)CALL p2_snes_stop
     $        ("Time-step crash. dt reduced below dtmin")
         IF(pt2%h /= pt2%h_old)THEN
            newjac=.TRUE.
            pc_reset=.TRUE.
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_step
c-----------------------------------------------------------------------
c     subprogram 4. p2_snes_function.
c     the nonlinear function to be solved.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_function(snes,x,f,ctx,ierr)

      SNES :: snes
      Vec :: x,f
      INTEGER, INTENT(IN) ::  ctx
      INTEGER, INTENT(OUT) :: ierr

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: np,iqty,jqty,ix,iy,jx,jy,jx0,jy0,ii
      REAL(r8), DIMENSION(pt2%nqty,(pt2%np+1)**2) :: cell
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: du,mdu,rhs
      REAL(r8), DIMENSION(pt2%nqty,pt2%couple_msize,0:pt2%np,0:pt2%np) 
     $     :: mmat_loc
c-----------------------------------------------------------------------
c     compute residual, mass terms.
c-----------------------------------------------------------------------
      mdu=0
      CALL VecWAXPY(stv%du,-one,stv%ug0,x,ierr)
      IF(pt2%couple_msize > 0)THEN
         CALL p2_sel_get_local(pt2%jscat,pt2%jvec,stv%du,du)
         np=pt2%np
         DO iy=pt2%iymin,pt2%iymax-1
            jy0=iy*(np+1)
            DO ix=pt2%ixmin,pt2%ixmax-1
               jx0=ix*(np+1)

               cell=RESHAPE(du(:,jx0:jx0+np,jy0:jy0+np),
     $              (/pt2%nqty,(np+1)**2/))
               DO jy=0,np
                  DO jx=0,np
                     mmat_loc(:,:,jx,jy)
     $                    = MATMUL(cell,pt2%mass(:,:,jx0+jx,jy0+jy))
                     DO ii=1,pt2%couple_msize
                        jqty=pt2%couple_mlist(1,ii)
                        iqty=pt2%couple_mlist(2,ii)
                        mdu(iqty,jx0+jx,jy0+jy)=mdu(iqty,jx0+jx,jy0+jy)
     $                       + mmat_loc(jqty,ii,jx,jy)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     compute residual, theta scheme.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%step_type)
      CASE("theta")
         CALL p2_rj_rhs(pt2%t,coord,x,rhs)
         stv%res = mdu - pt2%h*pt2%theta*rhs
         WHERE(pt2%theta_fac)
            stv%res = stv%res - pt2%h*(one-pt2%theta)*stv%rhs
         END WHERE
      CASE("theta2")
         CALL VecAYPX(stv%du,pt2%theta-one,x,ierr)
         CALL p2_rj_rhs(pt2%t+(pt2%theta-one)*pt2%h,coord,stv%du,du)
         CALL p2_rj_rhs(pt2%t,coord,x,rhs)
         stv%res = mdu - pt2%h*pt2%theta*rhs
         WHERE(pt2%theta_fac)
            stv%res = stv%res - pt2%h*(du - rhs*pt2%theta)
         ENDWHERE
c-----------------------------------------------------------------------
c     compute residual, bdf2.
c-----------------------------------------------------------------------
      CASE("bdf2")
         CALL p2_rj_rhs(pt2%t,coord,x,stv%rhs)
         stv%res = mdu - pt2%bdfc*stv%rhs
c-----------------------------------------------------------------------
c     default.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL p2_snes_stop
     $        ("p2_snes_function: cannot recognize step_type = 
     $        "//TRIM(pt2%step_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     copy residual into a PETSc vector.
c-----------------------------------------------------------------------
      CALL p2_sel_set(.FALSE.,stv%res,pt2%rscat,pt2%jvec,f)
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL p2_snes_diagnose1(x,du,mdu,rhs,stv%rhs,stv%res)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      ierr=0
      RETURN
      END SUBROUTINE p2_snes_function
c-----------------------------------------------------------------------
c     subprogram 5. p2_snes_norm.
c     computes error norm.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_norm(rhs,errnorm)

      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(IN) :: rhs
      REAL(r8), DIMENSION(pt2%nqty), INTENT(OUT) :: errnorm

      INTEGER :: iqty
      REAL(r8), DIMENSION(pt2%nqty) :: norm0
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: rhs0
c-----------------------------------------------------------------------
c     compute vector.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%step_type)
      CASE("theta","theta2")
         rhs0=pt2%h*rhs
      CASE("bdf2")
         rhs0=(2._r8/3._r8)*pt2%h*rhs
      CASE DEFAULT
         CALL p2_snes_stop("p2_snes_norm: cannot recognize step_type = 
     $        "//TRIM(pt2%step_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     collect values from adjacent grid cells.
c-----------------------------------------------------------------------
      CALL p2_sel_set(.FALSE.,rhs0,pt2%rscat,pt2%jvec,stv%f)
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,stv%f,rhs0)
c-----------------------------------------------------------------------
c     compute maximum norm.
c-----------------------------------------------------------------------
      SELECT CASE(norm_type)
      CASE("max")
         DO iqty=1,pt2%nqty
            norm0(iqty)=MAXVAL(ABS(rhs0(iqty,:,:)))
         ENDDO
         CALL MPI_Allreduce(norm0,errnorm,pt2%nqty,
     $        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
c-----------------------------------------------------------------------
c     compute L2 norm.
c-----------------------------------------------------------------------
      CASE("l2","l22")
         DO iqty=1,pt2%nqty
            norm0(iqty)=SUM(rhs0(iqty,:,:)**2)
         ENDDO
         CALL MPI_Allreduce(norm0,errnorm,pt2%nqty,
     $        MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
         IF(norm_type == "l2")errnorm=SQRT(ABS(errnorm))
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL p2_snes_stop("p2_snes_norm: "
     $        //"cannot recognize norm_type = "//TRIM(norm_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     set minimum value for static equations.
c-----------------------------------------------------------------------
      DO iqty=1,pt2%nqty
         IF(pt2%stat_qty(iqty) .AND. errnorm(iqty) < SQRT(min_eps))
     $        errnorm(iqty)=one
      ENDDO
      errnorm=MAX(errnorm,min_eps)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_norm
c-----------------------------------------------------------------------
c     subprogram 6. p2_snes_error.
c     computes error.     
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_error(snes,it,xnorm,gnorm,fnorm,reason,ctx,
     $     ierr)

      SNES :: snes
      INTEGER, INTENT(IN) :: it
      REAL(r8), INTENT(IN) :: xnorm,gnorm,fnorm
      SNESConvergedReason :: reason
      INTEGER, INTENT(IN) :: ctx
      INTEGER, INTENT(OUT) :: ierr

      INTEGER :: iqty
      REAL(r8) :: error0,error
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: res0
      Vec :: res
c-----------------------------------------------------------------------
c     get solution vector.
c-----------------------------------------------------------------------
      CALL SNESGetFunction(snes,res,
     $     PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,ierr)
c-----------------------------------------------------------------------
c     get local vector and normalize.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,res,res0)
      DO iqty=1,pt2%nqty
         res0(iqty,:,:)=res0(iqty,:,:)/stv%errnorm(iqty)
      ENDDO
c-----------------------------------------------------------------------
c     compute maximum error.
c-----------------------------------------------------------------------
      SELECT CASE(norm_type)
      CASE("max")
         error0=MAXVAL(ABS(res0))
         CALL MPI_Allreduce(error0,error,1,
     $        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
c-----------------------------------------------------------------------
c     compute L^2 error.
c-----------------------------------------------------------------------
      CASE("l2","l22")
         error0=SUM(res0**2)
         CALL MPI_Allreduce(error0,error,1,
     $        MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
         IF(norm_type == "l2")error=SQRT(ABS(error))
      END SELECT
c-----------------------------------------------------------------------
c     compute converged reason.
c-----------------------------------------------------------------------
      IF(mpi_rank==0 .AND. pt2%couple_msize <= 0)
     $        WRITE(*,'(a,e10.3)')"residual error =",error

      IF(error < errtol)THEN
         reason=SNES_CONVERGED_FNORM_RELATIVE
      ELSEIF(it >= pt2%itmax_decr)THEN
         reason=SNES_DIVERGED_MAX_IT
      ELSE
         reason=SNES_CONVERGED_ITERATING
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_error
c-----------------------------------------------------------------------
c     subprogram 7. p2_snes_jacobian.
c     jacobian formation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_jacobian(snes,x,a,b,ctx,ierr)

      SNES :: snes
      Vec :: x
      Mat :: a,b
      INTEGER, INTENT(IN) :: ctx
      INTEGER, INTENT(OUT) :: ierr

      REAL(r8), DIMENSION(pt2%elsize,pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: drdu
c-----------------------------------------------------------------------
c     assemble input matrix.
c-----------------------------------------------------------------------
      CALL MatAssemblyBegin(a,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd(a,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatSetOption(a,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE,ierr)
c-----------------------------------------------------------------------
c     update current solution vector and set return values.
c-----------------------------------------------------------------------
      ierr=0
      stv%jacs=stv%jacs+1
c-----------------------------------------------------------------------
c     diagnose and return if new jacobian not needed.
c-----------------------------------------------------------------------
      IF(.NOT. newjac)RETURN
c-----------------------------------------------------------------------
c     update count and reset newjac.
c-----------------------------------------------------------------------
      pt2%jacs=pt2%jacs+1
      stv%newjacs=stv%newjacs+1
      newjac=.FALSE.
c-----------------------------------------------------------------------
c     update jacobian
c-----------------------------------------------------------------------
      IF(pt2%couple_msize > 0)THEN
         CALL p2_rj_jac(pt2%t,coord,x,drdu,pt2%mass)
      ELSE
         CALL p2_rj_jac(pt2%t,coord,x,drdu)
      ENDIF
c-----------------------------------------------------------------------
c     transfer jacobian to appropriate location.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%solve_type)
      CASE("full")
         CALL p2_rj_j2petsc(drdu,pt2%jac)
      CASE("condense")
         CALL p2_condense_jac(ct_obj,drdu,ct_obj%jac)
      CASE("schur")
         CALL p2_schur_mat(pt2,scv,coord,drdu,x)
      END SELECT
c-----------------------------------------------------------------------
c     save full jacobian for diagnostics.
c-----------------------------------------------------------------------
      IF(du_diagnose .AND. pt2%solve_type /= "full")THEN
         CALL p2_sel_mat_alloc(pt2,pt2%jac)
         CALL p2_rj_j2petsc(drdu,pt2%jac)
      ENDIF
c-----------------------------------------------------------------------
c     reset preconditioner.
c-----------------------------------------------------------------------
      IF(pc_reset)THEN
         SELECT CASE(pt2%solve_type)
         CASE("full")
            CALL KSPSetOperators(pt2%ksp,pt2%jac,pt2%jac,ierr)
         CASE("condense")
            CALL KSPSetOperators(ct_obj%ksp,ct_obj%jac,ct_obj%jac,ierr)
         END SELECT
         IF(.NOT. always_pc_reset)pc_reset=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_jacobian
c-----------------------------------------------------------------------
c     subprogram 8. p2_snes_precon.
c     preconditioner.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE p2_snes_precon(matrix,x,y,ierr)

      Mat :: matrix
      Vec :: x,y
      INTEGER, INTENT(OUT) :: ierr

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER :: its
      REAL(r8) :: emax,emin
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: uu
c-----------------------------------------------------------------------
c     apply preconditioner.
c-----------------------------------------------------------------------
      SELECT CASE(pt2%solve_type)
      CASE("full")
         CALL KSPSolve(pt2%ksp,x,y,ierr)
         CALL KSPComputeExtremeSingularValues(pt2%ksp,emax,emin,ierr)
         CALL KSPGetIterationNumber(pt2%ksp,its,ierr)
      CASE("condense")
         CALL p2_condense_rhs(ct_obj,stv%res,ct_obj%cres)
         CALL KSPSolve(ct_obj%ksp,ct_obj%cres,ct_obj%cdu,ierr)
         CALL KSPComputeExtremeSingularValues(ct_obj%ksp,emax,emin,ierr)
         CALL KSPGetIterationNumber(ct_obj%ksp,its,ierr)
         CALL p2_condense_solve(ct_obj,ct_obj%cdu,uu)
         CALL p2_sel_set(.TRUE.,uu,pt2%jscat,pt2%jvec,y)
      CASE("schur")
         CALL p2_schur_solve(scv,stv%res,uu,emax,emin,its)
         CALL p2_sel_set(.TRUE.,uu,pt2%jscat,pt2%jvec,y)
      CASE DEFAULT
         CALL p2_snes_stop
     $        ("p2_snes_precon: cannot recognize solve_type =
     $        "//TRIM(pt2%solve_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     compute condition number.
c-----------------------------------------------------------------------
      IF(its >= ksp_restart)pc_reset=.TRUE.
      pt2%ksps=pt2%ksps+its
      ierr=0
      IF(emin == 0)THEN
         condno = 0
      ELSE
         condno=emax/emin
      ENDIF
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)CALL p2_snes_diagnose2(stv%u,x,y)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_precon
c-----------------------------------------------------------------------
c     subprogram 9. p2_snes_reduce_dt.
c     reduce time step.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_reduce_dt

      REAL(r8) :: hsum,denom,bdfa,bdfb
c-----------------------------------------------------------------------
c     general case.
c-----------------------------------------------------------------------
      pc_reset=.TRUE.
      pt2%t=pt2%t-pt2%h*(1-stv%dt_dfac)
      pt2%h=pt2%h*stv%dt_dfac
      stv%dt_dfac=stv%dt_dfac**2

      IF(pt2%h < pt2%dtmin)
     $     CALL p2_snes_stop("Time-step crash. dt reduced below dtmin")
c-----------------------------------------------------------------------
c     special treatment for bdf2.
c-----------------------------------------------------------------------
      IF(pt2%step_type == "bdf2")THEN
         hsum=pt2%h_old+pt2%h
         denom=pt2%h_old+2*pt2%h
         bdfa=hsum**2/(pt2%h_old*denom)
         bdfb=pt2%h**2/(pt2%h_old*denom)
         pt2%bdfc=pt2%h*hsum/denom
         CALL VecCopy(stv%ug00,stv%ug0,ierr)
         CALL VecScale(stv%ug0,bdfa,ierr)
         CALL VecAXPY(stv%ug0,-bdfb,stv%ug000,ierr)
      ENDIF
c-----------------------------------------------------------------------
c     reset counter.
c-----------------------------------------------------------------------
      stv%newjacs=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_reduce_dt
c-----------------------------------------------------------------------
c     subprogram 10. p2_snes_diagnose1.
c     monitor function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_diagnose1(u,du,mdu,rhs,rhs0,res)

      Vec, INTENT(IN) :: u
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax), INTENT(IN) :: du,mdu,rhs,rhs0,res

      CHARACTER(80) :: filename
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: uu
c-----------------------------------------------------------------------
c     get local solution vector.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,u,uu)
c-----------------------------------------------------------------------
c     open output file.
c-----------------------------------------------------------------------
      WRITE(filename,'(a,2i1,a)')"sdout",mpi_rank,mpi_size,".out"
      OPEN(UNIT=debug_unit,FILE=TRIM(filename),STATUS="UNKNOWN")
c-----------------------------------------------------------------------
c     write output.
c-----------------------------------------------------------------------
      CALL debug2(uu,"uu",debug_unit)
      CALL debug2(du,"du",debug_unit)
      CALL debug2(mdu,"mdu",debug_unit)
      CALL debug2(rhs,"rhs",debug_unit)
      CALL debug2(rhs0,"rhs0",debug_unit)
      CALL debug2(res,"res",debug_unit)
c-----------------------------------------------------------------------
c     close output unit and terminate.
c-----------------------------------------------------------------------
      CLOSE(UNIT=debug_unit)
      CALL p2_rj_stop("Termination by p2_snes_diagnose1.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      ierr=0
      RETURN
      END SUBROUTINE p2_snes_diagnose1
c-----------------------------------------------------------------------
c     subprogram 11. p2_snes_diagnose2.
c     diagnose preconditioned vector.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_snes_diagnose2(u,x,y)

      Vec :: u,x,y

      CHARACTER(80) :: filename
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: uu,rr,du,ratio
c-----------------------------------------------------------------------
c     open output file.
c-----------------------------------------------------------------------
      WRITE(filename,'(a,2i1,a)')
     $     "du",mpi_rank,mpi_size,".out."//TRIM(pt2%solve_type)
      OPEN(UNIT=debug_unit,FILE=TRIM(filename),status="UNKNOWN")
c-----------------------------------------------------------------------
c     get local vectors.
c-----------------------------------------------------------------------
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,u,uu)
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,x,rr)
      CALL p2_sel_get_local(pt2%jscat,pt2%jvec,y,du)
c-----------------------------------------------------------------------
c     compute ratio.
c-----------------------------------------------------------------------
      WHERE(uu /= 0)
         ratio=du/uu
      ELSEWHERE
         ratio=0
      ENDWHERE
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      CALL debug2(uu,"uu",debug_unit)
      CALL debug2(rr,"rr",debug_unit)
      CALL debug2(du,"du",debug_unit)
      CALL debug2(ratio,"ratio",debug_unit)
c-----------------------------------------------------------------------
c     close output file and terminate.
c-----------------------------------------------------------------------
      CLOSE(UNIT=debug_unit)
      CALL p2_rj_stop("Termination by p2_snes_diagnose2.")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_snes_diagnose2
      END MODULE p2_snes_mod
