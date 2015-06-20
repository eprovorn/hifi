c-----------------------------------------------------------------------
c     file p2_advance.F.
c     integrates PDEs with spectral element method.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_advance_mod.
c     1. p2_advance_run.
c     2. p2_advance_init.
c     3. p2_advance_restart.
c     4. p2_advance_switchsolve.
c     5. p2_advance_write_step.
c     6. p2_advance_write_header.
c-----------------------------------------------------------------------
c     subprogram 0. p2_advance_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_advance_mod
      USE fd_mod
      USE p2_grid_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      INTEGER, PRIVATE :: iout

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_advance_run.
c     sets up and integrates equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_advance_run(nx,ny,nbx,np,nq,nqty,nqty_schur,nstep,
     $     tmax,dmout,h,dtmax,dtmin,theta,nodal,xperiodic,yperiodic,
     $     restart_flag,restart_dir,restart_step,quad_type,
     $     step_type,itmax,itmax_incr,itmax_decr,dt_incr,dt_decr,
     $     solve_type,restart_time,restart_t0,outfile_type)

      INTEGER, INTENT(IN) :: nx,ny,nbx,np,nq,nqty,nqty_schur,
     $     restart_step,dmout,itmax,itmax_incr,itmax_decr,nstep
      REAL(r8), INTENT(IN) :: h,dtmax,dtmin,tmax,theta,dt_incr,dt_decr,
     $     restart_t0
      LOGICAL, INTENT(IN) :: nodal,xperiodic,yperiodic,restart_flag,
     $     restart_time
      CHARACTER(*), INTENT(IN) :: restart_dir,quad_type,step_type,
     $     solve_type,outfile_type

      LOGICAL :: change_grid=.FALSE.,bdf2first=.TRUE.
      CHARACTER(80) :: filename
      INTEGER :: iqty,adapt_start=0
      INTEGER, DIMENSION(:), ALLOCATABLE :: rnge
      REAL(r8) :: gridtol,gridtol1=1.,gridtol2=1.
c-----------------------------------------------------------------------
c     initial wallclock timer.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0)THEN
         CALL DATE_AND_TIME(date,time,zone,clock)
         wclock=(clock(5)*60+clock(6))*60+clock(7)+clock(8)/1000._r8
      ENDIF
c-----------------------------------------------------------------------
c     allocate and initialize "solver" objects.
c-----------------------------------------------------------------------
      CALL p2_sel_alloc(solver,nx,ny,nbx,np,nq,nqty,nqty_schur,nodal,
     $     xperiodic,yperiodic,quad_type,step_type,solve_type,
     $     restart_step)
c-----------------------------------------------------------------------
c     finish allocation and initialization of "solver" objects.
c-----------------------------------------------------------------------
      CALL p2_advance_init(restart_dir,outfile_type,restart_flag,
     $     restart_time,nodal,itmax,itmax_incr,itmax_decr,theta,h,
     $     dtmax,dtmin,tmax,dt_incr,dt_decr,restart_t0,filename)

      CALL p2_grid_error(solver)
c-----------------------------------------------------------------------
c     choose solver type.
c-----------------------------------------------------------------------
      SELECT CASE(solver%solve_type)
      CASE("full")
         CALL KSPGMRESSetRestart(solver%ksp,ksp_restart,ierr)
      CASE("condense")
         ALLOCATE(rnge(nqty))
         rnge=(/(iqty,iqty=1,nqty)/)
         CALL p2_condense_alloc(solver,ct_solver,rnge)
         DEALLOCATE(rnge)
         CALL KSPGMRESSetRestart(ct_solver%ksp,ksp_restart,ierr)
      CASE("schur")
         CALL p2_schur_alloc(solver,scv)
      CASE DEFAULT
         CALL program_stop("p2_advance_run: cannot recognize "
     $        //"solve_type = "//TRIM(solver%solve_type)//".")
      END SELECT
c-----------------------------------------------------------------------
c     initialize snes and switch to solver.
c-----------------------------------------------------------------------
      CALL p2_snes_init(solver,st_solver)
      CALL p2_advance_switchsolve("solver")
c-----------------------------------------------------------------------
c     prepare for adaptive grid
c-----------------------------------------------------------------------
      IF(adapt_grid)THEN
c-----------------------------------------------------------------------
c     allocate grid(p2_sel_type) variable
c-----------------------------------------------------------------------
         CALL p2_sel_alloc(grid,nx,ny,nbx,np,nq,ndim,0,nodal,
     $        .FALSE.,.FALSE.,quad_type,"theta",solve_type,0)
c-----------------------------------------------------------------------
c     alloc. and init. coord. transformation and "solver" mass matrix.
c-----------------------------------------------------------------------
         CALL p2_ct_alloc(coord,restart_dir,outfile_type)
         CALL p2_ct_calc(coord)
         IF(solver%couple_msize > 0)CALL p2_rj_mass
     $        (coord,solver%stat_qty,solver%theta_fac,solver%mass)
c-----------------------------------------------------------------------
c     initialize "grid" objects
c-----------------------------------------------------------------------
         CALL p2_grid_init(nodal,adapt_start,gridtol1,gridtol2)
         gridtol=gridtol1
         SELECT CASE(grid%solve_type)
         CASE("condense")
            ALLOCATE(rnge(ndim))
            rnge=(/(iqty,iqty=1,ndim)/)
            CALL p2_condense_alloc(grid,ct_grid,rnge)
            DEALLOCATE(rnge)
         END SELECT
         CALL p2_snes_init(grid,st_grid)
      ELSE
c-----------------------------------------------------------------------
c     alloc. and init. coord. transformation and "solver" mass matrix.
c-----------------------------------------------------------------------
         CALL p2_ct_alloc(coord,restart_dir,outfile_type)
         CALL p2_ct_calc(coord)
         IF(solver%couple_msize > 0)CALL p2_rj_mass
     $        (coord,solver%stat_qty,solver%theta_fac,solver%mass)
      ENDIF

      IF(monitor)CALL p2_diagnose_output(solver,outfile_type,"grid",
     $     dmout,solver%t,coord%crd_mat)
c-----------------------------------------------------------------------
c     set pointers to the "solver" object.
c-----------------------------------------------------------------------
      IF(adapt_grid .AND. adapt_start <= 0 
     $     .AND. solver%griderr > gridtol)THEN
         CALL p2_advance_switchsolve("grid")
         IF(mpi_rank == 0)CALL p2_grid_write_header
      ELSE
         CALL p2_advance_switchsolve("solver")
         IF(mpi_rank == 0)CALL p2_advance_write_header
      ENDIF
c-----------------------------------------------------------------------
c     initialize counters.
c-----------------------------------------------------------------------
      iout=0
      bdf2first=.TRUE.
      matzero=.FALSE.
      change_grid=.FALSE.
      solver%nstep=nstep+solver%istep
c-----------------------------------------------------------------------
c     monitor solution.
c-----------------------------------------------------------------------
      DO
         CALL p2_grid_error(pt2)
         IF((MOD(pt2%istep-restart_step,dmout) == 0 .OR. 
     $        pt2%istep >= pt2%nstep .OR. pt2%t >= pt2%tmax) .AND. 
     $        .NOT. grid_solver)THEN
            IF(monitor)CALL p2_diagnose_output_p(solver,outfile_type,
     $           "data",dmout,solver%t,solver%global)
            IF(mpi_rank==0)CALL p2_advance_write_step

            pt2%ksps=0
            pt2%jacs=0
            pt2%newts=0
            iout=iout+1
         ENDIF

         IF(solver%istep >= solver%nstep .OR. solver%t >= solver%tmax
     $        .OR. (solver%couple_msize <= 0 
     $        .AND. solver%istep == restart_step+1))EXIT
c-----------------------------------------------------------------------
c     use adaptive grid.
c-----------------------------------------------------------------------
         IF(adapt_grid .AND. solver%istep >= adapt_start)THEN
            IF(grid_solver)THEN
               CALL p2_grid_output
               IF(pt2%istep >= 1)THEN
c-----------------------------------------------------------------------
c     interpolate grid, switch pointers
c-----------------------------------------------------------------------
                  CALL p2_grid_interp(coord%crd_vec)
                  CALL p2_advance_switchsolve("solver")
c-----------------------------------------------------------------------
c     map solution onto new grid
c-----------------------------------------------------------------------
                  CALL p2_grid_interp(pt2%global)
                  IF(solver%step_type=="bdf2" .AND. .NOT. bdf2first)
     $                 CALL p2_grid_interp(stv%ug00)
                  CALL p2_grid_error(solver)
c-----------------------------------------------------------------------
c     if grid-resolution too low, increase gridtol or cut the run
c-----------------------------------------------------------------------
                  IF(solver%griderr > gridtol2)THEN
                     IF(mpi_rank==0)WRITE(*,*)"griderr=",solver%griderr
                     CALL p2_snes_stop("Out of resolution!")
                  ELSEIF(solver%griderr > gridtol/(1.+bel_phifac*.5)
     $                    .AND. solver%griderr > gridtol1
     $                    .AND. gridtol /= gridtol2)THEN
                     gridtol=MIN(MIN(solver%griderr*(1.+bel_phifac*.25),
     $                    gridtol*(1.+bel_phifac*.25)),gridtol2)
                     IF(mpi_rank==0)CALL p2_grid_write_gridtol(gridtol)
                  ENDIF
c-----------------------------------------------------------------------
c     output data
c-----------------------------------------------------------------------
                  IF(monitor)CALL p2_diagnose_output_p(grid,
     $                 outfile_type,"grid",dmout,solver%t,coord%crd_vec)
                  IF(mpi_rank == 0)CALL p2_advance_write_header
c-----------------------------------------------------------------------
c     allocate and compute new matrices.
c-----------------------------------------------------------------------
                  CALL p2_ct_calc(coord)
                  IF(solver%couple_msize > 0)CALL p2_rj_mass(coord,
     $                 solver%stat_qty,solver%theta_fac,solver%mass)
c-----------------------------------------------------------------------
c     produce output to header.txt
c-----------------------------------------------------------------------
                  IF(mpi_rank == 0 .AND. monitor)THEN
                     OPEN(UNIT=bin_unit,FILE=TRIM(filename),
     $                    ACTION="WRITE",STATUS="OLD",POSITION="APPEND")
                     WRITE(bin_unit,'(i8)')solver%istep-restart_step
                     CLOSE(UNIT=bin_unit)
                  ENDIF
                  CYCLE
               ENDIF
c-----------------------------------------------------------------------
c     switch pointers
c-----------------------------------------------------------------------
            ELSEIF(solver%griderr >= gridtol .AND. change_grid)THEN
               change_grid=.FALSE.
               CALL p2_advance_switchsolve("grid")
               CALL p2_grid_reinit
               CYCLE
            ELSEIF(solver%griderr >= gridtol)THEN
               change_grid=.TRUE.
c-----------------------------------------------------------------------
c     if necessary, reduce gridtol.
c-----------------------------------------------------------------------
            ELSEIF(solver%griderr < .2*gridtol 
     $              .AND. gridtol > gridtol1)THEN
               gridtol = MAX(.5*gridtol,gridtol1)
               IF(mpi_rank==0)CALL p2_grid_write_gridtol(gridtol)
            ENDIF
         ENDIF
c-----------------------------------------------------------------------
c     advance.
c-----------------------------------------------------------------------
         pt2%istep=pt2%istep+1
         SELECT CASE(pt2%step_type)
         CASE("theta","theta2")
            CALL p2_snes_step(coord,stv)
         CASE("bdf2")
            IF(bdf2first)THEN
               pt2%step_type="theta"
               pt2%theta=.5
               CALL p2_snes_step(coord,stv)
               CALL VecCopy(stv%ug0,stv%ug00,ierr)
               CALL p2_rj_rhs(pt2%t,coord,pt2%global,stv%rhs)
               pt2%step_type="bdf2"
               pc_reset=.TRUE.
               newjac=.TRUE.
               bdf2first=.FALSE.
            ELSE
               CALL p2_snes_step(coord,stv)
            ENDIF
         CASE DEFAULT
            CALL p2_snes_stop("p2_advance_run: "
     $           //"cannot recognize step_type = "
     $           //TRIM(pt2%step_type)//".")
         END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL p2_snes_stop("Normal termination.")
      END SUBROUTINE p2_advance_run
c-----------------------------------------------------------------------
c     subprogram 2. p2_advance_init.
c     initializes "solver" solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_advance_init(restart_dir,outfile_type,restart_flag,
     $     restart_time,nodal,itmax,itmax_incr,itmax_decr,theta,h,
     $     dtmax,dtmin,tmax,dt_incr,dt_decr,restart_t0,filename)

      CHARACTER(*), INTENT(IN) :: restart_dir,outfile_type
      LOGICAL, INTENT(IN) :: restart_flag,restart_time,nodal
      INTEGER, INTENT(IN) :: itmax,itmax_incr,itmax_decr
      REAL(r8), INTENT(IN) :: theta,h,dtmax,dtmin,tmax,dt_incr,dt_decr,
     $     restart_t0
      CHARACTER(*), INTENT(OUT) :: filename

      LOGICAL :: xper,yper,polar_axis
      INTEGER :: ix,iy,nx,ny,np,iqty,ixmin,ixmax,iymin,iymax
      REAL(r8), DIMENSION(solver%nqty,solver%jxmin:solver%jxmax,
     $     solver%jymin:solver%jymax) :: uu
c-----------------------------------------------------------------------
c     define local variables.
c-----------------------------------------------------------------------
      nx=solver%nx
      ny=solver%ny
      np=solver%np
      ixmin=solver%ixmin
      ixmax=solver%ixmax
      iymin=solver%iymin
      iymax=solver%iymax
      xper=solver%xperiodic
      yper=solver%yperiodic
c-----------------------------------------------------------------------
c     initialize solution and boundary conditions.
c-----------------------------------------------------------------------
      solver%t=0
      CALL job2_init("physics",xper,yper,ixmin,iymin,ixmax,iymax,nx,ny,
     $     np,solver%ground,solver%couple_msize,solver%couple_mlist,
     $     solver%edge_order,solver%edges,uu)
c-----------------------------------------------------------------------
c     set polar_axis flag.
c-----------------------------------------------------------------------
      polar_axis=.FALSE.
      DO iqty=1,solver%nqty
         solver%polar_axis(iqty)
     $        = (solver%edges(1)%bc_type(iqty) == "polar")
         IF(solver%polar_axis(iqty))polar_axis=.TRUE.
      ENDDO
      CALL MPI_Bcast(polar_axis,1,MPI_LOGICAL,0,comm,ierr)
      IF(.NOT.(polar .EQV. polar_axis))
     $     CALL program_stop("Inconsistent polar axis settings.")
c-----------------------------------------------------------------------
c     create native mass matrix and initialize time-step parameters.
c-----------------------------------------------------------------------
      IF(solver%couple_msize > 0)THEN
         ALLOCATE(solver%mass((np+1)**2,solver%couple_msize,
     $        solver%jxmin:solver%jxmax,solver%jymin:solver%jymax))
         solver%mass=0
         solver%h=h
         solver%h_old=h
         solver%theta=theta
      ELSE
         solver%h=1.
         solver%h_old=1.
         solver%theta=1.         
      ENDIF
c-----------------------------------------------------------------------
c     create scatter context for rhs calculation
c-----------------------------------------------------------------------
      CALL p2_sel_scatter(solver,.TRUE.,solver%rscat)
c-----------------------------------------------------------------------
c     store initial solution in PETSc framework.
c-----------------------------------------------------------------------
      CALL p2_sel_set(.TRUE.,uu,solver%jscat,solver%jvec,solver%global)
c-----------------------------------------------------------------------
c     initialize other parameters.
c-----------------------------------------------------------------------
      solver%dtmax=dtmax
      solver%dtmin=dtmin
      solver%tmax=tmax
      solver%itmax=itmax
      solver%itmax_incr=itmax_incr
      solver%itmax_decr=itmax_decr
      solver%dt_incr=dt_incr
      solver%dt_decr=dt_decr
      solver%curve_grid=.TRUE.
c-----------------------------------------------------------------------
c     re-initialize solution and compute grid increments.
c-----------------------------------------------------------------------
      IF(restart_flag)CALL p2_advance_restart(restart_dir,
     $     restart_time,restart_t0,outfile_type)
c-----------------------------------------------------------------------
c     check analytical expressions for jacobian.
c-----------------------------------------------------------------------
      IF(fd_test)THEN
         CALL fd_run(solver)
         CALL program_stop("Termination by fd_run.")
      ENDIF
      IF(fd_bound_test)THEN
         IF(mpi_size > 1)
     $        CALL program_stop("Run fd_bound_test only on one proc!")
         IF(xper .OR. yper)CALL program_stop("Run fd_bound_test only"
     $        //" with non-periodic boundaries!")
         CALL fd_bound_run(solver)
         CALL program_stop("Termination by fd_bound_run.")
      ENDIF
c-----------------------------------------------------------------------
c     write header data.
c-----------------------------------------------------------------------
      filename=TRIM(outdir)//"/header.txt"
      IF(mpi_rank==0 .AND. monitor)THEN
         OPEN(UNIT=bin_unit,FILE=TRIM(filename),ACTION="WRITE",
     $        STATUS="REPLACE")
         WRITE(bin_unit,'(i5,i5,i5,i5,l1)')solver%nqty,nx,ny,np,nodal
         CLOSE(UNIT=bin_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_advance_init
c-----------------------------------------------------------------------
c     subprogram 3. p2_advance_restart.
c     broadcasts input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_advance_restart(restart_dir,restart_time,restart_t0,
     $     outfile_type)

      CHARACTER(*), INTENT(IN) :: restart_dir,outfile_type
      LOGICAL, INTENT(IN) :: restart_time
      REAL(r8), INTENT(IN) :: restart_t0

      CHARACTER(80) :: filename
      LOGICAL :: nodal_old,stop_flag=.FALSE.
      INTEGER :: nqty_old,nx_old,ny_old,np_old
c-----------------------------------------------------------------------
c     read restart header file.
c-----------------------------------------------------------------------
      IF(mpi_rank==0)THEN
         filename=TRIM(restart_dir)//"/header.txt"
         OPEN(UNIT=bin_unit,FILE=TRIM(filename),ACTION="READ",
     $        STATUS="OLD")
         READ(bin_unit,'(i5,i5,i5,i5,l1)')nqty_old,nx_old,ny_old,np_old,
     $        nodal_old
         CLOSE(UNIT=bin_unit)
         IF(solver%nqty /= nqty_old .OR. solver%nx /= nx_old
     $        .OR. solver%ny /= ny_old .OR. solver%np /= np_old)
     $        stop_flag=.TRUE.
      ENDIF
      CALL MPI_Bcast(stop_flag,1,MPI_LOGICAL,0,comm,ierr)
      IF(stop_flag)
     $     CALL program_stop("New and old parameters don't match.")
c-----------------------------------------------------------------------
c     read restart solution.
c-----------------------------------------------------------------------
      CALL p2_diagnose_input_p(solver,restart_dir,outfile_type,"data",
     $     solver%restart_step,solver%t,solver%global)
      IF(restart_time)solver%t=restart_t0
      CALL MPI_Bcast(solver%t,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_advance_restart
c-----------------------------------------------------------------------
c     subprogram 4. p2_advance_switchsolve.
c     switches between main PDE solver and grid solver.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_advance_switchsolve(pde_type)

      CHARACTER(*), INTENT(IN) :: pde_type
c-----------------------------------------------------------------------
c     switch all pointers to the appropriate PDE problem.
c-----------------------------------------------------------------------
      SELECT CASE(pde_type)
      CASE("solver")
         pt2 => solver
         SELECT CASE(solver%solve_type)
         CASE("condense")
            ct_obj => ct_solver
         END SELECT
         stv => st_solver
         grid_solver=.FALSE.
         job_type="physics"
      CASE("grid")
         pt2 => grid
         SELECT CASE(grid%solve_type)
         CASE("condense")
            ct_obj => ct_grid
         END SELECT
         stv => st_grid
         grid_solver=.TRUE.
         job_type="beltrami"
      END SELECT
      
      newjac=.TRUE.
      pc_reset=.TRUE.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_advance_switchsolve
c-----------------------------------------------------------------------
c     subprogram 5. p2_advance_write_step.
c     writes time-step summary to sel.out and standard output.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_advance_write_step

c-----------------------------------------------------------------------
c     write time-step summary.
c-----------------------------------------------------------------------
 10   FORMAT(2i6,2i4,i6,1p,e10.3,3e9.2,e10.3)
      CALL DATE_AND_TIME(date,time,zone,clock)
      wclock0=wclock
      wclock=(clock(5)*60+clock(6))*60+clock(7)+clock(8)/1000._r8

      WRITE(*,10)iout,pt2%istep,pt2%newts,pt2%jacs,pt2%ksps,pt2%t,
     $     pt2%h_old,wclock-wclock0,pt2%griderr,condno
      
      IF(monitor)THEN
         OPEN(UNIT=out_unit,FILE=TRIM(out_name),ACTION="WRITE",
     $        STATUS="OLD",POSITION="APPEND")
         WRITE(out_unit,10)iout,pt2%istep,pt2%newts,pt2%jacs,pt2%ksps,
     $        pt2%t,pt2%h_old,wclock-wclock0,pt2%griderr,condno
         CLOSE(UNIT=out_unit)
      ENDIF

      IF(MOD(iout,10) == 0 .AND. iout > 0 .OR. 
     $     pt2%istep == pt2%nstep .OR. pt2%t >= pt2%tmax)
     $     CALL p2_advance_write_header
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_advance_write_step
c-----------------------------------------------------------------------
c     subprogram 6. p2_advance_write_header.
c     writes head_string to sel.out and standard output.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_advance_write_header

      CHARACTER(80) :: head_string
c-----------------------------------------------------------------------
c     write head_string.
c-----------------------------------------------------------------------
      WRITE(head_string,'(3x,a,4x,a,3(2x,a),5x,a,8x,a,5x,a,3x,a,3x,a)')
     $     "iout","m","it","jac","ksp","t","dt","wclock","griderr",
     $     "condno"

      WRITE(*,*)
      WRITE(*,'(a)')TRIM(head_string)
      WRITE(*,*)

      IF(monitor)THEN
         OPEN(UNIT=out_unit,FILE=TRIM(out_name),ACTION="WRITE",
     $        STATUS="OLD",POSITION="APPEND")
         WRITE(out_unit,*)
         WRITE(out_unit,'(a)')TRIM(head_string)
         WRITE(out_unit,*)
         CLOSE(UNIT=out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_advance_write_header
      END MODULE p2_advance_mod
