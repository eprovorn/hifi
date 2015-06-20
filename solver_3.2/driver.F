c-----------------------------------------------------------------------
c     file driver.F.
c     spectral element code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. driver_mod.
c     1. driver_run.
c-----------------------------------------------------------------------
c     subprogram 0. driver_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE driver_mod
      USE p2_advance_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. driver_run.
c     real main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE driver_run

      LOGICAL :: nodal=.FALSE.,restart_flag=.FALSE.,xperiodic=.FALSE.,
     $     yperiodic=.FALSE.,restart_time=.FALSE.
      CHARACTER(3) :: quad_type="gl0"
      CHARACTER(4) :: outfile_type="bin"
      CHARACTER(8) :: step_type="theta",solve_type="condense"
      CHARACTER(80) :: restart_dir=".",physics_type="physics"
      CHARACTER(90) :: message1=".",message2=".",message3="."
      INTEGER :: itmax=19,itmax_incr=4,itmax_decr=6,np=1,nq=0,nx=1,ny=1,
     $     nbx=1,nstep=0,nqty=1,restart_step=0,dmout=1,myios=0,
     $     exit_flag=0,i,nqty_schur=0
      REAL(r8) :: tmax=0.,dt=0.,theta=.5,dtmax=1.,dtmin=0.,
     $     restart_t0=0.,dt_incr=1.2,dt_decr=.8

      NAMELIST/algorithm_input/monitor,nodal,quad_type,grid_inv_type,
     $     solve_type,fd_test,polar,adapt_grid,adapt_dt,itmax,
     $     itmax_incr,itmax_decr,dt_incr,dt_decr,step_type,theta,
     $     errtol,ksp_restart,ksps_max,always_pc_reset,parallel_write,
     $     parallel_read,du_diagnose,grid_type,outfile_type,
     $     fd_bound_test,schur_solve_type,mass_solve_type
      NAMELIST/universal_input/nx,ny,np,nq,nbx,xperiodic,yperiodic,
     $     dt,dtmax,dtmin,tmax,nstep,dmout,outdir,restart_flag,
     $     restart_dir,restart_step,read_grid,grid_step,restart_time,
     $     restart_t0,cubit_file
c-----------------------------------------------------------------------
c     Initialize petsc and mpi.
c-----------------------------------------------------------------------
      CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      comm=PETSC_COMM_WORLD
      CALL MPI_Comm_size(comm,mpi_size,ierr)
      CALL MPI_Comm_rank(comm,mpi_rank,ierr)
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0)THEN
         myios=0
         OPEN(UNIT=in_unit,FILE="./sel.in",ACTION="READ",STATUS="OLD")
c-----------------------------------------------------------------------
c     set general input parameters.
c-----------------------------------------------------------------------
         READ(in_unit,NML=algorithm_input,IOSTAT=myios)
         IF(myios /= 0)exit_flag=97
         READ(in_unit,NML=universal_input,IOSTAT=myios)
         IF(myios /= 0)exit_flag=98
c-----------------------------------------------------------------------
c     set job-specific input parameters.
c-----------------------------------------------------------------------
         CALL job2_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $        nqty_schur,dt,dtmax,tmax,nstep,du_diagnose,
     $        physics_type,exit_flag)
         CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     preset some initial parameters.
c-----------------------------------------------------------------------
         nbx=MIN(MAX(nbx,1),mpi_size)
         np=MAX(np,1)
         IF(np < 2)solve_type="full"
         dmout=MAX(dmout,1)
         IF(quad_type=="gl0")THEN
            nq=MAX(nq,np)
            IF(adapt_grid)nq=MAX(nq,np+2)
         ELSEIF(quad_type=="gll")THEN
            nq=MAX(nq,np+1)
            IF(adapt_grid)nq=MAX(nq,np+3)
         ENDIF
         IF(.NOT. restart_flag)restart_step=0
         IF(restart_flag)grid_type="bin"
         IF(restart_flag .AND. outdir==restart_dir)exit_flag=95
         IF((parallel_write .OR. parallel_read) 
     $        .AND. outfile_type /= "hdf5")exit_flag=96
c-----------------------------------------------------------------------
c     write message to screen and output file.
c-----------------------------------------------------------------------
         WRITE(message1,'(a,i6,a,i4)')" physics = "//TRIM(physics_type)
     $        //", step_type = "//TRIM(step_type)
     $        //", nproc = ",mpi_size,", nbx = ",nbx
         WRITE(message2,'(2(a,l1),2(a,i4),2(a,i2))')
     $        " xper = ",xperiodic,", yper = ",yperiodic,
     $        ", nx = ",nx,", ny = ",ny,", np =",np,", nq =",nq
         IF(solve_type == "schur")THEN
            WRITE(message3,'(a)')
     $           " solve_type = "//TRIM(solve_type)
     $           //", schur_solve_type = "//TRIM(schur_solve_type)
     $           //", mass_solve_type = "//TRIM(mass_solve_type)
         ELSE
            WRITE(message3,'(a)')
     $           " solve_type = "//TRIM(solve_type)
         ENDIF
         WRITE(*,'(a)')TRIM(message1)
         WRITE(*,'(a)')TRIM(message2)
         WRITE(*,'(a)')TRIM(message3)
         CALL system("mkdir -p "//TRIM(outdir))
         IF(outdir /= ".")
     $        CALL system("cp -p go sel.in beltrami.in "//TRIM(outdir))
         out_name=TRIM(outdir)//"/sel.out"
         CALL timer(0,out_name,out_unit)
         IF(monitor)THEN
            OPEN(UNIT=out_unit,FILE=TRIM(out_name),ACTION="WRITE",
     $           STATUS="UNKNOWN",POSITION="APPEND")
            WRITE(out_unit,'(72(a1))')("-",i=1,72)
            WRITE(out_unit,'(a)')TRIM(message1)
            WRITE(out_unit,'(a)')TRIM(message2)
            WRITE(out_unit,'(a)')TRIM(message3)
            CLOSE(UNIT=out_unit)
         ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     exit if exit_flag is not zero.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(exit_flag,1,MPI_INTEGER,0,comm,ierr)
      SELECT CASE(exit_flag)
      CASE(1)
         CALL program_stop("No Multiblock Support Yet")
      CASE(2)
         CALL program_stop("2D Data Only")
      CASE(95)
         CALL program_stop("Not allowed to output data to restart 
     $        directory")
      CASE(96)
         CALL program_stop("parallel I/O works only with 
     $        outfile_type='hdf5'.")
      CASE(97)
         CALL program_stop("algorithm_input is incorrect.")
      CASE(98)
         CALL program_stop("universal_input is incorrect.")
      CASE(99)
         CALL program_stop("physics_list is incorrect.")
      END SELECT
c-----------------------------------------------------------------------
c     broadcast character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(quad_type,3,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(outfile_type,4,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(grid_inv_type,8,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(grid_type,8,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(step_type,8,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(solve_type,8,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(schur_solve_type,8,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(mass_solve_type,8,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(outdir,80,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(restart_dir,80,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast logical input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(parallel_write,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(parallel_read,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(nodal,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(xperiodic,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(yperiodic,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(polar,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(monitor,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(restart_flag,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(fd_test,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(fd_bound_test,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(adapt_dt,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(adapt_grid,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(read_grid,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(always_pc_reset,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(du_diagnose,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(restart_time,1,MPI_LOGICAL,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast integer input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(nx,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(ny,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(nbx,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(np,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(nq,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(nqty,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(nqty_schur,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(nstep,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(dmout,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(restart_step,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(grid_step,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(itmax,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(itmax_incr,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(itmax_decr,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(ksp_restart,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(ksps_max,1,MPI_INTEGER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(tmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(dtmax,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(dtmin,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(theta,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(errtol,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(dt_incr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(dt_decr,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(restart_t0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     allocate job2_mod module objects
c-----------------------------------------------------------------------
      CALL job2_alloc(nodal,np)
c-----------------------------------------------------------------------
c     run spectral element simulation.
c-----------------------------------------------------------------------
      CALL p2_advance_run(nx,ny,nbx,np,nq,nqty,nqty_schur,nstep,tmax,
     $     dmout,dt,dtmax,dtmin,theta,nodal,xperiodic,yperiodic,
     $     restart_flag,restart_dir,restart_step,quad_type,step_type,
     $     itmax,itmax_incr,itmax_decr,dt_incr,dt_decr,solve_type,
     $     restart_time,restart_t0,outfile_type)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE driver_run
      END MODULE driver_mod
