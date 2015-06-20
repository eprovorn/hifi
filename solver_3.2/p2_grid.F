c-----------------------------------------------------------------------
c     file p2_grid.F.
c     spatial discretization of 2D PDEs with spectral element method.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. p2_grid_mod.
c     1. p2_grid_error.
c     2. p2_grid_init.
c     3. p2_grid_output.
c     4. p2_grid_interp.
c     5. p2_grid_gmatfit.
c     6. p2_grid_reinit.
c     7. p2_grid_write_gridtol.
c     8. p2_grid_write_header.
c-----------------------------------------------------------------------
c     subprogram 0. p2_grid_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE p2_grid_mod
      USE p2_snes_mod
      IMPLICIT NONE

#include "finclude/petscdef.h"

      CHARACTER(80), PRIVATE :: grid_name,grid_message
      REAL(r8), PRIVATE :: g_wclock=0,g_wclock0=0
      
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. p2_grid_error
c     computes grid error.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_grid_error(psv)

      TYPE(p2_sel_type) :: psv

      INTEGER :: ix,jx,iy,jy,nqty,np,iqty,ixmin,ixmax,iymin,iymax
      REAL(r8) :: griderr,maxnorm
      REAL(r8), DIMENSION(0:psv%np,0:psv%np) :: u,delu,mass0
      REAL(r8), DIMENSION(psv%nqty,psv%jxmin:psv%jxmax,
     $     psv%jymin:psv%jymax) :: uu
      REAL(r8), DIMENSION(psv%nqty,psv%nx,psv%ny) :: norm,errx,erry
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      nqty=psv%nqty
      np=psv%np
      ixmin=psv%ixmin
      ixmax=psv%ixmax
      iymin=psv%iymin
      iymax=psv%iymax
      mass0=psv%mass0
      CALL p2_sel_get_local(psv%jscat,psv%jvec,psv%global,uu)
c-----------------------------------------------------------------------
c     compute norm and error.
c-----------------------------------------------------------------------
      norm=0
      errx=0
      erry=0
      jy=psv%jymin
      DO iy=iymin+1,iymax
         jx=psv%jxmin
         DO ix=ixmin+1,ixmax
            DO iqty=1,nqty
               IF(.NOT. adapt_qty(iqty))CYCLE
               u=uu(iqty,jx:jx+np,jy:jy+np)
               IF(psv%nodal)
     $              u=MATMUL(MATMUL(psv%tinv,u),TRANSPOSE(psv%tinv))
               norm(iqty,ix,iy)=SUM(u*MATMUL(MATMUL(mass0,u),mass0))
               delu=0
               delu(np-1,:)=u(np-1,:)
               errx(iqty,ix,iy)=SUM(delu*
     $              MATMUL(MATMUL(mass0,delu),mass0))
               delu=0
               delu(:,np-1)=u(:,np-1)
               erry(iqty,ix,iy)=SUM(delu*
     $              MATMUL(MATMUL(mass0,delu),mass0))
            ENDDO
            jx=jx+np+1
         ENDDO
         jy=jy+np+1
      ENDDO
c-----------------------------------------------------------------------
c     normalize and collect error over the domain.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         CALL MPI_Allreduce(MAXVAL(norm(iqty,:,:)),maxnorm,1,
     $        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)                  
         WHERE(norm(iqty,:,:) <= maxnorm*min_eps**.25)
            norm(iqty,:,:) = one/min_eps
         END WHERE
      ENDDO
      WHERE(norm <= min_eps)
         norm = one/min_eps
      END WHERE
      
      griderr=MAX(MAXVAL(SQRT(ABS(errx/norm))),
     $     MAXVAL(SQRT(ABS(erry/norm))))
      
      CALL MPI_Allreduce(griderr,psv%griderr,1,
     $     MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_grid_error
c-----------------------------------------------------------------------
c     subprogram 2. p2_grid_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_grid_init(nodal,adapt_start,gridtol1,gridtol2)

      LOGICAL, INTENT(IN) :: nodal 
      INTEGER, INTENT(OUT) :: adapt_start
      REAL(r8), INTENT(OUT) :: gridtol1,gridtol2

      CHARACTER(8) :: solve_type="condense"
      CHARACTER(80) :: dat_name
      INTEGER :: nx,ny,np,itmax=50,itmax_decr=10,ix,iy,myios
      REAL(r8), DIMENSION(grid%nqty,grid%jxmin:grid%jxmax,
     $     grid%jymin:grid%jymax) :: uu
      REAL(r8), DIMENSION(solver%nqty,0:solver%nxp,0:solver%nyp) :: 
     $     uglobal      

      NAMELIST/beltrami_input/bel_diagnose,bel_phifac,itmax,itmax_decr,
     $     solve_type,adapt_start,gridtol1,gridtol2
c-----------------------------------------------------------------------
c     read input.
c-----------------------------------------------------------------------
      myios=0
      adapt_start=0
      gridtol1=1.
      gridtol2=1.
      IF(mpi_rank == 0)THEN
         g_wclock=wclock
         OPEN(UNIT=belin_unit,FILE="beltrami.in",ACTION="READ",
     $        STATUS="OLD")
         READ(belin_unit,NML=beltrami_input,IOSTAT=myios)
         CLOSE(UNIT=belin_unit)
      ENDIF
      CALL MPI_Bcast(myios,1,MPI_INTEGER,0,comm,ierr)      
      IF(myios /= 0)CALL program_stop("beltrami_input incorrect")
      IF(grid%np < 2)solve_type="full"
c-----------------------------------------------------------------------
c     broadcast input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(solve_type,8,MPI_CHARACTER,0,comm,ierr)
      CALL MPI_Bcast(bel_diagnose,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(itmax,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(itmax_decr,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(adapt_start,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(bel_phifac,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(gridtol1,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(gridtol2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     transfer input values to appropriate global variables.
c-----------------------------------------------------------------------
      grid%itmax=itmax
      grid%itmax_decr=itmax_decr
      IF(grid%solve_type == "full" .AND. solve_type /= "full")
     $     CALL program_stop
     $     ("Cannot overwrite solve_type = full in p2_grid_init")
      grid%solve_type=solve_type
c-----------------------------------------------------------------------
c     initialize global grid size
c-----------------------------------------------------------------------
      nx=grid%nx
      ny=grid%ny
      np=grid%np
c-----------------------------------------------------------------------
c     initialize beltrami gmat tensor from the solution vector.
c-----------------------------------------------------------------------
      CALL p2_sel_get_global(solver,solver%global,uglobal)
      CALL beltrami_getdata(solver%xperiodic,solver%yperiodic,
     $     solver%nx,solver%ny,solver%np,solver%mass0,.TRUE.)
      IF(solver%nodal)THEN
         CALL beltrami_griderr(.TRUE.,uglobal,solver%tinv)
      ELSE
         CALL beltrami_griderr(.TRUE.,uglobal)
      ENDIF
      CALL p2_grid_gmatfit(.TRUE.)
c-----------------------------------------------------------------------
c     initialize solution and boundary conditions.
c-----------------------------------------------------------------------
      CALL job2_init("beltrami",grid%xperiodic,grid%yperiodic,
     $     grid%ixmin,grid%iymin,grid%ixmax,grid%iymax,nx,ny,np,
     $     grid%ground,grid%couple_msize,grid%couple_mlist,
     $     grid%edge_order,grid%edges,uu)
c-----------------------------------------------------------------------
c     initialize other scalars.
c-----------------------------------------------------------------------
      grid%t=0
      grid%nstep=1
      grid%itmax_incr=1
      grid%dt_incr=1.
      grid%dt_decr=1.
      grid%h=1.
      grid%h_old=1.
      grid%dtmax=1.
      grid%dtmin=0.
      grid%tmax=1.
      grid%istep=0
      grid%theta=half
      grid%curve_grid=.FALSE.
      grid%polar_axis=.FALSE.
      grid%stat_qty=.TRUE.
c-----------------------------------------------------------------------
c     create scatter context for rhs calculation
c-----------------------------------------------------------------------
      CALL p2_sel_scatter(grid,.TRUE.,grid%rscat)
c-----------------------------------------------------------------------
c     store initial solution in PETSc framework.
c-----------------------------------------------------------------------
      CALL p2_sel_set(.TRUE.,uu,grid%jscat,grid%jvec,grid%global)
c-----------------------------------------------------------------------
c     write header data.
c-----------------------------------------------------------------------
      dat_name=TRIM(outdir)//"/grid.txt"
      IF(mpi_rank==0 .AND. monitor)THEN
         OPEN(UNIT=belbin_unit,FILE=TRIM(dat_name),ACTION="WRITE",
     $        STATUS="REPLACE")
         WRITE(belbin_unit,'(i5,i5,i5,i5,l1)')ndim,nx,ny,np,nodal
         CLOSE(UNIT=belbin_unit)
      ENDIF
      grid_name=TRIM(outdir)//"/grid.out"
      WRITE(grid_message,'(2(a,l1),a)')
     $     "job_type = beltrami, xper = ",grid%xperiodic,
     $     ", yper = ",grid%yperiodic,
     $     ", solve_type = "//TRIM(grid%solve_type)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_grid_init
c-----------------------------------------------------------------------
c     subprogram 3. p2_grid_output.
c     produces text output as new grids are calculated.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_grid_output

c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/2x,"iout",1x,"it",2x,"jac",2x,"ksp",7x,"wclock"/)
 20   FORMAT(i5,2i4,i6,e11.2)
c-----------------------------------------------------------------------
c     monitor solution.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0)THEN
         CALL DATE_AND_TIME(date,time,zone,clock)
         g_wclock0=g_wclock
         g_wclock=(clock(5)*60+clock(6))*60+clock(7)+clock(8)/1000._r8

         IF(grid%istep==0)WRITE(*,10)
         WRITE(*,20)grid%istep,grid%newts,grid%jacs,grid%ksps,
     $        g_wclock-g_wclock0

         IF(monitor)THEN
            OPEN(UNIT=grid_unit,FILE=TRIM(grid_name),ACTION="WRITE",
     $           STATUS="OLD",POSITION="APPEND")
            IF(grid%istep==0)WRITE(grid_unit,10)
            WRITE(grid_unit,20)grid%istep,grid%newts,grid%jacs,
     $           grid%ksps,g_wclock-g_wclock0
            CLOSE(UNIT=grid_unit)
         ENDIF
      ENDIF

      grid%ksps=0
      grid%jacs=0
      grid%newts=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_grid_output
c-----------------------------------------------------------------------
c     subprogram 4. p2_grid_interp.
c     interpolate variables from old logical onto new logical grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_grid_interp(global)

      Vec :: global

      INTEGER :: info,iqty,ix,iy,iz1,iz2,jxx,jyy,np,
     $     ixmin,ixmax,iymin,iymax,nqty,j
      INTEGER, DIMENSION(0:pt2%np) :: bpiv
      REAL(r8), DIMENSION(0:pt2%np) :: xpos,ypos
      REAL(r8), DIMENSION(0:pt2%np,0:pt2%np) :: bmat
      REAL(r8), DIMENSION(ndim,0:grid%np,0:grid%np) :: uigrid
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%np,0:pt2%np) :: ui
      REAL(r8), DIMENSION(pt2%nqty,pt2%jxmin:pt2%jxmax,
     $     pt2%jymin:pt2%jymax) :: u
      REAL(r8), DIMENSION(ndim,0:pt2%nxp,0:pt2%nyp) :: gridglobal
      REAL(r8), DIMENSION(pt2%nqty,0:pt2%nxp,0:pt2%nyp) :: uglobal
c-----------------------------------------------------------------------
c     define local parameters.
c-----------------------------------------------------------------------
      np=pt2%np
      nqty=pt2%nqty
      ixmin=pt2%ixmin
      ixmax=pt2%ixmax
      iymin=pt2%iymin
      iymax=pt2%iymax
      bpiv = pt2%bpiv
      bmat = pt2%bmat
c-----------------------------------------------------------------------
c     get global solution.
c-----------------------------------------------------------------------
      CALL p2_sel_get_global(pt2,global,uglobal)
c-----------------------------------------------------------------------
c     get global new -> old grid mapping.
c-----------------------------------------------------------------------
      CALL p2_sel_get_global(grid,grid%global,gridglobal)
c-----------------------------------------------------------------------
c     interpolate solution to the new grid.
c-----------------------------------------------------------------------
      u=0
      jyy=pt2%jymin
      DO iy=iymin,iymax-1
         ypos = (iy + (one+pt2%basis%qzero)*half)/REAL(pt2%ny,r8)
         jxx=pt2%jxmin
         DO ix=ixmin,ixmax-1
            xpos = (ix + (one+pt2%basis%qzero)*half)/REAL(pt2%nx,r8)
            DO iz2=0,np
               DO iz1=0,np
                  CALL p2_sel_interp(xpos(iz1),ypos(iz2),
     $                 gridglobal,uigrid(:,iz1,iz2))
                  CALL p2_sel_interp(uigrid(1,iz1,iz2),
     $                 uigrid(2,iz1,iz2),uglobal,ui(:,iz1,iz2))
               ENDDO
            ENDDO
            DO iqty=1,nqty
               DO j=1,2
                  CALL dgetrs('N',np+1,np+1,bmat,np+1,bpiv,ui(iqty,:,:),
     $                 np+1,info)
                  ui(iqty,:,:)=TRANSPOSE(ui(iqty,:,:))
               ENDDO
            ENDDO
            u(:,jxx:jxx+np,jyy:jyy+np)=ui(:,0:np,0:np)
            jxx=jxx+np+1
         ENDDO
         jyy=jyy+np+1
      ENDDO
c-----------------------------------------------------------------------
c     restore global solution.
c-----------------------------------------------------------------------
      CALL p2_sel_set(.TRUE.,u,pt2%jscat,pt2%jvec,global)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_grid_interp
c-----------------------------------------------------------------------
c     subprogram 5. p2_grid_gmatfit.
c     construct the g tensor for beltrami equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_grid_gmatfit(alloc)

      LOGICAL, INTENT(IN) :: alloc

      INTEGER :: ix,iy,i,j,iqty,jx,jy,nx,ny,np,info
      INTEGER, DIMENSION(0:grid%np) :: bpiv
      REAL(r8), DIMENSION(0:grid%np) :: xi
      REAL(r8), DIMENSION(0:grid%np,0:grid%np) :: bmat
      REAL(r8), DIMENSION(ndim,0:grid%np,0:grid%np) :: u_xy
      REAL(r8), DIMENSION(ndim**2,0:grid%np,0:grid%np) :: ui
      REAL(r8), DIMENSION(ndim**2,0:grid%nxp,0:grid%nyp) :: gfit
c-----------------------------------------------------------------------
c     copy variables.
c-----------------------------------------------------------------------
      nx=grid%nx
      ny=grid%ny
      np=grid%np
      xi=grid%basis%qzero
      bpiv=grid%bpiv
      bmat=grid%bmat
c-----------------------------------------------------------------------
c     start loop over finite elements.
c-----------------------------------------------------------------------
      ui=0
      jy=0
      DO iy=0,ny-1
         jx=0
         DO ix=0,nx-1
            DO i=0,np
               u_xy(1,:,i) = (ix + (1+xi)*half)/REAL(nx,r8)
               u_xy(2,i,:) = (iy + (1+xi)*half)/REAL(ny,r8)
            ENDDO
c-----------------------------------------------------------------------
c     solve for basis amplitudes of gfit.
c-----------------------------------------------------------------------
            CALL beltrami_gmat1(u_xy(1,:,:),u_xy(2,:,:),ui)
            DO iqty=1,ndim**2
               IF(SUM(ui(iqty,:,:)**2)==0)CYCLE
               DO j=1,2
                  CALL dgetrs('N',np+1,np+1,bmat,np+1,bpiv,
     $                 ui(iqty,:,:),np+1,info)
                  ui(iqty,:,:)=TRANSPOSE(ui(iqty,:,:))
               ENDDO
            ENDDO
            gfit(:,jx:jx+np,jy:jy+np)=ui
c-----------------------------------------------------------------------
c     finish loop over finite elements.
c-----------------------------------------------------------------------
            jx=jx+np+1
         ENDDO
         jy=jy+np+1
      ENDDO
      CALL beltrami_gmat_fit(gfit,alloc)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_grid_gmatfit
c-----------------------------------------------------------------------
c     subprogram 6. p2_grid_reinit.
c     reinitializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_grid_reinit

      REAL(r8), DIMENSION(ndim,grid%jxmin:grid%jxmax,
     $     grid%jymin:grid%jymax) :: uu
      REAL(r8), DIMENSION(solver%nqty,0:solver%nxp,0:solver%nyp) :: 
     $     uglobal,uglobal0
c-----------------------------------------------------------------------
c     reinitialze timer.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0)THEN
         CALL DATE_AND_TIME(date,time,zone,clock)
         g_wclock=(clock(5)*60+clock(6))*60+clock(7)+clock(8)/1000._r8
      ENDIF
c-----------------------------------------------------------------------
c     initialize beltrami gmat tensor from the solution vector.
c-----------------------------------------------------------------------
      CALL p2_sel_get_global(solver,solver%global,uglobal)
      CALL beltrami_getdata(solver%xperiodic,solver%yperiodic,
     $     solver%nx,solver%ny,solver%np,solver%mass0,.FALSE.)

      SELECT CASE(solver%step_type)
      CASE("bdf2")
         CALL p2_sel_get_global(solver,st_solver%ug00,uglobal0)
      CASE("theta","theta2")
         CALL p2_sel_get_global(solver,st_solver%ug0,uglobal0)
      END SELECT
      uglobal0=uglobal0+uglobal
      IF(solver%nodal)THEN
         CALL beltrami_griderr(.FALSE.,uglobal0,solver%tinv)
      ELSE
         CALL beltrami_griderr(.FALSE.,uglobal0)
      ENDIF
      CALL p2_grid_gmatfit(.FALSE.)
c-----------------------------------------------------------------------
c     initialize solution and modify positions.
c-----------------------------------------------------------------------
      CALL job2_init("beltrami",grid%xperiodic,grid%yperiodic,
     $     grid%ixmin,grid%iymin,grid%ixmax,grid%iymax,grid%nx,grid%ny,
     $     grid%np,grid%ground,grid%couple_msize,grid%couple_mlist,
     $     grid%edge_order,grid%edges,uu)
c-----------------------------------------------------------------------
c     store initial solution in PETSc framework.
c-----------------------------------------------------------------------
      CALL p2_sel_set(.TRUE.,uu,grid%jscat,grid%jvec,grid%global)
c-----------------------------------------------------------------------
c     reinitialize scalars.
c-----------------------------------------------------------------------
      grid%istep=0
      grid%t=0
      grid%h=1.
      grid%h_old=1.
c-----------------------------------------------------------------------
c     write header for "grid" output.
c-----------------------------------------------------------------------
      IF(mpi_rank == 0)CALL p2_grid_write_header
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_grid_reinit
c-----------------------------------------------------------------------
c     subprogram 7. p2_grid_write_gridtol.
c     writes a message whenever gridtol is changed.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_grid_write_gridtol(gridtol)

      REAL(r8), INTENT(IN) :: gridtol
c-----------------------------------------------------------------------
c     write message.
c-----------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,'(a,1p,e9.2)')"Changed gridtol to ",gridtol
      WRITE(*,*)
      
      IF(monitor)THEN
         OPEN(UNIT=out_unit,FILE=TRIM(out_name),ACTION="WRITE",
     $        STATUS="OLD",POSITION="APPEND")
         WRITE(out_unit,*)
         WRITE(out_unit,'(a,1p,e9.2)')"Changed gridtol to ",gridtol
         WRITE(out_unit,*)
         CLOSE(UNIT=out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_grid_write_gridtol
c-----------------------------------------------------------------------
c     subprogram 8. p2_grid_write_header.
c     writes a header whenever switching to "grid" solver.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE p2_grid_write_header

c-----------------------------------------------------------------------
c     write grid header.
c-----------------------------------------------------------------------
      WRITE(*,*)
      WRITE(*,'(a)')TRIM(grid_message)
      
      IF(monitor)THEN
         OPEN(UNIT=grid_unit,FILE=TRIM(grid_name),ACTION=
     $        "WRITE",STATUS="UNKNOWN",POSITION="APPEND")
         WRITE(grid_unit,*)
         WRITE(grid_unit,'(a)')TRIM(grid_message)
         CLOSE(UNIT=grid_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p2_grid_write_header
      END MODULE p2_grid_mod

