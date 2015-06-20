c-----------------------------------------------------------------------
c     file extra.f.
c     reads external equilibrium data, coil configurations, etc..
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. extra_mod.
c     1. extra_read_barnes.
c     2. extra_read_marklin.
c     3. extra_interp.
c     4. extra_lininterp.
c     5. extra_xysearch.
c     6. extra_equil_dealloc.
c     7. extra_coilalloc.
c     8. extra_coildealloc.
c     9. extra_coileval.
c-----------------------------------------------------------------------
c     subprogram 0. extra_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE extra_mod
      USE bicube_mod
      IMPLICIT NONE

      TYPE :: coil_type
      INTEGER :: ncoils
      REAL(r8) :: lam
      REAL(r8), DIMENSION(:), POINTER :: zb,ze,tfire,tqtr,tcrow,vc0
      END TYPE coil_type

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. extra_read_barnes.
c     reads equilibrium data from an equilibrium input file.
c     input file is from a code by Dan Barnes.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extra_read_barnes(filename,interp,
     $     equil_bc,equil,equilxy)

      CHARACTER(*), INTENT(IN) :: filename,interp
      TYPE(bicube_type) :: equil_bc
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: equil,equilxy

      INTEGER, PARAMETER :: nqty=3
      INTEGER :: nr,nz,i,iqty,ierr,myios
      REAL(r8), PARAMETER :: prs=0.5_r8
      REAL(r8) :: rw,zl
c-----------------------------------------------------------------------
c     open data file and read, broadcast, and normalize scalars.
c-----------------------------------------------------------------------
      IF(mpi_rank==0)OPEN(UNIT=equil_unit,FILE=TRIM(filename),
     $     ACTION="READ",STATUS="OLD",IOSTAT=myios)
      CALL MPI_Bcast(myios,1,MPI_INTEGER,0,comm,ierr)
      IF(myios /= 0)
     $     CALL program_stop("Error reading Barnes equilibrium")
      IF(mpi_rank==0)THEN
         READ(equil_unit,*)nr,nz,rw,zl
c-----------------------------------------------------------------------
c        modify scalars (arrays are zero-based).
c-----------------------------------------------------------------------
         nr=nr-1
         nz=nz-1
      ENDIF
      CALL MPI_Bcast(nr,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(nz,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(rw,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(zl,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      
      SELECT CASE(interp)
      CASE("bicube")
         CALL bicube_alloc(equil_bc,nz,nr,nqty)
c-----------------------------------------------------------------------
c        read 2D equilibrium data.
c-----------------------------------------------------------------------
         IF(mpi_rank==0)THEN
            DO iqty=1,nqty
               READ(equil_unit,*)equil_bc%fs(:,:,iqty)
            ENDDO
            CLOSE(equil_unit)
c-----------------------------------------------------------------------
c           scale variables.
c-----------------------------------------------------------------------
            equil_bc%xs=(/(i,i=0,nz)/)*(zl/rw)/REAL(nz,r8)
            equil_bc%ys=(/(i,i=0,nr)/)/REAL(nr,r8)
            equil_bc%fs(:,:,2)=equil_bc%fs(:,:,2)*prs
            DO i=0,nr
               equil_bc%fs(:,i,3)=equil_bc%fs(:,i,3)
     $              *equil_bc%ys(i)*half
            ENDDO
         ENDIF
c-----------------------------------------------------------------------
c       broadcast equilibrium data.
c-----------------------------------------------------------------------
         CALL MPI_Bcast(equil_bc%fs,nqty*(nr+1)*(nz+1),
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
         CALL MPI_Bcast(equil_bc%xs,nz+1,
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
         CALL MPI_Bcast(equil_bc%ys,nr+1,
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c       fit to bicubic splines.
c-----------------------------------------------------------------------
         CALL bicube_fit(equil_bc,"extrap","extrap")
      CASE("linear")
         ALLOCATE(equil(nqty,0:nz,0:nr),equilxy(2,0:nz,0:nr))
c-----------------------------------------------------------------------
c        read 2D equilibrium data.
c-----------------------------------------------------------------------
         IF(mpi_rank==0)THEN
            DO iqty=1,nqty
               READ(equil_unit,*)equil(iqty,:,:)
            ENDDO
            CLOSE(equil_unit)
c-----------------------------------------------------------------------
c           compute physical coordinates.
c-----------------------------------------------------------------------
            DO i=0,nz
               equilxy(1,i,:)=i*(zl/rw)/REAL(nz,r8)
            ENDDO
            DO i=0,nr
               equilxy(2,:,i)=i/REAL(nr,r8)
            ENDDO
c-----------------------------------------------------------------------
c           scale variables.
c-----------------------------------------------------------------------
            equil(2,:,:)=equil(2,:,:)*prs
            equil(3,:,:)=equil(3,:,:)*equilxy(2,:,:)*half
         ENDIF
c-----------------------------------------------------------------------
c       broadcast equilibrium data.
c-----------------------------------------------------------------------
         CALL MPI_Bcast(equil,nqty*(nr+1)*(nz+1),
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
         CALL MPI_Bcast(equilxy,2*(nz+1)*(nr+1),
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
      END SELECT
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_read_barnes
c-----------------------------------------------------------------------
c     subprogram 2. extra_read_marklin.
c     reads equilibrium data from an equilibrium input file.
c     input file is from a code by George Marklin.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extra_read_marklin(filename,interp,pmax,
     $     equil_bc,equil,equilxy)

      CHARACTER(*), INTENT(IN) :: filename,interp
      REAL(r8), INTENT(IN) :: pmax
      TYPE(bicube_type) :: equil_bc
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: equil,equilxy

      CHARACTER(8) :: aux
      INTEGER, PARAMETER :: nqty=3
      INTEGER :: nr,nz,ir,iz,myios
      REAL(r8) :: rw,zl,flx,prs,cur,fac
c-----------------------------------------------------------------------
c     read data file header.
c-----------------------------------------------------------------------
      IF(mpi_rank==0)OPEN(UNIT=equil_unit,FILE=TRIM(filename),
     $     ACTION="READ",STATUS="OLD",IOSTAT=myios)
      CALL MPI_Bcast(myios,1,MPI_INTEGER,0,comm,ierr)
      IF(myios /= 0)
     $     CALL program_stop("Error reading Marklin equilibrium")
      IF(mpi_rank==0)THEN
         READ(equil_unit,*)
         READ(equil_unit,*)aux,nr,nz
         nr=nr-1
         nz=nz-1
      ENDIF
      CALL MPI_Bcast(nr,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(nz,1,MPI_INTEGER,0,comm,ierr)

      SELECT CASE(interp)
      CASE("bicube")
         CALL bicube_alloc(equil_bc,nz,nr,nqty)
c-----------------------------------------------------------------------
c        read 2D equilibrium data.
c-----------------------------------------------------------------------
         IF(mpi_rank==0)THEN
            DO iz=0,nz
               DO ir=0,nr
                  READ(equil_unit,*)equil_bc%ys(ir),equil_bc%xs(iz),
     $                 equil_bc%fs(iz,ir,1),equil_bc%fs(iz,ir,2),
     $                 equil_bc%fs(iz,ir,3)
               ENDDO
            ENDDO
            CLOSE(equil_unit)
c-----------------------------------------------------------------------
c     scale variables.
c-----------------------------------------------------------------------
            rw=MAXVAL(equil_bc%ys)
            fac=pmax/MAXVAL(equil_bc%fs(:,:,2))
            flx=SQRT(fac)
            prs=fac
            cur=SQRT(fac)
            equil_bc%xs=equil_bc%xs/rw
            equil_bc%ys=equil_bc%ys/rw
            equil_bc%fs(:,:,1)=equil_bc%fs(:,:,1)*flx
            equil_bc%fs(:,:,2)=equil_bc%fs(:,:,2)*prs
            equil_bc%fs(:,:,3)=equil_bc%fs(:,:,3)*cur
         ENDIF
c-----------------------------------------------------------------------
c       broadcast equilibrium data.
c-----------------------------------------------------------------------
         CALL MPI_Bcast(equil_bc%fs,nqty*(nr+1)*(nz+1),
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
         CALL MPI_Bcast(equil_bc%xs,nz+1,
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
         CALL MPI_Bcast(equil_bc%ys,nr+1,
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c       fit to bicubic splines.
c-----------------------------------------------------------------------
         CALL bicube_fit(equil_bc,"extrap","extrap")
      CASE("linear")
         ALLOCATE(equil(nqty,0:nz,0:nr),equilxy(2,0:nz,0:nr))
c-----------------------------------------------------------------------
c        read 2D equilibrium data.
c-----------------------------------------------------------------------
         IF(mpi_rank==0)THEN
            DO iz=0,nz
               DO ir=0,nr
                  READ(equil_unit,*) equilxy(2,iz,ir),equilxy(1,iz,ir),
     $                 equil(1,iz,ir),equil(2,iz,ir),equil(3,iz,ir)
               ENDDO
            ENDDO
            CLOSE(equil_unit)
            rw=MAXVAL(equilxy(2,:,:))
c-----------------------------------------------------------------------
c           scale variables.
c-----------------------------------------------------------------------
            fac=pmax/MAXVAL(equil(2,:,:))
            flx=SQRT(fac)
            prs=fac
            cur=SQRT(fac)
            equilxy=equilxy/rw
            equil(1,:,:)=equil(1,:,:)*flx
            equil(2,:,:)=equil(2,:,:)*prs
            equil(3,:,:)=equil(3,:,:)*cur
         ENDIF
c-----------------------------------------------------------------------
c       broadcast equilibrium data.
c-----------------------------------------------------------------------
         CALL MPI_Bcast(equil,nqty*(nr+1)*(nz+1),
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
         CALL MPI_Bcast(equilxy,2*(nz+1)*(nr+1),
     $        MPI_DOUBLE_PRECISION,0,comm,ierr)
      END SELECT
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_read_marklin
c-----------------------------------------------------------------------
c     subprogram 3. extra_interp.
c     interpolate the frc equilibrium data onto a grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extra_interp(interp,xpi,ypi,equil_bc,equil,equilxy,
     $     eq)

      CHARACTER(*), INTENT(IN) :: interp
      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      TYPE(bicube_type) :: equil_bc
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: equil,equilxy
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: eq

      INTEGER :: ix,iy
      REAL(r8), DIMENSION(SIZE(eq,1)) :: values
c-----------------------------------------------------------------------
c     loop over all points and get psi, p, and jphi.
c-----------------------------------------------------------------------
      eq=zero

      DO ix=0,SIZE(xpi,1)-1
         DO iy=0,SIZE(xpi,2)-1
            SELECT CASE(interp)
            CASE("bicube")
               CALL bicube_eval(equil_bc,ABS(xpi(ix,iy)),
     $              ABS(ypi(ix,iy)))
               eq(:,ix,iy)=equil_bc%f
            CASE("linear")
               CALL extra_lininterp(ABS(xpi(ix,iy)),ABS(ypi(ix,iy)),
     $              equil,equilxy,values)
               eq(:,ix,iy)=values
            END SELECT
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_interp
c-----------------------------------------------------------------------
c     subprogram 4. extra_lininterp.
c     find values at (x0,y0) using linear interpolation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extra_lininterp(x0,y0,data,xy,values)
      
      REAL(r8), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN):: data,xy
      REAL(r8), DIMENSION(:), INTENT(OUT) :: values

      LOGICAL :: err
      INTEGER :: ix,iy,counter
      REAL(r8), DIMENSION(5) :: dist
c-----------------------------------------------------------------------
c     find the location in the xy grid.
c-----------------------------------------------------------------------
      ix=0
      iy=0
      counter=0
      err=.FALSE.
      CALL extra_xysearch(x0,y0,xy,ix,iy,counter,err)
c-----------------------------------------------------------------------
c     find the values of data at (x0,y0).
c-----------------------------------------------------------------------
      dist(1)=SQRT((x0-xy(1,ix,iy))**2+(y0-xy(2,ix,iy))**2)
      dist(2)=SQRT((x0-xy(1,ix+1,iy))**2+(y0-xy(2,ix+1,iy))**2)
      dist(3)=SQRT((x0-xy(1,ix,iy+1))**2+(y0-xy(2,ix,iy+1))**2)
      dist(4)=SQRT((x0-xy(1,ix+1,iy+1))**2+(y0-xy(2,ix+1,iy+1))**2)
      WHERE(dist < 1e-12) dist=1e-12
      dist(5)=SUM(1._r8/dist(1:4))
      values=(data(:,ix,iy)/dist(1)+data(:,ix+1,iy)/dist(2)
     $     +data(:,ix,iy+1)/dist(3)+data(:,ix+1,iy+1)/dist(4))/dist(5)
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_lininterp
c-----------------------------------------------------------------------
c     subprogram 5. extra_xysearch.
c     finds indices of the point (x,y) in the array xy(:,:,:).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE extra_xysearch(x0,y0,xy,jx,jy,counter,err)

      REAL(r8), INTENT(IN) :: x0,y0
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xy
      INTEGER, INTENT(INOUT) :: jx,jy,counter
      LOGICAL, INTENT(INOUT) :: err

      INTEGER :: ix,iy
c-----------------------------------------------------------------------
c     find the location of (x0,y0) point in xy matrix.
c-----------------------------------------------------------------------
      IF(counter > 100)THEN
         jx=0
         jy=0
         err=.TRUE.
         IF(err)THEN
            write(*,*) "x0,y0",x0,y0
            stop "Error in xysearch. Check boundary limits."
         ENDIF
         RETURN
      ENDIF

      counter = counter + 1
      DO ix=0,SIZE(xy,2)-2
         IF(x0 <= xy(1,ix+1,jy))EXIT
      ENDDO
      jx=ix
      DO iy=0,SIZE(xy,3)-2
         IF(y0 <= xy(2,jx,iy+1))EXIT
      ENDDO
      jy=iy

      IF(xy(1,jx,jy) <= x0 .AND. xy(1,jx+1,jy) >= x0
     $     .AND. xy(2,jx,jy) <= y0 .AND. xy(2,jx,jy+1) >= y0)RETURN
      CALL extra_xysearch(x0,y0,xy,jx,jy,counter,err)
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_xysearch
c-----------------------------------------------------------------------
c     subprogram 6. extra_equil_dealloc.
c     deallocates objects used to store equilibrium data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extra_equil_dealloc(interp,equil_bc,equil,equilxy)

      CHARACTER(*), INTENT(IN) :: interp
      TYPE(bicube_type) :: equil_bc
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: equil,equilxy
c-----------------------------------------------------------------------
c     deallocate appropriate arrays.
c-----------------------------------------------------------------------
      SELECT CASE(interp)
      CASE("bicube")
         CALL bicube_dealloc(equil_bc)
      CASE("linear")
         DEALLOCATE(equil,equilxy)
      END SELECT
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_equil_dealloc
c-----------------------------------------------------------------------
c     subprogram 7. extra_coilalloc.
c     allocates and reads in data for a coil_type variable.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extra_coilalloc(coilfile,coils)

      CHARACTER(*), INTENT(IN) :: coilfile
      TYPE(coil_type) :: coils

      INTEGER :: i,n,myios
c-----------------------------------------------------------------------
c     open input coilfile with coil shape and location info.
c-----------------------------------------------------------------------
      IF(mpi_rank==0)OPEN(UNIT=coil_unit,FILE=TRIM(coilfile),
     $     ACTION="READ",STATUS="OLD",IOSTAT=myios)
      CALL MPI_Bcast(myios,1,MPI_INTEGER,0,comm,ierr)
      IF(myios /= 0)CALL program_stop("Error reading coilfile")
      IF(mpi_rank==0)THEN
         READ(coil_unit,*)
         READ(coil_unit,*)coils%ncoils
         READ(coil_unit,*)coils%lam
         READ(coil_unit,*)
      ENDIF
c-----------------------------------------------------------------------
c     broadcast scalar info and allocate arrays of appropriate size.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(coils%ncoils,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(coils%lam,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      n=coils%ncoils
      ALLOCATE(coils%zb(n),coils%ze(n),coils%vc0(n),coils%tfire(n),
     $     coils%tqtr(n),coils%tcrow(n))
c-----------------------------------------------------------------------
c     read coil array info and close input coilfile.
c
c     the following are the explanations for the array fields: 
c     zb and ze -- axial locations of the beginning/end of each coil.
c     tfire, tqtr, tcrow -- fire, quarter-cycle and crowbar time 
c     for each coil.
c     vc0 -- initial voltage applied to each coil at t=tfire.
c-----------------------------------------------------------------------
      IF(mpi_rank==0)THEN
         DO i=1,n
            READ(coil_unit,*)coils%zb(i),coils%ze(i),coils%vc0(i),
     $           coils%tfire(i),coils%tqtr(i),coils%tcrow(i)
         ENDDO
         CLOSE(coil_unit)
      ENDIF
c-----------------------------------------------------------------------
c     broadcast coil array info.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(coils%zb,n,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(coils%ze,n,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(coils%vc0,n,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(coils%tfire,n,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(coils%tqtr,n,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(coils%tcrow,n,MPI_DOUBLE_PRECISION,0,comm,ierr)
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_coilalloc
c-----------------------------------------------------------------------
c     subprogram 8. extra_coildealloc.
c     deallocates all fields of a coil_type variable.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extra_coildealloc(coils)

      TYPE(coil_type) :: coils
c-----------------------------------------------------------------------
c     deallocate.
c-----------------------------------------------------------------------
      DEALLOCATE(coils%zb,coils%ze,coils%vc0,coils%tfire,coils%tqtr,
     $     coils%tcrow)
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_coildealloc
c-----------------------------------------------------------------------
c     subprogram 9. extra_coileval.
c     computes voltages due to external coils
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extra_coileval(t,x,coils,volt)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x
      TYPE(coil_type) :: coils
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: volt

      INTEGER :: i
      REAL(r8), PARAMETER :: tft=2.e-1
      REAL(r8) :: tq,lamn,zn
      REAL(r8), DIMENSION(coils%ncoils) :: vc,ft
c-----------------------------------------------------------------------
c     account for time-dependence.
c-----------------------------------------------------------------------
      ft=one
      vc=zero
      volt=zero

      WHERE(t >= coils%tfire .AND. t <= coils%tfire + coils%tqtr*tft)
         ft = (one - COS(pi*(t-coils%tfire)/(coils%tqtr*tft)))*half
      ELSEWHERE(t <= coils%tcrow 
     $        .AND. t >= coils%tcrow - coils%tqtr*tft)
         ft = (one + COS(pi*(t-coils%tcrow+coils%tqtr*tft)
     $        /(coils%tqtr*tft)))*half
      END WHERE

      WHERE(t > coils%tfire .AND. t < coils%tcrow)
     $     vc = ft*coils%vc0*COS(pi*(t-coils%tfire)/(two*coils%tqtr))

      IF(SUM(vc)==0)THEN
         volt=zero
         RETURN
      ENDIF
c-----------------------------------------------------------------------
c     calculate voltage using tanh functions to smooth field.
c-----------------------------------------------------------------------
      volt=vc(1)
      DO i=2,coils%ncoils
         zn = half*(coils%zb(i) + coils%ze(i-1))
         lamn = (coils%zb(i) - coils%ze(i-1))/coils%lam
         volt = volt + half*(vc(i)-vc(i-1))*(one + TANH((x-zn)/lamn))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extra_coileval
      END MODULE extra_mod
