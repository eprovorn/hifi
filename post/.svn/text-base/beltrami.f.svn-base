c-----------------------------------------------------------------------
c     file beltrami.f.
c     post-processes output from sel code for solutions of Beltrami Eq..
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. beltrami_mod.
c     1. beltrami_read.
c     2. beltrami_grid.
c-----------------------------------------------------------------------
c     subprogram 0. beltrami_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE beltrami_mod
      USE slice_mod
      IMPLICIT NONE

      CHARACTER(8), PRIVATE :: bel_grid_type="align",solve_type="full"
      LOGICAL, PRIVATE :: bel_diagnose
      INTEGER, PRIVATE :: bel_npsi,nx_grid,ny_grid,npoints=64,nt_prev=0
      REAL(r8), PRIVATE :: bel_k0,bel_k1,bel_eps0,bel_phifac

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. beltrami_read.
c     read necessary post-processing parameters from  beltrami.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_read

      CHARACTER(80) :: filename
      INTEGER :: itmax,itmax_decr

      NAMELIST/beltrami_input/bel_diagnose,bel_eps0,bel_k0,bel_k1,
     $     bel_npsi,bel_grid_type,bel_phifac,itmax,itmax_decr,solve_type
      NAMELIST/output_control/nx_grid,ny_grid,npoints
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      filename=TRIM(indir)//"/beltrami.in"
      OPEN(UNIT=in_unit,FILE=filename,STATUS="OLD")
      READ(in_unit,NML=beltrami_input)
      READ(in_unit,NML=output_control)
      CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_read
c-----------------------------------------------------------------------
c     subprogram 2. beltrami_grid.
c     draws the resultant grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE beltrami_grid(basis,xx,yy,uu,polar)

      TYPE(jacobi_type) :: basis
      REAL(r8), DIMENSION(0:), INTENT(IN) :: xx,yy
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      LOGICAL, INTENT(IN) :: polar

      INTEGER :: ix,iy,ixx=0,iyy=0
      REAL(r8) :: x,y,dx,dy
      REAL(r8), DIMENSION(SIZE(uu,1)) :: u,ux,uy
c-----------------------------------------------------------------------
c     define increments and open file.
c-----------------------------------------------------------------------
      OPEN(UNIT=belgrid_unit,FILE="belgrid.bin",STATUS="UNKNOWN",
     $     FORM="UNFORMATTED")
c-----------------------------------------------------------------------
c     draw horizontal curves.
c-----------------------------------------------------------------------
      dx=(xx(SIZE(xx)-1)-xx(0))/npoints
      dy=(yy(SIZE(yy)-1)-yy(0))/ny_grid
      y=yy(0)
      DO iy=0,ny_grid
         x=xx(0)
         DO ix=0,npoints
            CALL slice_interp(basis,xx,yy,uu,x,y,ixx,iyy,u,ux,uy,0)
            IF(polar)THEN
               u(1)=u(1)*COS(u(2))
               u(2)=u(1)*TAN(u(2))
            ENDIF
            WRITE(belgrid_unit)REAL(u(1:2),4)
            x=x+dx
         ENDDO
         WRITE(belgrid_unit)
         y=y+dy
      ENDDO
c-----------------------------------------------------------------------
c     draw vertical curves.
c-----------------------------------------------------------------------
      dx=(xx(SIZE(xx)-1)-xx(0))/nx_grid
      dy=(yy(SIZE(yy)-1)-yy(0))/npoints
      x=xx(0)
      DO ix=0,nx_grid
         y=yy(0)
         DO iy=0,npoints
            CALL slice_interp(basis,xx,yy,uu,x,y,ixx,iyy,u,ux,uy,0)
            IF(polar)THEN
               u(1)=u(1)*COS(u(2))
               u(2)=u(1)*TAN(u(2))
            ENDIF
            WRITE(belgrid_unit)REAL(u(1:2),4)
            y=y+dy
         ENDDO
         WRITE(belgrid_unit)
         x=x+dx
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CLOSE(UNIT=belgrid_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE beltrami_grid
      END MODULE beltrami_mod
