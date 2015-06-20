c-----------------------------------------------------------------------
c     file slice.f.
c     creates slice plots.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. slice_mod.
c     1. slice_interp.
c     2. slice_xt.
c     3. slice_yt.
c     4. slice_xy.
c     5. slice_yx.
c-----------------------------------------------------------------------
c     subprogram 0. slice_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE slice_mod
      USE jacobi_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. slice_interp.
c     interpolate solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE slice_interp(basis,xx,yy,uu,x,y,ix,iy,u,ux,uy,mode)

      REAL(r8), DIMENSION(0:), INTENT(IN) :: xx,yy
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(r8), INTENT(IN) :: x,y
      INTEGER, INTENT(INOUT) :: ix,iy
      REAL(r8), DIMENSION(SIZE(uu,1)), INTENT(OUT) :: u,ux,uy
      TYPE(jacobi_type), INTENT(INOUT) :: basis
      INTEGER, INTENT(IN) :: mode

      INTEGER :: jx,jy,inx,iny,nx,ny,nqty,np
      REAL(r8) :: z
      REAL(r8), DIMENSION(0:basis%np) :: px,py,qx,qy
      REAL(r8), DIMENSION(SIZE(xx)-1) :: dx
      REAL(r8), DIMENSION(SIZE(yy)-1) :: dy
      REAL(r8) :: dxfac,dyfac
c-----------------------------------------------------------------------
c     define local sizes.
c-----------------------------------------------------------------------
      nx=SIZE(xx)-1
      ny=SIZE(yy)-1
      nqty=SIZE(uu,1)
      np=basis%np
c-----------------------------------------------------------------------
c     locate x interval, compute offset and basis.
c-----------------------------------------------------------------------
      ix=MAX(MIN(ix,nx),1)
      DO
         IF(x >= xx(ix-1) .OR. ix == 1)EXIT
         ix=ix-1
      ENDDO
      DO
         IF(x <= xx(ix) .OR. ix == nx)EXIT
         ix=ix+1
      ENDDO
      jx=(ix-1)*(np+1)
      z=2*(x-xx(ix-1))/(xx(ix)-xx(ix-1))-1
      CALL jacobi_basis(z,basis)
      px=basis%pb
      qx=basis%qb
c-----------------------------------------------------------------------
c     locate y interval and compute offset and basis
c-----------------------------------------------------------------------
      iy=MAX(MIN(iy,ny),1)
      DO
         IF(y >= yy(iy-1) .OR. iy == 1)EXIT
         iy=iy-1
      ENDDO
      DO
         IF(y <= yy(iy) .OR. iy == ny)EXIT
         iy=iy+1
      ENDDO
      jy=(iy-1)*(np+1)
      z=2*(y-yy(iy-1))/(yy(iy)-yy(iy-1))-1
      CALL jacobi_basis(z,basis)
      py=basis%pb
      qy=basis%qb
c-----------------------------------------------------------------------
c     interpolate solutions.
c-----------------------------------------------------------------------
      u=0
      DO inx=0,np
         DO iny=0,np
            u=u+uu(:,jx+inx,jy+iny)*px(inx)*py(iny)
         ENDDO
      ENDDO
      IF(mode == 0)RETURN
c-----------------------------------------------------------------------
c     evaluate derivatives.
c-----------------------------------------------------------------------
      ux=0
      uy=0
      DO inx=0,np
         DO iny=0,np
            ux=ux+uu(:,jx+inx,jy+iny)*qx(inx)*py(iny)
            uy=uy+uu(:,jx+inx,jy+iny)*px(inx)*qy(iny)
         ENDDO
      ENDDO
      dx=xx(1:nx)-xx(0:nx-1)
      dy=yy(1:ny)-yy(0:ny-1)
      dyfac=2/dy(iy)
      dxfac=2/dx(ix)
      ux=ux*dxfac
      uy=uy*dyfac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slice_interp
c-----------------------------------------------------------------------
c     subprogram 2. slice_xt.
c     draws a t family of curves with independent variable x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE slice_xt(basis,xx,yy,uu,t,ys,nx)

      TYPE(jacobi_type) :: basis
      REAL(r8), DIMENSION(0:), INTENT(IN) :: xx,yy
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(8), INTENT(IN) :: ys
      INTEGER, INTENT(IN) :: t,nx

      INTEGER :: ix,ixx=0,iyy=0
      REAL(r8) :: x,dx,y
      REAL(r8), DIMENSION(SIZE(uu,1)) :: u,ux,uy
c-----------------------------------------------------------------------
c     draw curve.
c-----------------------------------------------------------------------
      x=xx(0)
      dx=(xx(SIZE(xx)-1)-xx(0))/nx
      y=yy(0)+ys*(yy(SIZE(yy)-1)-yy(0))
      DO ix=0,nx
         CALL slice_interp(basis,xx,yy,uu,x,y,ixx,iyy,u,ux,uy,0)
         WRITE(xt_unit)REAL(x,4),REAL(t,4),REAL(u,4)
         x=x+dx
      ENDDO
      WRITE(xt_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slice_xt
c-----------------------------------------------------------------------
c     subprogram 3. slice_yt.
c     draws a t family of curves with independent variable y.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE slice_yt(basis,xx,yy,uu,t,xs,ny)

      TYPE(jacobi_type) :: basis
      REAL(r8), DIMENSION(0:), INTENT(IN) :: xx,yy
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      REAL(8), INTENT(IN) :: xs
      INTEGER, INTENT(IN) :: t,ny

      INTEGER :: iy,ixx=0,iyy=0
      REAL(r8) :: x,y,dy
      REAL(r8), DIMENSION(SIZE(uu,1)) :: u,ux,uy
c-----------------------------------------------------------------------
c     draw curve.
c-----------------------------------------------------------------------
      x=xx(0)+xs*(xx(SIZE(xx)-1)-xx(0))
      y=yy(0)
      dy=(yy(SIZE(yy)-1)-yy(0))/ny
      DO iy=0,ny
         CALL slice_interp(basis,xx,yy,uu,x,y,ixx,iyy,u,ux,uy,0)
         WRITE(yt_unit)REAL(y,4),REAL(t,4),REAL(u,4)
         y=y+dy
      ENDDO
      WRITE(yt_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slice_yt
c-----------------------------------------------------------------------
c     subprogram 4. slice_xy.
c     draws a y family of curves with independent variable x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE slice_xy(basis,xx,yy,uu,nx,ny,name)

      TYPE(jacobi_type) :: basis
      REAL(r8), DIMENSION(0:), INTENT(IN) :: xx,yy
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      INTEGER, INTENT(IN) :: nx,ny
      CHARACTER(*), INTENT(IN) :: name

      INTEGER :: ix,iy,ixx=0,iyy=0
      REAL(r8) :: x,y,dx,dy
      REAL(r8), DIMENSION(SIZE(uu,1)) :: u,ux,uy
c-----------------------------------------------------------------------
c     define increments and open file.
c-----------------------------------------------------------------------
      dx=(xx(SIZE(xx)-1)-xx(0))/nx
      dy=(yy(SIZE(yy)-1)-yy(0))/ny
      OPEN(UNIT=xy_unit,FILE=TRIM(name)//".bin",STATUS="UNKNOWN",
     $     FORM="UNFORMATTED")
c-----------------------------------------------------------------------
c     draw curves.
c-----------------------------------------------------------------------
      y=yy(0)
      DO iy=0,ny
         x=xx(0)
         DO ix=0,nx
            CALL slice_interp(basis,xx,yy,uu,x,y,ixx,iyy,u,ux,uy,0)
            WRITE(xy_unit)REAL(x,4),REAL(y,4),REAL(u,4)
            x=x+dx
         ENDDO
         WRITE(xy_unit)
         y=y+dy
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CLOSE(UNIT=xy_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slice_xy
c-----------------------------------------------------------------------
c     subprogram 5. slice_yx.
c     draws a y family of curves with independent variable x.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE slice_yx(basis,xx,yy,uu,nx,ny,name)

      TYPE(jacobi_type) :: basis
      REAL(r8), DIMENSION(0:), INTENT(IN) :: xx,yy
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: uu
      INTEGER, INTENT(IN) :: nx,ny
      CHARACTER(*), INTENT(IN) :: name

      INTEGER :: ix,iy,ixx=0,iyy=0
      REAL(r8) :: x,y,dx,dy
      REAL(r8), DIMENSION(SIZE(uu,1)) :: u,ux,uy
c-----------------------------------------------------------------------
c     define increments and open file.
c-----------------------------------------------------------------------
      dx=(xx(SIZE(xx)-1)-xx(0))/nx
      dy=(yy(SIZE(yy)-1)-yy(0))/ny
      OPEN(UNIT=xy_unit,FILE=TRIM(name)//".bin",STATUS="UNKNOWN",
     $     FORM="UNFORMATTED")
c-----------------------------------------------------------------------
c     draw curves.
c-----------------------------------------------------------------------
      x=xx(0)
      DO ix=0,nx
         y=yy(0)
         DO iy=0,ny
            CALL slice_interp(basis,xx,yy,uu,x,y,ixx,iyy,u,ux,uy,0)
            WRITE(xy_unit)REAL(x,4),REAL(y,4),REAL(u,4)
            y=y+dy
         ENDDO
         WRITE(xy_unit)
         x=x+dx
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CLOSE(UNIT=xy_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE slice_yx
      END MODULE slice_mod
