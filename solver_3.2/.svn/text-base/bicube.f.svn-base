c-----------------------------------------------------------------------
c     file bicube.f.
c     fits functions to bicubic splines.
c     Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
c     Translated from the German by W. D. Hoskins and H. W. Sager.
c     Utilitas Mathematica Publishing Inc., Winnepeg, 1974.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. bicube_type definition.
c      1. bicube_alloc.
c      2. bicube_dealloc.
c      3. bicube_lsfit.
c      4. bicube_lsfit_xp.
c      5. bicube_lsfit_yp.
c      6. bicube_eval.
c      7. bicube_getco.
c      8. bicube_fit.
c-----------------------------------------------------------------------
c     subprogram 0. bicube_type definition.
c     defines bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE bicube_mod
      USE spline_mod
      IMPLICIT NONE

      TYPE :: bicube_type
      INTEGER :: mx,my,nqty,ix,iy
      REAL(r8), DIMENSION(2) :: x0,y0
      REAL(r8), DIMENSION(:), POINTER :: xs,ys
      REAL(r8), DIMENSION(:,:), POINTER :: xpower,ypower
      REAL(r8), DIMENSION(:,:,:), POINTER :: fs,fsx,fsy,fsxy
      REAL(r8), DIMENSION(:), POINTER :: f
      LOGICAL, DIMENSION(2) :: periodic
      END TYPE bicube_type

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. bicube_alloc.
c     allocates space for bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_alloc(bcs,mx,my,nqty)

      INTEGER, INTENT(IN) :: mx,my,nqty
      TYPE(bicube_type), INTENT(OUT) :: bcs
c-----------------------------------------------------------------------
c     set scalars.
c-----------------------------------------------------------------------
      bcs%mx=mx
      bcs%my=my
      bcs%ix=0
      bcs%iy=0
      bcs%nqty=nqty
      bcs%periodic=.FALSE.
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(bcs%xs(0:mx))
      ALLOCATE(bcs%ys(0:my))
      ALLOCATE(bcs%fs(0:mx,0:my,nqty))
      ALLOCATE(bcs%fsx(0:mx,0:my,nqty))
      ALLOCATE(bcs%fsy(0:mx,0:my,nqty))
      ALLOCATE(bcs%fsxy(0:mx,0:my,nqty))
      ALLOCATE(bcs%f(nqty))
      ALLOCATE(bcs%xpower(2,nqty),bcs%ypower(2,nqty))
      bcs%xpower=0
      bcs%ypower=0
      bcs%x0=0
      bcs%y0=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_alloc
c-----------------------------------------------------------------------
c     subprogram 2. bicube_dealloc.
c     deallocates space for bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_dealloc(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(bcs%xs)
      DEALLOCATE(bcs%ys)
      DEALLOCATE(bcs%fs)
      DEALLOCATE(bcs%fsx)
      DEALLOCATE(bcs%fsy)
      DEALLOCATE(bcs%fsxy)
      DEALLOCATE(bcs%f)
      DEALLOCATE(bcs%xpower,bcs%ypower)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. bicube_lsfit.
c     least-square fit to cubic splines of piecewise-constant functions.
c     minimization with respect to f,fx,fy,fxy.
c     no imposed boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_lsfit(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs

      INTEGER :: ix,iy,nx,ny,mx,my,nrhs,n,info,i,j,k,l,kd,ldab
      INTEGER, DIMENSION(7*(bcs%mx+1)*(bcs%my+1)) :: ipiv
      REAL(r8), PARAMETER ::
     $     a1=169/1225._r8,a2=117/2450._r8,a3=81/4900._r8,
     $     a4=143/7350._r8,a5=169/14700._r8,a6=33/4900._r8,
     $     a7=39/9800._r8,a8=121/44100._r8,a9=143/88200._r8,
     $     a10=169/176400._r8,a11=13/3675._r8,a12=3/2450._r8,
     $     a13=13/4900._r8,a14=9/9800._r8,a15=11/22050._r8,
     $     a16=13/44100._r8,a17=11/29400._r8,a18=13/58800._r8,
     $     a19=1/4._r8,a20=1/24._r8,a21=1/11025._r8,
     $     a22=1/14700._r8,a23=1/19600._r8,a24=1/144._r8
      REAL(r8), DIMENSION(0:bcs%mx+1) :: dx
      REAL(r8), DIMENSION(0:bcs%my+1) :: dy
      REAL(r8), DIMENSION(7,0:bcs%mx,0:bcs%my,bcs%nqty) :: rhs
      REAL(r8), DIMENSION(0:bcs%mx+1,0:bcs%my+1,bcs%nqty) :: g
      REAL(r8), DIMENSION(7,7,-1:1,-1:1,0:bcs%mx,0:bcs%my) :: amat
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ab
c-----------------------------------------------------------------------
c     define sizes.
c-----------------------------------------------------------------------
      nx=bcs%mx
      ny=bcs%my
      nrhs=bcs%nqty
      kd=7*(nx+1)+13
      ldab=3*kd+1
      n=7*(nx+1)*(ny+1)
c-----------------------------------------------------------------------
c     initialize arrays.
c-----------------------------------------------------------------------
      amat=0
      rhs=0
      dx=0
      dy=0
      g=0
      dx(1:nx)=bcs%xs(1:nx)-bcs%xs(0:nx-1)
      dy(1:ny)=bcs%ys(1:ny)-bcs%ys(0:ny-1)
      g(1:nx,1:ny,:)=bcs%fs(1:nx,1:ny,:)
      g(0,:,:)=g(1,:,:)
      g(nx+1,:,:)=g(nx,:,:)
      g(:,0,:)=g(:,1,:)
      g(:,ny+1,:)=g(:,ny,:)
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,1,0,0,0:nx,iy)=+a1
     $        *(dx(0:nx)+dx(1:nx+1))*(dy(iy)+dy(iy+1))
         amat(1,1,-1,0,0:nx,iy)=+a2*dx(0:nx)*(dy(iy)+dy(iy+1))
         amat(1,1,+1,0,0:nx,iy)=+a2*dx(1:nx+1)*(dy(iy)+dy(iy+1))
         amat(1,1,0,-1,0:nx,iy)=+a2*(dx(0:nx)+dx(1:nx+1))*dy(iy)
         amat(1,1,0,+1,0:nx,iy)=+a2*(dx(0:nx)+dx(1:nx+1))*dy(iy+1)
         amat(1,1,-1,-1,0:nx,iy)=+a3*dx(0:nx)*dy(iy)
         amat(1,1,-1,+1,0:nx,iy)=+a3*dx(0:nx)*dy(iy+1)
         amat(1,1,+1,-1,0:nx,iy)=+a3*dx(1:nx+1)*dy(iy)
         amat(1,1,+1,+1,0:nx,iy)=+a3*dx(1:nx+1)*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,2,0,0,0:nx,iy)=+a4
     $        *(dx(1:nx+1)**2-dx(0:nx)**2)*(dy(iy)+dy(iy+1))
         amat(1,2,-1,0,0:nx,iy)=+a5*(dy(iy)+dy(iy+1))*dx(0:nx)**2
         amat(1,2,+1,0,0:nx,iy)=-a5*(dy(iy)+dy(iy+1))*dx(1:nx+1)**2
         amat(1,2,0,-1,0:nx,iy)=+a6*(dx(1:nx+1)**2-dx(0:nx)**2)*dy(iy)
         amat(1,2,0,+1,0:nx,iy)=+a6*(dx(1:nx+1)**2-dx(0:nx)**2)*
     $        dy(iy+1)
         amat(1,2,-1,-1,0:nx,iy)=+a7*dx(0:nx)**2*dy(iy)
         amat(1,2,-1,+1,0:nx,iy)=+a7*dx(0:nx)**2*dy(iy+1)
         amat(1,2,+1,-1,0:nx,iy)=-a7*dx(1:nx+1)**2*dy(iy)
         amat(1,2,+1,+1,0:nx,iy)=-a7*dx(1:nx+1)**2*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, y derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(1,3,0,0,ix,0:ny)=+a4
     $        *(dy(1:ny+1)**2-dy(0:ny)**2)*(dx(ix)+dx(ix+1))
         amat(1,3,-1,0,ix,0:ny)=+a6*dx(ix)*
     $        (dy(1:ny+1)**2-dy(0:ny)**2)
         amat(1,3,+1,0,ix,0:ny)=+a6*dx(ix+1)*
     $        (dy(1:ny+1)**2-dy(0:ny)**2)
         amat(1,3,0,-1,ix,0:ny)=+a5*(dx(ix)+dx(ix+1))*dy(0:ny)**2
         amat(1,3,0,+1,ix,0:ny)=-a5*(dx(ix)+dx(ix+1))*dy(1:ny+1)**2
         amat(1,3,-1,-1,ix,0:ny)=+a7*dy(0:ny)**2*dx(ix)
         amat(1,3,-1,+1,ix,0:ny)=-a7*dy(1:ny+1)**2*dx(ix)
         amat(1,3,+1,-1,ix,0:ny)=+a7*dy(0:ny)**2*dx(ix+1)
         amat(1,3,+1,+1,ix,0:ny)=-a7*dy(1:ny+1)**2*dx(ix+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,4,0,0,0:nx,iy)=+a8
     $        *(dx(0:nx)**2-dx(1:nx+1)**2)*(dy(iy)**2-dy(iy+1)**2)
         amat(1,4,-1,0,0:nx,iy)=+a9*dx(0:nx)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(1,4,+1,0,0:nx,iy)=-a9*dx(1:nx+1)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(1,4,0,-1,0:nx,iy)=+a9*(dx(1:nx+1)**2-dx(0:nx)**2)
     $        *dy(iy)**2
         amat(1,4,0,+1,0:nx,iy)=-a9*(dx(1:nx+1)**2-dx(0:nx)**2)
     $        *dy(iy+1)**2
         amat(1,4,-1,-1,0:nx,iy)=+a10*dx(0:nx)**2*dy(iy)**2
         amat(1,4,-1,+1,0:nx,iy)=-a10*dx(0:nx)**2*dy(iy+1)**2
         amat(1,4,+1,-1,0:nx,iy)=-a10*dx(1:nx+1)**2*dy(iy)**2
         amat(1,4,+1,+1,0:nx,iy)=+a10*dx(1:nx+1)**2*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, lambdas.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,5,-1,0,2:nx,iy)=-3*dx(1:nx-1)/dx(2:nx)
         amat(1,5,0,0,1:nx-1,iy)=3*
     $        (dx(1:nx-1)/dx(2:nx)-dx(2:nx)/dx(1:nx-1))
         amat(1,5,1,0,0:nx-2,iy)=3*dx(2:nx)/dx(1:nx-1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, mus.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(1,6,0,-1,ix,2:ny)=-3*dy(1:ny-1)/dy(2:ny)
         amat(1,6,0,0,ix,1:ny-1)=3*
     $        (dy(1:ny-1)/dy(2:ny)-dy(2:ny)/dy(1:ny-1))
         amat(1,6,0,1,ix,0:ny-2)=3*dy(2:ny)/dy(1:ny-1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx
            rhs(1,ix,iy,:)=a19*((dx(ix)*g(ix,iy,:)
     $           +dx(ix+1)*g(ix+1,iy,:))*dy(iy)
     $           +(dx(ix)*g(ix,iy+1,:)
     $           +dx(ix+1)*g(ix+1,iy+1,:))*dy(iy+1))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(2,1,0,0,0:nx,iy)=+a4
     $        *(dx(1:nx+1)**2-dx(0:nx)**2)*(dy(iy)+dy(iy+1))
         amat(2,1,-1,0,0:nx,iy)=-a5*dx(0:nx)**2*(dy(iy)+dy(iy+1))
         amat(2,1,+1,0,0:nx,iy)=+a5*dx(1:nx+1)**2*(dy(iy)+dy(iy+1))
         amat(2,1,0,-1,0:nx,iy)=+a6*(dx(1:nx+1)**2-dx(0:nx)**2)*dy(iy)
         amat(2,1,0,+1,0:nx,iy)=+a6*(dx(1:nx+1)**2-dx(0:nx)**2)*dy(iy+1)
         amat(2,1,-1,-1,0:nx,iy)=-a7*dx(0:nx)**2*dy(iy)
         amat(2,1,-1,+1,0:nx,iy)=-a7*dx(0:nx)**2*dy(iy+1)
         amat(2,1,+1,-1,0:nx,iy)=+a7*dx(1:nx+1)**2*dy(iy)
         amat(2,1,+1,+1,0:nx,iy)=+a7*dx(1:nx+1)**2*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(2,2,0,0,0:nx,iy)=+a11
     $        *(dx(1:nx+1)**3+dx(0:nx)**3)*(dy(iy)+dy(iy+1))
         amat(2,2,-1,0,0:nx,iy)=-a13*(dy(iy)+dy(iy+1))*dx(0:nx)**3
         amat(2,2,+1,0,0:nx,iy)=-a13*(dy(iy)+dy(iy+1))*dx(1:nx+1)**3
         amat(2,2,0,-1,0:nx,iy)=+a12*(dx(1:nx+1)**3+dx(0:nx)**3)*dy(iy)
         amat(2,2,0,+1,0:nx,iy)=+a12*(dx(1:nx+1)**3+dx(0:nx)**3)*
     $        dy(iy+1)
         amat(2,2,-1,-1,0:nx,iy)=-a14*dx(0:nx)**3*dy(iy)
         amat(2,2,-1,+1,0:nx,iy)=-a14*dx(0:nx)**3*dy(iy+1)
         amat(2,2,+1,-1,0:nx,iy)=-a14*dx(1:nx+1)**3*dy(iy)
         amat(2,2,+1,+1,0:nx,iy)=-a14*dx(1:nx+1)**3*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, y derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(2,3,0,0,ix,0:ny)=+a8
     $        *(dy(1:ny+1)**2-dy(0:ny)**2)*(dx(ix+1)**2-dx(ix)**2)
         amat(2,3,-1,0,ix,0:ny)=-a9*dx(ix)**2*
     $        (dy(1:ny+1)**2-dy(0:ny)**2)
         amat(2,3,+1,0,ix,0:ny)=+a9*dx(ix+1)**2*
     $        (dy(1:ny+1)**2-dy(0:ny)**2)
         amat(2,3,0,-1,ix,0:ny)=+a9*(dx(ix+1)**2-dx(ix)**2)*dy(0:ny)**2
         amat(2,3,0,+1,ix,0:ny)=-a9*(dx(ix+1)**2-dx(ix)**2)
     $        *dy(1:ny+1)**2
         amat(2,3,-1,-1,ix,0:ny)=-a10*dy(0:ny)**2*dx(ix)**2
         amat(2,3,-1,+1,ix,0:ny)=+a10*dy(1:ny+1)**2*dx(ix)**2
         amat(2,3,+1,-1,ix,0:ny)=+a10*dy(0:ny)**2*dx(ix+1)**2
         amat(2,3,+1,+1,ix,0:ny)=-a10*dy(1:ny+1)**2*dx(ix+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(2,4,0,0,0:nx,iy)=+a15
     $        *(dx(0:nx)**3+dx(1:nx+1)**3)*(dy(iy+1)**2-dy(iy)**2)
         amat(2,4,-1,0,0:nx,iy)=-a17*dx(0:nx)**3
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(2,4,+1,0,0:nx,iy)=-a17*dx(1:nx+1)**3
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(2,4,0,-1,0:nx,iy)=+a16*(dx(1:nx+1)**3+dx(0:nx)**3)
     $        *dy(iy)**2
         amat(2,4,0,+1,0:nx,iy)=-a16*(dx(1:nx+1)**3+dx(0:nx)**3)
     $        *dy(iy+1)**2
         amat(2,4,-1,-1,0:nx,iy)=-a18*dx(0:nx)**3*dy(iy)**2
         amat(2,4,-1,+1,0:nx,iy)=+a18*dx(0:nx)**3*dy(iy+1)**2
         amat(2,4,+1,-1,0:nx,iy)=-a18*dx(1:nx+1)**3*dy(iy)**2
         amat(2,4,+1,+1,0:nx,iy)=+a18*dx(1:nx+1)**3*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, lambdas.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(2,5,-1,0,2:nx,iy)=dx(1:nx-1)
         amat(2,5,0,0,1:nx-1,iy)=2*(dx(1:nx-1)+dx(2:nx))
         amat(2,5,1,0,0:nx-2,iy)=dx(2:nx)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, nus.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(2,7,0,-1,ix,2:ny)=-3*dy(1:ny-1)/dy(2:ny)
         amat(2,7,0,0,ix,1:ny-1)=3*
     $        (dy(1:ny-1)/dy(2:ny)-dy(2:ny)/dy(1:ny-1))
         amat(2,7,0,1,ix,0:ny-2)=3*dy(2:ny)/dy(1:ny-1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx
            rhs(2,ix,iy,:)=a20*((-dx(ix)**2*g(ix,iy,:)
     $           +dx(ix+1)**2*g(ix+1,iy,:))*dy(iy)
     $           +(-dx(ix)**2*g(ix,iy+1,:)
     $           +dx(ix+1)**2*g(ix+1,iy+1,:))*dy(iy+1))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(3,1,0,0,0:nx,iy)=+a4
     $        *(dx(0:nx)+dx(1:nx+1))*(dy(iy+1)**2-dy(iy)**2)
         amat(3,1,-1,0,0:nx,iy)=+a6*dx(0:nx)*(dy(iy+1)**2-dy(iy)**2)
         amat(3,1,+1,0,0:nx,iy)=+a6*dx(1:nx+1)*(dy(iy+1)**2-dy(iy)**2)
         amat(3,1,0,-1,0:nx,iy)=-a5*(dx(0:nx)+dx(1:nx+1))*dy(iy)**2
         amat(3,1,0,+1,0:nx,iy)=+a5*(dx(0:nx)+dx(1:nx+1))*dy(iy+1)**2
         amat(3,1,-1,-1,0:nx,iy)=-a7*dx(0:nx)*dy(iy)**2
         amat(3,1,-1,+1,0:nx,iy)=+a7*dx(0:nx)*dy(iy+1)**2
         amat(3,1,+1,-1,0:nx,iy)=-a7*dx(1:nx+1)*dy(iy)**2
         amat(3,1,+1,+1,0:nx,iy)=+a7*dx(1:nx+1)*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(3,2,0,0,0:nx,iy)=+a8
     $        *(dx(1:nx+1)**2-dx(0:nx)**2)*(dy(iy+1)**2-dy(iy)**2)
         amat(3,2,-1,0,0:nx,iy)=+a9*(dy(iy+1)**2-dy(iy)**2)*dx(0:nx)**2
         amat(3,2,+1,0,0:nx,iy)=-a9*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(1:nx+1)**2
         amat(3,2,0,-1,0:nx,iy)=-a9*(dx(1:nx+1)**2-dx(0:nx)**2)
     $        *dy(iy)**2
         amat(3,2,0,+1,0:nx,iy)=+a9*(dx(1:nx+1)**2-dx(0:nx)**2)*
     $        dy(iy+1)**2
         amat(3,2,-1,-1,0:nx,iy)=-a10*dx(0:nx)**2*dy(iy)**2
         amat(3,2,-1,+1,0:nx,iy)=+a10*dx(0:nx)**2*dy(iy+1)**2
         amat(3,2,+1,-1,0:nx,iy)=+a10*dx(1:nx+1)**2*dy(iy)**2
         amat(3,2,+1,+1,0:nx,iy)=-a10*dx(1:nx+1)**2*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, y derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(3,3,0,0,ix,0:ny)=+a11
     $        *(dy(1:ny+1)**3+dy(0:ny)**3)*(dx(ix)+dx(ix+1))
         amat(3,3,-1,0,ix,0:ny)=+a12*dx(ix)*
     $        (dy(1:ny+1)**3+dy(0:ny)**3)
         amat(3,3,+1,0,ix,0:ny)=+a12*dx(ix+1)*
     $        (dy(1:ny+1)**3+dy(0:ny)**3)
         amat(3,3,0,-1,ix,0:ny)=-a13*(dx(ix)+dx(ix+1))*dy(0:ny)**3
         amat(3,3,0,+1,ix,0:ny)=-a13*(dx(ix)+dx(ix+1))*dy(1:ny+1)**3
         amat(3,3,-1,-1,ix,0:ny)=-a14*dy(0:ny)**3*dx(ix)
         amat(3,3,-1,+1,ix,0:ny)=-a14*dy(1:ny+1)**3*dx(ix)
         amat(3,3,+1,-1,ix,0:ny)=-a14*dy(0:ny)**3*dx(ix+1)
         amat(3,3,+1,+1,ix,0:ny)=-a14*dy(1:ny+1)**3*dx(ix+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(3,4,0,0,0:nx,iy)=+a15
     $        *(dx(1:nx+1)**2-dx(0:nx)**2)*(dy(iy)**3+dy(iy+1)**3)
         amat(3,4,-1,0,0:nx,iy)=+a16*dx(0:nx)**2
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(3,4,+1,0,0:nx,iy)=-a16*dx(1:nx+1)**2
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(3,4,0,-1,0:nx,iy)=-a17*(dx(1:nx+1)**2-dx(0:nx)**2)
     $        *dy(iy)**3
         amat(3,4,0,+1,0:nx,iy)=-a17*(dx(1:nx+1)**2-dx(0:nx)**2)
     $        *dy(iy+1)**3
         amat(3,4,-1,-1,0:nx,iy)=-a18*dx(0:nx)**2*dy(iy)**3
         amat(3,4,-1,+1,0:nx,iy)=-a18*dx(0:nx)**2*dy(iy+1)**3
         amat(3,4,+1,-1,0:nx,iy)=+a18*dx(1:nx+1)**2*dy(iy)**3
         amat(3,4,+1,+1,0:nx,iy)=+a18*dx(1:nx+1)**2*dy(iy+1)**3
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, etas.
c-----------------------------------------------------------------------
      DO iy=0,ny,ny
         amat(3,7,-1,0,2:nx,iy)=-3*dx(1:nx-1)/dx(2:nx)
         amat(3,7,0,0,1:nx-1,iy)=3*
     $        (dx(1:nx-1)/dx(2:nx)-dx(2:nx)/dx(1:nx-1))
         amat(3,7,1,0,0:nx-2,iy)=3*dx(2:nx)/dx(1:nx-1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, mus.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(3,6,0,-1,ix,2:ny)=dy(1:ny-1)
         amat(3,6,0,0,ix,1:ny-1)=2*(dy(1:ny-1)+dy(2:ny))
         amat(3,6,0,1,ix,0:ny-2)=dy(2:ny)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx
            rhs(3,ix,iy,:)=a20*((-dx(ix)*g(ix,iy,:)
     $           -dx(ix+1)*g(ix+1,iy,:))*dy(iy)**2
     $           +(dx(ix)*g(ix,iy+1,:)
     $           +dx(ix+1)*g(ix+1,iy+1,:))*dy(iy+1)**2)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(4,1,0,0,0:nx,iy)=+a8
     $        *(dx(1:nx+1)**2-dx(0:nx)**2)*(dy(iy+1)**2-dy(iy)**2)
         amat(4,1,-1,0,0:nx,iy)=-a9*dx(0:nx)**2*(dy(iy+1)**2-dy(iy)**2)
         amat(4,1,+1,0,0:nx,iy)=+a9*dx(1:nx+1)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(4,1,0,-1,0:nx,iy)=-a9*(dx(1:nx+1)**2-dx(0:nx)**2)
     $        *dy(iy)**2
         amat(4,1,0,+1,0:nx,iy)=+a9*(dx(1:nx+1)**2-dx(0:nx)**2)
     $        *dy(iy+1)**2
         amat(4,1,-1,-1,0:nx,iy)=+a10*dx(0:nx)**2*dy(iy)**2
         amat(4,1,-1,+1,0:nx,iy)=-a10*dx(0:nx)**2*dy(iy+1)**2
         amat(4,1,+1,-1,0:nx,iy)=-a10*dx(1:nx+1)**2*dy(iy)**2
         amat(4,1,+1,+1,0:nx,iy)=+a10*dx(1:nx+1)**2*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(4,2,0,0,0:nx,iy)=+a15
     $        *(dx(1:nx+1)**3+dx(0:nx)**3)*(dy(iy+1)**2-dy(iy)**2)
         amat(4,2,-1,0,0:nx,iy)=-a17*(dy(iy+1)**2-dy(iy)**2)*dx(0:nx)**3
         amat(4,2,+1,0,0:nx,iy)=-a17*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(1:nx+1)**3
         amat(4,2,0,-1,0:nx,iy)=-a16*(dx(1:nx+1)**3+dx(0:nx)**3)
     $        *dy(iy)**2
         amat(4,2,0,+1,0:nx,iy)=+a16*(dx(1:nx+1)**3+dx(0:nx)**3)*
     $        dy(iy+1)**2
         amat(4,2,-1,-1,0:nx,iy)=+a18*dx(0:nx)**3*dy(iy)**2
         amat(4,2,-1,+1,0:nx,iy)=-a18*dx(0:nx)**3*dy(iy+1)**2
         amat(4,2,+1,-1,0:nx,iy)=+a18*dx(1:nx+1)**3*dy(iy)**2
         amat(4,2,+1,+1,0:nx,iy)=-a18*dx(1:nx+1)**3*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, y derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(4,3,0,0,ix,0:ny)=+a15
     $        *(dy(1:ny+1)**3+dy(0:ny)**3)*(dx(ix+1)**2-dx(ix)**2)
         amat(4,3,-1,0,ix,0:ny)=-a16*dx(ix)**2
     $        *(dy(1:ny+1)**3+dy(0:ny)**3)
         amat(4,3,+1,0,ix,0:ny)=+a16*dx(ix+1)**2
     $        *(dy(1:ny+1)**3+dy(0:ny)**3)
         amat(4,3,0,-1,ix,0:ny)=-a17*(dx(ix+1)**2-dx(ix)**2)
     $        *dy(0:ny)**3
         amat(4,3,0,+1,ix,0:ny)=-a17*(dx(ix+1)**2-dx(ix)**2)
     $        *dy(1:ny+1)**3
         amat(4,3,-1,-1,ix,0:ny)=+a18*dy(0:ny)**3*dx(ix)**2
         amat(4,3,-1,+1,ix,0:ny)=+a18*dy(1:ny+1)**3*dx(ix)**2
         amat(4,3,+1,-1,ix,0:ny)=-a18*dy(0:ny)**3*dx(ix+1)**2
         amat(4,3,+1,+1,ix,0:ny)=-a18*dy(1:ny+1)**3*dx(ix+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(4,4,0,0,0:nx,iy)=+a21
     $        *(dx(1:nx+1)**3+dx(0:nx)**3)*(dy(iy)**3+dy(iy+1)**3)
         amat(4,4,-1,0,0:nx,iy)=-a22*dx(0:nx)**3
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(4,4,+1,0,0:nx,iy)=-a22*dx(1:nx+1)**3
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(4,4,0,-1,0:nx,iy)=-a22*(dx(1:nx+1)**3+dx(0:nx)**3)
     $        *dy(iy)**3
         amat(4,4,0,+1,0:nx,iy)=-a22*(dx(1:nx+1)**3+dx(0:nx)**3)
     $        *dy(iy+1)**3
         amat(4,4,-1,-1,0:nx,iy)=+a23*dx(0:nx)**3*dy(iy)**3
         amat(4,4,-1,+1,0:nx,iy)=+a23*dx(0:nx)**3*dy(iy+1)**3
         amat(4,4,+1,-1,0:nx,iy)=+a23*dx(1:nx+1)**3*dy(iy)**3
         amat(4,4,+1,+1,0:nx,iy)=+a23*dx(1:nx+1)**3*dy(iy+1)**3
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, etas.
c-----------------------------------------------------------------------
      DO iy=0,ny,ny
         amat(4,7,-1,0,2:nx,iy)=dx(1:nx-1)
         amat(4,7,0,0,1:nx-1,iy)=2*(dx(1:nx-1)+dx(2:nx))
         amat(4,7,1,0,0:nx-2,iy)=dx(2:nx)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, nus.
c-----------------------------------------------------------------------
      DO ix=0,nx
         amat(4,7,0,-1,ix,2:ny)=dy(1:ny-1)
         amat(4,7,0,0,ix,1:ny-1)=2*(dy(1:ny-1)+dy(2:ny))
         amat(4,7,0,1,ix,0:ny-2)=dy(2:ny)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx
            rhs(4,ix,iy,:)=a24*((dx(ix)**2*g(ix,iy,:)
     $           -dx(ix+1)**2*g(ix+1,iy,:))*dy(iy)**2
     $           +(-dx(ix)**2*g(ix,iy+1,:)
     $           +dx(ix+1)**2*g(ix+1,iy+1,:))*dy(iy+1)**2)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     continuity of second x-derivatives.
c-----------------------------------------------------------------------
      DO ix=1,nx-1
         amat(5,1,-1,0,ix,0:ny)=3*dx(ix+1)/dx(ix)
         amat(5,1,0,0,ix,0:ny)=3*(dx(ix)/dx(ix+1)-dx(ix+1)/dx(ix))
         amat(5,1,1,0,ix,0:ny)=-3*dx(ix)/dx(ix+1)
         amat(5,2,-1,0,ix,0:ny)=dx(ix+1)
         amat(5,2,0,0,ix,0:ny)=2*(dx(ix)+dx(ix+1))
         amat(5,2,1,0,ix,0:ny)=dx(ix)
      ENDDO
      amat(5,5,0,0,0:nx:nx,0:ny)=1
c-----------------------------------------------------------------------
c     continuity of second y-derivatives.
c-----------------------------------------------------------------------
      DO iy=1,ny-1
         amat(6,1,0,-1,0:nx,iy)=3*dy(iy+1)/dy(iy)
         amat(6,1,0,0,0:nx,iy)=3*(dy(iy)/dy(iy+1)-dy(iy+1)/dy(iy))
         amat(6,1,0,1,0:nx,iy)=-3*dy(iy)/dy(iy+1)
         amat(6,3,0,-1,0:nx,iy)=dy(iy+1)
         amat(6,3,0,0,0:nx,iy)=2*(dy(iy)+dy(iy+1))
         amat(6,3,0,1,0:nx,iy)=dy(iy)
      ENDDO
      amat(6,6,0,0,0:nx,0:ny:ny)=1
c-----------------------------------------------------------------------
c     continuity of mixed second derivatives.
c-----------------------------------------------------------------------
      DO iy=1,ny-1
         amat(7,2,0,-1,0:nx,iy)=3*dy(iy+1)/dy(iy)
         amat(7,2,0,0,0:nx,iy)=3*(dy(iy)/dy(iy+1)-dy(iy+1)/dy(iy))
         amat(7,2,0,1,0:nx,iy)=-3*dy(iy)/dy(iy+1)
         amat(7,4,0,-1,0:nx,iy)=dy(iy+1)
         amat(7,4,0,0,0:nx,iy)=2*(dy(iy)+dy(iy+1))
         amat(7,4,0,1,0:nx,iy)=dy(iy)
      ENDDO

      DO ix=1,nx-1
         amat(7,3,-1,0,ix,0:ny:ny)=3*dx(ix+1)/dx(ix)
         amat(7,3,0,0,ix,0:ny:ny)=3*(dx(ix)/dx(ix+1)-dx(ix+1)/dx(ix))
         amat(7,3,1,0,ix,0:ny:ny)=-3*dx(ix)/dx(ix+1)
         amat(7,4,-1,0,ix,0:ny:ny)=dx(ix+1)
         amat(7,4,0,0,ix,0:ny:ny)=2*(dx(ix)+dx(ix+1))
         amat(7,4,1,0,ix,0:ny:ny)=dx(ix)
      ENDDO

      amat(7,7,0,0,0:nx:nx,0:ny:ny)=1
c-----------------------------------------------------------------------
c     transfer matrix to lapack band storage.
c-----------------------------------------------------------------------
      ALLOCATE(ab(ldab,n))
      ab=0
      DO ix=0,nx
         DO mx=MAX(-ix,-1),MIN(nx-ix,1)
            DO iy=0,ny
               DO my=MAX(-iy,-1),MIN(ny-iy,1)
                  DO k=1,7
                     i=7*(iy*(nx+1)+ix)+k
                     DO l=1,7
                        j=7*((iy+my)*(nx+1)+ix+mx)+l
                        ab(2*kd+1+i-j,j)=amat(k,l,mx,my,ix,iy)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     factor and solve.
c-----------------------------------------------------------------------
      CALL dgbtrf(n,n,kd,kd,ab,ldab,ipiv,info)
      CALL dgbtrs('N',n,kd,kd,nrhs,ab,ldab,ipiv,rhs,n,info)
      DEALLOCATE(ab)
c-----------------------------------------------------------------------
c     compute output.
c-----------------------------------------------------------------------
      bcs%fs=rhs(1,:,:,:)
      bcs%fsx=rhs(2,:,:,:)
      bcs%fsy=rhs(3,:,:,:)
      bcs%fsxy=rhs(4,:,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_lsfit
c-----------------------------------------------------------------------
c     subprogram 4. bicube_lsfit_xp.
c     least-square fit to cubic splines of piecewise-constant functions.
c     minimization with respect to f,fx,fy,fxy.
c     periodic in x boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_lsfit_xp(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs

      INTEGER :: ix,iy,nx,ny,mx,my,nrhs,n,info,i,j,k,l,kd,ldab
      INTEGER, DIMENSION(7*bcs%mx*(bcs%my+1)) :: ipiv
      REAL(r8), PARAMETER ::
     $     a1=169/1225._r8,a2=117/2450._r8,a3=81/4900._r8,
     $     a4=143/7350._r8,a5=169/14700._r8,a6=33/4900._r8,
     $     a7=39/9800._r8,a8=121/44100._r8,a9=143/88200._r8,
     $     a10=169/176400._r8,a11=13/3675._r8,a12=3/2450._r8,
     $     a13=13/4900._r8,a14=9/9800._r8,a15=11/22050._r8,
     $     a16=13/44100._r8,a17=11/29400._r8,a18=13/58800._r8,
     $     a19=1/4._r8,a20=1/24._r8,a21=1/11025._r8,
     $     a22=1/14700._r8,a23=1/19600._r8,a24=1/144._r8
      REAL(r8), DIMENSION(0:bcs%mx+1) :: dx
      REAL(r8), DIMENSION(0:bcs%my+1) :: dy
      REAL(r8), DIMENSION(7,0:bcs%mx-1,0:bcs%my,bcs%nqty) :: rhs
      REAL(r8), DIMENSION(0:bcs%mx,0:bcs%my+1,bcs%nqty) :: g
      REAL(r8), DIMENSION(7,7,1-bcs%mx:bcs%mx-1,-1:1,
     $     0:bcs%mx-1,0:bcs%my) :: amat
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ab
c-----------------------------------------------------------------------
c     define sizes.
c-----------------------------------------------------------------------
      nx=bcs%mx
      ny=bcs%my
      nrhs=bcs%nqty
      kd=14*nx-1
      ldab=3*kd+1
      n=7*nx*(ny+1)
c-----------------------------------------------------------------------
c     initialize arrays.
c-----------------------------------------------------------------------
      amat=0
      rhs=0
      dx=0
      dy=0
      g=0
      dx(1:nx)=bcs%xs(1:nx)-bcs%xs(0:nx-1)
      dx(0)=dx(nx)
      dx(nx+1)=dx(1)
      dy(1:ny)=bcs%ys(1:ny)-bcs%ys(0:ny-1)
      g(1:nx,1:ny,:)=bcs%fs(1:nx,1:ny,:)
      g(0,:,:)=g(nx,:,:)
      g(:,0,:)=g(:,1,:)
      g(:,ny+1,:)=g(:,ny,:)
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,1,0,0,0:nx-1,iy)=+a1
     $        *(dx(0:nx-1)+dx(1:nx))*(dy(iy)+dy(iy+1))
         amat(1,1,-1,0,1:nx-1,iy)=+a2*dx(1:nx-1)*(dy(iy)+dy(iy+1))
         amat(1,1,+1,0,0:nx-2,iy)=+a2*dx(1:nx-1)*(dy(iy)+dy(iy+1))
         amat(1,1,0,-1,0:nx-1,iy)=+a2*(dx(0:nx-1)+dx(1:nx))*dy(iy)
         amat(1,1,0,+1,0:nx-1,iy)=+a2*(dx(0:nx-1)+dx(1:nx))*dy(iy+1)
         amat(1,1,-1,-1,1:nx-1,iy)=+a3*dx(1:nx-1)*dy(iy)
         amat(1,1,-1,+1,1:nx-1,iy)=+a3*dx(1:nx-1)*dy(iy+1)
         amat(1,1,+1,-1,0:nx-2,iy)=+a3*dx(1:nx-1)*dy(iy)
         amat(1,1,+1,+1,0:nx-2,iy)=+a3*dx(1:nx-1)*dy(iy+1)

         amat(1,1,nx-1,0,0,iy)=+a2*dx(0)*(dy(iy)+dy(iy+1))
         amat(1,1,nx-1,-1,0,iy)=+a3*dx(0)*dy(iy)
         amat(1,1,nx-1,+1,0,iy)=+a3*dx(0)*dy(iy+1)

         amat(1,1,1-nx,0,nx-1,iy)=+a2*dx(nx)*(dy(iy)+dy(iy+1))
         amat(1,1,1-nx,-1,nx-1,iy)=+a3*dx(nx)*dy(iy)
         amat(1,1,1-nx,+1,nx-1,iy)=+a3*dx(nx)*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,2,0,0,0:nx-1,iy)=+a4
     $        *(dx(1:nx)**2-dx(0:nx-1)**2)*(dy(iy)+dy(iy+1))
         amat(1,2,-1,0,1:nx-1,iy)=+a5*(dy(iy)+dy(iy+1))*dx(1:nx-1)**2
         amat(1,2,+1,0,0:nx-2,iy)=-a5*(dy(iy)+dy(iy+1))*dx(1:nx-1)**2
         amat(1,2,0,-1,0:nx-1,iy)=+a6*(dx(1:nx)**2-dx(0:nx-1)**2)*dy(iy)
         amat(1,2,0,+1,0:nx-1,iy)=+a6*(dx(1:nx)**2-dx(0:nx-1)**2)*
     $        dy(iy+1)
         amat(1,2,-1,-1,1:nx-1,iy)=+a7*dx(1:nx-1)**2*dy(iy)
         amat(1,2,-1,+1,1:nx-1,iy)=+a7*dx(1:nx-1)**2*dy(iy+1)
         amat(1,2,+1,-1,0:nx-2,iy)=-a7*dx(1:nx-1)**2*dy(iy)
         amat(1,2,+1,+1,0:nx-2,iy)=-a7*dx(1:nx-1)**2*dy(iy+1)

         amat(1,2,nx-1,0,0,iy)=+a5*(dy(iy)+dy(iy+1))*dx(0)**2
         amat(1,2,nx-1,-1,0,iy)=+a7*dx(0)**2*dy(iy)
         amat(1,2,nx-1,+1,0,iy)=+a7*dx(0)**2*dy(iy+1)

         amat(1,2,1-nx,0,nx-1,iy)=-a5*(dy(iy)+dy(iy+1))*dx(nx)**2
         amat(1,2,1-nx,-1,nx-1,iy)=-a7*dx(nx)**2*dy(iy)
         amat(1,2,1-nx,+1,nx-1,iy)=-a7*dx(nx)**2*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, y derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         amat(1,3,0,0,ix,0:ny)=+a4
     $        *(dy(1:ny+1)**2-dy(0:ny)**2)*(dx(ix)+dx(ix+1))
         IF(ix>0)THEN
            amat(1,3,-1,0,ix,0:ny)=+a6*dx(ix)*
     $           (dy(1:ny+1)**2-dy(0:ny)**2)
            amat(1,3,-1,-1,ix,0:ny)=+a7*dy(0:ny)**2*dx(ix)
            amat(1,3,-1,+1,ix,0:ny)=-a7*dy(1:ny+1)**2*dx(ix)
         ENDIF
         amat(1,3,0,-1,ix,0:ny)=+a5*(dx(ix)+dx(ix+1))*dy(0:ny)**2
         amat(1,3,0,+1,ix,0:ny)=-a5*(dx(ix)+dx(ix+1))*dy(1:ny+1)**2
         IF(ix<nx-1)THEN
            amat(1,3,+1,0,ix,0:ny)=+a6*dx(ix+1)*
     $           (dy(1:ny+1)**2-dy(0:ny)**2)
            amat(1,3,+1,-1,ix,0:ny)=+a7*dy(0:ny)**2*dx(ix+1)
            amat(1,3,+1,+1,ix,0:ny)=-a7*dy(1:ny+1)**2*dx(ix+1)
         ENDIF
      ENDDO

         amat(1,3,nx-1,0,0,0:ny)=+a6*dx(0)*
     $        (dy(1:ny+1)**2-dy(0:ny)**2)
         amat(1,3,nx-1,-1,0,0:ny)=+a7*dy(0:ny)**2*dx(0)
         amat(1,3,nx-1,+1,0,0:ny)=-a7*dy(1:ny+1)**2*dx(0)

         amat(1,3,1-nx,0,nx-1,0:ny)=+a6*dx(nx)*
     $        (dy(1:ny+1)**2-dy(0:ny)**2)
         amat(1,3,1-nx,-1,nx-1,0:ny)=+a7*dy(0:ny)**2*dx(nx)
         amat(1,3,1-nx,+1,nx-1,0:ny)=-a7*dy(1:ny+1)**2*dx(nx)
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,4,0,0,0:nx-1,iy)=+a8
     $        *(dx(0:nx-1)**2-dx(1:nx)**2)*(dy(iy)**2-dy(iy+1)**2)
         amat(1,4,-1,0,1:nx-1,iy)=+a9*dx(1:nx-1)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(1,4,+1,0,0:nx-2,iy)=-a9*dx(1:nx-1)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(1,4,0,-1,0:nx-1,iy)=+a9*(dx(1:nx)**2-dx(0:nx-1)**2)
     $        *dy(iy)**2
         amat(1,4,0,+1,0:nx-1,iy)=-a9*(dx(1:nx)**2-dx(0:nx-1)**2)
     $        *dy(iy+1)**2
         amat(1,4,-1,-1,1:nx-1,iy)=+a10*dx(1:nx-1)**2*dy(iy)**2
         amat(1,4,-1,+1,1:nx-1,iy)=-a10*dx(1:nx-1)**2*dy(iy+1)**2
         amat(1,4,+1,-1,0:nx-2,iy)=-a10*dx(1:nx-1)**2*dy(iy)**2
         amat(1,4,+1,+1,0:nx-2,iy)=+a10*dx(1:nx-1)**2*dy(iy+1)**2

         amat(1,4,nx-1,0,0,iy)=+a9*dx(0)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(1,4,nx-1,-1,0,iy)=+a10*dx(0)**2*dy(iy)**2
         amat(1,4,nx-1,+1,0,iy)=-a10*dx(0)**2*dy(iy+1)**2

         amat(1,4,1-nx,0,nx-1,iy)=-a9*dx(nx)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(1,4,1-nx,-1,nx-1,iy)=-a10*dx(nx)**2*dy(iy)**2
         amat(1,4,1-nx,+1,nx-1,iy)=+a10*dx(nx)**2*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, lambdas.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(1,5,-1,0,1:nx-1,iy)=-3*dx(0:nx-2)/dx(1:nx-1)
         amat(1,5,0,0,0:nx-1,iy)=3*
     $        (dx(0:nx-1)/dx(1:nx)-dx(1:nx)/dx(0:nx-1))
         amat(1,5,1,0,0:nx-2,iy)=3*dx(2:nx)/dx(1:nx-1)

         amat(1,5,nx-1,0,0,iy)=-3*dx(nx-1)/dx(nx)
         amat(1,5,1-nx,0,nx-1,iy)=3*dx(nx+1)/dx(nx)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, mus.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         amat(1,6,0,-1,ix,2:ny)=-3*dy(1:ny-1)/dy(2:ny)
         amat(1,6,0,0,ix,1:ny-1)=3*
     $        (dy(1:ny-1)/dy(2:ny)-dy(2:ny)/dy(1:ny-1))
         amat(1,6,0,1,ix,0:ny-2)=3*dy(2:ny)/dy(1:ny-1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. f, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx-1
            rhs(1,ix,iy,:)=a19*((dx(ix)*g(ix,iy,:)
     $           +dx(ix+1)*g(ix+1,iy,:))*dy(iy)
     $           +(dx(ix)*g(ix,iy+1,:)
     $           +dx(ix+1)*g(ix+1,iy+1,:))*dy(iy+1))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(2,1,0,0,0:nx-1,iy)=+a4
     $        *(dx(1:nx)**2-dx(0:nx-1)**2)*(dy(iy)+dy(iy+1))
         amat(2,1,-1,0,1:nx-1,iy)=-a5*dx(1:nx-1)**2*(dy(iy)+dy(iy+1))
         amat(2,1,+1,0,0:nx-2,iy)=+a5*dx(1:nx-1)**2*(dy(iy)+dy(iy+1))
         amat(2,1,0,-1,0:nx-1,iy)=+a6*(dx(1:nx)**2-dx(0:nx-1)**2)*dy(iy)
         amat(2,1,0,+1,0:nx-1,iy)=+a6*(dx(1:nx)**2-dx(0:nx-1)**2)
     $        *dy(iy+1)
         amat(2,1,-1,-1,1:nx-1,iy)=-a7*dx(1:nx-1)**2*dy(iy)
         amat(2,1,-1,+1,1:nx-1,iy)=-a7*dx(1:nx-1)**2*dy(iy+1)
         amat(2,1,+1,-1,0:nx-2,iy)=+a7*dx(1:nx-1)**2*dy(iy)
         amat(2,1,+1,+1,0:nx-2,iy)=+a7*dx(1:nx-1)**2*dy(iy+1)

         amat(2,1,nx-1,0,0,iy)=-a5*dx(0)**2*(dy(iy)+dy(iy+1))
         amat(2,1,nx-1,-1,0,iy)=-a7*dx(0)**2*dy(iy)
         amat(2,1,nx-1,+1,0,iy)=-a7*dx(0)**2*dy(iy+1)

         amat(2,1,1-nx,0,nx-1,iy)=+a5*dx(nx)**2*(dy(iy)+dy(iy+1))
         amat(2,1,1-nx,-1,nx-1,iy)=+a7*dx(nx)**2*dy(iy)
         amat(2,1,1-nx,+1,nx-1,iy)=+a7*dx(nx)**2*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(2,2,0,0,0:nx-1,iy)=+a11
     $        *(dx(1:nx)**3+dx(0:nx-1)**3)*(dy(iy)+dy(iy+1))
         amat(2,2,-1,0,1:nx-1,iy)=-a13*(dy(iy)+dy(iy+1))*dx(1:nx-1)**3
         amat(2,2,+1,0,0:nx-2,iy)=-a13*(dy(iy)+dy(iy+1))*dx(1:nx-1)**3
         amat(2,2,0,-1,0:nx-1,iy)=+a12*(dx(1:nx)**3+dx(0:nx-1)**3)
     $        *dy(iy)
         amat(2,2,0,+1,0:nx-1,iy)=+a12*(dx(1:nx)**3+dx(0:nx-1)**3)*
     $        dy(iy+1)
         amat(2,2,-1,-1,1:nx-1,iy)=-a14*dx(1:nx-1)**3*dy(iy)
         amat(2,2,-1,+1,1:nx-1,iy)=-a14*dx(1:nx-1)**3*dy(iy+1)
         amat(2,2,+1,-1,0:nx-2,iy)=-a14*dx(1:nx-1)**3*dy(iy)
         amat(2,2,+1,+1,0:nx-2,iy)=-a14*dx(1:nx-1)**3*dy(iy+1)

         amat(2,2,nx-1,0,0,iy)=-a13*(dy(iy)+dy(iy+1))*dx(0)**3
         amat(2,2,nx-1,-1,0,iy)=-a14*dx(0)**3*dy(iy)
         amat(2,2,nx-1,+1,0,iy)=-a14*dx(0)**3*dy(iy+1)

         amat(2,2,1-nx,0,nx-1,iy)=-a13*(dy(iy)+dy(iy+1))*dx(nx)**3
         amat(2,2,1-nx,-1,nx-1,iy)=-a14*dx(nx)**3*dy(iy)
         amat(2,2,1-nx,+1,nx-1,iy)=-a14*dx(nx)**3*dy(iy+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, y derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         amat(2,3,0,0,ix,0:ny)=+a8
     $        *(dy(1:ny+1)**2-dy(0:ny)**2)*(dx(ix+1)**2-dx(ix)**2)
         IF(ix>0)THEN
            amat(2,3,-1,0,ix,0:ny)=-a9*dx(ix)**2*
     $           (dy(1:ny+1)**2-dy(0:ny)**2)
            amat(2,3,-1,-1,ix,0:ny)=-a10*dy(0:ny)**2*dx(ix)**2
            amat(2,3,-1,+1,ix,0:ny)=+a10*dy(1:ny+1)**2*dx(ix)**2
         ENDIF
         amat(2,3,0,-1,ix,0:ny)=+a9*(dx(ix+1)**2-dx(ix)**2)*dy(0:ny)**2
         amat(2,3,0,+1,ix,0:ny)=-a9*(dx(ix+1)**2-dx(ix)**2)
     $        *dy(1:ny+1)**2
         IF(ix<nx-1)THEN
            amat(2,3,+1,0,ix,0:ny)=+a9*dx(ix+1)**2*
     $           (dy(1:ny+1)**2-dy(0:ny)**2)
            amat(2,3,+1,-1,ix,0:ny)=+a10*dy(0:ny)**2*dx(ix+1)**2
            amat(2,3,+1,+1,ix,0:ny)=-a10*dy(1:ny+1)**2*dx(ix+1)**2
         ENDIF
      ENDDO

      amat(2,3,nx-1,0,0,0:ny)=-a9*dx(0)**2*
     $     (dy(1:ny+1)**2-dy(0:ny)**2)
      amat(2,3,nx-1,-1,0,0:ny)=-a10*dy(0:ny)**2*dx(0)**2
      amat(2,3,nx-1,+1,0,0:ny)=+a10*dy(1:ny+1)**2*dx(0)**2

      amat(2,3,1-nx,0,nx-1,0:ny)=+a9*dx(nx)**2*
     $     (dy(1:ny+1)**2-dy(0:ny)**2)
      amat(2,3,1-nx,-1,nx-1,0:ny)=+a10*dy(0:ny)**2*dx(nx)**2
      amat(2,3,1-nx,+1,nx-1,0:ny)=-a10*dy(1:ny+1)**2*dx(nx)**2
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(2,4,0,0,0:nx-1,iy)=+a15
     $        *(dx(0:nx-1)**3+dx(1:nx)**3)*(dy(iy+1)**2-dy(iy)**2)
         amat(2,4,-1,0,1:nx-1,iy)=-a17*dx(1:nx-1)**3
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(2,4,+1,0,0:nx-2,iy)=-a17*dx(1:nx-1)**3
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(2,4,0,-1,0:nx-1,iy)=+a16*(dx(1:nx)**3+dx(0:nx-1)**3)
     $        *dy(iy)**2
         amat(2,4,0,+1,0:nx-1,iy)=-a16*(dx(1:nx)**3+dx(0:nx-1)**3)
     $        *dy(iy+1)**2
         amat(2,4,-1,-1,1:nx-1,iy)=-a18*dx(1:nx-1)**3*dy(iy)**2
         amat(2,4,-1,+1,1:nx-1,iy)=+a18*dx(1:nx-1)**3*dy(iy+1)**2
         amat(2,4,+1,-1,0:nx-2,iy)=-a18*dx(1:nx-1)**3*dy(iy)**2
         amat(2,4,+1,+1,0:nx-2,iy)=+a18*dx(1:nx-1)**3*dy(iy+1)**2

         amat(2,4,nx-1,0,0,iy)=-a17*dx(0)**3
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(2,4,nx-1,-1,0,iy)=-a18*dx(0)**3*dy(iy)**2
         amat(2,4,nx-1,+1,0,iy)=+a18*dx(0)**3*dy(iy+1)**2

         amat(2,4,1-nx,0,nx-1,iy)=-a17*dx(nx)**3
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(2,4,1-nx,-1,nx-1,iy)=-a18*dx(nx)**3*dy(iy)**2
         amat(2,4,1-nx,+1,nx-1,iy)=+a18*dx(nx)**3*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, lambdas.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(2,5,-1,0,1:nx-1,iy)=dx(0:nx-2)
         amat(2,5,0,0,0:nx-1,iy)=2*(dx(0:nx-1)+dx(1:nx))
         amat(2,5,1,0,0:nx-2,iy)=dx(2:nx)

         amat(2,5,nx-1,0,0,iy)=dx(nx-1)
         amat(2,5,1-nx,0,nx-1,iy)=dx(nx+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, nus.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         amat(2,7,0,-1,ix,2:ny)=-3*dy(1:ny-1)/dy(2:ny)
         amat(2,7,0,0,ix,1:ny-1)=3*
     $        (dy(1:ny-1)/dy(2:ny)-dy(2:ny)/dy(1:ny-1))
         amat(2,7,0,1,ix,0:ny-2)=3*dy(2:ny)/dy(1:ny-1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fx, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx-1
            rhs(2,ix,iy,:)=a20*((-dx(ix)**2*g(ix,iy,:)
     $           +dx(ix+1)**2*g(ix+1,iy,:))*dy(iy)
     $           +(-dx(ix)**2*g(ix,iy+1,:)
     $           +dx(ix+1)**2*g(ix+1,iy+1,:))*dy(iy+1))
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(3,1,0,0,0:nx-1,iy)=+a4
     $        *(dx(0:nx-1)+dx(1:nx))*(dy(iy+1)**2-dy(iy)**2)
         amat(3,1,-1,0,1:nx-1,iy)=+a6*dx(1:nx-1)*(dy(iy+1)**2-dy(iy)**2)
         amat(3,1,+1,0,0:nx-2,iy)=+a6*dx(1:nx-1)*(dy(iy+1)**2-dy(iy)**2)
         amat(3,1,0,-1,0:nx-1,iy)=-a5*(dx(0:nx-1)+dx(1:nx))*dy(iy)**2
         amat(3,1,0,+1,0:nx-1,iy)=+a5*(dx(0:nx-1)+dx(1:nx))*dy(iy+1)**2
         amat(3,1,-1,-1,1:nx-1,iy)=-a7*dx(1:nx-1)*dy(iy)**2
         amat(3,1,-1,+1,1:nx-1,iy)=+a7*dx(1:nx-1)*dy(iy+1)**2
         amat(3,1,+1,-1,0:nx-2,iy)=-a7*dx(1:nx-1)*dy(iy)**2
         amat(3,1,+1,+1,0:nx-2,iy)=+a7*dx(1:nx-1)*dy(iy+1)**2

         amat(3,1,nx-1,0,0,iy)=+a6*dx(0)*(dy(iy+1)**2-dy(iy)**2)
         amat(3,1,nx-1,-1,0,iy)=-a7*dx(0)*dy(iy)**2
         amat(3,1,nx-1,+1,0,iy)=+a7*dx(0)*dy(iy+1)**2

         amat(3,1,1-nx,0,nx-1,iy)=+a6*dx(nx)*(dy(iy+1)**2-dy(iy)**2)
         amat(3,1,1-nx,-1,nx-1,iy)=-a7*dx(nx)*dy(iy)**2
         amat(3,1,1-nx,+1,nx-1,iy)=+a7*dx(nx)*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(3,2,0,0,0:nx-1,iy)=+a8
     $        *(dx(1:nx)**2-dx(0:nx-1)**2)*(dy(iy+1)**2-dy(iy)**2)
         amat(3,2,-1,0,1:nx-1,iy)=+a9*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(1:nx-1)**2
         amat(3,2,+1,0,0:nx-2,iy)=-a9*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(1:nx-1)**2
         amat(3,2,0,-1,0:nx-1,iy)=-a9*(dx(1:nx)**2-dx(0:nx-1)**2)
     $        *dy(iy)**2
         amat(3,2,0,+1,0:nx-1,iy)=+a9*(dx(1:nx)**2-dx(0:nx-1)**2)*
     $        dy(iy+1)**2
         amat(3,2,-1,-1,1:nx-1,iy)=-a10*dx(1:nx-1)**2*dy(iy)**2
         amat(3,2,-1,+1,1:nx-1,iy)=+a10*dx(1:nx-1)**2*dy(iy+1)**2
         amat(3,2,+1,-1,0:nx-2,iy)=+a10*dx(1:nx-1)**2*dy(iy)**2
         amat(3,2,+1,+1,0:nx-2,iy)=-a10*dx(1:nx-1)**2*dy(iy+1)**2

         amat(3,2,nx-1,0,0,iy)=+a9*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(0)**2
         amat(3,2,nx-1,-1,0,iy)=-a10*dx(0)**2*dy(iy)**2
         amat(3,2,nx-1,+1,0,iy)=+a10*dx(0)**2*dy(iy+1)**2

         amat(3,2,1-nx,0,nx-1,iy)=-a9*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(nx)**2
         amat(3,2,1-nx,-1,nx-1,iy)=+a10*dx(nx)**2*dy(iy)**2
         amat(3,2,1-nx,+1,nx-1,iy)=-a10*dx(nx)**2*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, y derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         amat(3,3,0,0,ix,0:ny)=+a11
     $        *(dy(1:ny+1)**3+dy(0:ny)**3)*(dx(ix)+dx(ix+1))
         IF(ix>0)THEN
            amat(3,3,-1,0,ix,0:ny)=+a12*dx(ix)*
     $           (dy(1:ny+1)**3+dy(0:ny)**3)
            amat(3,3,-1,-1,ix,0:ny)=-a14*dy(0:ny)**3*dx(ix)
            amat(3,3,-1,+1,ix,0:ny)=-a14*dy(1:ny+1)**3*dx(ix)
         ENDIF
         amat(3,3,0,-1,ix,0:ny)=-a13*(dx(ix)+dx(ix+1))*dy(0:ny)**3
         amat(3,3,0,+1,ix,0:ny)=-a13*(dx(ix)+dx(ix+1))*dy(1:ny+1)**3
         IF(ix<nx-1)THEN
            amat(3,3,+1,0,ix,0:ny)=+a12*dx(ix+1)*
     $           (dy(1:ny+1)**3+dy(0:ny)**3)
            amat(3,3,+1,-1,ix,0:ny)=-a14*dy(0:ny)**3*dx(ix+1)
            amat(3,3,+1,+1,ix,0:ny)=-a14*dy(1:ny+1)**3*dx(ix+1)
         ENDIF
      ENDDO

      amat(3,3,nx-1,0,0,0:ny)=+a12*dx(0)*
     $     (dy(1:ny+1)**3+dy(0:ny)**3)
      amat(3,3,nx-1,-1,0,0:ny)=-a14*dy(0:ny)**3*dx(0)
      amat(3,3,nx-1,+1,0,0:ny)=-a14*dy(1:ny+1)**3*dx(0)

      amat(3,3,1-nx,0,nx-1,0:ny)=+a12*dx(nx)*
     $     (dy(1:ny+1)**3+dy(0:ny)**3)
      amat(3,3,1-nx,-1,nx-1,0:ny)=-a14*dy(0:ny)**3*dx(nx)
      amat(3,3,1-nx,+1,nx-1,0:ny)=-a14*dy(1:ny+1)**3*dx(nx)
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(3,4,0,0,0:nx-1,iy)=+a15
     $        *(dx(1:nx)**2-dx(0:nx-1)**2)*(dy(iy)**3+dy(iy+1)**3)
         amat(3,4,-1,0,1:nx-1,iy)=+a16*dx(1:nx-1)**2
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(3,4,+1,0,0:nx-2,iy)=-a16*dx(1:nx-1)**2
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(3,4,0,-1,0:nx-1,iy)=-a17*(dx(1:nx)**2-dx(0:nx-1)**2)
     $        *dy(iy)**3
         amat(3,4,0,+1,0:nx-1,iy)=-a17*(dx(1:nx)**2-dx(0:nx-1)**2)
     $        *dy(iy+1)**3
         amat(3,4,-1,-1,1:nx-1,iy)=-a18*dx(1:nx-1)**2*dy(iy)**3
         amat(3,4,-1,+1,1:nx-1,iy)=-a18*dx(1:nx-1)**2*dy(iy+1)**3
         amat(3,4,+1,-1,0:nx-2,iy)=+a18*dx(1:nx-1)**2*dy(iy)**3
         amat(3,4,+1,+1,0:nx-2,iy)=+a18*dx(1:nx-1)**2*dy(iy+1)**3

         amat(3,4,nx-1,0,0,iy)=+a16*dx(0)**2
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(3,4,nx-1,-1,0,iy)=-a18*dx(0)**2*dy(iy)**3
         amat(3,4,nx-1,+1,0,iy)=-a18*dx(0)**2*dy(iy+1)**3

         amat(3,4,1-nx,0,nx-1,iy)=-a16*dx(nx)**2
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(3,4,1-nx,-1,nx-1,iy)=+a18*dx(nx)**2*dy(iy)**3
         amat(3,4,1-nx,+1,nx-1,iy)=+a18*dx(nx)**2*dy(iy+1)**3
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, etas.
c-----------------------------------------------------------------------
      DO iy=0,ny,ny
         amat(3,7,-1,0,1:nx-1,iy)=-3*dx(0:nx-2)/dx(1:nx-1)
         amat(3,7,0,0,0:nx-1,iy)=3*
     $        (dx(0:nx-1)/dx(1:nx)-dx(1:nx)/dx(0:nx-1))
         amat(3,7,1,0,0:nx-2,iy)=3*dx(2:nx)/dx(1:nx-1)

         amat(3,7,nx-1,0,0,iy)=-3*dx(nx-1)/dx(nx)
         amat(3,7,1-nx,0,nx-1,iy)=3*dx(nx+1)/dx(nx)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, mus.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         amat(3,6,0,-1,ix,2:ny)=dy(1:ny-1)
         amat(3,6,0,0,ix,1:ny-1)=2*(dy(1:ny-1)+dy(2:ny))
         amat(3,6,0,1,ix,0:ny-2)=dy(2:ny)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx-1
            rhs(3,ix,iy,:)=a20*((-dx(ix)*g(ix,iy,:)
     $           -dx(ix+1)*g(ix+1,iy,:))*dy(iy)**2
     $           +(dx(ix)*g(ix,iy+1,:)
     $           +dx(ix+1)*g(ix+1,iy+1,:))*dy(iy+1)**2)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, function values.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(4,1,0,0,0:nx-1,iy)=+a8
     $        *(dx(1:nx)**2-dx(0:nx-1)**2)*(dy(iy+1)**2-dy(iy)**2)
         amat(4,1,-1,0,1:nx-1,iy)=-a9*dx(1:nx-1)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(4,1,+1,0,0:nx-2,iy)=+a9*dx(1:nx-1)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(4,1,0,-1,0:nx-1,iy)=-a9*(dx(1:nx)**2-dx(0:nx-1)**2)
     $        *dy(iy)**2
         amat(4,1,0,+1,0:nx-1,iy)=+a9*(dx(1:nx)**2-dx(0:nx-1)**2)
     $        *dy(iy+1)**2
         amat(4,1,-1,-1,1:nx-1,iy)=+a10*dx(1:nx-1)**2*dy(iy)**2
         amat(4,1,-1,+1,1:nx-1,iy)=-a10*dx(1:nx-1)**2*dy(iy+1)**2
         amat(4,1,+1,-1,0:nx-2,iy)=-a10*dx(1:nx-1)**2*dy(iy)**2
         amat(4,1,+1,+1,0:nx-2,iy)=+a10*dx(1:nx-1)**2*dy(iy+1)**2

         amat(4,1,nx-1,0,0,iy)=-a9*dx(0)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(4,1,nx-1,-1,0,iy)=+a10*dx(0)**2*dy(iy)**2
         amat(4,1,nx-1,+1,0,iy)=-a10*dx(0)**2*dy(iy+1)**2

         amat(4,1,1-nx,0,nx-1,iy)=+a9*dx(nx)**2
     $        *(dy(iy+1)**2-dy(iy)**2)
         amat(4,1,1-nx,-1,nx-1,iy)=-a10*dx(nx)**2*dy(iy)**2
         amat(4,1,1-nx,+1,nx-1,iy)=+a10*dx(nx)**2*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, x derivatiives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(4,2,0,0,0:nx-1,iy)=+a15
     $        *(dx(1:nx)**3+dx(0:nx-1)**3)*(dy(iy+1)**2-dy(iy)**2)
         amat(4,2,-1,0,1:nx-1,iy)=-a17*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(1:nx-1)**3
         amat(4,2,+1,0,0:nx-2,iy)=-a17*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(1:nx-1)**3
         amat(4,2,0,-1,0:nx-1,iy)=-a16*(dx(1:nx)**3+dx(0:nx-1)**3)
     $        *dy(iy)**2
         amat(4,2,0,+1,0:nx-1,iy)=+a16*(dx(1:nx)**3+dx(0:nx-1)**3)*
     $        dy(iy+1)**2
         amat(4,2,-1,-1,1:nx-1,iy)=+a18*dx(1:nx-1)**3*dy(iy)**2
         amat(4,2,-1,+1,1:nx-1,iy)=-a18*dx(1:nx-1)**3*dy(iy+1)**2
         amat(4,2,+1,-1,0:nx-2,iy)=+a18*dx(1:nx-1)**3*dy(iy)**2
         amat(4,2,+1,+1,0:nx-2,iy)=-a18*dx(1:nx-1)**3*dy(iy+1)**2

         amat(4,2,nx-1,0,0,iy)=-a17*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(0)**3
         amat(4,2,nx-1,-1,0,iy)=+a18*dx(0)**3*dy(iy)**2
         amat(4,2,nx-1,+1,0,iy)=-a18*dx(0)**3*dy(iy+1)**2

         amat(4,2,1-nx,0,nx-1,iy)=-a17*(dy(iy+1)**2-dy(iy)**2)
     $        *dx(nx)**3
         amat(4,2,1-nx,-1,nx-1,iy)=+a18*dx(nx)**3*dy(iy)**2
         amat(4,2,1-nx,+1,nx-1,iy)=-a18*dx(nx)**3*dy(iy+1)**2
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, y derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         amat(4,3,0,0,ix,0:ny)=+a15
     $        *(dy(1:ny+1)**3+dy(0:ny)**3)*(dx(ix+1)**2-dx(ix)**2)
         IF(ix>0)THEN
            amat(4,3,-1,0,ix,0:ny)=-a16*dx(ix)**2
     $           *(dy(1:ny+1)**3+dy(0:ny)**3)
            amat(4,3,-1,-1,ix,0:ny)=+a18*dy(0:ny)**3*dx(ix)**2
            amat(4,3,-1,+1,ix,0:ny)=+a18*dy(1:ny+1)**3*dx(ix)**2
         ENDIF
         amat(4,3,0,-1,ix,0:ny)=-a17*(dx(ix+1)**2-dx(ix)**2)
     $        *dy(0:ny)**3
         amat(4,3,0,+1,ix,0:ny)=-a17*(dx(ix+1)**2-dx(ix)**2)
     $        *dy(1:ny+1)**3
         IF(ix<nx-1)THEN
            amat(4,3,+1,0,ix,0:ny)=+a16*dx(ix+1)**2
     $           *(dy(1:ny+1)**3+dy(0:ny)**3)
            amat(4,3,+1,-1,ix,0:ny)=-a18*dy(0:ny)**3*dx(ix+1)**2
            amat(4,3,+1,+1,ix,0:ny)=-a18*dy(1:ny+1)**3*dx(ix+1)**2
         ENDIF
      ENDDO

      amat(4,3,nx-1,0,0,0:ny)=-a16*dx(0)**2
     $     *(dy(1:ny+1)**3+dy(0:ny)**3)
      amat(4,3,nx-1,-1,0,0:ny)=+a18*dy(0:ny)**3*dx(0)**2
      amat(4,3,nx-1,+1,0,0:ny)=+a18*dy(1:ny+1)**3*dx(0)**2

      amat(4,3,1-nx,0,nx-1,0:ny)=+a16*dx(nx)**2
     $     *(dy(1:ny+1)**3+dy(0:ny)**3)
      amat(4,3,1-nx,-1,nx-1,0:ny)=-a18*dy(0:ny)**3*dx(nx)**2
      amat(4,3,1-nx,+1,nx-1,0:ny)=-a18*dy(1:ny+1)**3*dx(nx)**2
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, mixed derivatives.
c-----------------------------------------------------------------------
      DO iy=0,ny
         amat(4,4,0,0,0:nx-1,iy)=+a21
     $        *(dx(1:nx)**3+dx(0:nx-1)**3)*(dy(iy)**3+dy(iy+1)**3)
         amat(4,4,-1,0,1:nx-1,iy)=-a22*dx(1:nx-1)**3
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(4,4,+1,0,0:nx-2,iy)=-a22*dx(1:nx-1)**3
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(4,4,0,-1,0:nx-1,iy)=-a22*(dx(1:nx)**3+dx(0:nx-1)**3)
     $        *dy(iy)**3
         amat(4,4,0,+1,0:nx-1,iy)=-a22*(dx(1:nx)**3+dx(0:nx-1)**3)
     $        *dy(iy+1)**3
         amat(4,4,-1,-1,1:nx-1,iy)=+a23*dx(1:nx-1)**3*dy(iy)**3
         amat(4,4,-1,+1,1:nx-1,iy)=+a23*dx(1:nx-1)**3*dy(iy+1)**3
         amat(4,4,+1,-1,0:nx-2,iy)=+a23*dx(1:nx-1)**3*dy(iy)**3
         amat(4,4,+1,+1,0:nx-2,iy)=+a23*dx(1:nx-1)**3*dy(iy+1)**3

         amat(4,4,nx-1,0,0,iy)=-a22*dx(0)**3
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(4,4,nx-1,-1,0,iy)=+a23*dx(0)**3*dy(iy)**3
         amat(4,4,nx-1,+1,0,iy)=+a23*dx(0)**3*dy(iy+1)**3

         amat(4,4,1-nx,0,nx-1,iy)=-a22*dx(nx)**3
     $        *(dy(iy+1)**3+dy(iy)**3)
         amat(4,4,1-nx,-1,nx-1,iy)=+a23*dx(nx)**3*dy(iy)**3
         amat(4,4,1-nx,+1,nx-1,iy)=+a23*dx(nx)**3*dy(iy+1)**3
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, etas.
c-----------------------------------------------------------------------
      DO iy=0,ny,ny
         amat(4,7,-1,0,1:nx-1,iy)=dx(0:nx-2)
         amat(4,7,0,0,0:nx-1,iy)=2*(dx(0:nx-1)+dx(1:nx))
         amat(4,7,1,0,0:nx-2,iy)=dx(2:nx)

         amat(4,7,nx-1,0,0,iy)=dx(nx-1)
         amat(4,7,1-nx,0,nx-1,iy)=dx(nx+1)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fxy, nus.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         amat(4,7,0,-1,ix,2:ny)=dy(1:ny-1)
         amat(4,7,0,0,ix,1:ny-1)=2*(dy(1:ny-1)+dy(2:ny))
         amat(4,7,0,1,ix,0:ny-2)=dy(2:ny)
      ENDDO
c-----------------------------------------------------------------------
c     lsq, minimizing w.r.t. fy, rhs.
c-----------------------------------------------------------------------
      DO iy=0,ny
         DO ix=0,nx-1
            rhs(4,ix,iy,:)=a24*((dx(ix)**2*g(ix,iy,:)
     $           -dx(ix+1)**2*g(ix+1,iy,:))*dy(iy)**2
     $           +(-dx(ix)**2*g(ix,iy+1,:)
     $           +dx(ix+1)**2*g(ix+1,iy+1,:))*dy(iy+1)**2)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     continuity of second x-derivatives.
c-----------------------------------------------------------------------
      DO ix=0,nx-1
         IF(ix>0)THEN
            amat(5,1,-1,0,ix,0:ny)=3*dx(ix+1)/dx(ix)
            amat(5,2,-1,0,ix,0:ny)=dx(ix+1)
         ENDIF
         amat(5,1,0,0,ix,0:ny)=3*(dx(ix)/dx(ix+1)-dx(ix+1)/dx(ix))
         amat(5,2,0,0,ix,0:ny)=2*(dx(ix)+dx(ix+1))
         IF(ix<nx-1)THEN
            amat(5,1,1,0,ix,0:ny)=-3*dx(ix)/dx(ix+1)
            amat(5,2,1,0,ix,0:ny)=dx(ix)
         ENDIF
      ENDDO

      amat(5,1,nx-1,0,0,0:ny)=3*dx(1)/dx(0)
      amat(5,2,nx-1,0,0,0:ny)=dx(1)

      amat(5,1,1-nx,0,nx-1,0:ny)=-3*dx(nx-1)/dx(nx)
      amat(5,2,1-nx,0,nx-1,0:ny)=dx(nx-1)
c-----------------------------------------------------------------------
c     continuity of second y-derivatives.
c-----------------------------------------------------------------------
      DO iy=1,ny-1
         amat(6,1,0,-1,0:nx-1,iy)=3*dy(iy+1)/dy(iy)
         amat(6,1,0,0,0:nx-1,iy)=3*(dy(iy)/dy(iy+1)-dy(iy+1)/dy(iy))
         amat(6,1,0,1,0:nx-1,iy)=-3*dy(iy)/dy(iy+1)
         amat(6,3,0,-1,0:nx-1,iy)=dy(iy+1)
         amat(6,3,0,0,0:nx-1,iy)=2*(dy(iy)+dy(iy+1))
         amat(6,3,0,1,0:nx-1,iy)=dy(iy)
      ENDDO
      amat(6,6,0,0,0:nx-1,0:ny:ny)=1
c-----------------------------------------------------------------------
c     continuity of mixed second derivatives.
c-----------------------------------------------------------------------
      DO iy=1,ny-1
         amat(7,2,0,-1,0:nx-1,iy)=3*dy(iy+1)/dy(iy)
         amat(7,2,0,0,0:nx-1,iy)=3*(dy(iy)/dy(iy+1)-dy(iy+1)/dy(iy))
         amat(7,2,0,1,0:nx-1,iy)=-3*dy(iy)/dy(iy+1)
         amat(7,4,0,-1,0:nx-1,iy)=dy(iy+1)
         amat(7,4,0,0,0:nx-1,iy)=2*(dy(iy)+dy(iy+1))
         amat(7,4,0,1,0:nx-1,iy)=dy(iy)
      ENDDO

      DO ix=0,nx-1
         IF(ix>0)THEN
            amat(7,3,-1,0,ix,0:ny:ny)=3*dx(ix+1)/dx(ix)
            amat(7,4,-1,0,ix,0:ny:ny)=dx(ix+1)
         ENDIF
         amat(7,3,0,0,ix,0:ny:ny)=3*(dx(ix)/dx(ix+1)-dx(ix+1)/dx(ix))
         amat(7,4,0,0,ix,0:ny:ny)=2*(dx(ix)+dx(ix+1))
         IF(ix<nx-1)THEN
            amat(7,3,1,0,ix,0:ny:ny)=-3*dx(ix)/dx(ix+1)
            amat(7,4,1,0,ix,0:ny:ny)=dx(ix)
         ENDIF
      ENDDO

      amat(7,3,nx-1,0,0,0:ny:ny)=3*dx(1)/dx(0)
      amat(7,4,nx-1,0,0,0:ny:ny)=dx(1)

      amat(7,3,1-nx,0,nx-1,0:ny:ny)=-3*dx(nx-1)/dx(nx)
      amat(7,4,1-nx,0,nx-1,0:ny:ny)=dx(nx-1)
c-----------------------------------------------------------------------
c     transfer matrix to lapack band storage.
c-----------------------------------------------------------------------
      ALLOCATE(ab(ldab,n))
      ab=0
      DO ix=0,nx-1
         DO mx=MAX(-ix,-nx+1),MIN(nx-1-ix,nx-1)
            DO iy=0,ny
               DO my=MAX(-iy,-1),MIN(ny-iy,1)
                  DO k=1,7
                     i=7*(iy*nx+ix)+k
                     DO l=1,7
                        j=7*((iy+my)*nx+ix+mx)+l
                        ab(2*kd+1+i-j,j)=amat(k,l,mx,my,ix,iy)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     factor and solve.
c-----------------------------------------------------------------------
      CALL dgbtrf(n,n,kd,kd,ab,ldab,ipiv,info)
      CALL dgbtrs('N',n,kd,kd,nrhs,ab,ldab,ipiv,rhs,n,info)
      DEALLOCATE(ab)
c-----------------------------------------------------------------------
c     compute output.
c-----------------------------------------------------------------------
      bcs%fs(0:nx-1,:,:)=rhs(1,0:nx-1,:,:)
      bcs%fsx(0:nx-1,:,:)=rhs(2,0:nx-1,:,:)
      bcs%fsy(0:nx-1,:,:)=rhs(3,0:nx-1,:,:)
      bcs%fsxy(0:nx-1,:,:)=rhs(4,0:nx-1,:,:)
      bcs%fs(nx,:,:)=rhs(1,0,:,:)
      bcs%fsx(nx,:,:)=rhs(2,0,:,:)
      bcs%fsy(nx,:,:)=rhs(3,0,:,:)
      bcs%fsxy(nx,:,:)=rhs(4,0,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_lsfit_xp
c-----------------------------------------------------------------------
c     subprogram 5. bicube_lsfit_yp.
c     least-square fit to cubic splines of piecewise-constant functions.
c     minimization with respect to f,fx,fy,fxy.
c     periodic in y boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_lsfit_yp(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs

      TYPE(bicube_type) :: temp
      INTEGER :: iqty
c-----------------------------------------------------------------------
c     initiate temporary bicubic spline.
c-----------------------------------------------------------------------
      CALL bicube_alloc(temp,bcs%my,bcs%mx,bcs%nqty)
c-----------------------------------------------------------------------
c     switch x and y directions.
c-----------------------------------------------------------------------
      temp%xs=bcs%ys
      temp%ys=bcs%xs
      DO iqty=1,bcs%nqty
         temp%fs(:,:,iqty)=TRANSPOSE(bcs%fs(:,:,iqty))
      ENDDO
c-----------------------------------------------------------------------
c     fit temporary bicubic spline.
c-----------------------------------------------------------------------
      CALL bicube_lsfit_xp(temp)
c-----------------------------------------------------------------------
c     return to original bicubic spline.
c-----------------------------------------------------------------------
      DO iqty=1,bcs%nqty
         bcs%fs(:,:,iqty)=TRANSPOSE(temp%fs(:,:,iqty))
         bcs%fsx(:,:,iqty)=TRANSPOSE(temp%fsy(:,:,iqty))
         bcs%fsy(:,:,iqty)=TRANSPOSE(temp%fsx(:,:,iqty))
         bcs%fsxy(:,:,iqty)=TRANSPOSE(temp%fsxy(:,:,iqty))
      ENDDO
c-----------------------------------------------------------------------
c     deallocate temporary bicubic spline.
c-----------------------------------------------------------------------
      CALL bicube_dealloc(temp)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_lsfit_yp
c-----------------------------------------------------------------------
c     subprogram 6. bicube_eval.
c     evaluates bicubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_eval(bcs,x,y)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      REAL(r8), INTENT(IN) :: x,y

      INTEGER :: i,iqty,iside
      REAL(r8) :: dx,dy,xx,yy,g,xfac,yfac
      REAL(r8), DIMENSION (4,4,bcs%nqty) :: c

c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      bcs%ix=MAX(bcs%ix,0)
      bcs%ix=MIN(bcs%ix,bcs%mx-1)
      bcs%iy=MAX(bcs%iy,0)
      bcs%iy=MIN(bcs%iy,bcs%my-1)
      xx=x
      yy=y
c-----------------------------------------------------------------------
c     normalize x interval for periodic splines.
c-----------------------------------------------------------------------
      IF(bcs%periodic(1))THEN
         DO
            IF(xx < bcs%xs(bcs%mx))EXIT
            xx=xx-bcs%xs(bcs%mx)
         ENDDO
         DO
            IF(xx >= bcs%xs(0))EXIT
            xx=xx+bcs%xs(bcs%mx)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find x interval.
c-----------------------------------------------------------------------
      DO
         IF(xx >= bcs%xs(bcs%ix) .OR. bcs%ix <= 0)EXIT
         bcs%ix=bcs%ix-1
      ENDDO
      DO
         IF(xx < bcs%xs(bcs%ix+1) .OR. bcs%ix >= bcs%mx-1)EXIT
         bcs%ix=bcs%ix+1
      ENDDO
c-----------------------------------------------------------------------
c     normalize y interval for periodic splines.
c-----------------------------------------------------------------------
      IF(bcs%periodic(2))THEN
         DO
            IF(yy < bcs%ys(bcs%my))EXIT
            yy=yy-bcs%ys(bcs%my)
         ENDDO
         DO
            IF(yy >= bcs%ys(0))EXIT
            yy=yy+bcs%ys(bcs%my)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find y interval.
c-----------------------------------------------------------------------
      DO
         IF(yy >= bcs%ys(bcs%iy) .OR. bcs%iy <= 0)EXIT
         bcs%iy=bcs%iy-1
      ENDDO
      DO
         IF(yy < bcs%ys(bcs%iy+1) .OR. bcs%iy >= bcs%my-1)EXIT
         bcs%iy=bcs%iy+1
      ENDDO
c-----------------------------------------------------------------------
c     find offsets and compute local coefficients.
c-----------------------------------------------------------------------
      dx=xx-bcs%xs(bcs%ix)
      dy=yy-bcs%ys(bcs%iy)
      c=bicube_getco(bcs)
c-----------------------------------------------------------------------
c     evaluate f.
c-----------------------------------------------------------------------
      bcs%f=0
      DO i=4,1,-1
         bcs%f=bcs%f*dx+((c(i,4,:)*dy+c(i,3,:))*dy+c(i,2,:))*dy+c(i,1,:)
      ENDDO
c-----------------------------------------------------------------------
c     restore x powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=x-bcs%x0(iside)
         DO iqty=1,bcs%nqty
            IF(bcs%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**bcs%xpower(iside,iqty)
            g=bcs%f(iqty)*xfac
            bcs%f(iqty)=g
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     restore y powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dy=y-bcs%y0(iside)
         DO iqty=1,bcs%nqty
            IF(bcs%ypower(iside,iqty) == 0)CYCLE
            yfac=ABS(dy)**bcs%ypower(iside,iqty)
            g=bcs%f(iqty)*yfac
            bcs%f(iqty)=g
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_eval
c-----------------------------------------------------------------------
c     subprogram 7. bicube_getco.
c     computes coefficient matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION bicube_getco(bcs) RESULT(cmat)

      TYPE(bicube_type), INTENT(IN) :: bcs
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: cmat

      REAL(r8) :: hxfac,hyfac
      REAL(r8), DIMENSION(3:4,4) :: gxmat,gymat
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: temp
c-----------------------------------------------------------------------
c     compute gxmat.
c-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(bcs%ix+1)-bcs%xs(bcs%ix))
      gxmat(3,1)=-3*hxfac**2
      gxmat(3,2)=-2*hxfac
      gxmat(3,3)=3*hxfac**2
      gxmat(3,4)=-hxfac
      gxmat(4,1)=2*hxfac**3
      gxmat(4,2)=hxfac**2
      gxmat(4,3)=-2*hxfac**3
      gxmat(4,4)=hxfac**2
c-----------------------------------------------------------------------
c     compute gymat.
c-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(bcs%iy+1)-bcs%ys(bcs%iy))
      gymat(3,1)=-3*hyfac**2
      gymat(3,2)=-2*hyfac
      gymat(3,3)=3*hyfac**2
      gymat(3,4)=-hyfac
      gymat(4,1)=2*hyfac**3
      gymat(4,2)=hyfac**2
      gymat(4,3)=-2*hyfac**3
      gymat(4,4)=hyfac**2
c-----------------------------------------------------------------------
c     compute smat.
c-----------------------------------------------------------------------
      cmat(1,1,:)=bcs%fs(bcs%ix,bcs%iy,:)
      cmat(1,2,:)=bcs%fsy(bcs%ix,bcs%iy,:)
      cmat(1,3,:)=bcs%fs(bcs%ix,bcs%iy+1,:)
      cmat(1,4,:)=bcs%fsy(bcs%ix,bcs%iy+1,:)
      cmat(2,1,:)=bcs%fsx(bcs%ix,bcs%iy,:)
      cmat(2,2,:)=bcs%fsxy(bcs%ix,bcs%iy,:)
      cmat(2,3,:)=bcs%fsx(bcs%ix,bcs%iy+1,:)
      cmat(2,4,:)=bcs%fsxy(bcs%ix,bcs%iy+1,:)
      cmat(3,1,:)=bcs%fs(bcs%ix+1,bcs%iy,:)
      cmat(3,2,:)=bcs%fsy(bcs%ix+1,bcs%iy,:)
      cmat(3,3,:)=bcs%fs(bcs%ix+1,bcs%iy+1,:)
      cmat(3,4,:)=bcs%fsy(bcs%ix+1,bcs%iy+1,:)
      cmat(4,1,:)=bcs%fsx(bcs%ix+1,bcs%iy,:)
      cmat(4,2,:)=bcs%fsxy(bcs%ix+1,bcs%iy,:)
      cmat(4,3,:)=bcs%fsx(bcs%ix+1,bcs%iy+1,:)
      cmat(4,4,:)=bcs%fsxy(bcs%ix+1,bcs%iy+1,:)
c-----------------------------------------------------------------------
c     multiply by gymat^T.
c-----------------------------------------------------------------------
      temp(:,1:2,:)=cmat(:,1:2,:)
      temp(:,3,:) = cmat(:,1,:)*gymat(3,1) + cmat(:,2,:)*gymat(3,2)
     $     + cmat(:,3,:)*gymat(3,3) + cmat(:,4,:)*gymat(3,4)
      temp(:,4,:) = cmat(:,1,:)*gymat(4,1) + cmat(:,2,:)*gymat(4,2)
     $     + cmat(:,3,:)*gymat(4,3) + cmat(:,4,:)*gymat(4,4)
c-----------------------------------------------------------------------
c     multiply by gxmat.
c-----------------------------------------------------------------------
      cmat(1:2,:,:) = temp(1:2,:,:)
      cmat(3,:,:) = gxmat(3,1)*temp(1,:,:) + gxmat(3,2)*temp(2,:,:)
     $     + gxmat(3,3)*temp(3,:,:) + gxmat(3,4)*temp(4,:,:)
      cmat(4,:,:) = gxmat(4,1)*temp(1,:,:) + gxmat(4,2)*temp(2,:,:)
     $     + gxmat(4,3)*temp(3,:,:) + gxmat(4,4)*temp(4,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION bicube_getco
c-----------------------------------------------------------------------
c     subprogram 8. bicube_fit.
c     fits functions to bicubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_fit(bcs,endmode1,endmode2)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      CHARACTER(*), INTENT(IN) :: endmode1,endmode2

      INTEGER :: iqty,iside,ix,iy
      REAL(r8), DIMENSION(0:bcs%mx) :: xfac
      REAL(r8), DIMENSION(0:bcs%my) :: yfac
      TYPE(spline_type) :: spl
c-----------------------------------------------------------------------
c     extract x powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,bcs%nqty
            IF(bcs%xpower(iside,iqty) /= 0)THEN
               xfac=1/ABS(bcs%xs-bcs%x0(iside))**bcs%xpower(iside,iqty)
               DO iy=0,bcs%my
                  bcs%fs(:,iy,iqty)=bcs%fs(:,iy,iqty)*xfac
               ENDDO
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     extract y powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,bcs%nqty
            IF(bcs%ypower(iside,iqty) /= 0)THEN
               yfac=1/ABS(bcs%ys-bcs%y0(iside))**bcs%ypower(iside,iqty)
               DO ix=0,bcs%mx
                  bcs%fs(ix,:,iqty)=bcs%fs(ix,:,iqty)*yfac
               ENDDO
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set periodicity.
c-----------------------------------------------------------------------
      bcs%periodic=(/endmode1 == "periodic",endmode2 == "periodic"/)
      IF(bcs%periodic(1))bcs%fs(bcs%mx,:,:)=bcs%fs(0,:,:)
      IF(bcs%periodic(2))bcs%fs(:,bcs%my,:)=bcs%fs(:,0,:)
c-----------------------------------------------------------------------
c     evaluate y derivatives.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl,bcs%my,bcs%mx+1)
      spl%xs=bcs%ys
      DO iqty=1,bcs%nqty
         spl%fs=TRANSPOSE(bcs%fs(:,:,iqty))
         CALL spline_fit(spl,endmode2)
         bcs%fsy(:,:,iqty)=TRANSPOSE(spl%fs1)
      ENDDO
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     evaluate x derivatives.
c-----------------------------------------------------------------------
      spl%mx=bcs%mx
      spl%nqty=bcs%my+1
      CALL spline_alloc(spl,bcs%mx,bcs%my+1)
      spl%xs=bcs%xs
      DO iqty=1,bcs%nqty
         spl%fs=bcs%fs(:,:,iqty)
         CALL spline_fit(spl,endmode1)
         bcs%fsx(:,:,iqty)=spl%fs1
      ENDDO
c-----------------------------------------------------------------------
c     evaluate mixed derivatives.
c-----------------------------------------------------------------------
      DO iqty=1,bcs%nqty
         spl%fs=bcs%fsy(:,:,iqty)
         CALL spline_fit(spl,endmode1)
         bcs%fsxy(:,:,iqty)=spl%fs1
      ENDDO
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_fit
      END MODULE bicube_mod
