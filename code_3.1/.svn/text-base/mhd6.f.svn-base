c-----------------------------------------------------------------------
c     file mhd6.f.
c     contains specifications for Hall MHD model with 6 variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      module organization.
c-----------------------------------------------------------------------
c     0. mhd6_mod.
c     1. mhd6_bubble.
c     2. mhd6_cross.
c     3. mhd6_eig.
c-----------------------------------------------------------------------
c     external subprograms.
c-----------------------------------------------------------------------
c     a. physics_input.
c     b. physics_init_parameters.
c     c. physics_init.
c     d. physics_boundary.
c     e. physics_edge_rhs.
c     f. physics_edge_drdu.
c     g. physics_edge_mass.
c     j. physics_rhs.
c     k. physics_drdu.
c     l. physics_mass.
c     m. physics_equil.
c     n. physics_grid.
c     o. physics_schur.
c     p. physics_dealloc.
c     q. physics_main.
c-----------------------------------------------------------------------
c     subprogram 0. mhd6_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE mhd6_mod
      USE local_mod
      IMPLICIT NONE

      LOGICAL :: hall_flag=.TRUE.
      CHARACTER(5) :: init_type
      INTEGER :: knx=1,kny=1,lx=1,ly=1,hall_type=1
      REAL(r8), PARAMETER :: gamma_mhd=5._r8/3._r8
      REAL(r8) :: Bx,By,Bz,kx,ky,k0=1,thetaB=0,phiB=0,beta0=0,
     $     ca2=1,cs2,omci2,omci,di=1,nu=0
      COMPLEX(r8), DIMENSION(6) :: u0fac
      COMPLEX(r8), PARAMETER :: czero=0,ifac=(0,1)

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. mhd6_bubble.
c     performs a bubble sort in decreasing order of value.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE mhd6_bubble(key,index,mmin,mmax)

      REAL(r8), DIMENSION(:), INTENT(IN) :: key
      INTEGER, DIMENSION(:), INTENT(INOUT) :: index
      INTEGER :: mmin,mmax

      LOGICAL :: switch
      INTEGER :: i,temp
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      switch= .TRUE.
      DO while(switch)
         switch= .FALSE.
         DO i=mmin,mmax-1
            IF(key(index(i)) < key(index(i+1)))THEN
               temp=index(i)
               index(i)=index(i+1)
               index(i+1)=temp
               switch= .TRUE.
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mhd6_bubble
c-----------------------------------------------------------------------
c     subprogram 2. mhd6_cross.
c     computes vector product.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION mhd6_cross(a,b) RESULT(c)

      REAL(r8), DIMENSION(3), INTENT(IN) :: a,b
      REAL(r8), DIMENSION(3) :: c
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      c=(/a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)/)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION mhd6_cross
c-----------------------------------------------------------------------
c     subprogram 3. mhd6_eig.
c     computes eigenvalues and eigenvectors numerically.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE mhd6_eig(kvec,e,kpar,kprp,evals,omega,v,b,p,
     $     init_type,jeig,errmax)

      REAL(r8), DIMENSION(3), INTENT(IN) :: kvec
      REAL(r8), DIMENSION(3,3), INTENT(IN) :: e
      REAL(r8), INTENT(OUT) :: kpar,kprp,errmax
      REAL(r8), DIMENSION(3), INTENT(OUT) :: evals,omega
      COMPLEX(r8), DIMENSION(3,3), INTENT(OUT) :: v,b
      COMPLEX(r8), DIMENSION(3), INTENT(OUT) :: p
      CHARACTER(*), INTENT(IN) :: init_type
      INTEGER, INTENT(OUT) :: jeig

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER, PARAMETER :: nqty=7
      CHARACTER(2), DIMENSION(nqty) ::
     $     var=(/" p","b1","b2","b3","v1","v2","v3"/)
      INTEGER, PARAMETER :: lwork=4*nqty
      INTEGER :: info,i,j
      INTEGER, DIMENSION(nqty) :: index
      REAL(r8), PARAMETER :: eps=1e-13
      REAL(r8), DIMENSION(2*nqty) :: rwork
      REAL(r8), DIMENSION(nqty) :: key,error
      REAL(r8), DIMENSION(3,3) :: rmat,ident
      COMPLEX(r8) :: vl
      COMPLEX(r8), DIMENSION(nqty) :: wvec
      COMPLEX(r8), DIMENSION(lwork) :: work
      COMPLEX(r8), DIMENSION(nqty,nqty) :: cmat,cmat0,vr
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/5x,"i",3x,"re wvec",4x,"im wvec"/)
 20   FORMAT(i6,1p,2e11.3)
 30   FORMAT(/"var",3(2x,"re mode ",i1,2x,"im mode ",i1)/)
 40   FORMAT(a3,1p,6e11.3)
c-----------------------------------------------------------------------
c     compute k.B and related quantities.
c-----------------------------------------------------------------------
      kpar=SUM(kvec*e(:,1))
      kprp=SUM(kvec*e(:,3))
c-----------------------------------------------------------------------
c     compute matrix, ideal terms.
c-----------------------------------------------------------------------
      ident=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      cmat=0
      cmat(1,5:7)=-cs2*kvec
      cmat(5:7,1)=-kvec
      cmat(2:4,5:7)=kpar*ident
     $     -RESHAPE((/((e(i,1)*kvec(j),i=1,3),j=1,3)/),(/3,3/))
      cmat(5:7,2:4)=kpar*ident
     $     -RESHAPE((/((kvec(i)*e(j,1),i=1,3),j=1,3)/),(/3,3/))
c-----------------------------------------------------------------------
c     compute matrix, Hall terms.
c-----------------------------------------------------------------------
      IF(hall_flag)THEN
         rmat=cmat(5:7,2:4)/omci
         DO j=1,3
            rmat(:,j)=mhd6_cross(kvec,rmat(:,j))
         ENDDO
         cmat(2:4,2:4)=-ifac*rmat
      ENDIF
      cmat0=cmat
c-----------------------------------------------------------------------
c     compute eigenvalues and eigenvectors.
c-----------------------------------------------------------------------
      CALL zgeev('N','V',nqty,cmat,nqty,wvec,vl,1,vr,nqty,
     $     work,lwork,rwork,info)
c-----------------------------------------------------------------------
c     sort eigenvalues.
c-----------------------------------------------------------------------
      index=(/(i,i=1,nqty)/)
      key=wvec
      CALL mhd6_bubble(key,index,1,nqty)
c-----------------------------------------------------------------------
c     eliminate roundoff terms.
c-----------------------------------------------------------------------
      WHERE(ABS(wvec) < eps)
         wvec=0
      ENDWHERE
      WHERE(ABS(vr) < eps)
         vr=0
      ENDWHERE
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         OPEN(UNIT=dat_unit,FILE="eig.out",STATUS="UNKNOWN")
         WRITE(dat_unit,10)
         WRITE(dat_unit,20)(i,wvec(index(i)),i=1,nqty)
         WRITE(dat_unit,30)(i,i,i=1,3)
         WRITE(dat_unit,40)(var(i),(vr(i,index(j)),j=1,3),i=1,nqty)
         CLOSE(UNIT=dat_unit)
      ENDIF
c-----------------------------------------------------------------------
c     compute output.
c-----------------------------------------------------------------------
      DO i=1,3
         omega(i)=wvec(index(i))
         p(i)=vr(1,index(i))
         b(:,i)=vr(2:4,index(i))
         v(:,i)=vr(5:7,index(i))
      ENDDO
      evals=omega**2
c-----------------------------------------------------------------------
c     select eigenpair.
c-----------------------------------------------------------------------
      SELECT CASE(init_type)
      CASE("fast")
         jeig=1
      CASE("shear")
         jeig=2
      CASE("slow")
         jeig=3
      CASE DEFAULT
         CALL program_stop
     $        ("Cannot recognize init_type = "//TRIM(init_type))
      END SELECT
c-----------------------------------------------------------------------
c     compute error.
c-----------------------------------------------------------------------
      DO i=1,nqty
         error(i)=MAXVAL(ABS(MATMUL(cmat0,vr(:,i))-wvec(i)*vr(:,i)))
     $        /MAXVAL(ABS(vr(:,i)))
      ENDDO
      errmax=MAXVAL(error)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mhd6_eig
      END MODULE mhd6_mod
c-----------------------------------------------------------------------
c     subprogram a. physics_input.
c     sets up input constants.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_input(nx,ny,np,nq,nbx,xperiodic,yperiodic,nqty,
     $     nqty_schur,dt,dtmax,tmax,nstep,gr_curve_in,du_diagnose,
     $     physics_type,exit_flag)
      USE mhd6_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: physics_type
      LOGICAL, INTENT(INOUT) :: xperiodic,yperiodic,du_diagnose
      INTEGER, INTENT(OUT) :: nqty,nqty_schur
      INTEGER, INTENT(INOUT) :: nx,ny,np,nq,nbx,nstep,exit_flag
      REAL(r8), INTENT(INOUT) :: dt,dtmax,tmax,gr_curve_in

      INTEGER :: myios,i,jeig
      REAL(r8) :: kpar,kprp,k,kb_angle,tfac,errmax,ratio
      REAL(r8), DIMENSION(3) :: kvec,evals,omega
      REAL(r8), DIMENSION(3,3) :: e
      COMPLEX(r8), DIMENSION(3) :: p
      COMPLEX(r8), DIMENSION(3,3) :: v,b
c-----------------------------------------------------------------------
c     namelist.
c-----------------------------------------------------------------------
      NAMELIST/mhd_list/nx,ny,np,nq,nbx,xperiodic,yperiodic,dt,dtmax,
     $     tmax,nstep,gr_curve_in,du_diagnose,init_type,thetaB,phiB,
     $     beta0,knx,kny,lx,ly,k0,di,hall_flag,hall_type,nu
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,"knx",1x,"kny",2x,"hall",1x,"jeig",2x,"init",
     $     4x,"thetaB",6x,"phiB",6x,"beta0",8x,"di"
     $     //2i4,l5,i5,a7,1p,4e11.3)
 20   FORMAT(/5x,"kpar",7x,"kprp",5x,"kb_angle",4x,
     $     "omega1",5x,"omega2",5x,"omega3",5x,"ratio"//1p,7e11.3)
 30   FORMAT (/2x,"i",6x,"e1",9x,"e2",9x,"e3"/)
 40   FORMAT(i3,1p,3e11.3)
 50   FORMAT(a6,1p,3e11.3)
 60   FORMAT(/"var",3(2x,"re mode ",i1,2x,"im mode ",i1)/)
 70   FORMAT(a3,1p,6e11.3)
 80   FORMAT(/4x,"u0fac1",5x,"u0fac2",5x,"u0fac3",5x,"u0fac4",
     $     5x,"u0fac5",5x,"u0fac6"/)
 90   FORMAT(1p,6e11.3)
c-----------------------------------------------------------------------
c     read and set/reset input variables.
c-----------------------------------------------------------------------
      READ(in_unit,NML=mhd_list,IOSTAT=myios)
      IF(myios /= 0)exit_flag=99
      nqty=6
      nqty_schur=3
      physics_type="mhd6"
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      ca2=1
      cs2=gamma_mhd*beta0*half
      omci2=ca2/di**2
      omci=SQRT(omci2)
c-----------------------------------------------------------------------
c     compute equilibrium magnetic field.
c-----------------------------------------------------------------------
      IF(phib == 90)THEN
         Bx=0
      ELSE
         Bx=SIN(thetaB*dtor)*COS(phiB*dtor)
      ENDIF
      By=SIN(thetaB*dtor)*SIN(phiB*dtor)
      IF(thetaB == 90)THEN
         Bz=0
      ELSE
         Bz=COS(thetaB*dtor)
      ENDIF
c-----------------------------------------------------------------------
c     compute kx, ky, and kvec.
c-----------------------------------------------------------------------
      k=twopi*SQRT(REAL(knx**2+kny**2,8))
      kx=twopi*knx*k0/k
      ky=twopi*kny*k0/k
      kvec=(/kx,ky,zero/)
c-----------------------------------------------------------------------
c     compute orthonormal basis vectors.
c-----------------------------------------------------------------------
      e(:,1)=(/Bx,By,Bz/)
      e(:,2)=mhd6_cross(kvec,e(:,1))
      e(:,3)=mhd6_cross(e(:,1),e(:,2))
      DO i=1,3
         e(:,i)=e(:,i)/SQRT(SUM(e(:,i)**2))
      ENDDO
c-----------------------------------------------------------------------
c     compute eigenvalues and eigenvectors.
c-----------------------------------------------------------------------
      CALL mhd6_eig(kvec,e,kpar,kprp,evals,omega,v,b,p,init_type,
     $     jeig,errmax)
c-----------------------------------------------------------------------
c     compute u0fac.
c-----------------------------------------------------------------------
      u0fac(1)=p(jeig)
      u0fac(2)=b(2,jeig)/(ifac*kx)
      u0fac(3)=b(3,jeig)
      u0fac(4:6)=v(:,jeig)
c-----------------------------------------------------------------------
c     open file and write scalars.
c-----------------------------------------------------------------------
      kb_angle=ACOS(kpar/SQRT(SUM(kvec**2)))*rtod
      ratio=omega(1)/omega(3)
      OPEN(UNIT=debug_unit,FILE="mhd6.out",STATUS="UNKNOWN")
      WRITE(debug_unit,10)knx,kny,hall_flag,jeig,TRIM(init_type),
     $     thetaB,phiB,beta0,di
      WRITE(debug_unit,20)kpar,kprp,kb_angle,omega,ratio
      WRITE(debug_unit,30)
      WRITE(debug_unit,40)(i,e(i,:),i=1,3)
c-----------------------------------------------------------------------
c     diagnose eigenvectors.
c-----------------------------------------------------------------------
      WRITE(debug_unit,60)(i,i,i=1,3)
      WRITE(debug_unit,70)"p",p
      WRITE(debug_unit,70)"b1",b(1,:)
      WRITE(debug_unit,70)"b2",b(2,:)
      WRITE(debug_unit,70)"b3",b(3,:)
      WRITE(debug_unit,70)"v1",v(1,:)
      WRITE(debug_unit,70)"v2",v(2,:)
      WRITE(debug_unit,70)"v3",v(3,:)
c-----------------------------------------------------------------------
c     diagnose initial conditions.
c-----------------------------------------------------------------------
      WRITE(debug_unit,80)
      WRITE(debug_unit,90)REAL(u0fac)
      WRITE(debug_unit,90)IMAG(u0fac)
      CLOSE(UNIT=debug_unit)
c-----------------------------------------------------------------------
c     rescale dt-related input.
c-----------------------------------------------------------------------
      tfac=twopi/omega(jeig)
      dt=dt*tfac
      dtmax=dtmax*tfac
      tmax=tmax*tfac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_input
c-----------------------------------------------------------------------
c     subprogram b. physics_init_parameters.
c     special initial conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_init_parameters(static,ground,adapt_qty)
      USE mhd6_mod
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), INTENT(INOUT) :: static,ground,adapt_qty
c-----------------------------------------------------------------------
c     broadcast logical and character input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(hall_flag,1,MPI_LOGICAL,0,comm,ierr)
      CALL MPI_Bcast(init_type,16,MPI_CHARACTER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast integer input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(knx,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(kny,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(lx,1,MPI_INTEGER,0,comm,ierr)
      CALL MPI_Bcast(ly,1,MPI_INTEGER,0,comm,ierr)
c-----------------------------------------------------------------------
c     broadcast real input.
c-----------------------------------------------------------------------
      CALL MPI_Bcast(Bx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(By,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(Bz,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(kx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ky,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(thetaB,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(phiB,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(beta0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(ca2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(cs2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(omci2,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(omci,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(di,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
      CALL MPI_Bcast(u0fac,2*SIZE(u0fac),MPI_DOUBLE_PRECISION,0,comm,
     $     ierr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_init_parameters
c-----------------------------------------------------------------------
c     subprogram c. physics_init.
c     initializes solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_init(xpi,ypi,ui)
      USE mhd6_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: xpi,ypi
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: ui

      INTEGER :: iqty
      COMPLEX(r8), DIMENSION(0:SIZE(xpi,1)-1,0:SIZE(xpi,2)-1) :: expfac
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      expfac=EXP(ifac*(kx*xpi+ky*ypi))
      DO iqty=1,6
         ui(iqty,:,:)=u0fac(iqty)*expfac
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_init
c-----------------------------------------------------------------------
c     subprogram d. physics_boundary
c     initializes boundary conditions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_boundary(left,right,top,bottom,nqty,edge_order)
      USE mhd6_mod
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nqty
      INTEGER, DIMENSION(4), INTENT(INOUT) :: edge_order
      TYPE(edge_type) :: left,right,top,bottom

      INTEGER :: k,j
      REAL(r8), DIMENSION(nqty,nqty) :: ident
c-----------------------------------------------------------------------
c     initialize boundary conditions.
c-----------------------------------------------------------------------
      top%bc_type="robin"
      top%static=.FALSE.
      bottom%bc_type="robin"
      bottom%static=.FALSE.
      left%bc_type="robin"
      left%static=.FALSE.
      right%bc_type="robin"
      right%static=.FALSE.

      SELECT CASE(init_type)
      CASE("Harris")
         top%static(2:4)=.TRUE.
         bottom%static(2:4)=.TRUE.
      CASE("Harris_open")
         left%static=.TRUE.
         right%static=.TRUE.
         top%static=.TRUE.
         bottom%static=.TRUE.
      CASE("Harris_local")
         left%static=.TRUE.
         right%static=.TRUE.
         top%static(2:4)=.TRUE.
         bottom%static=.TRUE.
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_boundary
c-----------------------------------------------------------------------
c     subprogram e. physics_edge_rhs.
c     computes rhs for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_rhs(lrtb,t,x,y,nhat,u,ux,uy,uxx,uyy,uxy,c)
      USE mhd6_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: c
c-----------------------------------------------------------------------
c     set c(nqty,:,:) to the RHS for "robin" b.c. 
c     or to normal outward flux for "normflux" b.c.
c-----------------------------------------------------------------------
      c=0

      SELECT CASE(lrtb)
      CASE("left")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("right")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("top")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("bottom")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_edge_rhs
c-----------------------------------------------------------------------
c     subprogram f. physics_edge_drdu.
c     computes drdu for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_drdu(lrtb,t,x,y,nhat,
     $     u,ux,uy,uxx,uyy,uxy,c_u,c_ux,c_uy,c_uxx,c_uyy,c_uxy)
      USE mhd6_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat,u,ux,uy,uxx,uyy,uxy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: c_u,c_ux,c_uy,
     $     c_uxx,c_uyy,c_uxy
c-----------------------------------------------------------------------
c     set c_xxx(iqty,jqty), to the derivatives of c(iqty) 
c     given in physics_edge_rhs with respect to 
c     u(jqty),ux(jqty),uy(jqty),uxx(jqty),uyy(jqty),uxy(jqty).
c-----------------------------------------------------------------------
      c_u=0
      c_ux=0
      c_uy=0
      c_uxx=0
      c_uyy=0
      c_uxy=0

      SELECT CASE(lrtb)
      CASE("left")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("right")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("top")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("bottom")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_edge_drdu
c-----------------------------------------------------------------------
c     subprogram g. physics_edge_mass.
c     computes mass matrices for non-linear b.c.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_edge_mass(lrtb,x,y,nhat,mass,mass_x,mass_y)
      USE mhd6_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: lrtb
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: nhat
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     set mass,mass_x,mass_y(iqty,jqty,:,:), ... to coupling 
c     mass matrices for du/dt, d(du/dx)/dt, and d(du/dy)/dt terms
c     of the general "robin" b.c. equations
c-----------------------------------------------------------------------
      mass=0
      mass_x=0
      mass_y=0

      SELECT CASE(lrtb)
      CASE("left")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("right")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("top")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      CASE("bottom")
         SELECT CASE(init_type)
         CASE(" ")

         END SELECT
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_edge_mass
c-----------------------------------------------------------------------
c     subprogram j. physics_rhs.
c     computes fluxes and sources.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_rhs(t,x,y,u,ux,uy,fx,fy,s,first)
      USE mhd6_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: fx,fy,s
      LOGICAL, INTENT(INOUT) :: first
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx=0
      fy=0
      s=0
c-----------------------------------------------------------------------
c     pressure p.
c-----------------------------------------------------------------------
      fx(1,:,:)=cs2*u(4,:,:)
      fy(1,:,:)=cs2*u(5,:,:)
c-----------------------------------------------------------------------
c     magnetic flux psi.
c-----------------------------------------------------------------------
      s(2,:,:)=-u(4,:,:)*By+u(5,:,:)*Bx
c-----------------------------------------------------------------------
c     normal magnetic field bz.
c-----------------------------------------------------------------------
      fx(3,:,:)=Bz*u(4,:,:)-Bx*u(6,:,:)
      fy(3,:,:)=Bz*u(5,:,:)-By*u(6,:,:)
c-----------------------------------------------------------------------
c     vx.
c-----------------------------------------------------------------------
      fx(4,:,:)=u(1,:,:)+ux(2,:,:)*By+u(3,:,:)*Bz
      fy(4,:,:)=uy(2,:,:)*By
c-----------------------------------------------------------------------
c     vy.
c-----------------------------------------------------------------------
      fx(5,:,:)=-ux(2,:,:)*Bx
      fy(5,:,:)=u(1,:,:)-uy(2,:,:)*Bx+u(3,:,:)*Bz
c-----------------------------------------------------------------------
c     vz.
c-----------------------------------------------------------------------
      fx(6,:,:)=-u(3,:,:)*Bx
      fy(6,:,:)=-u(3,:,:)*By
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_rhs
c-----------------------------------------------------------------------
c     subprogram k. physics_drdu.
c     computes contributions to the jacobian.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_drdu(t,x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
      USE mhd6_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx_u=0
      fx_ux=0
      fx_uy=0
      fy_u=0
      fy_ux=0
      fy_uy=0
      s_u=0
      s_ux=0
      s_uy=0
c-----------------------------------------------------------------------
c     pressure p.
c-----------------------------------------------------------------------
      fx_u(1,4,:,:)=cs2
      fy_u(1,5,:,:)=cs2
c-----------------------------------------------------------------------
c     magnetic flux psi.
c-----------------------------------------------------------------------
      s_u(2,4,:,:)=-By
      s_u(2,5,:,:)=Bx
c-----------------------------------------------------------------------
c     normal magnetic field bz.
c-----------------------------------------------------------------------
      fx_u(3,4,:,:)=Bz
      fx_u(3,6,:,:)=-Bx
      fy_u(3,5,:,:)=Bz
      fy_u(3,6,:,:)=-By
c-----------------------------------------------------------------------
c     vx.
c-----------------------------------------------------------------------
      fx_u(4,1,:,:)=one
      fx_u(4,3,:,:)=Bz
      fx_ux(4,2,:,:)=By
      fy_uy(4,2,:,:)=By
c-----------------------------------------------------------------------
c     vy.
c-----------------------------------------------------------------------
      fx_ux(5,2,:,:)=-Bx
      fy_u(5,1,:,:)=one
      fy_u(5,3,:,:)=Bz
      fy_uy(5,2,:,:)=-Bx
c-----------------------------------------------------------------------
c     vz.
c-----------------------------------------------------------------------
      fx_u(6,3,:,:)=-Bx
      fy_u(6,3,:,:)=-By
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_drdu
c-----------------------------------------------------------------------
c     subprogram l. physics_mass.
c     computes contributions to the mass matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_mass(x,y,mass,mass_x,mass_y)
      USE mhd6_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mass,mass_x,mass_y
c-----------------------------------------------------------------------
c     modify mass matrices.
c-----------------------------------------------------------------------
      IF(hall_flag)THEN
         mass(2,6,:,:)=-1/omci
         mass_x(3,5,:,:)=1/omci
         mass_y(3,4,:,:)=-1/omci
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_mass
c-----------------------------------------------------------------------
c     subprogram m. physics_equil.
c     computes equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_equil(x,y,u,ux,uy,derivs)
      USE mhd6_mod
      IMPLICIT NONE
      
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: u,ux,uy
      LOGICAL, INTENT(IN) :: derivs
c-----------------------------------------------------------------------
c     compute equilibrium.
c-----------------------------------------------------------------------
      u=0
      ux=0
      uy=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_equil
c-----------------------------------------------------------------------
c     subprogram n. physics_grid.
c     computes initial physical 2D grid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_grid(x,y,ksi,eta)
      USE mhd6_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(0:,0:), INTENT(OUT) :: ksi,eta
c-----------------------------------------------------------------------
c     y positions.
c-----------------------------------------------------------------------
      ksi=twopi*x*lx/kx
      eta=twopi*y*ly/ky
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_grid
c-----------------------------------------------------------------------
c     subprogram o. physics_schur.
c     computes schur complement.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_schur(hfac,x,y,u,ux,uy,
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy)
      USE mhd6_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: hfac
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: u,ux,uy
      REAL(r8), DIMENSION(:,:,:,:), INTENT(OUT) ::
     $     fx_u,fx_ux,fx_uy,fy_u,fy_ux,fy_uy,s_u,s_ux,s_uy
 
      REAL(r8) :: fac
c-----------------------------------------------------------------------
c     zero output.
c-----------------------------------------------------------------------
      fx_u=0
      fx_ux=0
      fx_uy=0
      fy_u=0
      fy_ux=0
      fy_uy=0
      s_u=0
      s_ux=0
      s_uy=0
c-----------------------------------------------------------------------
c     derivatives of txx.
c-----------------------------------------------------------------------
      fx_ux(1,1,:,:)=Bx*Bx-(1+cs2)
      fx_ux(1,2,:,:)=Bx*By
      fx_uy(1,2,:,:)=Bx*Bx+By*By-(1+cs2)
      fx_ux(1,3,:,:)=Bx*Bz
      fx_uy(1,3,:,:)=By*Bz
c-----------------------------------------------------------------------
c     derivatives of tyy.
c-----------------------------------------------------------------------
      fy_ux(2,1,:,:)=Bx*Bx+By*By-(1+cs2)
      fy_uy(2,1,:,:)=Bx*By
      fy_uy(2,2,:,:)=By*By-(1+cs2)
      fy_ux(2,3,:,:)=Bx*Bz
      fy_uy(2,3,:,:)=By*Bz
c-----------------------------------------------------------------------
c     derivatives of txy.
c-----------------------------------------------------------------------
      fx_ux(2,1,:,:)=Bx*By
      fx_ux(2,2,:,:)=-Bx*Bx
      fx_ux(2,3,:,:)=0
      fx_uy(2,3,:,:)=0
c-----------------------------------------------------------------------
c     derivatives of tyx.
c-----------------------------------------------------------------------
      fy_uy(1,1,:,:)=-By*By
      fy_uy(1,2,:,:)=Bx*By
      fy_ux(1,3,:,:)=0
      fy_uy(1,3,:,:)=0
c-----------------------------------------------------------------------
c     derivatives of txz.
c-----------------------------------------------------------------------
      fx_ux(3,1,:,:)=Bx*Bz
      fx_ux(3,2,:,:)=0
      fx_uy(3,2,:,:)=Bx*Bz
      fx_ux(3,3,:,:)=-Bx*Bx
      fx_uy(3,3,:,:)=-Bx*By
c-----------------------------------------------------------------------
c     derivatives of tyz.
c-----------------------------------------------------------------------
      fy_ux(3,1,:,:)=By*Bz
      fy_uy(3,1,:,:)=0
      fy_uy(3,2,:,:)=By*Bz
      fy_ux(3,3,:,:)=-Bx*By
      fy_uy(3,3,:,:)=-By*By
c-----------------------------------------------------------------------
c     multiply by hfac**2.
c-----------------------------------------------------------------------
      fac=hfac*hfac
      fx_ux=fx_ux*fac
      fx_uy=fx_uy*fac
      fy_ux=fy_ux*fac
      fy_uy=fy_uy*fac
c-----------------------------------------------------------------------
c     Hall terms.
c-----------------------------------------------------------------------
      IF(hall_flag)THEN
         fac=-hfac/omci

         fx_uy(1,1,:,:)=fx_uy(1,1,:,:)-Bz*fac
         fx_ux(1,2,:,:)=fx_ux(1,2,:,:)+Bz*fac
         fx_ux(1,3,:,:)=fx_ux(1,3,:,:)-By*fac
         
         fy_uy(2,1,:,:)=fy_uy(2,1,:,:)-Bz*fac
         fy_ux(2,2,:,:)=fy_ux(2,2,:,:)+Bz*fac
         fy_uy(2,3,:,:)=fy_uy(2,3,:,:)+Bx*fac
         
         fx_ux(2,3,:,:)=fx_ux(2,3,:,:)+Bx*fac
         fy_uy(1,3,:,:)=fy_uy(1,3,:,:)-By*fac
         
         fx_uy(3,1,:,:)=fx_uy(3,1,:,:)+Bx*fac
         fx_ux(3,2,:,:)=fx_ux(3,2,:,:)-Bx*fac
         
         fy_uy(3,1,:,:)=fy_uy(3,1,:,:)+By*fac
         fy_ux(3,2,:,:)=fy_ux(3,2,:,:)-By*fac

      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_schur
c-----------------------------------------------------------------------
c     subprogram p. physics_dealloc.
c     deallocates private module objects allocated above.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE physics_dealloc
      USE mhd6_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physics_dealloc
c-----------------------------------------------------------------------
c     subprogram q. physics_main.
c     trivial main program.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM physics_main
      USE driver_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     do it.
c-----------------------------------------------------------------------
      CALL driver_run
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL program_stop("Normal termination.")
      END PROGRAM physics_main
