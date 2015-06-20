c-----------------------------------------------------------------------
c     file bessel.f.
c     evaluates functions Jm and Ym for integer m>0 with accuracy eps.
c     evaluates eigenvalue lambda_n^(m).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module bessel_evaluation.
c     1. besselJ
c     2. besselJp
c     3. besselY
c     4. besselYp
c     5. fun_for_lam
c     6. lambda
c     7. lambda0
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module declarations.
c-----------------------------------------------------------------------
      MODULE bessel_evaluation
      USE local_mod
      IMPLICIT NONE
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. besselJ.
c     evaluates the series expansion for the J Bessel function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION besselJ(m, z, eps) RESULT(bessJ)

      INTEGER, INTENT(IN) :: m
      REAL(r8), INTENT(IN) :: z, eps
      REAL(r8) :: bessJ, nextterm
      INTEGER :: j, fac
c-----------------------------------------------------------------------
c     evaluate m!
c-----------------------------------------------------------------------
      fac=1
      DO j=2,m
         fac=fac*j
      ENDDO
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      bessJ=1
      nextterm=1
      j=1   
c-----------------------------------------------------------------------
c     evaluate Jm(z)
c-----------------------------------------------------------------------
      DO
        nextterm=-nextterm*z*z/(4*j*(j+m)) 
        bessJ=bessJ+nextterm
        IF (ABS(nextterm/bessJ)<eps) EXIT
        j=j+1
      ENDDO
      bessJ=bessJ*(0.5*z)**m/fac
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION besselJ
c-----------------------------------------------------------------------
c     subprogram 2. besselJp.
c     evaluates derivative of the J Bessel function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION besselJp(m, z, eps) RESULT(bessJp)

      INTEGER, INTENT(IN) :: m
      REAL(r8), INTENT(IN) :: z, eps
      REAL(r8) :: bessJp
c-----------------------------------------------------------------------
c     evaluate dJm(z)/dz
c-----------------------------------------------------------------------
      bessJp=REAL(m,8)/z*besselJ(m,z,eps)-besselJ(m+1,z,eps)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION besselJp
c-----------------------------------------------------------------------
c     subprogram 3. besselY.
c     evaluates the series expansion for the Y Bessel function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------    
      FUNCTION besselY(m, z, eps) RESULT(bessY)

      INTEGER, INTENT(IN) :: m
      REAL(r8), INTENT(IN) :: z, eps
      REAL(r8) :: bessY, nextterm, w, intgr
      INTEGER :: j, fac
      REAL(r8), PARAMETER :: gamma12=0.5772156649015328606_8
      REAL(r8), PARAMETER :: two=2.0_8,four=4.0_8
c-----------------------------------------------------------------------
c     evaluate m! and sum (1+1/2+1/3+...+1/m).
c-----------------------------------------------------------------------
      fac=1
      w=zero
      DO j=1,m
         fac=fac*j
         w=w+one/j
      ENDDO 
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      bessY=w
      nextterm=w
      intgr=one
      j=1   
c-----------------------------------------------------------------------
c     evaluate Ym(z).
c-----------------------------------------------------------------------
      DO 
        intgr=-intgr*(z*z/four)/(j*(j+m))
        w=w+one/j+one/(j+m)
        nextterm=intgr*w
        bessY=bessY+nextterm
        IF (ABS(nextterm/bessY)<eps) EXIT
        j=j+1
      ENDDO
      bessY =(-one/pi/fac)*bessY*(z/two)**m
     $     + (two/pi)*(LOG(z/two)+gamma12)*besselJ(m,z,eps)
      IF (m>0) THEN
         nextterm=fac/m
         intgr=fac/m
         DO j=1,(m-1)
            nextterm=nextterm*(z*z/four)/(j*(m-j))
            intgr=intgr+nextterm
         ENDDO
         intgr=-intgr/pi*(two/z)**m
         bessY=bessY+intgr
      ENDIF   
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION besselY
c-----------------------------------------------------------------------
c     subprogram 4. besselYp.
c     evaluates derivative of the Y Bessel function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION besselYp(m, z, eps) RESULT(bessYp)

      INTEGER, INTENT(IN) :: m
      REAL(r8), INTENT(IN) :: z, eps
      REAL(r8) :: bessYp
c-----------------------------------------------------------------------
c     evaluate dYm(z)/dz
c-----------------------------------------------------------------------
      bessYp=REAL(m,8)/z*besselY(m,z,eps)-besselY(m+1,z,eps)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION besselYp
c-----------------------------------------------------------------------
c     subprogram 5. fun_for_lam.
c     evaluates composite function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION fun_for_lam(m, rmin, rmax, eps, lamf) RESULT(fun)

      INTEGER, INTENT(IN) :: m
      REAL(r8) :: rmin, rmax, eps, fun, lamf
c-----------------------------------------------------------------------
c     evaluate fun_for_lam.
c-----------------------------------------------------------------------
      fun=besselJ(m,SQRT(lamf)*rmin,eps)*besselY(m,SQRT(lamf)*rmax,eps)-
     $     besselJ(m,SQRT(lamf)*rmax,eps)*besselY(m,SQRT(lamf)*rmin,eps)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION fun_for_lam
c-----------------------------------------------------------------------
c     subprogram 6. lambda.
c     evaluates nth root of function fun_for_lam.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------    
      FUNCTION lambda(m,n,rmin,rmax,eps) RESULT(lam)

      INTEGER, INTENT(IN) :: m, n
      REAL(r8), INTENT(IN) :: rmin, rmax, eps
      REAL(r8) :: lam, lam_min, lam_max, lam_cur
      INTEGER :: n0
      REAL(r8), PARAMETER :: two=2.0_r8,hund=100.0_r8
      REAL(r8) :: delta
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      delta=REAL(one/hund)
      lam=zero
      n0=0
c-----------------------------------------------------------------------
c     braket the nth root.
c-----------------------------------------------------------------------
      DO
         lam=lam+delta
         IF (fun_for_lam(m,rmin,rmax,eps,lam)*
     $   fun_for_lam(m,rmin,rmax,eps,lam+delta)<0) n0=n0+1
         IF (n0 == n) EXIT
      ENDDO      
c-----------------------------------------------------------------------
c     evaluate the nth root with accuarcy eps using bisection.
c-----------------------------------------------------------------------
      lam_min=lam
      lam_max=lam+delta
      DO
         IF ((lam_max-lam_min)<eps) EXIT
         lam_cur=(lam_max+lam_min)/two
         IF (fun_for_lam(m,rmin,rmax,eps,lam_min)*
     $   fun_for_lam(m,rmin,rmax,eps,lam_cur)>0) THEN
            lam_min=lam_cur
         ELSE
            lam_max=lam_cur
         ENDIF   
      ENDDO
      lam=lam_cur
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION lambda
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 7. lambda0.
c     evaluates nth (not counting z=0) root of BesselJ_m.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------    
      FUNCTION lambda0(m,n,eps) RESULT(root)

      INTEGER, INTENT(IN) :: m, n
      REAL(r8), INTENT(IN) :: eps
      REAL(r8) :: root, root_min, root_max, root_cur
      INTEGER :: n0
      REAL(r8), PARAMETER :: two=2.0_8,hund=100.0_8
      REAL(r8) :: delta
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      delta=REAL(one/hund)
      root=zero
      n0=0
c-----------------------------------------------------------------------
c     braket the nth root.
c-----------------------------------------------------------------------
      DO
         root=root+delta
         IF (besselJ(m,root,eps)*
     $   besselJ(m,root+delta,eps)<0) n0=n0+1
         IF (n0 == n) EXIT
      ENDDO
c-----------------------------------------------------------------------
c     evaluate the nth root with accuarcy eps using bisection.
c-----------------------------------------------------------------------
      root_min=root
      root_max=root+delta
      DO
         IF ((root_max-root_min)<eps) EXIT
         root_cur=(root_max+root_min)/two
         IF (besselJ(m,root_min,eps)*
     $   besselJ(m,root_cur,eps)>0) THEN
            root_min=root_cur
         ELSE
            root_max=root_cur
         ENDIF   
      ENDDO
      root=root_cur
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END FUNCTION lambda0
c-----------------------------------------------------------------------
      END MODULE bessel_evaluation


