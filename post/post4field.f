c-----------------------------------------------------------------------
c     file post4field.f.
c     post-processes output from sel code for "four_field" equations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. post4field_mod.
c     1. post4field_read.
c     2. post4field_UxyT.
c     3. post4field_VecSqInt.
c-----------------------------------------------------------------------
c     subprogram 0. post4field_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE post4field_mod
      USE plotter_mod
      IMPLICIT NONE

      REAL(r8), PRIVATE :: x0,y0,skin,hall,eta,mu
      CHARACTER(20), PRIVATE :: UxyT_name,VecSqInt_name

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. post4field_read.
c     read necessary post-processing parameters from  post.in 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post4field_read

      NAMELIST/fourfield_input/skin,hall,eta,mu,x0,y0,UxyT_name,
     $     VecSqInt_name
c-----------------------------------------------------------------------
c     read control file.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE="post.in",STATUS="OLD")
      READ(in_unit,NML=fourfield_input)
      CLOSE(UNIT=in_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post4field_read
c-----------------------------------------------------------------------
c     subprogram 2. post4field_UxyT.
c     generate u(x,y) vs. time plot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post4field_UxyT(t,xyw,xyw_kt,uw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: xyw_kt

      LOGICAL :: err
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: utemp
      REAL(r8), DIMENSION(SIZE(uw,1)) :: value
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      IF(UxyT_name=="Current")THEN
         ALLOCATE(utemp(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         utemp(1,:,:)=uw(5,:,:)
         CALL  plotter_UxyT(t,x0,y0,xyw,xyw_kt,utemp,.TRUE.,value,err)
      ENDIF
      DEALLOCATE(utemp)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post4field_UxyT
c-----------------------------------------------------------------------
c     subprogram 3. post4field_VecSqInt.
c     generate a time plot for an integral of a square of a vector
c     quantity over the computational domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE post4field_VecSqInt(t,xyw,jac,uw,uxyw)

      REAL(r8), INTENT(IN) :: t
      REAL(r8), DIMENSION(0:,0:), INTENT(IN) :: jac
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: xyw,uw
      REAL(r8), DIMENSION(0:,0:,:,:), INTENT(IN) :: uxyw

      REAL(r8), DIMENSION(:), ALLOCATABLE :: value
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: uw1,uw2,uw3
c-----------------------------------------------------------------------
c     call appropriate subroutines for necessary physical quantities
c-----------------------------------------------------------------------
      IF(VecSqInt_name=="DissEnergy")THEN
         ALLOCATE(value(1))
         ALLOCATE(uw1(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         ALLOCATE(uw2(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         ALLOCATE(uw3(1,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         uw1(1,:,:)=SQRT(eta)*(uxyw(:,:,1,3)+skin*uxyw(:,:,1,6))
         uw2(1,:,:)=SQRT(eta)*(uxyw(:,:,2,3)+skin*uxyw(:,:,2,6))
         uw3(1,:,:)=SQRT(eta)*uw(5,:,:)
         CALL plotter_VecSqInt(t,xyw,uw1,uw2,uw3,.TRUE.,value)
      ELSEIF(VecSqInt_name=="Energy")THEN
         ALLOCATE(value(2))
         ALLOCATE(uw1(2,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         ALLOCATE(uw2(2,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         ALLOCATE(uw3(2,0:SIZE(uw,2)-1,0:SIZE(uw,3)-1))
         uw1(1,:,:)=uxyw(:,:,1,1)+skin*uxyw(:,:,1,5)
         uw2(1,:,:)=uxyw(:,:,2,1)+skin*uxyw(:,:,2,5)
         uw3(1,:,:)=uw(3,:,:)+skin*uw(6,:,:)
         uw1(2,:,:)=uxyw(:,:,1,7)
         uw2(2,:,:)=uxyw(:,:,2,7)
         uw3(2,:,:)=uw(4,:,:)
         CALL plotter_VecSqInt(t,xyw,uw1,uw2,uw3,.TRUE.,value)
      ENDIF
      DEALLOCATE(value,uw1,uw2,uw3)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE post4field_VecSqInt
      END MODULE post4field_mod
