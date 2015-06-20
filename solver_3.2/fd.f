c-----------------------------------------------------------------------
c     file fd.f.
c     compares jacobians evaluated by quadrature and finite differences.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. fd_mod.
c     1. fd_run.
c     2. fd_bound_run.
c     3. fd_bound_edge.
c-----------------------------------------------------------------------
c     subprogram 0. fd_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE fd_mod
      USE p2_sel_mod
      IMPLICIT NONE

      INTEGER, PARAMETER, PRIVATE :: num_max=9,seed=5

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fd_run.
c     evaluates jacobian by finite differences in a set of random
c     points and compares the result with analytical expression.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fd_run(psv)
      IMPLICIT NONE

      TYPE(p2_sel_type) :: psv

      INTEGER, PARAMETER :: num_step=100
      INTEGER :: itemp,jtemp,jac_err_out=666,jac_err_xdraw=667,
     $     i_step,i_max
      INTEGER, DIMENSION(num_max,4) :: coord_max  
      REAL(r8), PARAMETER :: deltau_min=1.0e-5,deltau_max=1.0e-1
      REAL(r8) :: deltau,u_step,temp_fx,temp_fy,temp_s,t=0
      REAL(r8), DIMENSION(num_max) :: vec_err
      REAL(r8), DIMENSION(1) :: xtemp,ytemp
      REAL(r8), DIMENSION(psv%nqty,1) :: utemp,uxtemp,uytemp,delu
      REAL(r8), DIMENSION(1,psv%nqty) :: fxtemp,fytemp,stemp,
     $     fxtemp_pl,fytemp_pl,stemp_pl,fxtemp_mn,fytemp_mn,stemp_mn
      REAL(r8), DIMENSION(1,psv%nqty,psv%nqty) :: fx_utemp,
     $     fx_uxtemp,fx_uytemp,fy_utemp,fy_uxtemp,fy_uytemp,s_utemp,
     $     s_uxtemp,s_uytemp
      REAL(r8), DIMENSION(psv%nqty,psv%nqty,3,3) :: rel_err
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"i",4x,"j",5x,"fx_u(i,j)",6x,"fy_u(i,j)",6x,
     $     "s_u(i,j)"/)
 20   FORMAT(2i5,1p,3e15.5)
 30   FORMAT(/4x,"i",4x,"j",5x,"fx_ux(i,j)",5x,"fy_ux(i,j)",5x,
     $     "s_ux(i,j)"/)
 40   FORMAT(/4x,"i",4x,"j",5x,"fx_uy(i,j)",5x,"fy_uy(i,j)",5x,
     $     "s_uy(i,j)"/)
 50   FORMAT(/4x,"delta u =",1p,e9.2/)
 60   FORMAT(/6x,"rel_err",6x,"i",4x,"j",4x,"u#",3x,"f#"/)
 70   FORMAT(1p,e15.5,4i5)
 80   FORMAT(2i5,1p,3e15.5/)
c-----------------------------------------------------------------------
c     initialize values of x, y, u, ux, uy with random numbers.
c-----------------------------------------------------------------------
      CALL put_seed(seed)
      CALL RANDOM_NUMBER(xtemp(1))
      CALL RANDOM_NUMBER(ytemp(1))
      CALL RANDOM_NUMBER(utemp(:,1))
      CALL RANDOM_NUMBER(uxtemp(:,1))
      CALL RANDOM_NUMBER(uytemp(:,1))
c-----------------------------------------------------------------------
c     evaluate analytical expressions for fx, fy, s, jacobian.
c-----------------------------------------------------------------------
      CALL job2_rhs(t,xtemp,ytemp,utemp,uxtemp,uytemp,fxtemp,fytemp,
     $     stemp)
      CALL job2_drdu(t,xtemp,ytemp,utemp,uxtemp,uytemp,fx_utemp,
     $     fx_uxtemp,fx_uytemp,fy_utemp,fy_uxtemp,fy_uytemp,s_utemp,
     $     s_uxtemp,s_uytemp)
c-----------------------------------------------------------------------
c     initialize deltau, step in deltau, "coordinates" of maximum
c     errors at deltau=deltau_min.
c-----------------------------------------------------------------------
      coord_max=0
      deltau=deltau_min
      u_step=(deltau_max/deltau_min)**(1.0_r8/REAL(num_step,r8))
c-----------------------------------------------------------------------
c     evaluate relative differences between analytical and finite 
c     difference expressions for jacobian.
c
c     open files for a table (jac_err.out) and for error plot
c     (jac_err.bin).
c-----------------------------------------------------------------------
      IF(mpi_rank==0)THEN
         OPEN(UNIT=jac_err_out,FILE="jac_err.out",ACTION="WRITE",
     $        STATUS="REPLACE")
         OPEN(UNIT=jac_err_xdraw,FILE="jac_err.bin",ACTION="WRITE",
     $        STATUS="REPLACE",FORM="UNFORMATTED")
         WRITE(jac_err_out,'(a)')
     $        "Relative errors in the analytical Jacobian (top entry)"
         WRITE(jac_err_out,'(a)')
     $        " vs. finite difference approximation (middle entry)"
         WRITE(jac_err_out,'(a)')" are given in the bottom entry."
         WRITE(jac_err_out,'(a,i2,a)')"Highest ",num_max,
     $        " errors are shown at the bottom of the file."
         WRITE(jac_err_out,*)
         WRITE(jac_err_out,'(a)')"Note that equations are indexed by i"
     $        //" and dependent variables by j."
         WRITE(jac_err_out,'(a)')"u# 1,2,3 correspond to u,ux,uy."
         WRITE(jac_err_out,'(a)')"f# 1,2,3 correspond to fx,fy,s."
         WRITE(jac_err_out,50)deltau
         WRITE(jac_err_out,10)
c-----------------------------------------------------------------------
c     start loop over deltau.
c-----------------------------------------------------------------------
         DO i_step=1,num_step+1
c-----------------------------------------------------------------------
c     evaluate errors in dfx/du, dfy/du, ds/du.
c-----------------------------------------------------------------------
            DO itemp=1,psv%nqty
               DO jtemp=1,psv%nqty
                  delu(:,1)=0.0
                  delu(jtemp,1)=deltau
                  CALL job2_rhs(t,xtemp,ytemp,utemp+delu,uxtemp,uytemp,
     $                 fxtemp_pl,fytemp_pl,stemp_pl)
                  CALL job2_rhs(t,xtemp,ytemp,utemp-delu,uxtemp,uytemp,
     $                 fxtemp_mn,fytemp_mn,stemp_mn)
                  temp_fx=(fxtemp_pl(1,itemp)
     $                 -fxtemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,1,1)=
     $                 fx_utemp(1,jtemp,itemp)-temp_fx
                  IF(ABS(fx_utemp(1,jtemp,itemp))>0 
     $                 .AND. ABS(temp_fx)>0)THEN
                     rel_err(itemp,jtemp,1,1)=rel_err(itemp,jtemp,1,1)
     $                    /fx_utemp(1,jtemp,itemp)
                  ENDIF
                  temp_fy=(fytemp_pl(1,itemp)
     $                 -fytemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,1,2)=
     $                 fy_utemp(1,jtemp,itemp)-temp_fy
                  IF(ABS(fy_utemp(1,jtemp,itemp))>0 
     $                 .AND. ABS(temp_fy)>0)THEN
                     rel_err(itemp,jtemp,1,2)=rel_err(itemp,jtemp,1,2)
     $                    /fy_utemp(1,jtemp,itemp)
                  ENDIF
                  temp_s=(stemp_pl(1,itemp)
     $                 -stemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,1,3)=
     $                 s_utemp(1,jtemp,itemp)-temp_s
                  IF(ABS(s_utemp(1,jtemp,itemp))>0
     $                 .AND. ABS(temp_s)>0)THEN
                     rel_err(itemp,jtemp,1,3)=rel_err(itemp,jtemp,1,3)
     $                    /s_utemp(1,jtemp,itemp)
                  ENDIF
c-----------------------------------------------------------------------
c     at the first step write results in jac_err.out.
c-----------------------------------------------------------------------
                  IF(i_step==1)THEN
                     WRITE(jac_err_out,20)itemp,jtemp,
     $                    fx_utemp(1,jtemp,itemp),
     $                    fy_utemp(1,jtemp,itemp),
     $                    s_utemp(1,jtemp,itemp)
                     WRITE(jac_err_out,20)itemp,jtemp,
     $                    temp_fx,temp_fy,temp_s
                     WRITE(jac_err_out,80)itemp,jtemp,
     $                    rel_err(itemp,jtemp,1,1),
     $                    rel_err(itemp,jtemp,1,2),
     $                    rel_err(itemp,jtemp,1,3)
                  ENDIF
               ENDDO
            ENDDO
            IF(i_step==1)WRITE(jac_err_out,10)
            IF(i_step==1)WRITE(jac_err_out,30)
c-----------------------------------------------------------------------
c     evaluate errors in dfx/dux, dfy/dux, ds/dux.
c-----------------------------------------------------------------------
            DO itemp=1,psv%nqty
               DO jtemp=1,psv%nqty
                  delu(:,1)=0.0_r8
                  delu(jtemp,1)=deltau
                  CALL job2_rhs(t,xtemp,ytemp,utemp,uxtemp+delu,uytemp,
     $                 fxtemp_pl,fytemp_pl,stemp_pl)
                  CALL job2_rhs(t,xtemp,ytemp,utemp,uxtemp-delu,uytemp,
     $                 fxtemp_mn,fytemp_mn,stemp_mn)
                  temp_fx=(fxtemp_pl(1,itemp)
     $                 -fxtemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,2,1)=
     $                 fx_uxtemp(1,jtemp,itemp)-temp_fx
                  IF(ABS(fx_uxtemp(1,jtemp,itemp))>0
     $                 .AND. ABS(temp_fx)>0)THEN
                     rel_err(itemp,jtemp,2,1)=rel_err(itemp,jtemp,2,1)
     $                    /fx_uxtemp(1,jtemp,itemp)
                  ENDIF
                  temp_fy=(fytemp_pl(1,itemp)
     $                 -fytemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,2,2)=
     $                 fy_uxtemp(1,jtemp,itemp)-temp_fy
                  IF(ABS(fy_uxtemp(1,jtemp,itemp))>0
     $                 .AND. ABS(temp_fy)>0)THEN
                     rel_err(itemp,jtemp,2,2)=rel_err(itemp,jtemp,2,2)
     $                    /fy_uxtemp(1,jtemp,itemp)
                  ENDIF
                  temp_s=(stemp_pl(1,itemp)
     $                 -stemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,2,3)=
     $                 s_uxtemp(1,jtemp,itemp)-temp_s
                  IF(ABS(s_uxtemp(1,jtemp,itemp))>0
     $                 .AND. ABS(temp_s)>0)THEN
                     rel_err(itemp,jtemp,2,3)=rel_err(itemp,jtemp,2,3)
     $                    /s_uxtemp(1,jtemp,itemp)
                  ENDIF
c-----------------------------------------------------------------------
c     at the first step write results in jac_err.out.
c-----------------------------------------------------------------------
                  IF(i_step==1)THEN
                     WRITE(jac_err_out,20)itemp,jtemp,
     $                    fx_uxtemp(1,jtemp,itemp),
     $                    fy_uxtemp(1,jtemp,itemp),
     $                    s_uxtemp(1,jtemp,itemp)
                     WRITE(jac_err_out,20)itemp,jtemp,
     $                    temp_fx,temp_fy,temp_s
                     WRITE(jac_err_out,80)itemp,jtemp,
     $                    rel_err(itemp,jtemp,2,1),
     $                    rel_err(itemp,jtemp,2,2),
     $                    rel_err(itemp,jtemp,2,3)
                  ENDIF
               ENDDO
            ENDDO
            IF(i_step==1)WRITE(jac_err_out,30)
            IF(i_step==1)WRITE(jac_err_out,40)
c-----------------------------------------------------------------------
c     evaluate errors in dfx/duy, dfy/duy, ds/duy.
c-----------------------------------------------------------------------
            DO itemp=1,psv%nqty
               DO jtemp=1,psv%nqty
                  delu(:,1)=0.0_r8
                  delu(jtemp,1)=deltau
                  CALL job2_rhs(t,xtemp,ytemp,utemp,uxtemp,uytemp+delu,
     $                 fxtemp_pl,fytemp_pl,stemp_pl)
                  CALL job2_rhs(t,xtemp,ytemp,utemp,uxtemp,uytemp-delu,
     $                 fxtemp_mn,fytemp_mn,stemp_mn)
                  temp_fx=(fxtemp_pl(1,itemp)
     $                 -fxtemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,3,1)=
     $                 fx_uytemp(1,jtemp,itemp)-temp_fx
                  IF(ABS(fx_uytemp(1,jtemp,itemp))>0
     $                 .AND. ABS(temp_fx)>0)THEN
                     rel_err(itemp,jtemp,3,1)=rel_err(itemp,jtemp,3,1)
     $                    /fx_uytemp(1,jtemp,itemp)
                  ENDIF
                  temp_fy=(fytemp_pl(1,itemp)
     $                 -fytemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,3,2)=
     $                 fy_uytemp(1,jtemp,itemp)-temp_fy
                  IF(ABS(fy_uytemp(1,jtemp,itemp))>0
     $                 .AND. ABS(temp_fy)>0)THEN
                     rel_err(itemp,jtemp,3,2)=rel_err(itemp,jtemp,3,2)
     $                    /fy_uytemp(1,jtemp,itemp)
                  ENDIF
                  temp_s=(stemp_pl(1,itemp)
     $                 -stemp_mn(1,itemp))/(2.0_r8*deltau)
                  rel_err(itemp,jtemp,3,3)=
     $                 s_uytemp(1,jtemp,itemp)-temp_s
                  IF(ABS(s_uytemp(1,jtemp,itemp))>0
     $                 .AND. ABS(temp_s)>0)THEN
                     rel_err(itemp,jtemp,3,3)=rel_err(itemp,jtemp,3,3)
     $                    /s_uytemp(1,jtemp,itemp)
                  ENDIF
c-----------------------------------------------------------------------
c     at the first step write results in jac_err.out.
c-----------------------------------------------------------------------
                  IF(i_step==1)THEN
                     WRITE(jac_err_out,20)itemp,jtemp,
     $                    fx_uytemp(1,jtemp,itemp),
     $                    fy_uytemp(1,jtemp,itemp),
     $                    s_uytemp(1,jtemp,itemp)
                     WRITE(jac_err_out,20)itemp,jtemp,
     $                    temp_fx,temp_fy,temp_s
                     WRITE(jac_err_out,80)itemp,jtemp,
     $                    rel_err(itemp,jtemp,3,1),
     $                    rel_err(itemp,jtemp,3,2),
     $                    rel_err(itemp,jtemp,3,3)
                  ENDIF
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     at the first step find num_max largest errors and write them
c     in jac_err.out.
c-----------------------------------------------------------------------
            IF(i_step==1)THEN
               WRITE(jac_err_out,40)
               WRITE(jac_err_out,60)
               DO i_max=1,num_max
                  coord_max(i_max,:)=MAXLOC(ABS(rel_err))
                  vec_err(i_max)=rel_err(coord_max(i_max,1),
     $                 coord_max(i_max,2),coord_max(i_max,3),
     $                 coord_max(i_max,4))
                  WRITE(jac_err_out,70)vec_err(i_max),
     $                 coord_max(i_max,1),
     $                 coord_max(i_max,2),
     $                 coord_max(i_max,3),
     $                 coord_max(i_max,4)
                  rel_err(coord_max(i_max,1),coord_max(i_max,2),
     $                 coord_max(i_max,3),coord_max(i_max,4))=0.0_r8
               ENDDO
               WRITE(jac_err_out,60)
               CLOSE(UNIT=jac_err_out)
            ENDIF
c-----------------------------------------------------------------------
c     write num_max largest errors vs. deltau in jac_err.bin.
c-----------------------------------------------------------------------
            IF(i_step>1)THEN
               DO i_max=1,num_max
                  vec_err(i_max)=rel_err(coord_max(i_max,1),
     $                 coord_max(i_max,2),coord_max(i_max,3),
     $                 coord_max(i_max,4))
               ENDDO
            ENDIF
            WRITE(jac_err_xdraw)REAL(LOG10(deltau),4),
     $           REAL(LOG10(ABS(vec_err)),4)
            deltau=deltau*u_step
c-----------------------------------------------------------------------
c     finish loop over deltau.
c-----------------------------------------------------------------------
         ENDDO
         CLOSE(UNIT=jac_err_xdraw)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fd_run
c-----------------------------------------------------------------------
c     subprogram 2. fd_bound_run.
c     run fd test for each edge and write summary of test.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fd_bound_run(psv)
      IMPLICIT NONE

      TYPE(p2_sel_type) :: psv

      INTEGER :: iedge,imax
      INTEGER, PARAMETER :: jac_err_out=666
      INTEGER, DIMENSION(4,num_max,3) :: coord
      REAL(r8), DIMENSION(4,num_max) :: vec_err
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 60   FORMAT(/6x,"rel_err",6x,"i",4x,"j",4x,"u#"/)
 70   FORMAT(1p,e15.5,3i5)
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      coord=0
      vec_err=0
c-----------------------------------------------------------------------
c     run fd test for each edge.
c-----------------------------------------------------------------------
      DO iedge=1,4
         CALL fd_bound_edge(psv,iedge,vec_err(iedge,:),coord(iedge,:,:))
      ENDDO
c-----------------------------------------------------------------------
c     write fd results for each edge into a file.
c-----------------------------------------------------------------------
      IF(mpi_rank==0)THEN
         OPEN(UNIT=jac_err_out,FILE="jac_bound_errsum.out",
     $        ACTION="WRITE",STATUS="REPLACE")
         WRITE(jac_err_out,'(a)')
     $        "A summary of highest relative errors in the analytical"
         WRITE(jac_err_out,'(a)')
     $        "Jacobian vs. finite difference approximation are given"
         WRITE(jac_err_out,'(a)')
     $        "for the boundary equation right-hand sides."
         WRITE(jac_err_out,'(a,i2,a)')"Highest ",num_max,
     $        " errors are shown."
         WRITE(jac_err_out,*)
         WRITE(jac_err_out,'(a)')
     $        "Note that equations are indexed by i and dependent"
         WRITE(jac_err_out,'(a)')
     $        "variables by j."
         WRITE(jac_err_out,'(a)')"u# 1,2,3,4,5,6 correspond to u,ux,"
     $        //"uy,uxx,uyy,uxy."
         WRITE(jac_err_out,*)
         WRITE(jac_err_out,'(a)') "Edges 1,2,3,4 are l,t,r,b."
         DO iedge=1,4
            WRITE(jac_err_out,*)
            WRITE(jac_err_out,*)"Edge #",iedge
            WRITE(jac_err_out,60)
            DO imax=1,num_max
               WRITE(jac_err_out,70)vec_err(iedge,imax),
     $              coord(iedge,imax,1),
     $              coord(iedge,imax,2),
     $              coord(iedge,imax,3)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fd_bound_run
c-----------------------------------------------------------------------
c     subprogram 3. fd_bound_edge.
c     evaluates jacobian at boundary by finite differences in a set of
c     random points and compares the result with analytical expression.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fd_bound_edge(psv,ie,vec_err,coord_max_ret)
      IMPLICIT NONE

      TYPE(p2_sel_type) :: psv
      INTEGER, INTENT(IN) :: ie
      INTEGER, DIMENSION(num_max,3), INTENT(INOUT) :: coord_max_ret
      REAL(r8), DIMENSION(num_max), INTENT(INOUT) :: vec_err

      INTEGER :: itemp,jtemp,i_step,imax
      INTEGER, DIMENSION(num_max,3) :: coord_max
      REAL(r8), PARAMETER :: deltau=1.e-5,t=0.
      REAL(r8) :: temp_c
      REAL(r8), DIMENSION(1) :: xtemp,ytemp
      REAL(r8), DIMENSION(2,1) :: dum
      REAL(r8), DIMENSION(psv%nqty,1) :: utemp,uxtemp,uytemp,uxxtemp,
     $     uyytemp,uxytemp,delu
      REAL(r8), DIMENSION(1,psv%nqty) :: ctemp,ctemp_pl,ctemp_mn
      REAL(r8), DIMENSION(1,psv%nqty,psv%nqty) :: c_utemp,
     $     c_uxtemp,c_uytemp,c_uxxtemp,c_uyytemp,c_uxytemp
      REAL(r8), DIMENSION(psv%nqty,psv%nqty,6) :: rel_err
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/4x,"i",4x,"j",5x,"c_u(i,j)"/)
 20   FORMAT(2i5,1p,1e15.5)
 30   FORMAT(/4x,"i",4x,"j",5x,"c_ux(i,j)"/)
 40   FORMAT(/4x,"i",4x,"j",5x,"c_uy(i,j)"/)
 42   FORMAT(/4x,"i",4x,"j",5x,"c_uxx(i,j)"/)
 44   FORMAT(/4x,"i",4x,"j",5x,"c_uyy(i,j)"/)
 46   FORMAT(/4x,"i",4x,"j",5x,"c_uxy(i,j)"/)
 50   FORMAT(/4x,"delta u =",1p,e9.2/)
 60   FORMAT(/6x,"rel_err",6x,"i",4x,"j",4x,"u#"/)
 70   FORMAT(1p,e15.5,3i5)
 80   FORMAT(2i5,1p,1e15.5/)
c-----------------------------------------------------------------------
c     initialize values of x, y, u, ux, uy with random numbers.
c-----------------------------------------------------------------------
      CALL put_seed(seed)
      CALL RANDOM_NUMBER(xtemp)
      CALL RANDOM_NUMBER(ytemp)
      CALL RANDOM_NUMBER(utemp(:,1))
      CALL RANDOM_NUMBER(dum(:,1))
      CALL RANDOM_NUMBER(uxtemp(:,1))
      CALL RANDOM_NUMBER(uytemp(:,1))
      CALL RANDOM_NUMBER(uxxtemp(:,1))
      CALL RANDOM_NUMBER(uyytemp(:,1))
      CALL RANDOM_NUMBER(uxytemp(:,1))
c-----------------------------------------------------------------------
c     evaluate analytical expressions for fx, fy, s, jacobian.
c-----------------------------------------------------------------------
      CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,utemp,uxtemp,
     $     uytemp,uxxtemp,uyytemp,uxytemp,ctemp)
      CALL job2_edge_drdu(psv%edges(ie),t,xtemp,ytemp,dum,utemp,uxtemp,
     $     uytemp,uxxtemp,uyytemp,uxytemp,c_utemp,c_uxtemp,c_uytemp,
     $     c_uxxtemp,c_uyytemp,c_uxytemp)
c-----------------------------------------------------------------------
c     initialize "coordinates" of maximum errors.
c-----------------------------------------------------------------------
      coord_max=0
c-----------------------------------------------------------------------
c     evaluate relative differences between analytical and finite 
c     difference expressions for jacobian.
c-----------------------------------------------------------------------
      IF(mpi_rank==0)THEN
         WRITE(*,*)"performing fd test for edge #",ie
c-----------------------------------------------------------------------
c     evaluate errors in dc/du.
c-----------------------------------------------------------------------
         DO itemp=1,psv%nqty
            DO jtemp=1,psv%nqty
               delu(:,1)=0.0
               delu(jtemp,1)=deltau
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp+delu,uxtemp,uytemp,uxxtemp,uyytemp,uxytemp,
     $              ctemp_pl)
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp-delu,uxtemp,uytemp,uxxtemp,uyytemp,uxytemp,
     $              ctemp_mn)
               temp_c=(ctemp_pl(1,itemp)
     $              -ctemp_mn(1,itemp))/(2.0_r8*deltau)
               rel_err(itemp,jtemp,1)=
     $              c_utemp(1,jtemp,itemp)-temp_c
               IF(ABS(c_utemp(1,jtemp,itemp))>0 
     $              .AND. ABS(temp_c)>0)THEN
                  rel_err(itemp,jtemp,1)=rel_err(itemp,jtemp,1)
     $                 /c_utemp(1,jtemp,itemp)
               ENDIF
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     evaluate errors in dc/dux
c-----------------------------------------------------------------------
         DO itemp=1,psv%nqty
            DO jtemp=1,psv%nqty
               delu(:,1)=0.0_r8
               delu(jtemp,1)=deltau
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp+delu,uytemp,uxxtemp,uyytemp,uxytemp,
     $              ctemp_pl)
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp-delu,uytemp,uxxtemp,uyytemp,uxytemp,
     $              ctemp_mn)
               temp_c=(ctemp_pl(1,itemp)
     $              -ctemp_mn(1,itemp))/(2.0_r8*deltau)
               rel_err(itemp,jtemp,2)=
     $              c_uxtemp(1,jtemp,itemp)-temp_c
               IF(ABS(c_uxtemp(1,jtemp,itemp))>0
     $              .AND. ABS(temp_c)>0)THEN
                  rel_err(itemp,jtemp,2)=rel_err(itemp,jtemp,2)
     $                 /c_uxtemp(1,jtemp,itemp)
               ENDIF
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     evaluate errors in dc/duy
c-----------------------------------------------------------------------
         DO itemp=1,psv%nqty
            DO jtemp=1,psv%nqty
               delu(:,1)=0.0_r8
               delu(jtemp,1)=deltau
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp,uytemp+delu,uxxtemp,uyytemp,uxytemp,
     $              ctemp_pl)
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp,uytemp-delu,uxxtemp,uyytemp,uxytemp,
     $              ctemp_mn)
               temp_c=(ctemp_pl(1,itemp)
     $              -ctemp_mn(1,itemp))/(2.0_r8*deltau)
               rel_err(itemp,jtemp,3)=
     $              c_uytemp(1,jtemp,itemp)-temp_c
               IF(ABS(c_uytemp(1,jtemp,itemp))>0
     $              .AND. ABS(temp_c)>0)THEN
                  rel_err(itemp,jtemp,3)=rel_err(itemp,jtemp,3)
     $                 /c_uytemp(1,jtemp,itemp)
               ENDIF
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     evaluate errors in dc/duxx
c-----------------------------------------------------------------------
         DO itemp=1,psv%nqty
            DO jtemp=1,psv%nqty
               delu(:,1)=0.0_r8
               delu(jtemp,1)=deltau
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp,uytemp,uxxtemp+delu,uyytemp,uxytemp,
     $              ctemp_pl)
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp,uytemp,uxxtemp-delu,uyytemp,uxytemp,
     $              ctemp_mn)
               temp_c=(ctemp_pl(1,itemp)
     $              -ctemp_mn(1,itemp))/(2.0_r8*deltau)
               rel_err(itemp,jtemp,4)=
     $              c_uxxtemp(1,jtemp,itemp)-temp_c
               IF(ABS(c_uxxtemp(1,jtemp,itemp))>0
     $              .AND. ABS(temp_c)>0)THEN
                  rel_err(itemp,jtemp,4)=rel_err(itemp,jtemp,4)
     $                 /c_uxxtemp(1,jtemp,itemp)
               ENDIF
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     evaluate errors in dc/duyy
c-----------------------------------------------------------------------
         DO itemp=1,psv%nqty
            DO jtemp=1,psv%nqty
               delu(:,1)=0.0_r8
               delu(jtemp,1)=deltau
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp,uytemp,uxxtemp,uyytemp+delu,uxytemp,
     $              ctemp_pl)
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp,uytemp,uxxtemp,uyytemp-delu,uxytemp,
     $              ctemp_mn)
               temp_c=(ctemp_pl(1,itemp)
     $              -ctemp_mn(1,itemp))/(2.0_r8*deltau)
               rel_err(itemp,jtemp,5)=
     $              c_uyytemp(1,jtemp,itemp)-temp_c
               IF(ABS(c_uyytemp(1,jtemp,itemp))>0
     $              .AND. ABS(temp_c)>0)THEN
                  rel_err(itemp,jtemp,5)=rel_err(itemp,jtemp,5)
     $                 /c_uyytemp(1,jtemp,itemp)
               ENDIF
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     evaluate errors in dc/duxy
c-----------------------------------------------------------------------
         DO itemp=1,psv%nqty
            DO jtemp=1,psv%nqty
               delu(:,1)=0.0_r8
               delu(jtemp,1)=deltau
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp,uytemp,uxxtemp,uyytemp,uxytemp+delu,
     $              ctemp_pl)
               CALL job2_edge_rhs(psv%edges(ie),t,xtemp,ytemp,dum,
     $              utemp,uxtemp,uytemp,uxxtemp,uyytemp,uxytemp-delu,
     $              ctemp_mn)
               temp_c=(ctemp_pl(1,itemp)
     $              -ctemp_mn(1,itemp))/(2.0_r8*deltau)
               rel_err(itemp,jtemp,6)=
     $              c_uxytemp(1,jtemp,itemp)-temp_c
               IF(ABS(c_uxytemp(1,jtemp,itemp))>0
     $              .AND. ABS(temp_c)>0)THEN
                  rel_err(itemp,jtemp,6)=rel_err(itemp,jtemp,6)
     $                 /c_uxytemp(1,jtemp,itemp)
               ENDIF
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     write largest errors in returned arrays.
c-----------------------------------------------------------------------
         DO imax=1,num_max
            coord_max(imax,:)=MAXLOC(ABS(rel_err))
            vec_err(imax)=rel_err(coord_max(imax,1),
     $           coord_max(imax,2),coord_max(imax,3))
            rel_err(coord_max(imax,1),coord_max(imax,2),
     $           coord_max(imax,3))=zero
            coord_max_ret(imax,:)=coord_max(imax,:)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fd_bound_edge
      END MODULE fd_mod
