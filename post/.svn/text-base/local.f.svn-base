c-----------------------------------------------------------------------
c     file local.f
c     local defintions for most computers.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. local.
c     1. timer.
c     2. bin_open.
c     3. bin_close.
c     4. ascii_open.
c     5. ascii_close.
c     6. program_stop.
c-----------------------------------------------------------------------
c     subprogram 0. local.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE local_mod
      USE io_mod
      IMPLICIT NONE

      INTEGER, PARAMETER :: zero=0,one=1
      INTEGER, PARAMETER ::
     $     r4=SELECTED_REAL_KIND(6,37),
     $     r8=SELECTED_REAL_KIND(13,307)
      REAL(r8), PARAMETER :: pi=3.1415926535897932385_r8,
     $     twopi=2*pi,pisq=pi*pi,mu0=4e-7_r8*pi,
     $     rtod=180/pi,dtor=pi/180,alog10=2.302585093_r8,two=2._r8,
     $     third=1._r8/3._r8,half=.5_r8,min_eps=1.e-8

      COMPLEX(r8), PARAMETER :: ifac=(0,1)

      LOGICAL :: flag1d=.FALSE.,flag2d=.FALSE.,polar=.FALSE.,
     $     xperiodic,yperiodic
      CHARACTER(4) :: out_type="bin",outfile_type="bin"
      CHARACTER(16) :: job_type
      CHARACTER(80) :: indir=".",compdir=".",postout="."
      INTEGER :: dmout=1,nx,ny,np,nq,nbx,nstep
      REAL(r8) :: dt,dtmax,tmax,gr_curve

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. timer.
c     handles machine-dependent timing statistics.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     delcarations.
c-----------------------------------------------------------------------
      SUBROUTINE timer(mode,unit)
      
      INTEGER, INTENT(IN) :: mode,unit

      CHARACTER(10) :: date,time,zone
      INTEGER, DIMENSION(8) :: values
      REAL(4), SAVE :: seconds
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(1x,a,1p,e10.3,a)
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(mode == 0)THEN
         CALL DATE_AND_TIME(date,time,zone,values)
         seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
      ELSE
         CALL DATE_AND_TIME(date,time,zone,values)
         seconds=(values(5)*60+values(6))*60+values(7)+values(8)*1e-3
     $        -seconds
         WRITE(unit,10)"Wallclock time = ",seconds," seconds"
         WRITE(*,10)"Wallclock time = ",seconds," seconds"//CHAR(7)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     subprogram 2. bin_open.
c     opens a binary input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bin_open(unit,name,stat,pos)

      CHARACTER(*), INTENT(IN) :: name,stat,pos
      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     open file.
c-----------------------------------------------------------------------
      OPEN(UNIT=unit,FILE=name,STATUS=stat,POSITION=pos,
     $     FORM="UNFORMATTED")  !,CONVERT="BIG_ENDIAN")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bin_open
c-----------------------------------------------------------------------
c     subprogram 3. ascii_open.
c     opens a ascii input or output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE ascii_open(unit,name,stat)

      CHARACTER(*), INTENT(IN) :: name,stat
      INTEGER, INTENT(IN) :: unit
c-----------------------------------------------------------------------
c     open file.
c-----------------------------------------------------------------------
      OPEN(UNIT=unit,FILE=name,STATUS=stat)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ascii_open
c-----------------------------------------------------------------------
c     subprogram 4. program_stop.
c     terminates program with message, calls timer, closes output file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE program_stop(message)

      CHARACTER(*), INTENT(IN) :: message
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      CALL timer(1,out_unit)
      CLOSE(UNIT=out_unit)
      WRITE(*,'(1x,2a)') 'PROGRAM_STOP => ', TRIM(message)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      STOP
      END SUBROUTINE program_stop
      END MODULE local_mod
