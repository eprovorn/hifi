c-----------------------------------------------------------------------
c     file io.f.
c     input and output unit declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module declarations.
c-----------------------------------------------------------------------
      MODULE io_mod
      IMPLICIT NONE

      INTEGER :: in_unit=1
      INTEGER :: out_unit=2
      INTEGER :: bin_unit=3
      INTEGER :: header_unit=4
      INTEGER :: contour_unit=12
      INTEGER :: grid_unit=13
      INTEGER :: xdmf_unit=14

      INTEGER :: xt_unit=21
      INTEGER :: yt_unit=22
      INTEGER :: xy_unit=23

      INTEGER :: UxyT_unit=31
      INTEGER :: VecSqInt_unit=32
      INTEGER :: Ucontour_unit=33
      INTEGER :: Uprofile_unit=34
      INTEGER :: dUdt_unit=35
      INTEGER :: belgrid_unit=36
      INTEGER :: Uamp_unit=37
      INTEGER :: maxUvsT_unit=38
      INTEGER :: width_unit=39
      INTEGER :: belt_unit=40
      INTEGER :: Bmodes_unit=41
      INTEGER :: qprof_unit=42
      INTEGER :: swten_unit=43
      INTEGER :: oned_unit=44
      INTEGER :: Uprofile_unit_asc=45

      INTEGER :: file_unit=99

      END MODULE io_mod
