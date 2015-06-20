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
      INTEGER :: dat_unit=4
      INTEGER :: grid_unit=8
      INTEGER :: jac_out_unit=11
      INTEGER :: jac_bin_unit=12
      INTEGER :: belin_unit=13
      INTEGER :: belbin_unit=14
      INTEGER :: beldat_unit=15
      INTEGER :: coil_unit=16
      INTEGER :: equil_unit=17
      INTEGER :: debug_unit=98
      INTEGER :: debug2_unit=99

      END MODULE io_mod
