!------------------------------------------------------------------------------
! # Fortran90 Common Routines Set #
!------------------------------------------------------------------------------
! ## MODULE: Common
!
!> @author
!> Luca Malavolta, Department of Physics & Astronomy, UniPD
!
! DESCRIPTION:
!> This module contains common parameters to all subroutines and main programs
!> to avoid conflicts (i.e. by fixing the length of strings and the numerical
!> precision) and general physics parameters
!
! REVISION HISTORY:
!> @var 1.00  01 Apr 2013 - Initial Version
!> @var 1.01  09 May 2013 - Comments and Intents fixed
!> @var 1.02  10 May 2013 - Bug in sgm2fwhm corrected
!------------------------------------------------------------------------------
! For more information on the commenting standard and additional
! Special Commands see http://www.doxygen.nl/manual/commands.html#cmd_intro
!------------------------------------------------------------------------------

module Common
  implicit none

  !> PR : precision used by computational subroutines (LMDIF, etc)
  !! PD : precision used to store the data (it may be lower and fixed by the storing system)
  !! Right now, only PD=8 mode in reading subroutine is supported
  integer, parameter :: PR = 8, PD = 8

  !> These are the standard sizes for the string variables
  !! using personalized values may affect the functioning of the programs
  integer, parameter :: nch_file = 256, nch_rad = 80, nch_object = 14
  integer, parameter :: nch_syn = 63, nch_presyn = 18
  integer, parameter :: nch_key = 16, nch_com = 64, nch_system = 600

  real , parameter :: nullval_r = -99.999_4
  real (PR), parameter :: nullval_d = -99.99999_PR

  real (PR), parameter :: sgm2fwhm= 2._PR*sqrt(2._PR*log(2._PR)) !> @var FWHM/SIGMA ratio for a gaussian
  real (PR), parameter :: PId=3.1415926535898_PR   !> @var Pi
  real (PR), parameter :: Euler=2.71828183_PR      !> @var Euler number
  real (PR), parameter :: cc=299792.4580_PR        !> @var Speed of Light in Km/s
  real (PR), parameter :: G_grav=6.67398e-11_PR    !> @var Gravitational Constant in m^3 kg^-1 s^-2
  real (PR), parameter :: M_sun =1.98892e30_PR     !> @var Mass of the sun in Kg
  real (PR), parameter :: M_jup = M_sun/1047.56_PR, M_earth=M_sun/333000._PR !> @var Masses of other bodies
  real (PR), parameter :: AU_m = 149597870700._PR   !> @var Astronomical Unit, in meters

  !  integer, parameter :: answer=42

contains

  !---------------------------------------------------------------------------
  ! ### Subroutine: get_lun(unit)
  !---------------------------------------------------------------------------
  !> @author
  !> Luca Malavolta, Department of Physics & Astronomy, UniPD
  !
  ! DESCRIPTION:
  !> This routine will give back a free logical unit number by using the
  !> FITSIO subroutine FTGIOU. A successive call of this subroutine will give
  !> back a different number, thus avoiding conflicts
  !> @brief
  !> A free logical unit is given back
  !
  !> @param[in] none
  !> @param[out] unit
  !> @return unit
  !---------------------------------------------------------------------------

  subroutine get_lun(unit)
    integer, intent(out) :: unit
    integer :: status

    status = 0
    call FTGIOU(unit,status)
    !we just want a free unit number without allocating it

    return
  end subroutine get_lun

  !---------------------------------------------------------------------------
  !> @author
  !> Luca Malavolta, Department of Physics & Astronomy, UniPD
  !
  ! DESCRIPTION:
  !> The logical unit "occupied" before is being freed. This subroutine
  !> should be run at the end of the main program
  !> @brief
  !> The given free logical unit is freed
  !
  !> @param[in] unit
  !> @param[out] none
  !> @return none
  !---------------------------------------------------------------------------

  subroutine free_lun(unit)
    integer, intent(in) :: unit
    integer :: status

    status = 0
    call FTFIOU(unit,status)
    !we just want a free unit number without allocating it

    return
  end subroutine free_lun


  subroutine int2chr4(int,chr4)
    integer :: int
    character (len=4) :: chr4

    chr4 = '----'
    if (int.ge.0    .and.int.lt.10   ) write(chr4,'(A3,I1)') '000', int
    if (int.ge.10   .and.int.lt.100  ) write(chr4,'(A2,I2)') '00' , int
    if (int.ge.100  .and.int.lt.1000 ) write(chr4,'(A1,I3)') '0'  , int
    if (int.ge.1000 .and.int.lt.1000 ) write(chr4,'(A1,I3)') '0'  , int
    if (int.gt.1000                  ) write(chr4,'(   I4)')        int

    if (int.lt.0    .and.int.gt.-10  ) write(chr4,'(A3,I1)') '-00', abs(int)
    if (int.le.-10  .and.int.gt.-100 ) write(chr4,'(A2,I2)') '-0' , abs(int)
    if (int.le.-100 .and.int.gt.-1000) write(chr4,'(A1,I3)') '-'  , abs(int)
    return
  end subroutine int2chr4


end module common
