!!******************************************************************************
!!
!! module: params - subroutines to read parameters file.
!!
!! Copyright (C) 2007 Grzegorz Kowal <kowal@astro.wisc.edu>
!!
!!******************************************************************************
!!
!!  This file is part of Spectrum.
!!
!!  Godunov-MHD is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  Foobar is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!
module params

  implicit none

! input parameters
!
  character(len = 128), save :: idir    = "./"         ! input directory
  character(len = 128), save :: odir    = "./"         ! output directory
  character(len =   4), save :: field   = "velo"       ! spectrum for which field
  character(len =   3), save :: frame   = "glo"        ! reference frame 'glo' or 'loc'
  character(len =   1), save :: dtype   = "W"          ! method of decomposition:
                                                       ! 'F' - Fourier
                                                       ! 'W' - Wavelet
  character(len =   4), save :: fformat = 'fits'       ! file format:
                                                       ! 'fits' - FITS
                                                       ! 'hdf5' - HDF5
  character(len =   1), save :: ftype   = 'r'          ! file type: 'r', 'p', 'f'
  character(len =   1), save :: sdir    = 'r'          ! direction along which spectrum is taken
                                                       ! 'R' - radial
                                                       ! 'X', 'Y', 'Z' - x, y, z direction
  character(len=32)   , save :: ffmt    = "a1,i6.6,'_',i5.5,'.h5'"
  character(len= 1)   , save :: type    = "I"          ! simulation type:
                                                       ! 'I' for incompressible,
                                                       ! 'C' for compressible
  integer             , save :: fnumber = 0            ! file number
  integer             , save :: maxlen  = 1            ! maximum separation length
  integer             , save :: maxexp  = 1            ! maximum exponent
  integer             , save :: nvecs   = 1            ! number of vectors
  integer             , save :: nshots  = 1            ! number of shots
  integer             , save :: kmax    = 512          ! maximum wavelength
  integer             , save :: wvlen   = 12           ! length of wavelet
  real                , save :: xc      = 0.0          ! position of P spectrum
  real                , save :: yc      = 0.0
  real                , save :: zc      = 0.0
  real                , save :: rd      = 1.0          ! width of window
  real                , save :: kd      = 1.0          ! width of window in Fourier space
  real                , save :: xlen    = 1.0
  real                , save :: ylen    = 1.0
  real                , save :: zlen    = 1.0
  real                , save :: dens    = 1.0          ! the initial density
  real                , save :: csnd    = 1.0          ! the sound of speed value
  real                , save :: thres   = 1.0          ! wavelet amplitude threshold

! common variables
!
  character(len =   4), save :: vars(7)
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! read_config: subroutine to read config file and to fill proper parameters
!
!===============================================================================
!
  subroutine read_params
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    character(len =   *), parameter :: configfile = './params.in' ! config file with parameters
    character(len=255) :: line, name, value
    integer            :: l, i, ios
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    open(unit=1, file=configfile, err=100)

10  read(unit=1, fmt="(a)", end=20, err=200 ) line

    call parse_line(line, name, value)

    select case(name)
      case ('idir')
        l = len_trim(value)
        write(idir   , "(a)" ) value(2:l-1)
      case ('odir')
        l = len_trim(value)
        write(odir   , "(a)" ) value(2:l-1)
      case ('dtype')
        l = len_trim(value)
        write(dtype  , "(a)" ) value(2:l-1)
      case ('frame')
        l = len_trim(value)
        write(frame  , "(a)" ) value(2:l-1)
      case ('field')
        l = len_trim(value)
        write(field  , "(a)" ) value(2:l-1)
      case ('fformat')
        l = len_trim(value)
        write(fformat, "(a)" ) value(2:l-1)
      case ('ftype')
        l = len_trim(value)
        write(ftype  , "(a)" ) value(2:l-1)
      case ('sdir')
        l = len_trim(value)
        write(sdir   , "(a)" ) value(2:l-1)
      case ('ffmt')
        l = len_trim(value)
        write(ffmt   , "(a)" ) value(2:l-1)
      case ('type')
        l = len_trim(value)
        write(type   , "(a)" ) value(2:l-1)
      case ('fnumber')
        read (value  , "(i6)") fnumber
      case ('maxlen')
        read (value  , "(i6)") maxlen
      case ('maxexp')
        read (value  , "(i6)") maxexp
      case ('nvecs')
        read (value  , "(i6)") nvecs
      case ('nshots')
        read (value  , "(i6)") nshots
      case ('kmax')
        read (value  , "(i6)") kmax
      case ('wvlen')
        read (value  , "(i6)") wvlen
      case ('xc')
        read (value  , *) xc
      case ('yc')
        read (value  , *) yc
      case ('zc')
        read (value  , *) zc
      case ('rd')
        read (value  , *) rd
      case ('kd')
        read (value  , *) kd
      case ('xlen')
        read (value  , *) xlen
      case ('ylen')
        read (value  , *) ylen
      case ('zlen')
        read (value  , *) zlen
      case ('dens')
        read (value  , *) dens
      case ('csnd')
        read (value  , *) csnd
      case ('thres')
        read (value  , *) thres
      case default
    end select

    go to 10

20  close(1)
!
    return
!
100 print *, 'Error opening file ', configfile
    stop
200 print *, 'Error reading file ', configfile
    stop
!
    return
!
  end subroutine read_params
!
!===============================================================================
!
! parse_line: subroutine to parse line into name and value of parameter
!
!===============================================================================
!
  subroutine parse_line(line, name, value)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: name, value

    integer :: l, i
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    l = len_trim(line)

    i = index( line, '=' )

    name  = trim(adjustl(line(1:i-1)))
    value = trim(adjustl(line(i+1:l)))
!
  end subroutine parse_line

end module params
