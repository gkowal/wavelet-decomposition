module wavelets

  implicit none

  integer, save ::  ncof = 0, ioff, joff
  real, dimension(:), allocatable, save :: cc,cr
  real, dimension(:), allocatable, save :: mc
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! DWT: subroutine for the discrete wavelet transform
!
!===============================================================================
!
  subroutine dwt(a, nn, order)

    implicit none

    real   , dimension(:,:,:), intent(inout) :: a
    integer, dimension(3)    , intent(in)    :: nn
    integer                  , intent(in)    :: order

! local variables
!
    integer                         :: i, j, k
!
!------------------------------------------------------------------------------
!
    if (order .ne. ncof) then
      call pwtset(order)
    endif

! transform along x-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do k = 1, nn(3)
      do j = 1, nn(2)
        call dwtstep(nn(1), a(1:nn(1),j,k))
      enddo
    enddo
!$omp end parallel do

! transform along y-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do k = 1, nn(3)
      do i = 1, nn(1)
        call dwtstep(nn(2), a(i,1:nn(2),k))
      enddo
    enddo
!$omp end parallel do

! transform along z-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do j = 1, nn(2)
      do i = 1, nn(1)
        call dwtstep(nn(3), a(i,j,1:nn(3)))
      enddo
    enddo
!$omp end parallel do

  end subroutine dwt
!
!===============================================================================
!
! IDWT: subroutine for the inverse discrete wavelet transform
!
!===============================================================================
!
  subroutine idwt(a, nn, order)

    implicit none

    real   , dimension(:,:,:), intent(inout) :: a
    integer, dimension(3)    , intent(in)    :: nn
    integer                  , intent(in)    :: order

! local variables
!
    integer                         :: i, j, k
!
!------------------------------------------------------------------------------
!
    if (order .ne. ncof) then
      call pwtset(order)
    endif

! transform along x-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do k = 1, nn(3)
      do j = 1, nn(2)
        call idwtstep(nn(1), a(1:nn(1),j,k))
      enddo
    enddo
!$omp end parallel do

! transform along y-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do k = 1, nn(3)
      do i = 1, nn(1)
        call idwtstep(nn(2), a(i,1:nn(2),k))
      enddo
    enddo
!$omp end parallel do

! transform along z-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do j = 1, nn(2)
      do i = 1, nn(1)
        call idwtstep(nn(3), a(i,j,1:nn(3)))
      enddo
    enddo
!$omp end parallel do

  end subroutine idwt
!
!===============================================================================
!
! DWM: subroutine for discrete wavelet averaging
!
!===============================================================================
!
  subroutine dwm(a, nn, order)

    implicit none

    real   , dimension(:,:,:), intent(inout) :: a
    integer, dimension(3)    , intent(in)    :: nn
    integer                  , intent(in)    :: order

! local variables
!
    integer                         :: i, j, k
!
!------------------------------------------------------------------------------
!
    if (order .ne. ncof) then
      call pwtset(order)
    endif

! transform along x-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do k = 1, nn(3)
      do j = 1, nn(2)
        call dwmean(nn(1), a(1:nn(1),j,k))
      enddo
    enddo
!$omp end parallel do

! transform along y-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do k = 1, nn(3)
      do i = 1, nn(1)
        call dwmean(nn(2), a(i,1:nn(2),k))
      enddo
    enddo
!$omp end parallel do

! transform along z-ddirection
!
!$omp parallel do default(private) shared(nn,a)
    do j = 1, nn(2)
      do i = 1, nn(1)
        call dwmean(nn(3), a(i,j,1:nn(3)))
      enddo
    enddo
!$omp end parallel do

  end subroutine dwm
!
!===============================================================================
!
! DWTSTEP: subroutine for the 1D discrete wavelet transform
!
!===============================================================================
!
  subroutine dwtstep(n, a)

    implicit none

    integer              , intent(in)    :: n
    real   , dimension(n), intent(inout) :: a

! local variables
!
    integer            :: nt
!
!------------------------------------------------------------------------------
!
    nt = n
    do while(nt >= 4)
      call pwt(a(1:nt),1)
      nt = nt / 2
    end do
  end subroutine dwtstep
!
!===============================================================================
!
! IDWTSTEP: subroutine for the 1D discrete wavelet transform
!
!===============================================================================
!
  subroutine idwtstep(n, a)

    implicit none

    integer              , intent(in)    :: n
    real   , dimension(n), intent(inout) :: a

! local variables
!
    integer            :: nt
!
!------------------------------------------------------------------------------
!
    nt = 4
    do while(nt <= n)
      call pwt(a(1:nt),-1)
      nt = nt * 2
    end do
  end subroutine idwtstep
!
!===============================================================================
!
! DWMEAN: subroutine for the 1D discrete wavelet averaging
!
!===============================================================================
!
  subroutine dwmean(n, a)

    implicit none

    integer              , intent(in)    :: n
    real   , dimension(n), intent(inout) :: a

! local variables
!
    integer            :: nt
!
!------------------------------------------------------------------------------
!
    nt = n
    do while(nt >= 4)
      call pwm(a(1:nt))
      nt = nt / 2
    end do
  end subroutine dwmean
!
!===============================================================================
!
! PWT: single run of the 1D discrete wavelet transform
!
!===============================================================================
!
  subroutine pwt(a, isign)

    implicit none

    real, dimension(:), intent(inout) :: a
    integer           , intent(in)    :: isign

    real   , dimension(size(a))   :: wksp
    integer, dimension(size(a)/2) :: jf, jr
    integer                       :: k, n, nh, nmod
!
!------------------------------------------------------------------------------
!
    n = size(a)

    if (n < 2) return

    nmod = ncof*n
    nh   = n/2

    wksp(:) = 0.0

    jf = iand(n-1,arth(2+nmod+ioff,2,nh))
    jr = iand(n-1,arth(2+nmod+joff,2,nh))

    if (isign >= 0) then
      do k = 1, ncof
        wksp(1:nh)   = wksp(1:nh)+cc(k)*a(jf+1)
        wksp(nh+1:n) = wksp(nh+1:n)+cr(k)*a(jr+1)
        if (k == ncof) exit
        jf = iand(n-1,jf+1)
        jr = iand(n-1,jr+1)
      end do
    else
      do k = 1, ncof
        wksp(jf+1)   = wksp(jf+1)+cc(k)*a(1:nh)
        wksp(jr+1)   = wksp(jr+1)+cr(k)*a(nh+1:n)
        if (k == ncof) exit
        jf = iand(n-1,jf+1)
        jr = iand(n-1,jr+1)
      end do
    end if

    a(:) = wksp(:)

  end subroutine pwt
!
!===============================================================================
!
! PWM: single run of the 1D discrete wavelet averaging
!
!===============================================================================
!
  subroutine pwm(a)

    implicit none

    real, dimension(:), intent(inout) :: a

    real   , dimension(size(a))   :: wksp
    integer, dimension(size(a)/2) :: jf, jr
    integer                       :: k, n, nh, nmod
    real, dimension(4)            :: c
!
!------------------------------------------------------------------------------
!
    n = size(a)

    if (n < 2) return

    nmod = ncof*n
    nh   = n/2

    wksp(:) = 0.0

    jf = iand(n-1,arth(2+nmod-ioff,2,nh))
    jr = iand(n-1,arth(2+nmod-joff,2,nh))

    do k = 1, ncof
      wksp(1:nh)   = wksp(1:nh)+mc(k)*a(jf+1)
      wksp(nh+1:n) = wksp(nh+1:n)+mc(k)*a(jr+1)
      if (k == ncof) exit
      jf = iand(n-1,jf+1)
      jr = iand(n-1,jr+1)
    end do

    a(:) = wksp(:)

  end subroutine pwm
!
!===============================================================================
!
! PWSET: subroutine to initialize coefficients for transforms
!
!===============================================================================
!
  subroutine pwtset(n)

    implicit none

    integer, intent(in) :: n

    real            :: sig
    real, parameter :: &
      c4(4)=(/&
                0.4829629131445341, 0.8365163037378079, &
                0.2241438680420134,-0.1294095225512604 /), &
      c12(12)=(/&
                0.111540743350, 0.494623890398, 0.751133908021, &
                0.315250351709,-0.226264693965,-0.129766867567, &
                0.097501605587, 0.027522865530,-0.031582039318, &
                0.000553842201, 0.004777257511,-0.001077301085 /), &
      c20(20)=(/&
                0.026670057901, 0.188176800078, 0.527201188932, &
                0.688459039454, 0.281172343661,-0.249846424327, &
               -0.195946274377, 0.127369340336, 0.093057364604, &
               -0.071394147166,-0.029457536822, 0.033212674059, &
                0.003606553567,-0.010733175483, 0.001395351747, &
                0.001992405295,-0.000685856695,-0.000116466855, &
                0.000093588670,-0.000013264203 /)
!
!------------------------------------------------------------------------------
!
    if (allocated(cc)) deallocate(cc)
    if (allocated(cr)) deallocate(cr)
    if (allocated(mc)) deallocate(mc)

    allocate(cc(n),cr(n),mc(n))

    ncof =  n
    ioff = -n/2
    joff = -n/2
    sig  = -1.0

    select case(n)
      case(4)
        cc=c4
      case(12)
        cc=c12
      case(20)
        cc=c20
      case default
        cc=c12
    end select

    cr(n:1:-1) =  cc
    cr(n:1:-2) = -cr(n:1:-2)

    mc(:) = 1. / ncof

  end subroutine pwtset
!
!===============================================================================
!
! ARTH: function to retrieve indexes
!
!===============================================================================
!
  function arth(first,increment,n)
    integer, intent(in) :: first, increment
    integer, intent(in) :: n
    integer, dimension(n) :: arth
    integer :: k, k2, temp

    if (n > 0) arth(1)=first
    if (n <= 16) then
      do k=2,n
        arth(k)=arth(k-1)+increment
      end do
    else
      do k=2,8
        arth(k)=arth(k-1)+increment
      end do
      temp=increment*8
      k=8
      do
        if (k >= n) exit
        k2=k+k
        arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
        temp=temp+temp
        k=k2
      end do
    end if
  end function arth

end module wavelets
