program wtmodes

  use params  , only : read_params, idir, odir, fformat, dtype, wvlen, thres   &
                     , vars, csnd
  use mod_hdf5, only : hdf5_init, hdf5_get_dims, hdf5_read_var
  use mod_fits, only : fits_init, fits_get_dims, fits_read_var                 &
                     , fits_put_data_1d, fits_put_data_2d, fits_put_data
  use wavelets, only : dwt, idwt, dwm

  implicit none

! local variables
!
  character(len = 255) :: fl
  integer              :: info, nlen, dims(3)
  integer              :: lb, le, ll, i, j, k, p
  real                 :: timer, factor, ct, aa, dd, fp, fr, pr
  real                 :: uu, ux, uy, uz, px, py, pz, rx, ry, rz
  real                 :: kpx, kpy, kpz, krx, kry, krz

! local allocatable arrays
!
  real, dimension(:)    , allocatable :: kx, ky, kz
  real, dimension(:,:,:), allocatable :: tmp, al, dn, vx, vy, vz, bx, by, bz
  real, dimension(:,:,:), allocatable :: ax, ay, az
  real, dimension(:,:,:), allocatable :: sx, sy, sz
  real, dimension(:,:,:), allocatable :: fx, fy, fz
!
!-------------------------------------------------------------------------------
!
! read parameters from config.in
!
  write( *, "('TASK          : ',a)" ) "MHD modes decomposition of a vector field"
  write( *, "('INFO          : ',a)" ) "reading parameters"
  call read_params

! print some information
!
  write( *, "('INDIR         : ',a)"    ) trim(idir)
  write( *, "('OUDIR         : ',a)"    ) trim(odir)
  write( *, "('FORMAT        : ',a)"    ) trim(fformat)
  write( *, "('DECOMPOSITION : ',a)"    ) trim(dtype)
  write( *, "('WAVELET       : ',i3)"   ) wvlen
  write( *, "('THRESHOLD     : ',f10.6)") thres

! prepare variable names to read
!
  write(vars(1), '(a)') 'dens'
  write(vars(2), '(a)') 'velx'
  write(vars(3), '(a)') 'vely'
  write(vars(4), '(a)') 'velz'
  write(vars(5), '(a)') 'magx'
  write(vars(6), '(a)') 'magy'
  write(vars(7), '(a)') 'magz'

! get data dimensions
!
  select case (trim(fformat))
  case('hdf5')
    call hdf5_init()
    call hdf5_get_dims(dims)
  case default
    call fits_init()
    call fits_get_dims(dims)
  endselect

  write( *, "('CSND          : ',f10.6)") csnd

! allocate variables
!
  write( *, "('INFO          : ',a)" ) "allocating variables"

! allocate vectors
!
  allocate(kx(dims(1)))
  allocate(ky(dims(2)))
  allocate(kz(dims(3)))

! allocate arrays
!
  allocate(tmp( dims(1),  dims(2),  dims(3)))
  allocate(al ( dims(1),  dims(2),  dims(3)))
  allocate(dn ( dims(1),  dims(2),  dims(3)))
  allocate(vx ( dims(1),  dims(2),  dims(3)))
  allocate(vy ( dims(1),  dims(2),  dims(3)))
  allocate(vz ( dims(1),  dims(2),  dims(3)))
  allocate(bx ( dims(1),  dims(2),  dims(3)))
  allocate(by ( dims(1),  dims(2),  dims(3)))
  allocate(bz ( dims(1),  dims(2),  dims(3)))
  allocate(ax ( dims(1),  dims(2),  dims(3)))
  allocate(ay ( dims(1),  dims(2),  dims(3)))
  allocate(az ( dims(1),  dims(2),  dims(3)))
  allocate(sx ( dims(1),  dims(2),  dims(3)))
  allocate(sy ( dims(1),  dims(2),  dims(3)))
  allocate(sz ( dims(1),  dims(2),  dims(3)))
  allocate(fx ( dims(1),  dims(2),  dims(3)))
  allocate(fy ( dims(1),  dims(2),  dims(3)))
  allocate(fz ( dims(1),  dims(2),  dims(3)))
!
! read density, magnetic field, calculate mean field and local directions
!
  write( *, "('INFO          : ',a)" ) "reading variables"

  if (fformat .eq. 'hdf5') then
    call hdf5_read_var('dens', dn)
    call hdf5_read_var('velx', vx)
    call hdf5_read_var('vely', vy)
    call hdf5_read_var('velz', vz)
    call hdf5_read_var('magx', bx)
    call hdf5_read_var('magy', by)
    call hdf5_read_var('magz', bz)
  end if
  if (fformat .eq. 'fits') then
    call fits_read_var('dens', dn)
    call fits_read_var('velx', vx)
    call fits_read_var('vely', vy)
    call fits_read_var('velz', vz)
    call fits_read_var('magx', bx)
    call fits_read_var('magy', by)
    call fits_read_var('magz', bz)
  end if

  timer = secnds(0.0)

  write( *, "('TASK          : ',a)" ) "calculating local mean magnetic field"

  call dwm(bx, dims, wvlen)
  call dwm(by, dims, wvlen)
  call dwm(bz, dims, wvlen)

  write( *, "('TASK          : ',a)" ) "calculating direction of local mean field"

  al  = bx**2 + by**2 + bz**2
  tmp = sqrt(al)
  where(tmp .eq. 0.0)
    tmp = 1.0
  end where
  bx = bx / tmp
  by = by / tmp
  bz = bz / tmp

! prepare array of alfa values
!
  write( *, "('TASK          : ',a)" ) "calculating local mean density"
  call dwm(dn, dims, wvlen)

  write( *, "('TASK          : ',a)" ) "calculating alpha parameter"
  al = csnd**2 * dn / al

! deallocate unnecessary variables
!
  if (allocated(dn) )  deallocate(dn)

! calculate wavelet transform of the velocity components
!
  write( *, "('TASK          : ',a)" ) "calculating DWT of velocity field"
  call dwt(vx, dims, wvlen)
  call dwt(vy, dims, wvlen)
  call dwt(vz, dims, wvlen)

! preparing the wave vectors corresponding to the wavelet coefficients
!
  write( *, "('TASK          : ',a)" ) "preparing wave vectors"

  le = dims(1)
  ll = dims(1)
  do while(le .gt. 2)
    le = ll
    ll = ll / 2
    lb = le - ll + 1
    kx(lb:le) = ll
  end do
  kx(1) = 1.0

  le = dims(2)
  ll = dims(2)
  do while(le .gt. 2)
    le = ll
    ll = ll / 2
    lb = le - ll + 1
    ky(lb:le) = ll
  end do
  ky(1) = 1.0

  le = dims(3)
  ll = dims(3)
  do while(le .gt. 2)
    le = ll
    ll = ll / 2
    lb = le - ll + 1
    kz(lb:le) = ll
  end do
  kz(1) = 1.0

  if (wvlen .eq. 4) then
    kx = 0.726562 * kx
    ky = 0.726562 * ky
    kz = 0.726562 * kz
  end if
  if (wvlen .eq. 12) then
    kx = 0.687500 * kx
    ky = 0.687500 * ky
    kz = 0.687500 * kz
  end if
  if (wvlen .eq. 20) then
    kx = 0.679688 * kx
    ky = 0.679688 * ky
    kz = 0.679688 * kz
  end if

! progress variables
!
  p  = 0
  pr = 100.0 / (dims(1) * dims(2))

! decompose Alfven mode
!
  write( *, "('TASK          : ',a)" ) "decomposing MHD modes"

  do k = 1, dims(3)
    do j = 1, dims(2)
!$omp parallel do
      do i = 1, dims(1)

! obtain the direction of the local mean magnetic field
!
        ux = bx(i,j,k)
        uy = by(i,j,k)
        uz = bz(i,j,k)

! obtain direction parallel to the local field
!
        uu = kx(i) * ux + ky(j) * uy + kz(k) * uz
        kpx = uu * ux
        kpy = uu * uy
        kpz = uu * uz

! obtain direction perpendicular to the local field
!
        krx = kx(i) - kpx
        kry = ky(j) - kpy
        krz = kz(k) - kpz

! normalize parallel and perpendicular versors
!
        uu = sqrt(kpx * kpx + kpy * kpy + kpz * kpz)
        if (uu .eq. 0.0) uu = 1.0
        px = kpx / uu
        py = kpy / uu
        pz = kpz / uu
        uu = sqrt(krx * krx + kry * kry + krz * krz)
        if (uu .eq. 0.0) uu = 1.0
        rx = krx / uu
        ry = kry / uu
        rz = krz / uu

! calculate xi_alphen = k_par x k_per
!
        ux = py * rz - ry * pz
        uy = pz * rx - rz * px
        uz = px * ry - rx * py
        uu = sqrt(ux * ux + uy * uy + uz * uz)
        if (uu .eq. 0.0) uu = 1.0
        ux = ux / uu
        uy = uy / uu
        uz = uz / uu

! project velocity vector
!
        uu = ux * vx(i,j,k) + uy * vy(i,j,k) + uz * vz(i,j,k)

        ax(i,j,k) = uu * ux
        ay(i,j,k) = uu * uy
        az(i,j,k) = uu * uz

! calculate the angle between B and k
!
        uu = sqrt(kx(i) * kx(i) + ky(j) * ky(j) + kz(k) * kz(k))
        if (uu .eq. 0.0) uu = 1.0
        ct = ((bx(i,j,k) * kx(i) + by(i,j,k) * ky(j) + bz(i,j,k) * kz(k)) / uu)**2

        aa = al(i,j,k)
        dd = sqrt((1.0 + aa)**2 - 4.0 * aa * ct)

! prepare slow mode coefficients
!
        fp = - 1.0 + aa - dd
        fr =   1.0 + aa - dd

! prepare slow mode unit vector
!
        ux = fp * kpx + fr * krx
        uy = fp * kpy + fr * kry
        uz = fp * kpz + fr * krz
        uu = sqrt(ux * ux + uy * uy + uz * uz)
        if (uu .eq. 0.0) uu = 1.0
        ux = ux / uu
        uy = uy / uu
        uz = uz / uu

! project velocity vector on the slow mode versor
!
        uu = ux * vx(i,j,k) + uy * vy(i,j,k) + uz * vz(i,j,k)

        sx(i,j,k) = uu * ux
        sy(i,j,k) = uu * uy
        sz(i,j,k) = uu * uz

! prepare fast mode coefficients
!
        fp = - 1.0 + aa + dd
        fr =   1.0 + aa + dd

! prepare fast mode unit vector
!
        ux = fp * kpx + fr * krx
        uy = fp * kpy + fr * kry
        uz = fp * kpz + fr * krz
        uu = sqrt(ux * ux + uy * uy + uz * uz)
        if (uu .eq. 0.0) uu = 1.0
        ux = ux / uu
        uy = uy / uu
        uz = uz / uu

! project velocity vector on the fast mode versor
!
        uu = ux * vx(i,j,k) + uy * vy(i,j,k) + uz * vz(i,j,k)

        fx(i,j,k) = uu * ux
        fy(i,j,k) = uu * uy
        fz(i,j,k) = uu * uz

      end do
!$omp end parallel do

      p = p + 1
      write(*,"('PROGRESS      : ',f6.2,' % done',a1,$)") pr*p, char(13)
    end do
  end do

  write(*,'(a)') ''

! deallocate local variables
!
  write( *, "('INFO          : ',a)" ) "deallocating local variables"

  if (allocated(kx) )  deallocate(kx)
  if (allocated(ky) )  deallocate(ky)
  if (allocated(kz) )  deallocate(kz)

  if (allocated(tmp))  deallocate(tmp)
  if (allocated(al) )  deallocate(al)
  if (allocated(vx) )  deallocate(vx)
  if (allocated(vy) )  deallocate(vy)
  if (allocated(vz) )  deallocate(vz)
  if (allocated(bx) )  deallocate(bx)
  if (allocated(by) )  deallocate(by)
  if (allocated(bz) )  deallocate(bz)

! calculate the inverse wavelet transform of the velocity mode
!
  write( *, "('TASK          : ',a)" ) "calculating the inverse DWT of"
  write( *, "('TASK          : ',a)" ) "  the Alfven mode"
  call idwt(ax, dims, wvlen)
  call idwt(ay, dims, wvlen)
  call idwt(az, dims, wvlen)
  write( *, "('TASK          : ',a)" ) "  the slow mode"
  call idwt(sx, dims, wvlen)
  call idwt(sy, dims, wvlen)
  call idwt(sz, dims, wvlen)
  write( *, "('TASK          : ',a)" ) "  the fast mode"
  call idwt(fx, dims, wvlen)
  call idwt(fy, dims, wvlen)
  call idwt(fz, dims, wvlen)

  timer = secnds(timer)
  write( *, "('COMPUTED      : ',a,1pe12.5,a)" ) 'computing done in ', timer, ' seconds'

! write Alfven mode components
!
  write( *, "('INFO          : ',a)" ) 'writing mode components of'
  write( *, "('INFO          : ',a)" ) "  the Alfven mode"

  write(fl, '("!",a)') trim(odir) // 'alfx.fits.gz'
  call fits_put_data(fl, ax, info)

  write(fl, '("!",a)') trim(odir) // 'alfy.fits.gz'
  call fits_put_data(fl, ay, info)

  write(fl, '("!",a)') trim(odir) // 'alfz.fits.gz'
  call fits_put_data(fl, az, info)

  write( *, "('INFO          : ',a)" ) "  the slow mode"

  write(fl, '("!",a)') trim(odir) // 'slox.fits.gz'
  call fits_put_data(fl, sx, info)

  write(fl, '("!",a)') trim(odir) // 'sloy.fits.gz'
  call fits_put_data(fl, sy, info)

  write(fl, '("!",a)') trim(odir) // 'sloz.fits.gz'
  call fits_put_data(fl, sz, info)

  write( *, "('INFO          : ',a)" ) "  the fast mode"

  write(fl, '("!",a)') trim(odir) // 'fasx.fits.gz'
  call fits_put_data(fl, fx, info)

  write(fl, '("!",a)') trim(odir) // 'fasy.fits.gz'
  call fits_put_data(fl, fy, info)

  write(fl, '("!",a)') trim(odir) // 'fasz.fits.gz'
  call fits_put_data(fl, fz, info)

! deallocate local variables
!
  write( *, "('INFO          : ',a)" ) "deallocating MHD wave components"

  if (allocated(ax) )  deallocate(ax)
  if (allocated(ay) )  deallocate(ay)
  if (allocated(az) )  deallocate(az)
  if (allocated(sx) )  deallocate(sx)
  if (allocated(sy) )  deallocate(sy)
  if (allocated(sz) )  deallocate(sz)
  if (allocated(fx) )  deallocate(fx)
  if (allocated(fy) )  deallocate(fy)
  if (allocated(fz) )  deallocate(fz)

end program
