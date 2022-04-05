module fstpack
  use fftpack
  use hilbrt
  implicit none
  private
  public freqdm
  public rdst2f
  public rfst1f
  public rfst1b

contains
  subroutine freqdm(x, y, s, vm)
    integer, intent(in) :: x, y
    complex, intent(in) :: s(:, :)
    complex, intent(out), allocatable :: vm(:, :)
    integer :: i, j, k, m, n, px, py
    real :: nx, ny, tx, ty, vx, vy, xs, ys

    k = size(s, 1)
    n = int(log(real(k))/log(2.)) - 1
    allocate(vm(-n:n, -n:n))
    vm(0, 0) = 0
    m = k / 2
    xs = x / float(k)
    ys = y / float(k)
    nx = 0
    ny = 0
    do px = 0, n
      tx = 2**(px - 1)
      vx = tx + nx
      nx = tx
      i = nint(vx - xs * nx)
      do py = 0, n
        if(px == 0 .and. py == 0) continue
        ty = 2**(py - 1)
        vy = ty + ny
        ny = ty
        j = nint(vy - ys * ny)
        vm(px, py) = s(m+i, m+j)
        vm(-px, -py) = s(m-i, m-j)
        vm(px, -py) = s(m+i, m-j)
        vm(-px, py) = s(m-i, m+j)
      end do
    end do
  end subroutine

  subroutine rfst1b(w, s, err)
    complex, intent(in) :: w(0:, 0:)
    real, intent(out) :: s(0:)
    integer, intent(out) :: err
    complex, allocatable :: h(:)
    integer :: i, l, l2, n

    l = size(s)
    allocate(h(0:l-1))
    h = 0
    l2 = size(w, 1)
    do concurrent(n = 0:l2-1, i = 0:l-1)
      h(n) = h(n) + w(n, i)
    end do

    call cht1b(h)
    call cfft1(h, 'b', err)
    if(err /= 0) return
    s = real(h) / l
  end subroutine

  subroutine rfst1f(s, w, err)
    real, intent(in) :: s(0:)
    complex, intent(out) :: w(0:, 0:)
    integer, intent(out) :: err
    complex, allocatable :: h(:)
    real, allocatable :: g(:)
    integer :: i, l, l2, n

    l = size(s)
    allocate(h(0:l-1))
    h = cmplx(s)
    call cfft1(h, 'f', err)
    if(err /= 0) return
    call cht1f(h)

    allocate(g(0:l-1))
    w(0, :) = cmplx(sum(s)/l)
    l2 = l / 2 + 1
    do n = 1, l2-1
      g(0) = gauss(n, 0)
      do concurrent(i = 1:l2-1)
        g(i) = gauss(n, i)
        g(l - i) = g(i)
      end do

      do concurrent(i = 0:l-1)
        w(n, i) = h(mod(n + i, l)) * g(i)
      end do

      call cfft1(w(n, :), 'b', err)
      if(err /= 0) return
    end do
  end subroutine

  subroutine rdst2b(w, s, err)
    complex, intent(in) :: w(0:, 0:)
    real, intent(out) :: s(0:, 0:)
    integer, intent(out) :: err
  end subroutine

  subroutine rdst2f(s, w, err)
    real, intent(in) :: s(0:, 0:)
    complex, intent(out) :: w(0:, 0:)
    integer, intent(out) :: err
    complex :: h(0:size(s, 1)-1, 0:size(s, 2)-1)
    integer :: rx, ry, k, m, n, nx, nx2,&
               ny, ny2, px, py, ty, ty2
    real :: sy, syx

    k = size(s, 1)
    if(k /= size(s, 2) .or. iand(k, k - 1) /= 0) then
      err = 1
      return
    end if
    h = cmplx(s)
    call cfft2(h, 'f', err)
    if(err /= 0) return
    m = k / 2
    w = 0
    w(0, 0) = h(0, 0)
    w(0, m) = h(0, m)
    w(m, 0) = h(m, 0)
    w(1, m) = h(1, m)
    w(m, 1) = h(m, 1)
    w(m, m) = h(m, m)

    n = int(log(real(k))/log(2.)) - 1
    do py = 1, n
      ny = 2**(py - 1)
      ny2 = ny * 2 - 1
      ry = floor(-ny / 2.)
      ty = k - ny
      ty2 = k - ny * 2 + 1
      sy = sqrt(real(ny))

      w(ny:ny2, 0) = shifft(h(ny:ny2, 0), ry) * sy
      w(0, ny:ny2) = shifft(h(0, ny:ny2), ry) * sy
      w(0, ty:ty2:-1) = shifft(h(0, ty2:ty), ry) * sy
      w(ny:ny2, m) = shifft(h(ny:ny2, m), ry) * sy
      w(m, ny:ny2) = shifft(h(m, ny:ny2), ry) * sy
      w(m, ty:ty2:-1) = shifft(h(m, ty2:ty), ry) * sy

      do px = 1, n
        nx = 2**(px - 1)
        nx2 = nx * 2 - 1
        rx = floor(-nx / 2.)
        syx = sqrt(real(ny*nx))

        w(nx:nx2, ny:ny2) = shifft2(h(nx:nx2, ny:ny2), [rx, ry]) * syx
        w(nx:nx2, ty:ty2:-1) = shifft2(h(nx:nx2, ty2:ty), [rx, ry]) * syx
      end do
    end do

    w(m+1:, 1:m-1) = conjg(w(m-1:1:-1, k-1:m+1:-1))
    w(m+1:, m+1:) = conjg(w(m-1:1:-1, m-1:1:-1))
    w(m+1:, 0) = conjg(w(m-1:1:-1, 0))
    w(m+1:, m) = conjg(w(m-1:1:-1, m))
  end subroutine

  pure function gauss(n, m)
    real, parameter :: pi = 4. * atan(1.)
    integer, intent(in) :: n, m
    real :: gauss

    gauss = exp(-2. * pi**2 * m**2 / n**2)
  end function

  function shifft(array, shift) result(a)
    complex, intent(in) :: array(:)
    integer, intent(in) :: shift
    complex :: a(size(array))
    integer :: err

    a = cshift(array, shift)
    call cfft1(a, 'b', err)
    if(err /= 0) error stop
  end function

  function shifft2(array, shift) result(a)
    complex, intent(in) :: array(:, :)
    integer, intent(in) :: shift(2)
    integer :: err
    complex :: a(size(array, 1), size(array, 2))

    a = cshift(array, shift(1), dim=1)
    a = cshift(a, shift(2), dim=2)
    call cfft2(a, 'b', err)
    if(err /= 0) error stop
  end function
end module

