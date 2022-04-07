module fstpack
  use fftpack
  use hilbrt
  implicit none
  private
  public lfrqdm
  public cdst2f
  public cfst1f
  public cfst1b

contains
  pure function lfrqdm(s, x, y) result(h)
    complex, intent(in) :: s(:, :)
    integer, intent(in) :: x, y
    complex, allocatable :: h(:, :)
    integer :: i, j, k, m, n, px, py
    real :: nx, ny, tx, ty, vx, vy, xs, ys

    k = size(s, 1)
    n = int(log(real(k))/log(2.)) - 1
    allocate(h(-n:n, -n:n))
    m = k / 2
    xs = 1. - x / float(k)
    ys = 1. - y / float(k)
    nx = 0
    ny = 0
    do px = 0, n
      tx = 2**(px - 1)
      vx = tx + nx
      nx = tx
      i = nint(vx - xs * nx)

      do py = 0, n
        ty = 2**(py - 1)
        vy = ty + ny
        ny = ty
        j = nint(vy - ys * ny)

        h(px, py) = s(m+i, m+j)
        h(px, -py) = s(m+i, m-j)
        h(-px, py) = s(m-i, m+j)
        h(-px, -py) = s(m-i, m-j)
      end do
    end do
  end function

  pure function cfst1b(s) result(h)
    complex, intent(in) :: s(0:, 0:)
    complex, allocatable :: h(:)
    integer :: i, err, l, n

    l = size(s, 1)
    allocate(h(0:l-1))
    h = 0
    do concurrent(n = 0:l-1, i = 0:l-1)
      h(n) = h(n) + s(n, i)
    end do

    call cht1b(h)
    call cfft1_('b', h, err)
    if(err /= 0) error stop
    h = h / l
  end function

  pure function cfst1f(h) result(s)
    complex, intent(in) :: h(0:)
    complex, allocatable :: s(:, :)
    complex, allocatable :: work(:)
    real, allocatable :: g(:)
    integer :: i, err, l, l2, n

    l = size(h)
    allocate(work(0:l-1))
    work = h
    call cfft1_('f', work, err)
    if(err /= 0) error stop
    call cht1f(work)

    allocate(g(0:l-1))
    allocate(s(0:l-1, 0:l-1))
    s(0, :) = sum(h) / l
    l2 = l / 2 + 1
    do concurrent(n = 1:l2-1)
      g(0) = gauss(n, 0)
      do concurrent(i = 1:l2-1)
        g(i) = gauss(n, i)
        g(l - i) = g(i)
      end do

      do concurrent(i = 0:l-1)
        s(n, i) = work(mod(n + i, l)) * g(i)
      end do

      call cfft1_('b', s(n, :), err)
      if(err /= 0) error stop
    end do
    deallocate(work)
  end function

  pure subroutine cdst2f(c)
    complex, intent(inout) :: c(0:, 0:)
    complex, allocatable :: work(:, :)
    integer :: err, k, m, n, nx, nx2, ny, ny2,&
               px, py, rx, ry, ty, ty2
    real :: sy, syx

    k = size(c, 1)
    if(k /= size(c, 2) .or. iand(k, k - 1) /= 0)&
      error stop
    allocate(work(0:k-1, 0:k-1))
    work = c
    call cfft2_('f', work, err)
    if(err /= 0) error stop
    m = k / 2
    c = 0
    c(0, 0) = work(0, 0)
    c(0, m) = work(0, m)
    c(m, 0) = work(m, 0)
    c(1, m) = work(1, m)
    c(m, 1) = work(m, 1)
    c(m, m) = work(m, m)

    n = int(log(real(k))/log(2.)) - 1
    do concurrent(py = 1:n)
      ny = 2**(py - 1)
      ny2 = ny * 2 - 1
      ry = floor(-ny / 2.)
      ty = k - ny
      ty2 = k - ny * 2 + 1
      sy = sqrt(real(ny))

      c(ny:ny2, 0) = shifft(work(ny:ny2, 0), ry) * sy
      c(0, ny:ny2) = shifft(work(0, ny:ny2), ry) * sy
      c(0, ty:ty2:-1) = shifft(work(0, ty2:ty), ry) * sy
      c(ny:ny2, m) = shifft(work(ny:ny2, m), ry) * sy
      c(m, ny:ny2) = shifft(work(m, ny:ny2), ry) * sy
      c(m, ty:ty2:-1) = shifft(work(m, ty2:ty), ry) * sy

      do concurrent(px = 1:n)
        nx = 2**(px - 1)
        nx2 = nx * 2 - 1
        rx = floor(-nx / 2.)
        syx = sqrt(real(ny*nx))

        c(nx:nx2, ny:ny2) = shifft2(work(nx:nx2, ny:ny2), [rx, ry]) * syx
        c(nx:nx2, ty:ty2:-1) = shifft2(work(nx:nx2, ty2:ty), [rx, ry]) * syx
      end do
    end do
    deallocate(work)

    c(m+1:, 1:m-1) = conjg(c(m-1:1:-1, k-1:m+1:-1))
    c(m+1:, m+1:) = conjg(c(m-1:1:-1, m-1:1:-1))
    c(m+1:, 0) = conjg(c(m-1:1:-1, 0))
    c(m+1:, m) = conjg(c(m-1:1:-1, m))
  end subroutine

  pure function gauss(n, m)
    real, parameter :: pi = 4. * atan(1.)
    integer, intent(in) :: n, m
    real :: gauss

    gauss = exp(-2. * pi**2 * m**2 / n**2)
  end function

  pure function shifft(array, shift) result(a)
    complex, intent(in) :: array(:)
    integer, intent(in) :: shift
    complex :: a(size(array))
    integer :: err

    a = cshift(array, shift)
    call cfft1_('b', a, err)
    if(err /= 0) error stop
  end function

  pure function shifft2(array, shift) result(a)
    complex, intent(in) :: array(:, :)
    integer, intent(in) :: shift(2)
    integer :: err
    complex :: a(size(array, 1), size(array, 2))

    a = cshift(array, shift(1), dim=1)
    a = cshift(a, shift(2), dim=2)
    call cfft2_('b', a, err)
    if(err /= 0) error stop
  end function
end module

