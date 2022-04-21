module fstpack
  use fftpack
  use hilbrt
  use mutl
  implicit none
  private
  public cdst2f
  public cfst1f
  public cfst1b
  public lfrqdm

contains
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
    c(m, 0) = work(m, 0)
    c(0, m) = work(0, m)
    c(m, 1) = work(m, 1)
    c(1, m) = work(1, m)
    c(m, m) = work(m, m)

    n = ilog2(k) - 1
    do py = 1, n
      ny  = 2**(py - 1)
      ny2 = ny * 2 - 1
      ry  = floor(-ny / 2.)
      ty  = k - ny
      ty2 = k - ny * 2 + 1
      sy  = sqrt(real(ny))

      c(0, ny:ny2)    = shifft(work(0, ny:ny2), ry) * sy
      c(ny:ny2, 0)    = shifft(work(ny:ny2, 0), ry) * sy
      c(ty:ty2:-1, 0) = shifft(work(ty2:ty, 0), ry) * sy
      c(m, ny:ny2)    = shifft(work(m, ny:ny2), ry) * sy
      c(ny:ny2, m)    = shifft(work(ny:ny2, m), ry) * sy
      c(ty:ty2:-1, m) = shifft(work(ty2:ty, m), ry) * sy

      do px = 1, n
        nx  = 2**(px - 1)
        nx2 = nx * 2 - 1
        rx  = floor(-nx / 2.)
        syx = sqrt(real(ny*nx))

        c(ny:ny2, nx:nx2)    = shifft2(work(ny:ny2, nx:nx2), [ry, rx]) * syx
        c(ty:ty2:-1, nx:nx2) = shifft2(work(ty2:ty, nx:nx2), [ry, rx]) * syx
      end do
    end do
    deallocate(work)

    c(1:m-1, m+1:) = conjg(c(k-1:m+1:-1, m-1:1:-1))
    c(m+1:, m+1:)  = conjg(c(m-1:1:-1, m-1:1:-1))
    c(0, m+1:)     = conjg(c(0, m-1:1:-1))
    c(m, m+1:)     = conjg(c(m, m-1:1:-1))
  end subroutine

  pure function cfst1b(s) result(h)
    complex, intent(in) :: s(:, :)
    complex, allocatable :: h(:)
    integer :: i, err, l, l2, n

    l2 = size(s, 1)
    l  = size(s, 2)
    allocate(h(l))
    h = 0
    do concurrent(n = 1:l2, i = 1:l)
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

    l2 = l / 2 + 1
    allocate(s(0:l2-1, 0:l-1))
    s(0, :) = sum(h) / l
    allocate(g(0:l-1))
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

  pure function lfrqdm(s, x, y) result(h)
    complex, intent(in) :: s(:, :)
    integer, intent(in) :: x, y
    complex, allocatable :: h(:, :)
    integer :: i, j, k, m, n, px, py
    real :: nx, ny, tx, ty, vx, vy, xs, ys

    k = size(s, 1)
    n = ilog2(k) - 1
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
      i  = nint(vx - xs * nx)
      do py = 0, n
        ty = 2**(py - 1)
        vy = ty + ny
        ny = ty
        j  = nint(vy - ys * ny)

        h(px, py)   = s(m+i, m+j)
        h(px, -py)  = s(m+i, m-j)
        h(-px, py)  = s(m-i, m+j)
        h(-px, -py) = s(m-i, m-j)
      end do
    end do
  end function

  pure function gauss(n, m)
    integer, intent(in) :: n, m
    real :: gauss

    gauss = exp(-2. * pi**2 * m**2 / n**2)
  end function

  pure function shifft(a, n) result(h)
    complex, intent(in) :: a(:)
    integer, intent(in) :: n
    complex :: h(size(a))
    integer :: err

    h = cshift(a, n)
    call cfft1_('b', h, err)
    if(err /= 0) error stop
  end function

  pure function shifft2(a, n) result(h)
    complex, intent(in) :: a(:, :)
    integer, intent(in) :: n(2)
    complex :: h(size(a, 1), size(a, 2))
    integer :: err

    h = cshift(a, n(1), dim=1)
    h = cshift(h, n(2), dim=2)
    call cfft2_('b', h, err)
    if(err /= 0) error stop
  end function
end module

