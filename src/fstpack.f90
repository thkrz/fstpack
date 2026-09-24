module fstpack
  use fftpack
  use hilbrt
  use mutl
  implicit none
  private
  public cdst2b
  public cdst2f
  public cfst1f
  public cfst1b
  public lspec2

contains
  pure subroutine cdst2b(c)
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
    call diagi(work, c)

    m = k / 2
    n = ilog2(k) - 1
    do py = 1, n
      ny  = 2**(py - 1)
      ny2 = ny * 2 - 1
      ry  = -floor(-ny / 2.)
      ty  = k - ny
      ty2 = k - ny * 2 + 1
      sy  = sqrt(real(ny))

      c(0, ny:ny2) = shifft(work(0, ny:ny2) / sy, ry)
      c(ny:ny2, 0) = shifft(work(ny:ny2, 0) / sy, ry)
      c(ty2:ty, 0) = shifft(work(ty:ty2:-1, 0) / sy, ry)
      c(m, ny:ny2) = shifft(work(m, ny:ny2) / sy, ry)
      c(ny:ny2, m) = shifft(work(ny:ny2, m) / sy, ry)
      c(ty2:ty, m) = shifft(work(ty:ty2:-1, m) / sy, ry)

      do px = 1, n
        nx  = 2**(px - 1)
        nx2 = nx * 2 - 1
        rx  = -floor(-nx / 2.)
        syx = sqrt(real(ny*nx))

        c(ny:ny2, nx:nx2) = shifft2(work(ny:ny2, nx:nx2) / syx, [ry, rx])
        c(ty2:ty, nx:nx2) = shifft2(work(ty:ty2:-1, nx:nx2) / syx, [ry, rx])
      end do
    end do
    deallocate(work)
    call hsymm2(c)
    call cfft2_('b', c, err)
    if(err /= 0) error stop
  end subroutine

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
    call diagi(work, c)

    m = k / 2
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
    call hsymm2(c)
  end subroutine

  pure function cfst1b(s) result(h)
    complex, intent(in) :: s(:, :)
    complex, allocatable :: h(:)
    integer :: i, err, l, l2, n

    l2 = size(s, 1)
    l  = size(s, 2)
    allocate(h(l))
    h = 0
    do n = 1, l2
      h(n) = sum(s(n, :))
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
    do n = 1, l2-1
      g(0) = gauss(n, 0)
      do i = 1, l2-1
        g(i) = gauss(n, i)
        g(l - i) = g(i)
      end do

      do i = 0, l-1
        s(n, i) = work(mod(n + i, l)) * g(i)
      end do

      call cfft1_('b', s(n, :), err)
      if(err /= 0) error stop
    end do
    deallocate(work)
  end function

  pure function lspec2(s, x, y) result(h)
    complex, intent(in) :: s(0:, 0:)
    integer, intent(in) :: x, y
    complex, allocatable :: h(:, :)
    integer :: k, m, n, px, py, bx, by, tx, ty, ix, iy

    k = size(s, 1)
    if (k < 2 .or. size(s, 2) /= k .or. iand(k, k - 1) /= 0) error stop
    if (x < 0 .or. x >= k .or. y < 0 .or. y >= k) error stop

    m = k / 2
    n = ilog2(k) - 1
    allocate(h(2*n + 2, 2*n + 2))

    do py = -n, n + 1
      if (py == 0) then
        iy = 0
      else if (py == n + 1) then
        iy = m
      else
        by = 2**(abs(py) - 1)
        ty = y * by / k
        if (py > 0) then
          iy = by + ty
        else
          iy = k - by - ty
        end if
      end if

      do px = -n, n + 1
        if (px == 0) then
          ix = 0
        else if (px == n + 1) then
          ix = m
        else
          bx = 2**(abs(px) - 1)
          tx = x * bx / k
          if (px > 0) then
            ix = bx + tx
          else
            ix = k - bx - tx
          end if
        end if

        h(px + n + 1, py + n + 1) = s(ix, iy)
      end do
    end do
  end function lspec2

  pure subroutine diagi(a, b)
    complex, intent(in) :: a(0:, 0:)
    complex, intent(out) :: b(0:, 0:)
    integer :: m

    m = size(a, 1) / 2
    b = 0
    b(0, 0) = a(0, 0)
    b(m, 0) = a(m, 0)
    b(0, m) = a(0, m)
    b(m, 1) = a(m, 1)
    b(1, m) = a(1, m)
    b(m, m) = a(m, m)
  end subroutine

  pure function gauss(n, m)
    integer, intent(in) :: n, m
    real :: gauss

    gauss = exp(-2. * pi**2 * m**2 / n**2)
  end function

  pure subroutine hsymm2(c)
    complex, intent(inout) :: c(0:, 0:)
    integer :: l, m

    l = size(c, 1)
    m = l / 2
    c(1:m-1, m+1:) = conjg(c(l-1:m+1:-1, m-1:1:-1))
    c(m+1:, m+1:)  = conjg(c(m-1:1:-1, m-1:1:-1))
    c(0, m+1:)     = conjg(c(0, m-1:1:-1))
    c(m, m+1:)     = conjg(c(m, m-1:1:-1))
  end subroutine

  pure function shifft(a, n) result(h)
    complex, intent(in) :: a(:)
    integer, intent(in) :: n
    complex :: h(size(a))
    integer :: err

    if(n < 0) then
      h = cshift(a, n)
      call cfft1_('b', h, err)
    else
      h = a
      call cfft1_('f', h, err)
      h = cshift(h, n)
    end if
    if(err /= 0) error stop
  end function

  pure function shifft2(a, n) result(h)
    complex, intent(in) :: a(:, :)
    integer, intent(in) :: n(2)
    complex :: h(size(a, 1), size(a, 2))
    integer :: err

    if(n(1) < 0 .or. n(2) < 0) then
      h = cshift(a, n(1), dim=1)
      h = cshift(h, n(2), dim=2)
      call cfft2_('b', h, err)
    else
      h = a
      call cfft2_('f', h, err)
      h = cshift(h, n(2), dim=2)
      h = cshift(h, n(1), dim=1)
    end if
    if(err /= 0) error stop
  end function
end module

