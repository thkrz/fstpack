program tfst2
  use fstpack, only: cdst2f
  implicit none

  integer, parameter :: n = 64
  real :: cr, eps
  complex, dimension(n, n) :: c, t
  integer :: c1, c2

  c = cmplx(chirp())
  t = c
  call system_clock(count_rate=cr)
  call system_clock(c1)
  call cdst2f(c)
  call system_clock(c2)
  print *, 'system_clock: ', (c2 - c1) / cr

  ! call cdst2b(c)
  eps = sqrt(epsilon(eps))
  print *, all(abs(c - t) < eps)

contains
  pure function chirp() result(h)
    real, parameter :: pi = 4.0 * atan(1.0)
    real :: f, h(0:n-1, 0:n-1), xx, yy
    integer :: n2, x, y

    n2 = floor(n / 2.)
    do concurrent(x = 0:n-1, y = 0:n-1)
      f = 10. * cos(2. * pi * (.15 * x) * x / n)
      xx = real(x - n2)**2
      yy = real(y - n2)**2
      h(x, y) = f * exp(-(xx / n + yy / n2))
    end do
  end function
end program
