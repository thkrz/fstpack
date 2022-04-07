program tfst2
  use fstpack, only: cdst2f, cdst2b
  implicit none

  integer, parameter :: n = 8
  real :: cr, eps
  complex, dimension(n, n) :: c, t
  integer :: c1, c2

  c = cmplx(chirp())
  print *, c
  print *, "--"
  t = c
  call system_clock(count_rate=cr)
  call system_clock(c1)
  call cdst2f(c)
  call system_clock(c2)
  print *, 'system_clock: ', (c2 - c1) / cr

  call cdst2b(c)
  print *, c
  eps = sqrt(epsilon(eps))
  print *, all(abs(c - t) < eps)

contains
  pure function chirp() result(h)
    real, parameter :: pi = 4.0 * atan(1.0)
    real :: f(0:n-1, 0:n-1), h(0:n-1, 0:n-1)
    integer :: n2, x, y, xx, yy

    do concurrent(x = 0:n-1)
      f(x, :) = 10. * cos(2. * pi * (.15 * x) * x / n)
    end do
    n2 = floor(n / 2.)
    do concurrent(x = 0:n-1, y = 0:n-1)
      xx = (x - n2)**2
      yy = (y - n2)**2
      h(x, y) = exp(-real(xx + yy) / 2.)
    end do
    h = matmul(f, h)
  end function
end program
