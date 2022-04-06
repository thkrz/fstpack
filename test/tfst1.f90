program tfst1
  use fstpack, only: cfst1f, cfst1b
  implicit none

  integer, parameter :: l = 128
  real :: cr, eps
  complex :: h(l), t(l)
  complex, allocatable :: s(:, :)
  integer :: c1, c2

  h = cmplx(chirp())
  call system_clock(count_rate=cr)
  call system_clock(c1)
  s = cfst1f(h)
  call system_clock(c2)
  print *, 'system_clock: ', (c2 - c1) / cr
  t = cfst1b(s)
  eps = sqrt(epsilon(eps))
  print *, all(abs(h - t) < eps)

contains
  pure function chirp() result(h)
    real, parameter :: pi = 4.0 * atan(1.0)
    real :: h(0:l-1)
    integer :: t

    do concurrent(t = 0:63)
      h(t) = cos(2. * pi * t * 6. / 128.)
    end do
    do concurrent(t = 64:127)
      h(t) = cos(2. * pi * t * 25. / 128.)
    end do
    do concurrent(t = 20:29)
      h(t) = h(t) + .5 * cos(2. * pi * t * 52. / 128.)
    end do
  end function
end program
