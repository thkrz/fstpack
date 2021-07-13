program trfst1
  use fstpack
  implicit none

  integer, parameter :: l = 128
  real :: cr, r(l), s(l)
  complex, allocatable :: w(:, :)
  integer :: c1, c2, err, n

  s = chirp()
  allocate(w(l/2+1, l))
  call system_clock(count_rate=cr)

  call system_clock(c1)
  call rfst1f(s, w, err)
  call system_clock(c2)
  print *, 'system_clock: ', (c2 - c1) / cr
  if(err /= 0) error stop
  call rfst1b(w, r, err)
  if(err /= 0) error stop
  do n = 1, l
    print '(F11.6,1X,F11.6,1X,F11.6)', s(n), r(n), abs(s(n) - r(n))
  end do

contains
  function chirp() result(s)
    real, parameter :: pi = 4.0 * atan(1.0)
    real :: s(0:l-1)
    integer :: t

    do concurrent(t = 0:63)
      s(t) = cos(2. * pi * t * 6. / 128.)
    end do
    do concurrent(t = 64:127)
      s(t) = cos(2. * pi * t * 25. / 128.)
    end do
    do concurrent(t = 20:29)
      s(t) = s(t) + .5 * cos(2. * pi * t * 52. / 128.)
    end do
  end function
end program
