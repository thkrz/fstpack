program trfst1
  use fftpack
  implicit none

  integer, parameter :: l = 128
  real :: s(l)
  complex :: h(l), t(l)
  integer :: err

  s = chirp()
  h = cmplx(s)
  t = h
  call cfft1(h, 'f', err)
  if(err /= 0) error stop
  call cfft1(h, 'b', err)
  if(err /= 0) error stop
  print *, all(abs(t-h) < 1e-6)

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
