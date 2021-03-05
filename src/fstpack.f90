module fstpack
  use fftpack
  use hilbrt
  implicit none
  private
  public rfst1f
  public rfst1b

contains
  subroutine rfst1b(w, s, err)
    complex, intent(in)  :: w(0:, 0:)
    real,    intent(out) :: s(0:)
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
    real,    intent(in)  :: s(0:)
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

  pure function gauss(n, m)
    real, parameter :: pi = 4.0 * atan(1.0)
    integer, intent(in) :: n, m
    real :: gauss

    gauss = exp(-2.0 * pi**2 * m**2 / n**2)
  end function
end module
