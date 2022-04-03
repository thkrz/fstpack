module fstpack
  use fftpack
  use hilbrt
  implicit none
  private
  public rdst2f
  public rfst1f
  public rfst1b

contains
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
    integer :: cx2, cy2, k, k2, n, nny, nny2,&
               nx, nx2, ny, ny2, px, py
    real :: sqny, sqnyx

    k = size(s, 1)
    if(k /= size(s, 2) .or. iand(k, k - 1) /= 0) then
      err = 1
      return
    end if
    h = cmplx(s)
    call cfft2(h, 'f', err)
    if(err /= 0) return
    k2 = k / 2
    w = 0
    w(0, 0) = h(0, 0)
    w(0, k2) = h(0, k2)
    w(k2, 0 ) = h(k2, 0)
    w(1, k2) = h(1, k2)
    w(k2, 1 ) = h(k2, 1)
    w(k2, k2) = h(k2, k2)

    n = int(log(real(k))/log(2.)) - 1
    do py = 1, n
      ny = 2**(py - 1)
      ny2 = ny * 2 - 1
      cy2 = floor(-ny / 2.)
      nny = k - ny
      nny2 = k - ny * 2 + 1
      sqny = sqrt(real(ny))

      w(ny:ny2, 0) = shifft(h(ny:ny2, 0), cy2) * sqny
      w(0, ny:ny2) = shifft(h(0, ny:ny2), cy2) * sqny
      w(0, nny:nny2:-1) = shifft(h(0, nny2:nny), cy2) * sqny
      w(ny:ny2, k2) = shifft(h(ny:ny2, k2), cy2) * sqny
      w(k2, ny:ny2) = shifft(h(k2, ny:ny2), cy2) * sqny
      w(k2, nny:nny2:-1) = shifft(h(k2, nny2:nny), cy2) * sqny

      do px = 1, n
        nx = 2**(px - 1)
        nx2 = nx * 2 - 1
        cx2 = floor(-nx / 2.)
        sqnyx = sqrt(real(ny*nx))

        w(nx:nx2, ny:ny2) = shifft2(h(nx:nx2, ny:ny2), [cx2, cy2]) * sqnyx
        w(nx:nx2, nny:nny2:-1) = shifft2(h(nx:nx2, nny2:nny), [cx2, cy2]) * sqnyx
      end do
    end do

    w(k2+1:, 1:k2-1) = conjg(w(k2-1:1:-1, k-1:k2+1:-1))
    w(k2+1:, k2+1:) = conjg(w(k2-1:1:-1, k2-1:1:-1))
    w(k2+1:, 0) = conjg(w(k2-1:1:-1, 0))
    w(k2+1:, k2) = conjg(w(k2-1:1:-1, k2))
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

