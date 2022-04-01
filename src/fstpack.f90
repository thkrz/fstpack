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

  subroutine rdst2f(s, w, err)
    real, intent(in) :: s(0:, 0:)
    complex, intent(out) :: w(0:, 0:)
    integer, intent(out) :: err
    complex, allocatable :: h(:, :)
    integer :: cx2, cy2, k, k2, n, nny, nny2,&
               nx, nx2, ny, ny2, px, py
    real :: sqny, sqnyx

    k = size(s, 1)
    if(k /= size(s, 2) .or. iand(k, k - 1) /= 0) then
      err = 1
      return
    end if
    allocate(h(0:k-1, 0:k-1))
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
      w(0, nny2:nny) = rev(shifft(h(0, nny2:nny), cy2)) * sqny
      w(ny:ny2, k2) = shifft(h(ny:ny2, k2), cy2) * sqny
      w(k2, ny:ny2) = shifft(h(k2, ny:ny2), cy2) * sqny
      w(k2, nny2:nny) = rev(shifft(h(k2, nny2:nny), cy2)) * sqny

      do px = 1, n
        nx = 2**(px - 1)
        nx2 = nx * 2 - 1
        cx2 = floor(-nx / 2.)
        sqnyx = sqrt(real(ny*nx))

        w(nx:nx2, ny:ny2) = shifft2(h(nx:nx2, ny:ny2), [cx2, cy2]) * sqnyx
        w(nx:nx2, nny2:nny) = rev2(shifft2(h(nx:nx2, nny2:nny), [cx2, cy2]), axis=2) * sqnyx
      end do
    end do
    return

    w(k2+1:, 1:k2-1) = conjg(rev2(w(1:k2-1, k2+1:)))
    w(k2+1:, k2+1:) = conjg(rev2(w(1:k2-1, 1:k2-1)))
    w(k2+1:, 0) = conjg(rev(w(1:k2-1, 0)))
    w(k2+1:, k2) = conjg(rev(w(k2+1:, k2)))
  end subroutine

  pure function gauss(n, m)
    real, parameter :: pi = 4. * atan(1.)
    integer, intent(in) :: n, m
    real :: gauss

    gauss = exp(-2. * pi**2 * m**2 / n**2)
  end function

  pure function rev(arr)
    complex, intent(in) :: arr(:)
    complex :: rev(size(arr))

    rev = arr(size(arr):1:-1)
  end function

  pure function rev2(arr, axis)
    complex, intent(in) :: arr(:, :)
    integer, optional, intent(in) :: axis
    complex :: rev2(size(arr, 1), size(arr, 2))
    integer :: a

    a = 0
    if(present(axis)) a = axis
    select case(a)
    case(0)
      rev2 = arr(size(arr, 1):1:-1, size(arr, 2):1:-1)
    case(1)
      rev2 = arr(size(arr, 1):1:-1, :)
    case(2)
      rev2 = arr(:, size(arr, 2):1:-1)
    end select
  end function

  function shifft(q, n) result(p)
    complex, intent(in) :: q(:)
    integer, intent(in) :: n
    complex :: p(size(q))
    integer :: err

    p = cshift(q, n)
    call cfft1(p, 'b', err)
    if(err /= 0) error stop
  end function

  function shifft2(q, n) result(p)
    complex, intent(in) :: q(:, :)
    integer, intent(in) :: n(2)
    integer :: err
    complex :: p(size(q, 1), size(q, 2))

    p(1, :) = cshift(q(2, :), n(1))
    p(2, :) = cshift(q(1, :), n(2))
    call cfft2(p, 'b', err)
    if(err /= 0) error stop
  end function
end module

