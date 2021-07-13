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
    w(k2, 0 ) = h(k2, 0)
    w(0, k2) = h(0, k2)
    w(k2, 1 ) = h(k2, 1)
    w(1, k2) = h(1, k2)
    w(k2, k2) = h(k2, k2)

    n = int(log(real(k))/log(2.))
    do py = 1, n
      do px = 1, n
        nx = 2**(px - 1)
        nx2 = nx * 2
        cx2 = -nx / 2
        ny = 2**(py - 1)
        ny2 = ny * 2
        cy2 = -ny / 2
        nny = k - ny + 1
        nny2 = k - ny2 + 1
        sqny = sqrt(real(ny))
        sqnyx = sqrt(real(ny*nx))

        call shifft(w(0, ny:ny2), h(0, ny:ny2), cy2, sqny, err)
        call shifft(w(ny:ny2, 0), h(ny:ny2, 0), cy2, sqny, err)
        call shifft(w(k2, ny:ny2), h(k2, ny:ny2), cy2, sqny, err)
        call shifft(w(ny:ny2, k2), h(ny:ny2, k2), cy2, sqny, err)
        call shifft(w(nny2:nny, k2), h(nny2:nny, k2), cy2, sqny, err, rev=.true.)

        call shifft2(w(ny:ny2, nx:nx2), h(ny:ny2, nx:nx2), [cy2, cx2], sqnyx)
        call shifft2(w(nny2:nny, nx:nx2), h(nny2:nny, nx:nx2), [cy2, cx2], sqnyx, rev=.true.)
      end do
    end do

    w(1:k2, k2+1:) = conjg(w(k-1:k2+1:-1, k2:1:-1))
    w(k2+1:, k2+1:) = conjg(w(k2:1:-1, k2:1:-1))
    w(0, k2+1:) = conjg(w(0, k2:1:-1))
    w(k2, k2+1:) = conjg(w(k2, k2:1:-1))
  end subroutine

  pure function gauss(n, m)
    real, parameter :: pi = 4. * atan(1.)
    integer, intent(in) :: n, m
    real :: gauss

    gauss = exp(-2. * pi**2 * m**2 / n**2)
  end function

  pure subroutine shifft(p, q, n, s, err, rev)
    complex, intent(out) :: p(:)
    complex, intent(in) :: q(:)
    integer, intent(in) :: n
    integer, intent(out) :: err
    logical, intent(in), optional :: rev

    p = cshift(q, n)
    call cfft1(p, 'b', err)
    if(present(rev) && rev == .true.) then
      p = p(size(p):1:-1) * s
    else
      p = p * s
    end if
  end subroutine

  pure subroutine shifft2(p, q, n, s, err, rev)
    complex, intent(inout) :: p(:, :)
    complex, intent(in) :: q(:, :)
    integer, intent(in) :: n(2)
    integer, intent(out) :: err
    logical, intent(in), optional :: rev

    p = cshift(q, n(1), 1)
    p = cshift(p, n(2), 2)
    call cfft2(p, 'b', err)
    if(present(rev) && rev == .true.) then
      p = p(size(p):1:-1, :) * s
    else
      p = p * s
    end if
  end subroutine
end module

