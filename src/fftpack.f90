module fftpack
  implicit none
  private
  public cfft1
  public cfft2

  real, allocatable :: wsave(:)
  real, allocatable :: work(:)
  integer :: lensav, lenwrk, n = 0

contains
  subroutine cfft1(s, norm, err)
    complex, intent(inout) :: s(:)
    character(1), intent(in) :: norm
    integer, intent(out) :: err
    integer, save :: l = 0

    if(n /= 1 .or. l /= size(s)) then
      n = 1
      if(allocated(wsave)) deallocate(wsave)
      l = size(s)
      lensav = 2 * l + int(log(real(l)) / log(2.)) + 4
      allocate(wsave(lensav))
      call cfft1i(l, wsave, lensav, err)
      if(err /= 0) return
    end if
    lenwrk = 2 * l
    allocate(work(lenwrk))
    if(norm == 'b') then
      call cfft1b(l, 1, s, l, wsave, lensav, work, lenwrk, err)
    else if(norm == 'f') then
      call cfft1f(l, 1, s, l, wsave, lensav, work, lenwrk, err)
    else
      err = 1
    end if
    deallocate(work)
    s = s * l
  end subroutine

  subroutine cfft2(s, norm, err)
    complex, intent(inout) :: s(:, :)
    character(1), intent(in) :: norm
    integer, intent(out) :: err
    integer, save :: l = 0, m = 0

    if(n /= 2 .or. l /= size(s, 1) .or. m /= size(s, 2)) then
      n = 2
      if(allocated(wsave)) deallocate(wsave)
      l = size(s, 1)
      m = size(s, 2)
      lensav = 2 * (l + m)&
        + int(log(real(l)) / log(2.))&
        + int(log(real(m)) / log(2.))&
        + 8
      allocate(wsave(lensav))
      call cfft2i(l, m, wsave, lensav, err)
      if(err /= 0) return
    end if
    lenwrk = 2 * l * m
    allocate(work(lenwrk))
    if(norm == 'b') then
      call cfft2b(l, l, m, s, wsave, lensav, work, lenwrk, err)
    else if(norm == 'f') then
      call cfft2f(l, l, m, s, wsave, lensav, work, lenwrk, err)
    else
      err = 1
    end if
    deallocate(work)
    s = s * l * m
  end subroutine
end module
