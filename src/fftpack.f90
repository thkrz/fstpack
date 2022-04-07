module fftpack
  implicit none
  private
  public cfft1_
  public cfft2_

interface
  pure subroutine cfft1i(n, wsave, lensav, ier)
    integer, intent(in) :: n, lensav
    integer, intent(out) :: ier
    real, intent(out) :: wsave(lensav)
  end subroutine

  pure subroutine cfft1b(n, inc, c, lenc, wsave,&
      lensav, work, lenwrk, ier)
    integer, intent(in) :: n, inc, lenc, lensav, lenwrk
    real, intent(in) :: wsave(lensav), work(lenwrk)
    complex, intent(inout) :: c(lenc)
    integer, intent(out) :: ier
  end subroutine

  pure subroutine cfft1f(n, inc, c, lenc, wsave,&
      lensav, work, lenwrk, ier)
    integer, intent(in) :: n, inc, lenc, lensav, lenwrk
    real, intent(in) :: wsave(lensav), work(lenwrk)
    complex, intent(inout) :: c(lenc)
    integer, intent(out) :: ier
  end subroutine

  pure subroutine cfft2i(l, m, wsave, lensav, ier)
    integer, intent(in) :: l, m, lensav
    integer, intent(out) :: ier
    real, intent(out) :: wsave(lensav)
  end subroutine

  pure subroutine cfft2b(ldim, l, m, c, wsave,&
      lensav, work, lenwrk, ier)
    integer, intent(in) :: ldim, l, m, lensav, lenwrk
    real, intent(in) :: wsave(lensav), work(lenwrk)
    complex, intent(inout) :: c(ldim, m)
    integer, intent(out) :: ier
  end subroutine

  pure subroutine cfft2f(ldim, l, m, c, wsave,&
      lensav, work, lenwrk, ier)
    integer, intent(in) :: ldim, l, m, lensav, lenwrk
    real, intent(in) :: wsave(lensav), work(lenwrk)
    complex, intent(inout) :: c(ldim, m)
    integer, intent(out) :: ier
  end subroutine
end interface

contains
  pure subroutine cfft1_(dir, c, err)
    character(1), intent(in) :: dir
    complex, intent(inout) :: c(:)
    integer, intent(out) :: err
    integer :: l, lensav, lenwrk
    real, allocatable :: wsave(:), work(:)

    l = size(c)
    lensav = 2 * l + int(log(real(l)) / log(2.)) + 4
    allocate(wsave(lensav))
    call cfft1i(l, wsave, lensav, err)
    if(err /= 0) return
    lenwrk = 2 * l
    allocate(work(lenwrk))
    if(dir == 'b') then
      call cfft1b(l, 1, c, l, wsave, lensav, work, lenwrk, err)
    else if(dir == 'f') then
      call cfft1f(l, 1, c, l, wsave, lensav, work, lenwrk, err)
    else
      err = 1
    end if
    deallocate(work)
    deallocate(wsave)
  end subroutine

  pure subroutine cfft2_(dir, c, err)
    character(1), intent(in) :: dir
    complex, intent(inout) :: c(:, :)
    integer, intent(out) :: err
    integer :: l, lensav, lenwrk, m
    real, allocatable :: wsave(:), work(:)

    l = size(c, 1)
    m = size(c, 2)
    lensav = 2 * (l + m)&
      + int(log(real(l)) / log(2.))&
      + int(log(real(m)) / log(2.))&
      + 8
    allocate(wsave(lensav))
    call cfft2i(l, m, wsave, lensav, err)
    if(err /= 0) return
    lenwrk = 2 * l * m
    allocate(work(lenwrk))
    if(dir == 'b') then
      call cfft2b(l, l, m, c, wsave, lensav, work, lenwrk, err)
    else if(dir == 'f') then
      call cfft2f(l, l, m, c, wsave, lensav, work, lenwrk, err)
    else
      err = 1
    end if
    deallocate(work)
    deallocate(wsave)
  end subroutine
end module
