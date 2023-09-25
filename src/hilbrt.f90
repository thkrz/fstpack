module hilbrt
  implicit none
  private
  public cht1f
  public cht1b

contains
  pure subroutine cht1b(h)
    complex, intent(inout) :: h(0:)
    integer :: i, l, m, n

    n = size(h)
    l = (n + 1) / 2
    m = n / 2 + 1
    h(1:l-1) = .5 * h(1:l-1)
    do concurrent(i = m:n-1)
      h(i) = conjg(h(n - i))
    end do
  end subroutine

  pure subroutine cht1f(h)
    complex, intent(inout) :: h(0:)
    integer :: l, m, n

    n = size(h)
    l = (n + 1) / 2
    m = n / 2 + 1
    h(1:l-1) = 2. * h(1:l-1)
    h(m:) = 0
  end subroutine
end module
