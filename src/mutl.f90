module mutl
  implicit none
  private
  public ilog2
  public pi

  real, parameter :: pi = 4. * atan(1.)

contains
  elemental function ilog2(x) result(y)
    integer, intent(in) :: x
    integer :: y

    y = int(log(real(x)) / log(2.))
  end function
end module
