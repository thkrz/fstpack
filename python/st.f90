subroutine cmsht2(l, m, r, x, y, h, coarse, fine, kernel_size)
  use hilbrt, only: cmsht2_ => cmsht2
  implicit none
  integer, intent(in) :: l, m, x, y
  real, intent(in) :: r(l, m)
  real, intent(out) :: h(4)
  real, intent(in), optional :: coarse, fine
  integer, intent(in), optional :: kernel_size

  call cmsht2_(r, x, y, coarse, fine, kernel_size, h)
end subroutine

subroutine imfreq(l, m, s, x, y, h)
  use fstpack, only: lfrqdm
  implicit none
  integer, intent(in)  :: l
  integer, intent(in)  :: m
  complex, intent(in)  :: s(l, l)
  integer, intent(in)  :: x
  integer, intent(in)  :: y
  complex, intent(out) :: h(m, m)

  h = lfrqdm(s, x, y)
end subroutine

subroutine idst2(n, s, h)
  use fstpack, only: cdst2b
  implicit none
  integer, intent(in)  :: n
  complex, intent(in)  :: s(n, n)
  complex, intent(out) :: h(n, n)

  h = s
  call cdst2b(h)
end subroutine

subroutine dst2(n, h, s)
  use fstpack, only: cdst2f
  implicit none
  integer, intent(in)  :: n
  complex, intent(in)  :: h(n, n)
  complex, intent(out) :: s(n, n)

  s = h
  call cdst2f(s)
end subroutine

subroutine ifst(n, s, h)
  use fstpack, only: cfst1b
  implicit none
  integer, intent(in)  :: n
  complex, intent(in)  :: s(n/2+1, n)
  complex, intent(out) :: h(n)

  h = cfst1b(s)
end subroutine

subroutine fst(n, h, s)
  use fstpack, only: cfst1f
  implicit none
  integer, intent(in)  :: n
  complex, intent(in)  :: h(n)
  complex, intent(out) :: s(n/2+1, n)

  s = cfst1f(h)
end subroutine
