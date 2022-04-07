subroutine freqdomain(l, m, s, x, y, h)
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
