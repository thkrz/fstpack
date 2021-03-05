subroutine fst(n, s, w)
  use fstpack, only: rfst1f
  implicit none
  integer, intent(in)  :: n
  real,    intent(in)  :: s(n)
  complex, intent(out) :: w(n/2+1, n)
  integer :: err

  call rfst1f(s, w, err)
  if(err /= 0) error stop
end subroutine

subroutine ifst(n, w, s)
  use fstpack, only: rfst1b
  implicit none
  integer, intent(in)  :: n
  complex, intent(in)  :: w(n/2+1, n)
  real,    intent(out) :: s(n)
  integer :: err

  call rfst1b(w, s, err)
  if(err /= 0) error stop
end subroutine

! subroutine fst2(m, n, s, w)
!   use fstpack, only: fst2_ => fst2
!   implicit none
!   integer, intent(in)  :: m, n
!   real,    intent(in)  :: s(m, n)
!   complex, intent(out) :: w(m/2+1, n/2+1, m, n)

!   call fst2_(s, w)
! end subroutine

! subroutine ifst2(m, n, w, s)
!   use fstpack, only: ifst2_ => ifst2
!   implicit none
!   integer, intent(in)  :: m, n
!   complex, intent(in)  :: w(m/2+1, n/2+1, m, n)
!   real,    intent(out) :: s(m, n)

!   call ifst2_(s, w)
! end subroutine
