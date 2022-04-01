program tfst2
  use fstpack, only: rdst2f
  implicit none

  integer, parameter :: l = 7
  real, allocatable :: r(:, :)
  complex, allocatable, dimension(:, :) :: s, t
  integer :: err, fid, i, j

  allocate(r(0:l, 0:l))
  open(newunit=fid, status='old', file='test/1.1.04.raw')
  do i = 0, l
    read(fid, *) (r(j, i), j = 0, l)
  end do
  close(fid)

  allocate(t(0:l, 0:l))
  ! open(newunit=fid, status='old', file='test/1.1.04.dost')
  ! do i = 0, l
  !   read(fid, *) (t(j, i), j = 0, l)
  ! end do
  ! close(fid)

  allocate(s(0:l, 0:l))
  call rdst2f(r, s, err)

  do j = 0, l
    print *, (s(i, j), i = 0, l)
  end do
  ! print *, all(abs(s - t) < 1e-6)
end program
