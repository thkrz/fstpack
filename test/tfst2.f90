program tfst2
  use fstpack, only: rdst2f
  implicit none

  integer, parameter :: l = 512
  real, allocatable :: r(:, :)
  complex, allocatable, dimension(:, :) :: s, t
  integer :: err, fid, i, j

  allocate(r(l, l))
  open(newunit=fid, status='old', file='test/1.1.03.raw')
  do i = 1, l
    read(fid, *) (r(j, i), j = 1, l)
  end do
  close(fid)

  allocate(t(l, l))
  open(newunit=fid, status='old', file='test/1.1.03.dost')
  do i = 1, l
    read(fid, *) (t(j, i), j = 1, l)
  end do
  close(fid)

  allocate(s(l, l))
  call rdst2f(r, s, err)
  if(err /= 0) error stop

  print *, all(abs(s - t) < 1e-4)
  do i = 1, l
    print *, s(:, i)
  end do
end program
