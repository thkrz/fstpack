program tfst2
  use fstpack, only: rdst2f
  implicit none

  real, allocatable :: r(:, :)
  complex, allocatable :: s(:, :)
  integer :: n, err, fid, i, j
  character(256) :: name
  namelist /tfst2ex/ n, name

  read(nml=tfst2ex
  allocate(r(0:l, 0:l))
  open(newunit=fid, status='old', file='test/1.1.04.raw')
  do i = 0, l
    read(fid, *) (r(j, i), j = 0, l)
  end do
  close(fid)

  allocate(t(0:l, 0:l))
  open(newunit=fid, status='old', file='test/1.1.04.dost')
  do i = 0, l
    read(fid, *) (t(j, i), j = 0, l)
  end do
  close(fid)

  allocate(s(0:l, 0:l))
  call rdst2f(r, s, err)

  print *, all(abs(s - t) < 1e-6)
end program
