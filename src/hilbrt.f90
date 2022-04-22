module hilbrt
  implicit none
  private
  public cmsht2
  public cht1f
  public cht1b

contains
  pure subroutine cmsht2(r, x, y, coarse, fine, ksize, h)
    real, intent(in) :: r(0:, 0:), coarse, fine
    integer, intent(in) :: x, y, ksize
    real, intent(out) :: h(4)
    real :: c, d, f, pf, pc, rp, rx, ry, rz, uvw, u, v, w
    integer :: cx, cy, l, m, i, j

    l = size(r, 1) - 1
    m = size(r, 2) - 1
    rp = 0
    rx = 0
    ry = 0
    rz = 0
    do concurrent(cx = -ksize:ksize, cy = -ksize:ksize)
      d = cx**2 + cy**2 + 1
      u = cx / d
      v = cy / d
      w = (d - 1) / d
      uvw = u**2 + v**2 + w**2
      pf = sqrt(fine**2 + uvw)
      pc = sqrt(coarse**2 + uvw)
      i = x + cx
      j = y + cy
      if(i * j < 0 .or. i > l .or. j > m) then
        f = 0
      else
        f = r(i, j)
      end if
      c = f * (pf - pc)
      rp = rp + f * (fine * pf - coarse * pc)
      rx = rx + u * c
      ry = ry + v * c
      rz = rz + w * c
    end do
    h(1) = hypot(rx, ry) / rz
    h(2) = atan2(ry, rx)
    h(3) = atan2(norm2([rx, ry, rz]), rp)
    h(4) = rp**2 + rx**2 + ry**2 + rz**2
  end subroutine

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
