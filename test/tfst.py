import numpy as np
import pyst
import unittest

# from numpy.fft import fft, fft2, ifft, ifftn
def fft(c):
    return np.fft.fft(c) / c.size


def fft2(c):
    return np.fft.fft2(c) / c.size


def ifft(c):
    return np.fft.ifft(c) * c.size


def ifftn(c):
    return np.fft.ifftn(c) * c.size


def chirp(order=8, dtype=np.float32):
    n = 2**order
    h = np.zeros((n, n), dtype=dtype)
    n2 = n // 2
    for i in range(n):
        f = 10.0 * np.cos(2.0 * np.pi * (0.15 * i) * i / 64.0)
        for j in range(n):
            x = i - n2
            y = j - n2
            h[i, j] = f * np.exp(-0.5 * (x**2 / n + y**2 / n2))
    return h


def shift(a, n):
    if len(a.shape) == 1:
        b = np.concatenate((a[n:], a[:n]))
    elif len(a.shape) == 2:
        c = np.concatenate((a[n[0] :, :], a[: n[0], :]))
        b = np.concatenate((c[:, n[1] :], c[:, : n[1]]), 1)
    return b


def idst2(S):
    N = len(S)
    N2 = N // 2
    IM = np.zeros((N, N), dtype=np.csingle)

    n = int(np.log2(N))
    for py in range(1, n):
        ny = 2 ** (py - 1)
        ny2 = ny * 2
        ty = N - ny2 + 1
        ty2 = N - ny + 1
        ry = ny // 2
        sy = np.sqrt(ny)

        IM[0, ny:ny2] = shift(fft(S[0, ny:ny2] / sy), ry)
        IM[ny:ny2, 0] = shift(fft(S[ny:ny2, 0] / sy), ry)
        IM[ty:ty2, 0] = shift(fft(S[ty:ty2, 0][::-1] / sy), ry)
        IM[N2, ny:ny2] = shift(fft(S[N2, ny:ny2] / sy), ry)
        IM[ny:ny2, N2] = shift(fft(S[ny:ny2, N2] / sy), ry)
        IM[ty:ty2, N2] = shift(fft(S[ty:ty2, N2][::-1] / sy), ry)

        for px in range(1, n):
            nx = 2 ** (px - 1)
            nx2 = nx * 2
            rx = nx // 2
            syx = np.sqrt(ny * nx)

            IM[ny:ny2, nx:nx2] = shift(fft2(S[ny:ny2, nx:nx2] / syx), (ry, rx))
            IM[ty:ty2, nx:nx2] = shift(
                fft2(S[ty:ty2, nx:nx2][::-1] / syx),
                (ry, rx),
            )

    IM[0, 0] = S[0, 0]
    IM[N2, 0] = S[N2, 0]
    IM[0, N2] = S[0, N2]
    IM[N2, 1] = S[N2, 1]
    IM[1, N2] = S[1, N2]
    IM[N2, N2] = S[N2, N2]

    IM[1:N2, N2 + 1 :] = np.conj(IM[N2 + 1 :, 1:N2][::-1, ::-1])
    IM[N2 + 1 :, N2 + 1 :] = np.conj(IM[1:N2, 1:N2][::-1, ::-1])
    IM[0, N2 + 1 :] = np.conj(IM[0, 1:N2][::-1])
    IM[N2, N2 + 1 :] = np.conj(IM[N2, 1:N2][::-1])

    im = ifftn(IM)
    return im.real


def dst2(im):
    N = len(im)
    IM = fft2(im)

    S = np.zeros((N, N), dtype=np.csingle)
    S[0, 0] = IM[0, 0]
    S[N // 2, 0] = IM[N // 2, 0]
    S[0, N // 2] = IM[0, N // 2]
    S[N // 2, 1] = IM[N // 2, 1]
    S[1, N // 2] = IM[1, N // 2]
    S[N // 2, N // 2] = IM[N // 2, N // 2]
    n = int(np.log2(N))

    for py in range(1, n):
        ny = 2 ** (py - 1)
        S[0, ny : ny * 2] = ifft(shift(IM[0, ny : ny * 2], -ny // 2)) * np.sqrt(ny)
        S[ny : ny * 2, 0] = ifft(shift(IM[ny : ny * 2, 0], -ny // 2)) * np.sqrt(ny)
        S[N - ny * 2 + 1 : N - ny + 1, 0] = ifft(
            shift(IM[N - ny * 2 + 1 : N - ny + 1, 0], -ny // 2)
        )[::-1] * np.sqrt(ny)
        S[N // 2, ny : ny * 2] = ifft(
            shift(IM[N // 2, ny : ny * 2], -ny // 2)
        ) * np.sqrt(ny)
        S[ny : ny * 2, N // 2] = ifft(
            shift(IM[ny : ny * 2, N // 2], -ny // 2)
        ) * np.sqrt(ny)
        S[N - ny * 2 + 1 : N - ny + 1, N // 2] = ifft(
            shift(IM[N - ny * 2 + 1 : N - ny + 1, N // 2], -ny // 2)
        )[::-1] * np.sqrt(ny)
        for px in range(1, n):
            nx = 2 ** (px - 1)
            S[ny : ny * 2, nx : nx * 2] = ifftn(
                shift(IM[ny : ny * 2, nx : nx * 2], (-ny // 2, -nx // 2))
            ) * np.sqrt(ny * nx)
            S[N - ny * 2 + 1 : N - ny + 1, nx : nx * 2] = ifftn(
                shift(
                    IM[N - ny * 2 + 1 : N - ny + 1, nx : nx * 2], (-ny // 2, -nx // 2)
                )
            )[::-1] * np.sqrt(ny * nx)

    S[1 : N // 2, N // 2 + 1 :] = np.conj(S[N // 2 + 1 :, 1 : N // 2][::-1, ::-1])
    S[N // 2 + 1 :, N // 2 + 1 :] = np.conj(S[1 : N // 2, 1 : N // 2][::-1, ::-1])
    S[0, N // 2 + 1 :] = np.conj(S[0, 1 : N // 2][::-1])
    S[N // 2, N // 2 + 1 :] = np.conj(S[N // 2, 1 : N // 2][::-1])
    return S


class Test(unittest.TestCase):
    def __init__(self, methodName="runTest"):
        super(Test, self).__init__(methodName=methodName)
        self.image = chirp()
        i = np.finfo(self.image.dtype)
        self.eps = np.power(1.0, -i.precision + 1)

    def test_dst2(self):
        S = dst2(self.image)
        s = pyst.dst2(self.image)
        self.assertTrue(np.all(np.abs(S - s) < self.eps))

    def test_inverse(self):
        t = idst2(dst2(self.image))
        self.assertTrue(np.all(np.abs(self.image - t) < self.eps))


if __name__ == "__main__":
    unittest.main()
