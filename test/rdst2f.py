# Source code from:
# Drabycz, S., Stockwell, R.G. & Mitchell, J.R. Image Texture Characterization
# Using the Discrete Orthonormal S-Transform. J Digit Imaging 22, 696 (2009).
# https://doi.org/10.1007/s10278-008-9138-8
import numpy as np
import re


def shift(a, n):
    if len(np.shape(a)) == 1:
        b = np.concatenate((a[n:], a[0:n]))
    elif len(np.shape(a)) == 2:
        c = np.concatenate((a[n[0] :, :], a[0 : n[0], :]))
        b = np.concatenate((c[:, n[1] :], c[:, 0 : n[1]]), 1)
    return b


# Image to be processed size N x N (2^n x 2^n)
def rdst2f(im):
    N = len(im)
    IM = np.fft.fft2(im) / (N*N)
    S = np.zeros((N, N), complex)
    S[0, 0] = IM[0, 0]
    S[N // 2, 0] = IM[N // 2, 0]
    S[0, N // 2] = IM[0, N // 2]
    S[N // 2, 1] = IM[N // 2, 1]
    S[1, N // 2] = IM[1, N // 2]
    S[N // 2, N // 2] = IM[N // 2, N // 2]
    n = int(np.log2(N))

    for py in range(1, n):
        ny = 2 ** (py - 1)
        S[0, ny : ny * 2] = np.fft.ifft(shift(IM[0, ny : ny * 2], -ny // 2)) * np.sqrt(
            ny
        )
        S[ny : ny * 2, 0] = np.fft.ifft(shift(IM[ny : ny * 2, 0], -ny // 2)) * np.sqrt(
            ny
        )
        S[N - ny * 2 + 1 : N - ny + 1, 0] = np.fft.ifft(
            shift(IM[N - ny * 2 + 1 : N - ny + 1, 0], -ny // 2)
        )[::-1] * np.sqrt(ny)
        S[N // 2, ny : ny * 2] = np.fft.ifft(
            shift(IM[N // 2, ny : ny * 2], -ny // 2)
        ) * np.sqrt(ny)
        S[ny : ny * 2, N // 2] = np.fft.ifft(
            shift(IM[ny : ny * 2, N // 2], -ny // 2)
        ) * np.sqrt(ny)
        S[N - ny * 2 + 1 : N - ny + 1, N // 2] = np.fft.ifft(
            shift(IM[N - ny * 2 + 1 : N - ny + 1, N // 2], -ny // 2)
        )[::-1] * np.sqrt(ny)
        for px in range(1, n):
            nx = 2 ** (px - 1)
            S[ny : ny * 2, nx : nx * 2] = np.fft.ifftn(
                shift(IM[ny : ny * 2, nx : nx * 2], (-ny // 2, -nx // 2))
            ) * np.sqrt(ny * nx)
            S[N - ny * 2 + 1 : N - ny + 1, nx : nx * 2] = np.fft.ifftn(
                shift(
                    IM[N - ny * 2 + 1 : N - ny + 1, nx : nx * 2], (-ny // 2, -nx // 2)
                )
            )[::-1] * np.sqrt(ny * nx)

    S[1 : N // 2, N // 2 + 1 :] = np.conjugate(S[N // 2 + 1 :, 1 : N // 2][::-1, ::-1])
    S[N // 2 + 1 :, N // 2 + 1 :] = np.conjugate(
        S[1 : N // 2 :, 1 : N // 2][::-1, ::-1]
    )
    S[0, N // 2 + 1 :] = np.conjugate(S[0, 1 : N // 2][::-1])
    S[N // 2, N // 2 + 1 :] = np.conjugate(S[N // 2, 1 : N // 2][::-1])
    return S


if __name__ == "__main__":
    r = np.loadtxt("1.1.03.raw")
    t = np.zeros(r.shape, complex)
    pat = re.compile(r"\(\s*(-?\d+\.\d+),\s*(-?\d+\.\d+)\)")
    with open("1.1.03.dost") as f:
        for i, s in enumerate(f):
            for j, m in enumerate(pat.finditer(s)):
                t[i, j] = complex(float(m.group(1)), float(m.group(2)))
    s = rdst2f(r)
    assert np.all(np.abs(s - t) < 1e-4)
