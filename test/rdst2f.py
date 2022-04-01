# Source code from:
# Drabycz, S., Stockwell, R.G. & Mitchell, J.R. Image Texture Characterization
# Using the Discrete Orthonormal S-Transform. J Digit Imaging 22, 696 (2009).
# https://doi.org/10.1007/s10278-008-9138-8
import numpy as np
import sys


def fft2(c):
    return np.fft.fft2(c) / c.size


def ifft(c):
    return np.fft.ifft(c) * c.size


def ifftn(c):
    return np.fft.ifft2(c) * c.size


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
    IM = fft2(im)
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

    return S
    S[1 : N // 2, N // 2 + 1 :] = np.conj(S[N // 2 + 1 :, 1 : N // 2][::-1, ::-1])
    S[N // 2 + 1 :, N // 2 + 1 :] = np.conj(S[1 : N // 2 :, 1 : N // 2][::-1, ::-1])
    S[0, N // 2 + 1 :] = np.conj(S[0, 1 : N // 2][::-1])
    S[N // 2, N // 2 + 1 :] = np.conj(S[N // 2, 1 : N // 2][::-1])
    return S


if __name__ == "__main__":
    # import fpy

    r = np.loadtxt(sys.argv[-1])
    s = rdst2f(r)
    print(s)
    # fpy.dumptxt(s, sys.stdout)
