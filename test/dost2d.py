import matplotlib.pyplot as plt
import numpy as np

from numpy.fft import fft, fft2, ifft, ifftn


def chirp(n=64):
    h = np.zeros((n, n))
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

    IM = np.zeros((N, N), dtype=complex)

    n = int(np.log2(N))
    for py in range(1, n):
        ny = 2 ** (py - 1)
        IM[0, ny : ny * 2] = shift(fft(S[0, ny : ny * 2] / np.sqrt(ny)), ny // 2)
        IM[ny : ny * 2, 0] = shift(fft(S[ny : ny * 2, 0] / np.sqrt(ny)), ny // 2)
        IM[N - ny * 2 + 1 : N - ny + 1, 0] = shift(
            fft(S[N - ny * 2 + 1 : N - ny + 1, 0][::-1] / np.sqrt(ny)), ny // 2
        )
        IM[N // 2, ny : ny * 2] = shift(
            fft(S[N // 2, ny : ny * 2] / np.sqrt(ny)), ny // 2
        )
        IM[ny : ny * 2, N // 2] = shift(
            fft(S[ny : ny * 2, N // 2] / np.sqrt(ny)), ny // 2
        )
        IM[N - ny * 2 + 1 : N - ny + 1, N // 2] = shift(
            fft(S[N - ny * 2 + 1 : N - ny + 1, N // 2][::-1] / np.sqrt(ny)), ny // 2
        )
        for px in range(1, n):
            nx = 2 ** (px - 1)
            IM[ny : ny * 2, nx : nx * 2] = shift(
                fft2(S[ny : ny * 2, nx : nx * 2] / np.sqrt(ny * nx)), (ny // 2, nx // 2)
            )
            IM[N - ny * 2 + 1 : N - ny + 1, nx : nx * 2] = shift(
                fft2(
                    S[N - ny * 2 + 1 : N - ny + 1, nx : nx * 2][::-1] / np.sqrt(ny * nx)
                ),
                (ny // 2, nx // 2),
            )

    IM[0, 0] = S[0, 0]
    IM[N // 2, 0] = S[N // 2, 0]
    IM[0, N // 2] = S[0, N // 2]
    IM[N // 2, 1] = S[N // 2, 1]
    IM[1, N // 2] = S[1, N // 2]
    IM[N // 2, N // 2] = S[N // 2, N // 2]

    im = ifftn(IM)
    return np.real(im)


def dst2(im):
    N = len(im)
    IM = fft2(im)

    S = np.zeros((N, N), dtype=complex)
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


if __name__ == "__main__":
    h = chirp()
    S = dst2(h)
    t = idst2(S)

    if np.all(np.abs(h - t) < 1e-6):
        print("T")
        exit(0)

    fig = plt.figure()
    ax = fig.add_subplot(211, projection="3d")
    ax.set_title("Synthetic signal")
    X = np.arange(h.shape[1])
    Y = np.arange(h.shape[0])
    X, Y = np.meshgrid(X, Y)
    ax.plot_surface(X, Y, h)

    ax = fig.add_subplot(212, projection="3d")
    ax.set_title("Reverted signal")
    ax.plot_surface(X, Y, t)

    plt.show()
