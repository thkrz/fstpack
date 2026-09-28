# fstpack

Fast 1-D Stockwell transforms and 2-D discrete orthonormal Stockwell
transforms (DOST), with Fortran and Python interfaces.

## Installation

Building requires `gfortran`, `make`, a C compiler, Python 3 with
development headers, and NumPy with `f2py`.

    make
    make tests
    sudo make install

`make install` installs the Python extension, but not `libfstpack.a`.
By default it installs to
`/usr/local/lib/pythonX.Y/dist-packages`, using the version of
`python3` found at build time. If that directory is not on your
Python import path, set `PYDIST` to the appropriate package directory
when running `make install`.

Use `make uninstall` with the same installation settings to remove
the extension.

## Quick start

    import numpy as np
    import fstpack

    signal = np.arange(8, dtype=np.complex64)
    spectrum = fstpack.fst(signal)
    restored_signal = fstpack.ifst(spectrum)

    image = np.arange(64, dtype=np.complex64).reshape(8, 8)
    coefficients = fstpack.dost(image)
    restored_image = fstpack.idost(coefficients)
    local_spectrum = fstpack.voices(coefficients)[:, :, 2, 3]

These round trips are intended for real-valued data represented as
complex arrays. Independent negative-frequency content in a complex
input is not preserved.

## Python API

| Function | Input | Output |
| --- | --- | --- |
| `fst(h)` | Series `(N,)` | Voices `(N//2+1, N)` |
| `ifst(s)` | Voices `(N//2+1, N)` | Series `(N,)` |
| `dost(h)` | Square image `(N, N)` | DOST coefficients `(N, N)` |
| `idost(s)` | DOST coefficients `(N, N)` | Image `(N, N)` |
| `voices(s)` | DOST coefficients `(N, N)` | Local spectra `(M, M, N, N)` |

For the 2-D functions, `N` must be a power of two and at least 2;
`M = 2*log2(N)`. For 1-D transforms, `N` must be at least 1 and
accepted by the linked FFT implementation.

The wrappers use default-kind Fortran `complex`, normally exposed as
`numpy.complex64`. Inputs are not modified. `dost` and `idost` copy
their inputs before calling the in-place Fortran routines. Output axes
retain Fortran index order: `s[f, t]` is frequency voice `f` at time
`t`, and `voices(s)[:, :, x, y]` is the local spectrum at `(x, y)`.
Multidimensional outputs are Fortran-contiguous.

Array extents are normally inferred. For odd-length `fst` output, pass
the length explicitly when calling `ifst`:

    restored = fstpack.ifst(s, n=s.shape[1])

Without `n`, the `ifst` wrapper infers an even length from the number
of voices and rejects odd-length output. `voices` allocates a spectrum
at every image position; use Fortran's `cvoc2c` when only one position
is needed.

## Fortran API

All public procedures are in module `fstpack` and are `pure`.

    call cdst2f(c)
    call cdst2b(c)
    s = cfst1f(h)
    h = cfst1b(s)
    local = cvoc2c(s, x, y)
    all_local = cvoc2a(s)

    pure subroutine cdst2f(c)
      complex, intent(inout) :: c(0:, 0:)

    pure subroutine cdst2b(c)
      complex, intent(inout) :: c(0:, 0:)

    pure function cfst1f(h) result(s)
      complex, intent(in) :: h(0:)
      complex, allocatable :: s(:, :)

    pure function cfst1b(s) result(h)
      complex, intent(in) :: s(:, :)
      complex, allocatable :: h(:)

    pure function cvoc2c(s, x, y) result(h)
      complex, intent(in) :: s(0:, 0:)
      integer, intent(in) :: x, y
      complex, allocatable :: h(:, :)

    pure function cvoc2a(s) result(h)
      complex, intent(in) :: s(0:, 0:)
      complex, allocatable :: h(:, :, :, :)

Assumed-shape dummies declared with lower bound zero index from the
*first element* of the actual argument. A caller may pass a 1-based
array; `x = 0, y = 0` still selects its first element.

### 1-D fast S-transform

`cfst1f(h)` returns an array allocated with bounds
`(0:N/2, 0:N-1)`. Its first index is frequency voice `f`; its second
is time `t`. Voice 0 is `sum(h)/N` at every time. For even `N`, voice
`N/2` is Nyquist; for odd `N`, the final voice is the highest positive
frequency.

The transform uses the analytic spectrum of `h`: negative-frequency
bins are cleared. For each `f = 1 .. N/2`, it inverse-transforms the
Fourier samples

    H[(f + m) mod N] * exp(-2*pi**2*m**2/f**2),

with the Gaussian mirrored about `m = 0`. The Gaussian has the fixed
Stockwell width; there is no adjustable α parameter.

`cfst1b(s)` expects shape `(N/2+1, N)` in voice-then-time order. It
returns `h(1:N)`, with `h(1)` the first sample. It sums each voice over
time, reconstructs the two-sided spectrum by conjugate symmetry, and
inverse-transforms it. Because negative frequencies are not retained
independently, this is not an inverse for arbitrary complex signals.

### 2-D DOST

`cdst2f(c)` transforms a square `N`-by-`N` array in place, where
`N = 2**p` and `p >= 1`. `cdst2b(c)` transforms DOST coefficients back
in place. The first array index is x and the second is y.

The Fourier grid is divided into dyadic bands. With
`n = log2(N)-1`, positive band `v = 1 .. n` occupies

    [2**(v-1), 2**v - 1]

and has width `w = 2**(v-1)`. Index 0 is DC and index `N/2` is
Nyquist. Interior tiles combine one band from each axis; axis tiles
use one band. Each tile is circularly shifted by `floor(-w/2)` on
each band axis, inverse-transformed, and scaled by `sqrt(w)` for an
axis tile or `sqrt(wx*wy)` for an interior tile. Its spatial samples
remain in the tile's index rectangle. The four DC/Nyquist
intersections are unscaled.

`cdst2f` fills the half-plane `y > N/2` by conjugate reflection:

    c(x, y) = conjg(c((-x) mod N, N-y))

For real-valued images this is the expected Fourier symmetry. For
arbitrary complex images it discards independent content in that
half-plane; `cdst2b` is an inverse only on the represented subspace.

### Local DOST spectra

`cvoc2c(s, x, y)` samples a DOST array at zero-based offsets
`0 <= x,y < N`. It returns a 1-based `(M, M)` array, where
`M = 2*log2(N)`. `cvoc2a(s)` returns the same spectra for every
position, with shape `(M, M, N, N)`. In both results the first two
indices are x-voice and y-voice.

For `n = log2(N)-1`, each voice axis runs from `-n` through `n+1`.
Voice `0` is DC, voices `1 .. n` are positive bands, voices
`-n .. -1` are negative bands, and voice `n+1` is Nyquist. A voice
`v` is at Fortran index `v + log2(N)` or Python index
`v + log2(N) - 1`.

For a band of width `b`, position `x` selects spatial sample
`x*b/N` using integer division; y works the same way. DC and
Nyquist ignore position. Negative voices read coefficients from the
conjugate-mirrored tiles without conjugating them again. These
functions sample existing coefficients; they do not perform another
transform.

## Errors

The Fortran routines use `error stop` for checked invalid shapes or
coordinates and for FFT failures; they have no status argument. A
Fortran `error stop` also terminates a Python process. f2py may reject
invalid wrapper shapes with a Python exception before calling
Fortran.

## References

1. Drabycz, S., Stockwell, R. G., & Mitchell, J. R. (2009). Image
   texture characterization using the discrete orthonormal
   S-transform. *Journal of Digital Imaging*, 22, 696–708.
2. Brown, R. A., & Frayne, R. (2008). A fast discrete S-transform for
   biomedical signal processing. *30th Annual International Conference
   of the IEEE Engineering in Medicine and Biology Society*, 2586–2589.
3. Mansinha, L., Stockwell, R. G., & Lowe, R. P. (1997). Pattern
   analysis with two-dimensional spectral localisation: Applications
   of two-dimensional S transforms. *Physica A: Statistical Mechanics
   and its Applications*, 239(1–3), 286–295.
4. Stockwell, R. G., & Mansinha, L. (1996). Localization of the complex
   spectrum: The S transform. *IEEE Transactions on Signal
   Processing*, 44(4), 998–1001.
