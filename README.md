# FSTPACK(3)

## NAME

fstpack — 1-D fast Stockwell transform and 2-D DOST

## SYNOPSIS

Fortran — module `fstpack`:

    call cdst2f(c)
    call cdst2b(c)
    s = cfst1f(h)
    h = cfst1b(s)
    p = lspec2(s, x, y)

Python — module `fstpack`:

    s = fst(h)
    h = ifst(s)
    s = dost(h)
    h = idost(s)
    p = local_spectrum(s, x, y)

## DESCRIPTION

`fstpack` computes Stockwell transforms. Values are default-kind
`complex` (single precision unless the compiler promotes default
real). The f2py wrappers expose that kind as `numpy.complex64` under
the same condition.

The 1-D pair is the fast S-transform: redundant, positive frequencies
only, Gaussian window fixed at the Stockwell width (no α parameter).
The 2-D pair is the discrete orthonormal Stockwell transform (DOST)
on a square power-of-two grid, dyadic partition fixed. `lspec2` reads
one DOST coefficient from each dyadic voice at one spatial point.

Negative frequencies are not an independent part of the result. The
1-D transform keeps the analytic spectrum. The 2-D transform writes
the half-plane `y > N/2` as the conjugate mirror of `y < N/2`, which
is the real-image case and a projection otherwise. `cdst2b` inverts
`cdst2f` on that subspace. `cfst1b` inverts `cfst1f`.

Every public procedure is `pure`. A failed precondition is
`error stop`. There is no status argument.

Names: `c` complex, `f` forward, `b` backward, `1`/`2` the dimension.
`dst` is the DOST, `fst` the fast S-transform, `lspec` the local
spectrum.

Fortran dummies declared `(0:,0:)` are indexed from the first element
of the actual argument. A 1-based caller array is legal; offsets below
are from element 1 of that array, not from its declared lower bound.
2-D Python axes follow the Fortran axes: index 0 is x, index 1 is y.
f2py keeps that index order and returns Fortran-contiguous arrays.

## FORTRAN

    subroutine cdst2f(c)
      complex, intent(inout) :: c(0:, 0:)

    subroutine cdst2b(c)
      complex, intent(inout) :: c(0:, 0:)

    function cfst1f(h) result(s)
      complex, intent(in) :: h(0:)
      complex, allocatable :: s(:, :)

    function cfst1b(s) result(h)
      complex, intent(in) :: s(:, :)
      complex, allocatable :: h(:)

    function lspec2(s, x, y) result(h)
      complex, intent(in) :: s(0:, 0:)
      integer, intent(in) :: x, y
      complex, allocatable :: h(:, :)

### cdst2f(c)

Forward 2-D DOST, in place.

`c` is square of order `N = 2**p`, `p >= 1`. On return the same array
holds DOST coefficients.

The Fourier grid is split into dyadic tiles. For voice `v = 1 .. n`,
`n = log2(N)-1`, the positive band is the index range

    [2**(v-1), 2**v - 1]

of width `w = 2**(v-1)`. DC is index 0. Nyquist is index `N/2`.
Interior tiles are the Cartesian product of two bands; axis tiles are
1-D bands. Each tile is circularly shifted by `floor(-w/2)` along each
of its band axes (`-1` when `w = 1`), inverse-transformed, and scaled
by `sqrt(w)` (by `sqrt(wx*wy)` for an interior tile). Samples inside a
tile are then spatial, still stored in that index rectangle. The
half-plane `y > N/2` is filled by

    c(x, y) = conjg(c(N-x, N-y))

with the axis rows handled the same way. `(0,0)`, `(N/2,0)`,
`(0,N/2)`, and `(N/2,N/2)` are the unscaled DC and Nyquist bins.

### cdst2b(c)

Inverse of `cdst2f`, in place. Same argument and constraints. Undoes
the shifts and the square-root scales, restores conjugate symmetry,
then inverse-transforms the full array.

### cfst1f(h)

Forward 1-D fast S-transform.

`h` has length `N >= 1`. `N` need not be a power of two; it must be a
length the linked 1-D FFT accepts. The result is allocated with bounds

    s(0:N/2, 0:N-1)

so `s(f, t)` is voice `f` at time `t`. Voice 0 is DC, voice `N/2` is
Nyquist.

Voice 0 is `sum(h)/N` at every `t`. For `f = 1 .. N/2`, let `H` be
the Fourier transform of `h` after the frequency-domain Hilbert
transform (negative bins cleared). Voice `f` is the inverse FFT of

    H[(f + m) mod N] * exp(-2 * pi**2 * m**2 / f**2),

`m = 0 .. N-1`, Gaussian mirrored about 0.

### cfst1b(s)

Inverse of `cfst1f`.

`s` must have shape `(N/2+1, N)` in storage order (voice, then time),
as returned by `cfst1f`. The result is allocated `h(1:N)`: `h(1)` is
the first sample of the reconstructed series. Each voice is summed
over time, the two-sided spectrum is restored, and the inverse FFT is
scaled by `1/N`.

### lspec2(s, x, y)

Local spectrum of a 2-D DOST array.

`s` is `N` by `N`, `N = 2**p`, `p >= 1`, in the layout `cdst2f` writes.
`x` and `y` are offsets from the first sample, `0 <= x,y < N`. The
result is allocated 1-based, shape `(M, M)`, `M = 2*log2(N)`.

Let `n = log2(N)-1`. Both axes run over voices `v = -n .. n+1`, stored
at 1-based index `v + log2(N)`. The first index is the x-voice, the
second the y-voice.

    v          what                         width
    -n .. -1   negative octave |v|          2**(|v|-1)
     0         DC                           1, ignores x and y
     1 .. n    positive band [2**(v-1), 2**v-1]
                                           2**(v-1)
     n+1       Nyquist                      1, ignores x and y

A band of width `b` holds `b` spatial samples. The sample taken for
position `x` is `x*b/N` (truncating division); likewise for `y`.
Positions in a block of length `N/b` share a coefficient. Negative
voices are read from the conjugate-mirrored tile, not conjugated
again. On an array just returned by `cdst2f`, octave voice `(-a,-b)`
is the conjugate of voice `(a,b)`.

The result is a sample of `s`. It is not a transform and has no
inverse here.

## PYTHON

    import fstpack

Thin f2py wrappers. Extent arguments are inferred and must not be
passed. Inputs are not modified. 2-D transforms are in place in
Fortran and returning in Python: the wrapper copies, then calls.

    fst(h) -> s
        cfst1f. h shape (N,). s shape (N/2+1, N), s[f, t].

    ifst(s) -> h
        cfst1b. s shape (N/2+1, N). h shape (N,).

    dost(h) -> s
        cdst2f. h and s shape (N, N), N = 2**p, p >= 1.

    idost(s) -> h
        cdst2b. Same shapes.

    local_spectrum(s, x, y) -> p
        lspec2. s shape (N, N). x, y integers, 0 <= x,y < N.
        p shape (M, M), M = 2*log2(N). 0-based voice index is
        v + log2(N) - 1, with v as in lspec2. See BUGS.

A failed precondition terminates the process. It is not raised as a
Python exception.

## DIAGNOSTICS

`error stop` if a 2-D array is not square, `N` is not a power of two,
`(x, y)` lies outside the image, or an FFT reports failure.

## SEE ALSO

1. Drabycz, S., Stockwell, R.G. & Mitchell, J.R. (2009). Image Texture
   Characterization Using the Discrete Orthonormal S-Transform. *Journal of
   Digital Imaging*, 22, 696-708.
2. Brown, R. A., &amp; Frayne, R. (2008). A Fast Discrete S-Transform for Biomedical Signal Processing.
   *30th Annual International Conference of the IEEE Engineering in Medicine and Biology Society*, 2586-2589.
3. Mansinha, L., Stockwell, R.G., &amp; Lowe, R.P. (1997).
   Pattern analysis with two-dimensional spectral localisation: Applications of two-dimensional S transforms.
   *Physica A: Statistical Mechanics and its Applications, 239*(**1-3**), 286-295.
4. Stockwell, R.G., &amp; Mansinha, L. (1996). Localization of the complex spectrum: the S transform.
   *IEEE Transactions on Signal Processing*, 44(**4**), 998-1001.
