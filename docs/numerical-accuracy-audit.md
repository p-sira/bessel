# Numerical accuracy audit

Audited baseline: `59f0a49`.

The audit combined the existing test suites with a deterministic boundary grid
(4,090 finite, nonzero SciPy references), a seeded random sweep (11,098 such
references), and high-precision checks of disputed values. The sweeps covered
real cylindrical J/Y/I/K and scaled I/K, spherical J/Y/I/K and their first
derivatives, and Struve H/L. Random inputs used NumPy seed 9817, 2,000 draws
per API, integer spherical orders in `[0, 1200)`, other orders uniform in
`[-499, 499)`, and positive arguments log-uniform in `[1e-3, 10^3.2)`.
Nonfinite and zero SciPy outputs were excluded from the comparison, so the
sweep alone does not certify overflow or underflow behavior.

SciPy 1.18.1 / NumPy 2.5.4 supplied the broad comparisons. Confirmed failures
and disputed SciPy values were evaluated with mpmath 1.4.1. The committed
178 regression references use exact binary64 inputs, evaluated at both 90 and
130 decimal digits; both evaluations round to the same binary64 result.
75 of these cases fail on the baseline and all pass with the fixes. Tests use
`2e-11 * abs(expected) + 4 * f64::from_bits(1)`; the absolute allowance is
needed for subnormal derivatives. This is a regression bound, not a global
accuracy guarantee.

## Confirmed issues and fixes

| Issue | Example before the fix | Fix |
| --- | --- | --- |
| Tiny-argument J loses the argument or parity | `J_0.5(2^-1074)` returns zero; `J_3(-1e-10)` has the wrong sign | Evaluate the leading term using log-gamma, preserve parity, and avoid forming `x/2` |
| J series cancellation | `J_283.8627249559925(135.55386283898832)` has relative error about `1.4e-8` | Retry cancellation-heavy sums using the existing double-double arithmetic for terms and accumulation |
| Fractional-order J/Y reflection overflows or loses precision | `Y_295.60954413659204(65.40689047242715)` returns NaN instead of about `-3.25e154`; negative J has similar failures | Use stable positive-order Y recurrence above order 50 and reflect afterward; include the reflection coefficient in recurrence seeds when unmultiplied Y overflows |
| Spherical I loses a finite result through scaled underflow | `i_500(100)` returns zero instead of about `1.41e-284`; its derivative also returns zero | Prefer the native unscaled I kernel, retaining the scaled fallback when cylindrical I overflows |
| Spherical derivatives underflow before multiplication by the order | `i_68'(0.0011338333600859067)` loses subnormal precision | Include the order in the leading derivative coefficient before the product |
| Half-order Struve H cancels near its double zeros | `H_0.5(710)` returns NaN; rounded multiples of `2*pi` can have large relative errors | Use `sqrt(2/pi)/sqrt(x) * 2*sin(x/2)^2`, from [DLMF 11.4.5](https://dlmf.nist.gov/11.4.E5) |

The spherical reduction and derivative identities follow
[DLMF 10.47](https://dlmf.nist.gov/10.47) and
[DLMF 10.51](https://dlmf.nist.gov/10.51).

## Reference disagreements

High precision confirmed that the crate was already correct for several
extreme cases where SciPy disagreed: `Y_-0.5(1e-310)`, spherical J derivatives
at `(n,x)=(2,1e-200)`, `(103,0.08776684995814026)` and
`(219,6.65254869143548)`, and the modified spherical derivative at
`(1,1e-200)`. Those results are retained as high-precision regressions.
SciPy's `H_0.5(710)` also differs from the high-precision value by about `5e-8`
relative; the new elementary evaluation agrees with high precision.

After fixing the confirmed issues, the boundary and random sweeps have no
unexplained disagreements exceeding `1e-8` relative and `1e-322` absolute.
The remaining disagreements are reference errors confirmed with mpmath.

## Reproduction and limits

```sh
uv run --script tests/fixtures/accuracy_audit_generate.py
cargo test --test accuracy_audit
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --all -- --check
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps
```

The existing static suites additionally cover complex branch cuts and Hankel
functions, Airy, Kelvin, integrals, zeros, ratios, and binary32 conversion.
This audit adds targeted coverage to those suites; it does not establish a
uniform error bound for every public API over all binary64 inputs. In
particular, the documented cancellation limitation of `besselpoly` at large
arguments remains outside the fixes here. The J series retries in double-double precision only when the largest term
exceeds the final sum's magnitude by more than `1e5`. This keeps ordinary
inputs on the binary64 path while retaining the accuracy regressions.

## Performance investigation

The initial PR revision (`dcab5e9`) used double-double arithmetic for every J
series. CodSpeed reported seven regressions of 42–51% efficiency, along with
a warning that base and head ran on different hardware. Local measurements
on the same machine confirmed a real slowdown across all seven paths.

After making extended precision conditional on cancellation, median local
wall times changed as follows (microseconds per benchmark sweep):

| Benchmark | Main (`f8ddffb`) | Initial PR | Conditional retry |
| --- | ---: | ---: | ---: |
| `jv` | 2.909 | 10.43 | 3.029 |
| `jv_f32` | 2.898 | 10.40 | 3.002 |
| `yv` | 4.873 | 17.34 | 4.995 |
| `jvp[1]` | 5.665 | 21.59 | 6.014 |
| `jvp[4]` | 17.00 | 56.08 | 17.59 |
| `hankel1` | 8.826 | 28.86 | 9.114 |
| `hankel2e` | 9.585 | 29.54 | 9.838 |

These local timings are diagnostic measurements, not CodSpeed simulation
results. The remaining local overhead is about 2–6%, including the check for
cancellation. Reproduce the benchmark selection with:

```sh
cargo bench --bench bessel -- 'cylindrical_real_order::(jv|yv)' \
  'derivatives::jvp' 'hankel::' --sample-count 300 --min-time 0.2
```

## Half-integer derivative performance (0.1.3 candidate)

The Sci `yvp`, `h1vp`, and `h2vp` regressions exercise 32 uniformly spaced
arguments in `[0.5, 20]`, with order `1.5` and derivative order `2`.
Half-integer J/Y have elementary initial values at `±1/2`; an upward recurrence
avoids repeated cancellation-sensitive J series for `|v| <= 64.5` when
`x >= |v|`. Other inputs retain the general kernels and their extended-precision
safeguards. Independent 100-digit mpmath references cover 496 J/Y and derivative
values, including nonbenchmark orders, switch boundaries, near zeros, large
arguments, and the full Sci grid.

Matched same-core release measurements (microseconds per 32-point sweep):

| Sci benchmark | 0.1.1 baseline | 0.1.2 | 0.1.3 candidate |
| --- | ---: | ---: | ---: |
| `yvp` | 19.68 | 23.23 | 4.67 |
| `h1vp` | 39.53 | 45.93 | 10.77 |
| `h2vp` | 36.99 | 43.88 | 8.49 |

These are local diagnostic timings, not CodSpeed results. The new Bessel
benchmarks `derivatives::{yvp,h1vp,h2vp}_sci_grid` retain this workload.
