# bessel

Standalone pure Rust cylindrical, spherical, Hankel, Kelvin, and Struve functions, extracted from the Bessel numerical family in [Sci](https://github.com/p-sira/sci). The library has no dependency on Sci, sci-special, or sibling crates. Arithmetic uses f64 internally; some checked APIs accept f32 through `num_traits::Float` and convert to/from f64.

```rust
use bessel::{jv, iv, spherical_jn, hankel1, Complex64};
assert!((jv(0.0, 0.0) - 1.0).abs() < 1e-14);
assert!((iv(0.0, 0.0) - 1.0).abs() < 1e-14);
assert!((spherical_jn(0, 1.0_f64, false) - 1.0_f64.sin()).abs() < 1e-14);
assert!((hankel1(0.0, 1.0) - Complex64::new(bessel::j0(1.0), bessel::y0(1.0))).norm() < 1e-13);
```

## Public functions and conventions

- **Cylindrical:** `j0`, `j1`, `jv`, `jn`, `y0`, `y1`, `yv`, `yn`, `i0`, `i1`, `iv`, `k0`, `k1`, `kv`, `kn`, with scaled `i0e`, `i1e`, `ive`, `k0e`, `k1e`, `kve`, `jve`, `yve` (the last two equal J and Y on real arguments). `iv_ratio`, `iv_ratio_c` and `iv_ratioinv` use a scaled continued fraction and bounded bisection.
- **Complex Bessel / Hankel:** `complex::{jv,yv,iv,kv,jve,yve,ive,kve,hankel1,hankel1e,hankel2,hankel2e}` take `v: f64` and `z: Complex64`, returning `Result<Complex64, complex_bessel::Error>`. Root-level `hankel1`, `hankel1e`, `hankel2`, `hankel2e` take real `x` and return `Complex64`. `checked_complex` exposes the original generic complex Hankel kernels.
- **Derivatives:** `jvp`, `yvp`, `ivp`, `kvp`, `h1vp`, and `h2vp` use adjacent-order derivative identities with order `n <= 1024`.
- **Spherical:** `spherical_jn`, `spherical_yn`, `spherical_in`, `spherical_kn` take `(n: i64, x: T, derivative: bool)` with integer `0 <= n <= 4096`, using stable recurrences and a half-integer AMOS fallback. The modified `spherical_kn` follows SciPy's `sqrt(pi/(2x)) K_(n+1/2)(x)` normalization. `riccati_jn` and `riccati_yn` return parallel vectors of `x*j_n(x)` or `x*y_n(x)` and their derivatives.
- **Struve:** `struve`, `modstruve` are H_v and L_v for real arguments; `checked::{struve_h,struve_l}` return explicit errors for invalid inputs or estimated precision loss.
- **Kelvin:** `ber`, `bei`, `ker`, `kei` and derivatives `berp`, `beip`, `kerp`, `keip`; `kelvin` returns four complex pairs (values and derivatives). Kelvin roots use `ber_zeros`, `bei_zeros`, `ker_zeros`, `kei_zeros`, `berp_zeros`, `beip_zeros`, `kerp_zeros`, `keip_zeros` or `kelvin_zeros(KelvinFunction, count)`. The `checked` module provides generic, fallible kernels.
- **Cylindrical real zeros:** `jn_zeros`, `jnp_zeros`, `yn_zeros`, `ynp_zeros`, `jnyn_zeros`, `jnjnp_zeros`, `y0_zeros`, `y1_zeros`, `y1p_zeros`; roots exclude the origin except the merged `jnjnp_zeros` convention.
- **Integrals:** `itj0y0`, `it2j0y0`, `iti0k0`, `it2i0k0`, `itstruve0`, `it2struve0`, `itmodstruve0` calculate the Bessel/Struve definite integral pairs and tails with bounded legacy series, corrected exponentially scaled K₀ tails and numerical quadrature. `itstruve0` is validated for `|x| <= 512`; larger arguments return NaN pending a reliable asymptotic.
- **Related:** `besselpoly(a, lambda, v)` computes the weighted integral of `J_v(2at)` on `[0,1]`; `airy`, `airye` expose the Airy quartet used by large-order Bessel expansions.

The real J and Y equations satisfy `x² y''+x y'+(x²-v²)y=0`; the modified solutions satisfy `x² y''+x y'-(x²+v²)y=0`. Hankel functions are `H_v^(1,2) = J_v ± iY_v`. The Kelvin convention is `ber+i bei = J_0(x e^(3πi/4))` and `ker+i kei = K_0(x e^(πi/4))`. Arguments, outputs, orders, and normalization are dimensionless. Negative real arguments require integer order for J/I, and Y/K use a positive real domain; complex functions use AMOS principal branches. Floating-point NaN represents invalid real domains; some singular limits produce infinities. Error-aware interfaces (`checked`, `checked_complex`, `complex`) can distinguish invalid inputs, nonconvergence and overflow. All numerical iterations have finite work limits. Relative accuracy may deteriorate near zeros, branch cuts, turning points and cancellation-sensitive high derivatives; this is a binary64 library, even when converting f32.

**Current scope:** all 79 original Sci Bessel namespace names are now implemented. The original Sci `kelvin_zeros(nt)` returns eight parallel sequences; here `kelvin_zeros(KelvinFunction, nt)` selects one sequence, and the eight individual convenience functions are available. Functions returning an error type in Sci may use a `&'static str` error type here. Full signature parity with SciPy is not intended.

## Verification and fixture regeneration

`tests/scipy/families_inputs.csv` contains 1,192 independent cases; `tests/fixtures/families_expected.csv` is generated with pinned **SciPy 1.18.1 / NumPy 2.5.4**, and tested alongside 1,259 real I/IE regression rows from Sci (`tests/scipy/bessel_iv_cases.json`, generated using the same pinned versions). The SciPy 1.18.1 implementation of `it2i0k0` loses accuracy in both components for large x; its `itstruve0` results also fail the defining integral for x >= 20. Those fixture rows use independent SciPy special-function quadrature with scaled kernels, as documented in the Python generator. Additional checked-in complex branch-cut, Hankel derivative, tiny-tail, and binary32-subnormal fixtures exercise near-zero and complex domains. Inputs and expected values are committed separately, and normal Rust tests do not run Python.

```bash
uv run --script tests/scipy/generate_families.py
uv run --script tests/scipy/generate_bessel_iv.py
uv run --script tests/scipy/generate_tails.py
cargo fmt --all -- --check
cargo test
cargo clippy --all-targets -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps
```

## Licensing and references

The extracted source comes from Sci (BSD-3-Clause), Cephes code by Stephen L. Moshier, public XSF/SciPy implementations, portions licensed under the Boost Software License 1.0, and the MIT-licensed `spec_math` translation. Retained copyright and distribution notices are in `LICENSE` and `LICENSES/`; dependency licensing follows each package's own terms. Mathematical and numerical references: [NIST DLMF Chapter 10 (Bessel)](https://dlmf.nist.gov/10), [NIST DLMF Chapter 11 (Struve)](https://dlmf.nist.gov/11), SciPy `special` documentation, and Gautschi and Slavik (1978) for modified Bessel ratios.
