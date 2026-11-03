# bessel

<p align="center">
    <a href="LICENSE" style="text-decoration:none">
        <img src="https://img.shields.io/badge/License-See%20LICENSE-brightgreen.svg" alt="License">
    </a>
    <a href="https://crates.io/crates/bessel" style="text-decoration:none">
        <img src="https://img.shields.io/crates/v/bessel" alt="Crate">
    </a>
    <a href="https://docs.rs/bessel" style="text-decoration:none">
        <img src="https://img.shields.io/badge/Docs-docs.rs-blue" alt="Documentation">
    </a>
    <a href="https://app.codspeed.io/p-sira/bessel?utm_source=badge" style="text-decoration:none">
        <img src="https://img.shields.io/endpoint?url=https://codspeed.io/badge.json" alt="CodSpeed"/>
    </a>
</p>

Pure Rust library for Bessel and related special functions, including cylindrical, spherical, Hankel, Kelvin, Struve, and Airy functions.

- Evaluate Bessel J, Y, I, and K functions, including scaled forms, derivatives, ratios, and zeros.
- Compute complex Bessel and Hankel functions, spherical Bessel functions, and Riccati–Bessel functions.
- Evaluate the real Airy quartet (Ai, Ai′, Bi, Bi′) with `airy(x)` and its scaled form `airye(x)` (for `x >= 0`).
- Evaluate Kelvin and Struve functions, their zeros, and related integrals.

```rust
use bessel::{jv, iv, hankel1, spherical_jn, Complex64};

assert_eq!(jv(0.0, 0.0), 1.0);
assert_eq!(iv(0.0, 0.0), 1.0);
assert!((spherical_jn(0, 1.0_f64, false) - 1.0_f64.sin()).abs() < 1e-14);
assert!((hankel1(0.0, 1.0) - Complex64::new(bessel::j0(1.0), bessel::y0(1.0))).norm() < 1e-13);

let (ai, aip, bi, bip) = bessel::airy(0.0_f64).unwrap();
assert!(ai > 0.0 && aip < 0.0 && bi > 0.0 && bip > 0.0);
let z = Complex64::new(1.0, 2.0);
let j = bessel::complex::jv(0.5_f64, z).unwrap();
assert!(j.re.is_finite() && j.im.is_finite());
```

The native `backend::complex` evaluator implements J, Y, I, K, and both Hankel kinds, including scaled variants, exposed through `bessel::complex::{jv, yv, iv, kv, jve, yve, ive, kve, hankel1, hankel2, hankel1e, hankel2e}`. The native `backend::airy` handles the real Airy quartet and internal complex Airy calculations; the public `airy(x)` and `airye(x)` return `Result<(Ai, Ai′, Bi, Bi′)>` for real `x`. Positive-real `kv` and `kve` use the independent `backend::real_k` kernel.

Real-valued scalar functions support `f32` and `f64` through `num_traits::Float` (evaluated internally in `f64`). Infallible scalar functions return `NaN` for invalid inputs; `airy`, `airye`, `checked`, `checked_complex`, and `complex` provide fallible interfaces. The `kelvin_zeros(nt)` API returns all eight zero sequences; `kelvin_zeros_for(family, nt)` selects one family. See [test fixture generation](tests/scipy/README.md) for independent SciPy/mpmath numerical references.
