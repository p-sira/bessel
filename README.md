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

Pure Rust library for cylindrical, spherical, Hankel, Kelvin, and Struve functions.

- Evaluate Bessel J, Y, I, and K functions, including scaled forms, derivatives, ratios, and zeros.
- Compute complex Bessel and Hankel functions, spherical Bessel functions, and Riccati–Bessel functions.
- Evaluate Kelvin and Struve functions, their zeros, and related integrals.

```rust
use bessel::{jv, iv, hankel1, spherical_jn, Complex64};

assert_eq!(jv(0.0, 0.0), 1.0);
assert_eq!(iv(0.0, 0.0), 1.0);
assert!((spherical_jn(0, 1.0_f64, false) - 1.0_f64.sin()).abs() < 1e-14);
assert!((hankel1(0.0, 1.0) - Complex64::new(bessel::j0(1.0), bessel::y0(1.0))).norm() < 1e-13);
```

Real-valued scalar functions support `f32` and `f64` through `num_traits::Float` (evaluated internally in `f64`). They return `NaN` for invalid inputs; the `checked`, `checked_complex`, and `complex` modules expose fallible alternatives for their respective functions. The `kelvin_zeros(nt)` API returns all eight zero sequences; `kelvin_zeros_for(family, nt)` selects one family. See [test fixture generation](tests/scipy/README.md) for numerical validation and reproducibility.

## Benchmarks

Benchmarks live in [`benches/`](benches) and use [divan](https://github.com/nvzqz/divan) through the CodSpeed compatibility layer. Run them locally with `cargo bench`; they are tracked continuously on [CodSpeed](https://app.codspeed.io/p-sira/bessel).
