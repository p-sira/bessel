//! Benchmarks for the public Bessel, Hankel, spherical, Kelvin, and Struve APIs.
//!
//! Scalar functions are evaluated over a fixed set of arguments spanning the
//! small-argument series, transition, and large-argument asymptotic regimes so
//! that each benchmark exercises every code path a typical caller hits.

use divan::black_box;

fn main() {
    divan::main();
}

/// Arguments covering small, moderate, and large `x`.
const XS: [f64; 8] = [0.05, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 100.0];
/// Representative real orders, including integer and fractional values.
const ORDERS: [f64; 4] = [0.0, 0.5, 2.0, 7.25];

#[inline(always)]
fn sweep(f: impl Fn(f64) -> f64) -> f64 {
    XS.iter().map(|&x| f(black_box(x))).sum()
}

#[inline(always)]
fn sweep_order(f: impl Fn(f64, f64) -> f64) -> f64 {
    let mut acc = 0.0;
    for &v in &ORDERS {
        for &x in &XS {
            acc += f(black_box(v), black_box(x));
        }
    }
    acc
}

mod cylindrical_integer_order {
    use super::*;

    #[divan::bench]
    fn j0() -> f64 {
        sweep(bessel::j0)
    }

    #[divan::bench]
    fn j1() -> f64 {
        sweep(bessel::j1)
    }

    #[divan::bench]
    fn y0() -> f64 {
        sweep(bessel::y0)
    }

    #[divan::bench]
    fn y1() -> f64 {
        sweep(bessel::y1)
    }

    #[divan::bench]
    fn i0() -> f64 {
        sweep(bessel::i0)
    }

    #[divan::bench]
    fn i1() -> f64 {
        sweep(bessel::i1)
    }

    #[divan::bench]
    fn k0() -> f64 {
        sweep(bessel::k0)
    }

    #[divan::bench]
    fn k1() -> f64 {
        sweep(bessel::k1)
    }

    #[divan::bench]
    fn jn_order_5() -> f64 {
        sweep(|x| bessel::jn(5, x))
    }

    #[divan::bench]
    fn yn_order_5() -> f64 {
        sweep(|x| bessel::yn(5, x))
    }
}

mod cylindrical_real_order {
    use super::*;

    #[divan::bench]
    fn jv() -> f64 {
        sweep_order(bessel::jv)
    }

    #[divan::bench]
    fn yv() -> f64 {
        sweep_order(bessel::yv)
    }

    #[divan::bench]
    fn iv() -> f64 {
        sweep_order(bessel::iv)
    }

    #[divan::bench]
    fn ive() -> f64 {
        sweep_order(bessel::ive)
    }

    #[divan::bench]
    fn kv() -> f64 {
        sweep_order(bessel::kv)
    }

    #[divan::bench]
    fn kve() -> f64 {
        sweep_order(bessel::kve)
    }

    #[divan::bench]
    fn jv_f32() -> f32 {
        let mut acc = 0.0_f32;
        for &v in &ORDERS {
            for &x in &XS {
                acc += bessel::jv(black_box(v as f32), black_box(x as f32));
            }
        }
        acc
    }
}

mod derivatives {
    use super::*;

    #[divan::bench(args = [1, 4])]
    fn jvp(n: u32) -> f64 {
        sweep_order(|v, x| bessel::jvp(v, x, n))
    }

    #[divan::bench(args = [1, 4])]
    fn kvp(n: u32) -> f64 {
        sweep_order(|v, x| bessel::kvp(v, x, n))
    }
}

mod hankel {
    use super::*;
    use bessel::Complex64;

    #[divan::bench]
    fn hankel1() -> Complex64 {
        let mut acc = Complex64::new(0.0, 0.0);
        for &v in &ORDERS {
            for &x in &XS {
                acc += bessel::hankel1(black_box(v), black_box(x));
            }
        }
        acc
    }

    #[divan::bench]
    fn hankel2e() -> Complex64 {
        let mut acc = Complex64::new(0.0, 0.0);
        for &v in &ORDERS {
            for &x in &XS {
                acc += bessel::hankel2e(black_box(v), black_box(x));
            }
        }
        acc
    }
}

mod complex_argument {
    use super::*;
    use bessel::Complex64;

    const ZS: [Complex64; 5] = [
        Complex64::new(0.5, 0.5),
        Complex64::new(2.0, -1.0),
        Complex64::new(-3.0, 4.0),
        Complex64::new(10.0, 2.0),
        Complex64::new(0.1, 25.0),
    ];

    fn sweep_complex(f: impl Fn(f64, Complex64) -> bessel::Result<Complex64>) -> Complex64 {
        let mut acc = Complex64::new(0.0, 0.0);
        for &v in &ORDERS {
            for &z in &ZS {
                // A failed evaluation must fail the benchmark, not silently
                // make an expensive numerical path look successful.
                let r = f(black_box(v), black_box(z))
                    .unwrap_or_else(|error| panic!("complex benchmark ({v}, {z:?}): {error:?}"));
                acc += r;
            }
        }
        acc
    }

    #[divan::bench]
    fn jv() -> Complex64 {
        sweep_complex(bessel::complex::jv)
    }

    #[divan::bench]
    fn kv() -> Complex64 {
        sweep_complex(bessel::complex::kv)
    }

    #[divan::bench]
    fn hankel1() -> Complex64 {
        sweep_complex(bessel::complex::hankel1)
    }

    // The nearly imaginary K sector used to enter an oscillatory integral
    // with tens of thousands of expensive complex exponentials.
    #[divan::bench]
    fn kv_near_imag_axis() -> Complex64 {
        let z = Complex64::new(0.11, 19.9);
        [2.25, 7.25, 7.99]
            .iter()
            .map(|&v| bessel::complex::kv(black_box(v), black_box(z)).unwrap())
            .sum()
    }

    #[divan::bench]
    fn hankel1_near_real_axis() -> Complex64 {
        let z = Complex64::new(19.9, 0.11);
        [2.25, 7.25, 7.99]
            .iter()
            .map(|&v| bessel::complex::hankel1(black_box(v), black_box(z)).unwrap())
            .sum()
    }
}

mod spherical {
    use super::*;

    #[divan::bench(args = [0, 5, 50])]
    fn spherical_jn(n: i64) -> f64 {
        sweep(|x| bessel::spherical_jn(n, x, false))
    }

    #[divan::bench(args = [0, 5, 50])]
    fn spherical_yn(n: i64) -> f64 {
        sweep(|x| bessel::spherical_yn(n, x, false))
    }

    #[divan::bench(args = [0, 5, 50])]
    fn spherical_in(n: i64) -> f64 {
        sweep(|x| bessel::spherical_in(n, x, false))
    }

    #[divan::bench(args = [0, 5, 50])]
    fn spherical_kn(n: i64) -> f64 {
        sweep(|x| bessel::spherical_kn(n, x, false))
    }

    #[divan::bench]
    fn spherical_jn_derivative() -> f64 {
        sweep(|x| bessel::spherical_jn(5, x, true))
    }

    #[divan::bench(args = [10, 100])]
    fn riccati_jn(n: usize) -> f64 {
        sweep(|x| {
            bessel::riccati_jn(n, x)
                .map(|(v, _)| v.iter().sum())
                .unwrap_or(0.0)
        })
    }
}

mod kelvin {
    use super::*;

    #[divan::bench]
    fn ber() -> f64 {
        sweep(bessel::ber)
    }

    #[divan::bench]
    fn kei() -> f64 {
        sweep(bessel::kei)
    }

    #[divan::bench]
    fn kelvin_all() -> f64 {
        sweep(|x| bessel::kelvin(x).iter().map(|c| c.re + c.im).sum())
    }

    #[divan::bench(args = [5, 20])]
    fn kelvin_zeros(nt: usize) -> usize {
        bessel::kelvin_zeros(black_box(nt))
            .map(|z| z.iter().map(Vec::len).sum())
            .unwrap_or(0)
    }
}

mod struve {
    use super::*;

    #[divan::bench]
    fn struve() -> f64 {
        sweep_order(bessel::struve)
    }

    #[divan::bench]
    fn modstruve() -> f64 {
        sweep_order(bessel::modstruve)
    }
}

mod zeros {
    use super::*;

    #[divan::bench(args = [10, 100])]
    fn jn_zeros(nt: usize) -> usize {
        bessel::jn_zeros(black_box(3), black_box(nt))
            .map(|z| z.len())
            .unwrap_or(0)
    }

    #[divan::bench(args = [10, 100])]
    fn jnyn_zeros(nt: usize) -> usize {
        bessel::jnyn_zeros(black_box(3), black_box(nt))
            .map(|z| z.iter().map(Vec::len).sum())
            .unwrap_or(0)
    }

    #[divan::bench(args = [1, 10, 64, 128])]
    fn jnjnp_zeros(nt: usize) -> usize {
        bessel::jnjnp_zeros(black_box(nt))
            .map(|z| z[0].len())
            .unwrap_or(0)
    }

    #[divan::bench(args = [10, 50])]
    fn y0_zeros(nt: usize) -> usize {
        bessel::y0_zeros(black_box(nt))
            .map(|(z, _)| z.len())
            .unwrap_or(0)
    }
}

mod ratios {
    use super::*;

    #[divan::bench]
    fn iv_ratio() -> f64 {
        sweep_order(|v, x| bessel::iv_ratio(v + 1.0, x))
    }

    #[divan::bench]
    fn iv_ratioinv() -> f64 {
        let mut acc = 0.0;
        for &v in &ORDERS {
            for &r in &[0.05, 0.3, 0.6, 0.9, 0.99] {
                acc += bessel::iv_ratioinv(black_box(v + 1.0), black_box(r));
            }
        }
        acc
    }
}

mod integrals {
    use super::*;

    #[divan::bench]
    fn itj0y0() -> f64 {
        sweep(|x| {
            let (a, b) = bessel::itj0y0(x);
            a + b
        })
    }

    #[divan::bench]
    fn it2i0k0() -> f64 {
        sweep(|x| {
            let (a, b) = bessel::it2i0k0(x);
            a + b
        })
    }

    #[divan::bench]
    fn itstruve0() -> f64 {
        sweep(bessel::itstruve0)
    }

    #[divan::bench]
    fn itmodstruve0() -> f64 {
        sweep(bessel::itmodstruve0)
    }

    #[divan::bench]
    fn besselpoly() -> f64 {
        sweep(|x| bessel::besselpoly(1.0, x, 2.0))
    }
}

mod airy {
    use super::*;

    #[divan::bench]
    fn airy() -> f64 {
        [-20.0, -5.0, -1.0, 0.0, 0.5, 2.0, 5.0, 20.0]
            .iter()
            .map(|&x| {
                bessel::airy(black_box(x))
                    .map(|a: bessel::Airy<f64>| a.0)
                    .unwrap_or(0.0)
            })
            .sum()
    }
}
