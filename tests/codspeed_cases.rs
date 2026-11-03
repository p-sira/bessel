//! Ensure benchmark sweeps exercise successful, finite numerical evaluations.
//! In particular, a failed Hankel evaluation must never be silently omitted
//! from a performance measurement.

use bessel::{self as b, Complex64};

const XS: [f64; 8] = [0.05, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 100.0];
const ORDERS: [f64; 4] = [0.0, 0.5, 2.0, 7.25];
const ZS: [Complex64; 5] = [
    Complex64::new(0.5, 0.5),
    Complex64::new(2.0, -1.0),
    Complex64::new(-3.0, 4.0),
    Complex64::new(10.0, 2.0),
    Complex64::new(0.1, 25.0),
];

type ComplexKernel = fn(f64, Complex64) -> b::Result<Complex64>;

#[test]
fn all_codspeed_complex_inputs_must_succeed() {
    let functions: [(&str, ComplexKernel); 3] = [
        ("jv", b::complex::jv),
        ("kv", b::complex::kv),
        ("hankel1", b::complex::hankel1),
    ];
    for (name, f) in functions {
        for v in ORDERS {
            for z in ZS {
                let value = f(v, z).unwrap_or_else(|e| panic!("{name}({v}, {z:?}) returned {e:?}"));
                assert!(
                    value.re.is_finite() && value.im.is_finite(),
                    "{name}({v}, {z:?}) was nonfinite: {value:?}"
                );
            }
        }
    }
}

#[test]
fn codspeed_hankel_turning_case_matches_mpmath() {
    // Independently evaluated with mpmath 1.4.1 at 75 decimal digits.
    // Previously this input triggered thousands of quadrature panels and
    // ultimately returned a convergence error, which the benchmark ignored.
    let v = 7.25;
    let z = Complex64::new(10.0, 2.0);
    let expected = Complex64::new(0.066_079_669_196_886_18, 0.016_689_765_055_237_67);
    let actual = b::complex::hankel1(v, z).unwrap();
    assert!(
        (actual - expected).norm() <= 3e-11 * expected.norm(),
        "hankel1({v}, {z:?}): actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn codspeed_complex_k_oscillatory_case_matches_mpmath() {
    // This sector is particularly sensitive to numerical integration error:
    // the original integration path erred by about 5e-7 relative, even
    // though it reported successful convergence.
    // Reference: mpmath 1.4.1 at 85 decimal digits.
    let v = 7.25;
    let z = Complex64::new(0.1, 25.0);
    let expected = Complex64::new(-0.030_745_806_021_712_767, -0.230_755_436_733_858_06);
    let actual = b::complex::kv(v, z).unwrap();
    assert!(
        (actual - expected).norm() <= 5e-10 * expected.norm(),
        "kv({v}, {z:?}): actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn codspeed_real_hankel_and_k_are_finite() {
    for v in ORDERS {
        for x in XS {
            for (name, value) in [
                ("hankel1", b::hankel1(v, x)),
                ("hankel2e", b::hankel2e(v, x)),
            ] {
                assert!(
                    value.re.is_finite() && value.im.is_finite(),
                    "{name}({v}, {x}) = {value:?}"
                );
            }
            for (name, value) in [("kv", b::kv(v, x)), ("kve", b::kve(v, x))] {
                assert!(value.is_finite(), "{name}({v}, {x}) = {value:?}");
            }
        }
    }
}
