//! First derivatives checked against elementary closed forms and the
//! adjacent-order identity over both Miller and upward-recurrence regimes.
use bessel::{spherical_in, spherical_jn};

#[test]
fn first_order_derivatives_match_elementary_formulas() {
    // j_1(x) = sin(x)/x^2 - cos(x)/x and
    // i_1(x) = cosh(x)/x - sinh(x)/x^2.
    for x in [0.75_f64, 1.0, 2.0, 5.0, 10.0, 20.0] {
        let (sin, cos) = x.sin_cos();
        let (sinh, cosh) = (x.sinh(), x.cosh());
        let j_prime = sin / x - 2.0 * sin / x.powi(3) + 2.0 * cos / (x * x);
        let i_prime = sinh / x - 2.0 * cosh / (x * x) + 2.0 * sinh / x.powi(3);
        let actual_j = spherical_jn(1, x, true);
        let actual_i = spherical_in(1, x, true);
        assert!(
            (actual_j - j_prime).abs() <= 2e-13 * j_prime.abs().max(1e-15),
            "j_1'({x}) = {actual_j:e}, closed form = {j_prime:e}"
        );
        assert!(
            (actual_i - i_prime).abs() <= 2e-13 * i_prime.abs().max(1e-15),
            "i_1'({x}) = {actual_i:e}, closed form = {i_prime:e}"
        );
    }
}

#[test]
fn higher_derivatives_match_adjacent_order_identity() {
    // Cross-check all optimized regimes, including x≈n and the n=256
    // Miller limit. Values use the public value-only kernels.
    for n in [1_i64, 2, 5, 20, 50, 128, 256] {
        for x in [0.5_f64, 0.75, 1.0, 5.0, 20.0, 50.0, 100.0, 256.0] {
            for (name, kernel) in [
                ("j", spherical_jn::<f64> as fn(i64, f64, bool) -> f64),
                ("i", spherical_in::<f64> as fn(i64, f64, bool) -> f64),
            ] {
                let previous = kernel(n - 1, x, false);
                let current = kernel(n, x, false);
                let derivative = kernel(n, x, true);
                let correction = (n as f64 + 1.0) / x * current;
                let expected = previous - correction;
                if expected.is_finite() && derivative.is_finite() {
                    let scale = previous.abs().max(correction.abs()).max(1e-290);
                    assert!(
                        (derivative - expected).abs() <= 5e-13 * scale,
                        "{name}_{n}'({x})={derivative:e}, adjacent orders={expected:e}"
                    );
                }
            }
        }
    }
}
