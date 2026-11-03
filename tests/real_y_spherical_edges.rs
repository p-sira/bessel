//! Targeted near-integer Neumann and large-argument spherical-I regressions.

use bessel as b;

#[test]
fn neumann_orders_close_to_integer_are_continuous() {
    for (order, x, tolerance) in [
        (0.0, 0.1, 2e-11),
        (0.0, 1.0, 2e-11),
        (1.0, 1.0, 2e-11),
        (2.0, 10.0, 2e-11),
        (7.0, 0.1, 2e-10),
    ] {
        let expected: f64 = b::yv(order, x);
        for delta in [-1e-14, 1e-14] {
            let actual = b::yv(order + delta, x);
            assert!(
                (actual - expected).abs() < tolerance * expected.abs().max(1e-100),
                "Y_{}({x})={actual}, integer limit {expected}",
                order + delta,
            );
        }
    }
}

#[test]
fn neumann_small_argument_is_stable_down_to_subnormals() {
    let euler = 0.577_215_664_901_532_9;
    for x in [f64::from_bits(1), 1e-320, 1e-308, 1e-100, 1e-6] {
        let y0 = (2.0 / core::f64::consts::PI) * (x.ln() - core::f64::consts::LN_2 + euler);
        let near_zero = b::yv(1e-16, x);
        assert!(
            (near_zero / y0 - 1.0).abs() < 2e-10,
            "Y_1e-16({x})={near_zero}, small-x Y0={y0}"
        );
    }
    // Exact half-order Neumann functions remain finite at the least
    // positive subnormal x: the earlier (x/2)^v path rounded x/2 to 0.
    let x = f64::from_bits(1);
    let actual = b::yv(0.5, x);
    let expected = -(2.0 / core::f64::consts::PI).sqrt() / x.sqrt();
    assert!(actual.is_finite());
    assert!((actual / expected - 1.0).abs() < 2e-13);
}

#[test]
fn large_half_integer_neumann_uses_finite_recurrence() {
    // Independent elementary half-order Y formulas, followed by the
    // three-term recurrence; this exercises the Cephes J_-v cutoff at 500.
    let x: f64 = 1000.0;
    let scale = (2.0 / (core::f64::consts::PI * x)).sqrt();
    let (sin, cos) = x.sin_cos();
    let mut previous = -scale * cos; // Y_(1/2)
    let mut current = -scale * (sin + cos / x); // Y_(3/2)
    for k in 1..500 {
        let next = 2.0 * (k as f64 + 0.5) / x * current - previous;
        previous = current;
        current = next;
    }
    let actual = b::yv(500.5, x);
    assert!(actual.is_finite());
    assert!(
        (actual / current - 1.0).abs() < 5e-10,
        "{actual} vs {current}"
    );
}

#[test]
fn negative_large_half_integer_j_matches_elementary_recurrence() {
    // Independent reference from J_{-1/2}(x)=sqrt(2/(pi*x))*cos(x) and
    // J_{-3/2}(x)=-sqrt(2/(pi*x))*(sin(x)+cos(x)/x).
    // Recur by J_{nu-1}=2*nu/x*J_nu-J_{nu+1} to J_{-500.5}.
    // This does not invoke Y or the reflected real-J implementation.
    for x in [100.0_f64, 500.0, 1000.0, 1500.0, 2000.0] {
        let scale = (2.0 / (core::f64::consts::PI * x)).sqrt();
        let (sin, cos) = x.sin_cos();
        let mut upper = scale * cos;
        let mut lower = -scale * (sin + cos / x);
        for k in 1..500 {
            let next = -2.0 * (k as f64 + 0.5) / x * lower - upper;
            upper = lower;
            lower = next;
        }
        let actual = b::jv(-500.5, x);
        assert!(actual.is_finite(), "J_-500.5({x})={actual}");
        assert!(
            (actual / lower - 1.0).abs() < 5e-10,
            "J_-500.5({x})={actual}, independent half-order recurrence={lower}",
        );
    }
}

#[test]
fn spherical_i0_remains_finite_until_true_overflow() {
    for x in [716.5_f64, 717.0] {
        let actual = b::spherical_in(0, x, false);
        let expected = (x - x.ln() - core::f64::consts::LN_2).exp();
        assert!(actual.is_finite());
        assert!(
            (actual / expected - 1.0).abs() < 2e-13,
            "i_0({x})={actual}, reference={expected}"
        );
    }
}
