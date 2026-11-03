//! Numerical transition tests at large cylindrical orders.
//! Values were independently evaluated with SciPy 1.18.1 (AMOS).
use bessel::{self as b, Complex64};

fn reference(got: Complex64, expected: Complex64, rtol: f64, context: &str) {
    assert!(
        got.re.is_finite() && got.im.is_finite(),
        "{context}: nonfinite {got}"
    );
    let error = (got - expected).norm() / expected.norm();
    assert!(
        error < rtol,
        "{context}: got={got}, expected={expected}, relative error={error:e}"
    );
}

#[test]
fn complex_j_and_i_at_integer_orders_and_zero() {
    for order in [-1000.0, -2.0, -1.0, 0.0, 1.0, 2.0, 1000.0] {
        for z in [
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, -0.0),
            Complex64::new(-0.0, 0.0),
        ] {
            let expected = Complex64::new(if order == 0.0 { 1.0 } else { 0.0 }, 0.0);
            for kernel in [
                b::complex::jv,
                b::complex::jve,
                b::complex::iv,
                b::complex::ive,
            ] {
                assert_eq!(kernel(order, z).unwrap(), expected, "order={order}, z={z}");
            }
        }
    }
}

#[test]
fn complex_bessel_large_order_turning_points() {
    for (order, argument, j, y, i, k, ie, ke) in [
        (
            800.0,
            800.0,
            0.048184634724495204,
            -0.08345853039915453,
            1.5910520281239164e183,
            2.7776697308505844e-187,
            5.835779296059744e-165,
            7.57296812391822e160,
        ),
        (
            1000.0,
            1000.0,
            0.044730672947964054,
            -0.07747600152072076,
            2.7234536469112603e229,
            1.2981802514665018e-233,
            1.3824138771101943e-205,
            2.557507414201687e201,
        ),
        (
            1000.5,
            1000.0,
            0.04267543252886601,
            -0.0810319100061619,
            1.752425070166447e229,
            2.0170049079328558e-233,
            8.895237627199093e-206,
            3.973643105949584e201,
        ),
        (
            1200.0,
            1200.0,
            0.04209318681656379,
            -0.07290768825886515,
            4.7579569529615025e275,
            6.192317812380819e-280,
            3.3422753395638886e-246,
            8.815186840414254e241,
        ),
    ] {
        let z = Complex64::new(argument, 0.0);
        for (name, kernel, expected) in [
            ("jv", b::complex::jv as fn(f64, Complex64) -> _, j),
            ("yv", b::complex::yv, y),
            ("iv", b::complex::iv, i),
            ("kv", b::complex::kv, k),
            ("ive", b::complex::ive, ie),
            ("kve", b::complex::kve, ke),
        ] {
            reference(
                kernel(order, z).unwrap(),
                Complex64::new(expected, 0.0),
                2e-8,
                &format!("{name}({order},{argument})"),
            );
        }
        reference(
            b::complex::hankel1(order, z).unwrap(),
            Complex64::new(j, y),
            2e-8,
            &format!("H1({order},{argument})"),
        );
    }
    reference(
        Complex64::new(b::spherical_jn(1000, 1000.0, false), 0.0),
        Complex64::new(0.001691367066787982, 0.0),
        2e-8,
        "spherical_jn(1000,1000)",
    );
}

#[test]
fn cylindrical_order_above_argument_stays_finite() {
    let z = Complex64::new(1000.0, 0.0);
    reference(
        b::complex::jv(1100.0, z).unwrap(),
        Complex64::new(2.42614418358924e-15, 0.0),
        2e-8,
        "J_1100(1000)",
    );
    reference(
        b::complex::yv(1100.0, z).unwrap(),
        Complex64::new(-286324273667.7516, 0.0),
        2e-8,
        "Y_1100(1000)",
    );
}

#[test]
fn scaled_real_k_resolves_the_narrow_large_x_saddle() {
    // Leading uniform saddle approximation with 70-digit arithmetic:
    // sqrt(pi/(2*hypot(x,v))) *
    // exp(v*asinh(v/x) - v²/(hypot(x,v)+x)).
    // Its relative first correction is O(1/x). At x >= 1e8 the
    // approximation is accurate enough to expose quadrature underresolution.
    for (v, x, expected, rtol) in [
        (10_000.0, 1e10, 1.259596400571861e-5, 2e-9),
        (1e7, 1e15, 4.1665314334512956e-8, 2e-9),
        (10_000.0, 1e8, 2.0663656710343466e-4, 5e-9),
        (100_000.0, 1e8, 6.498036037082404e17, 5e-9),
        (1000.0, 1e7, 4.1665314228613617e-4, 5e-8),
    ] {
        let got: f64 = b::kve(v, x);
        let relative = (got / expected - 1.0).abs();
        assert!(
            got.is_finite() && relative <= rtol,
            "kve({v},{x})={got:e}, saddle reference={expected:e}, rel={relative:e}"
        );
    }
}

#[test]
fn real_k_large_x_transition_matches_independent_scipy() {
    // Check both sides of the saddle-expansion transition at q ~= 1e5.
    for (v, x, expected) in [
        (100.0, 99_999.0, 0.004166548082339202),
        (100.0, 100_001.0, 0.004166502250709191),
        (1000.0, 99_500.0, 0.6046483170575085),
        (1000.0, 100_000.0, 0.588169978527613),
        (1000.0, 101_000.0, 0.5569846045495711),
    ] {
        let actual: f64 = b::kve(v, x);
        assert!(
            (actual / expected - 1.0).abs() <= 2e-9,
            "kve({v},{x}) = {actual:e}, SciPy reference = {expected:e}"
        );
    }
}
