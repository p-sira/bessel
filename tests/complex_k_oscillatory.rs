//! Independent 90-digit mpmath 1.4.1 references for the near-imaginary-axis
//! sector where the original real-axis K quadrature failed to converge.
#![allow(clippy::excessive_precision)]

use bessel::{self as b, Complex64};

fn close(name: &str, actual: Complex64, expected: Complex64, tolerance: f64) {
    let relative = (actual - expected).norm() / expected.norm();
    assert!(
        relative <= tolerance,
        "{name}: actual={actual:?}, expected={expected:?}, relative={relative:e}"
    );
}

#[test]
fn oscillatory_complex_k_matches_high_precision_references() {
    for (v, z, reference) in [
        (
            0.0,
            Complex64::new(0.11, 15.0),
            Complex64::new(-0.2891737404160174, 0.018955944796848414),
        ),
        (
            2.25,
            Complex64::new(0.11, 15.0),
            Complex64::new(-0.28384067278879496, 0.06766497390725665),
        ),
        (
            7.25,
            Complex64::new(0.11, 15.0),
            Complex64::new(0.08563178723759583, 0.3017588562838012),
        ),
        (
            7.99,
            Complex64::new(0.11, 15.0),
            Complex64::new(0.1983193038910039, 0.2510133035641146),
        ),
        (
            2.25,
            Complex64::new(0.11, 19.9),
            Complex64::new(-0.09443639925686441, -0.2343060042078653),
        ),
        (
            7.25,
            Complex64::new(0.11, 19.9),
            Complex64::new(-0.26262461143971035, 0.004856561289979073),
        ),
        (
            7.99,
            Complex64::new(0.11, 19.9),
            Complex64::new(-0.2527666959047173, 0.08057767251111397),
        ),
        (
            2.25,
            Complex64::new(0.09, 19.9),
            Complex64::new(-0.09645454364802216, -0.23896175639131854),
        ),
        (
            7.99,
            Complex64::new(0.11, 17.0),
            Complex64::new(0.19702730501186165, -0.21728323235888325),
        ),
        (
            7.25,
            Complex64::new(0.11, 13.0),
            Complex64::new(-0.34428743834823766, 0.045454728609264294),
        ),
        (
            7.99,
            Complex64::new(0.11, 16.0),
            Complex64::new(0.30533648183130513, 0.013032928316713277),
        ),
        (
            8.25,
            Complex64::new(1.0, 8.0),
            Complex64::new(0.43111502338681396, -0.2204503118023863),
        ),
        (
            8.25,
            Complex64::new(3.0, 10.0),
            Complex64::new(0.027093662614978503, -0.05881766838411253),
        ),
        (
            8.25,
            Complex64::new(0.11, 10.0),
            Complex64::new(-0.10307389492061148, -0.4666605229158861),
        ),
    ] {
        close(
            &format!("K_{v}({z:?})"),
            b::complex::kv(v, z).unwrap(),
            reference,
            3e-11,
        );
        close(
            &format!("K_{v}({z:?}) lower half-plane"),
            b::complex::kv(v, z.conj()).unwrap(),
            reference.conj(),
            3e-11,
        );
    }
}

#[test]
fn scaled_k_and_hankel_recover_near_real_axis() {
    for (v, z, expected) in [
        (
            2.25,
            Complex64::new(0.11, 15.0),
            Complex64::new(0.1915852631604792, -0.2634219790918179),
        ),
        (
            7.99,
            Complex64::new(0.11, 15.0),
            Complex64::new(-0.3503903956217087, -0.06890465373396824),
        ),
        (
            2.25,
            Complex64::new(0.11, 19.9),
            Complex64::new(0.17452090728249298, -0.22150395238127494),
        ),
    ] {
        close("scaled K", b::complex::kve(v, z).unwrap(), expected, 3e-11);
        close(
            "scaled K conjugate",
            b::complex::kve(v, z.conj()).unwrap(),
            expected.conj(),
            3e-11,
        );
    }

    for (v, z, h1, h1e) in [
        (
            2.25,
            Complex64::new(19.9, 0.11),
            Complex64::new(-0.1608163723908116, 0.0015388179073896974),
            Complex64::new(-0.08776230616630475, 0.15660985629833099),
        ),
        (
            7.25,
            Complex64::new(19.9, 0.11),
            Complex64::new(-0.15564845971195257, 0.061125181232764035),
            Complex64::new(-0.027182822473096233, 0.1846748556931634),
        ),
        (
            2.25,
            Complex64::new(15.0, 0.11),
            Complex64::new(-0.02935252499659557, -0.18342852452261052),
            Complex64::new(-0.10825952792155332, 0.17685865808641623),
        ),
    ] {
        close("H1", b::complex::hankel1(v, z).unwrap(), h1, 3e-11);
        close("H1e", b::complex::hankel1e(v, z).unwrap(), h1e, 3e-11);
        close(
            "H2 conjugate",
            b::complex::hankel2(v, z.conj()).unwrap(),
            h1.conj(),
            3e-11,
        );
        close(
            "H2e conjugate",
            b::complex::hankel2e(v, z.conj()).unwrap(),
            h1e.conj(),
            3e-11,
        );
    }
}

#[test]
fn tiny_complex_half_order_k_stays_representable() {
    for (z, expected) in [
        (
            Complex64::new(1e-200, 1e-200),
            Complex64::new(9.736834439223663e99, -4.033128879307871e99),
        ),
        (
            Complex64::new(1e-200, 2e-200),
            Complex64::new(7.129658962611928e99, -4.406371567089487e99),
        ),
    ] {
        for (name, kernel) in [
            ("K_1/2", b::complex::kv as fn(f64, Complex64) -> _),
            ("Kve_1/2", b::complex::kve),
        ] {
            close(name, kernel(0.5, z).unwrap(), expected, 8e-15);
            close(name, kernel(0.5, z.conj()).unwrap(), expected.conj(), 8e-15);
        }
    }
    let subnormal = Complex64::new(1e-320, 1e-320);
    for kernel in [b::complex::kv as fn(f64, Complex64) -> _, b::complex::kve] {
        let result = kernel(0.5, subnormal).unwrap();
        assert!(result.re.is_finite() && result.im.is_finite());
        assert!(result.norm() > 1e159);
    }
}
