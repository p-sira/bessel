//! High-precision references for the moderate-to-high-order real K quadrature.
use bessel::kve;

#[test]
fn fractional_order_quadrature_matches_mpmath() {
    // 85-digit mpmath values from tests/fixtures/native_reference.csv.
    // These inputs bypass the reduced-order recurrence and uniform asymptotic
    // shortcuts, exercising both broad and narrow quadrature peaks.
    let cases: [(f64, f64, f64); 4] = [
        (75.25, 75.25, 224_677_821_374_127.03),
        (100.25, 705.0, 57.817_055_297_498_68),
        (150.25, 10_000.0, 0.038_746_460_620_395_454),
        (75.25, 2831.28125, 0.064_009_073_969_170_15),
    ];
    for (v, x, expected) in cases {
        let actual: f64 = kve(v, x);
        assert!(
            actual.is_finite() && (actual / expected - 1.0).abs() <= 2e-10,
            "Kve({v}, {x}) = {actual:.17e}; expected {expected:.17e}"
        );
    }
}
