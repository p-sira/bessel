//! Real K for small fractional orders, checked against 90-digit mpmath.
use bessel::{kv, kve};

#[test]
fn reduced_order_k_matches_independent_high_precision_references() {
    let data = include_str!("fixtures/real_k_small_order.csv");
    for row in data.lines().skip(1) {
        let fields: Vec<_> = row.split(',').collect();
        let order: f64 = fields[0].parse().unwrap();
        let x: f64 = fields[1].parse().unwrap();
        let expected_scaled: f64 = fields[2].parse().unwrap();
        for v in [order, -order] {
            let got_scaled = kve(v, x);
            assert!(
                got_scaled.is_finite() && (got_scaled / expected_scaled - 1.0).abs() <= 5e-11,
                "kve({v}, {x}) = {got_scaled:.17e}, reference={expected_scaled:.17e}"
            );
            let expected_unscaled = expected_scaled * (-x).exp();
            let got_unscaled = kv(v, x);
            assert!(
                got_unscaled.is_finite() && (got_unscaled / expected_unscaled - 1.0).abs() <= 5e-11,
                "kv({v}, {x}) = {got_unscaled:.17e}, reference={expected_unscaled:.17e}"
            );
        }
    }
}
