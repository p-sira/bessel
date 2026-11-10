//! Independent audit regressions; regenerate with
//! `uv run --script tests/fixtures/accuracy_audit_generate.py`.
//! References use exact binary64 inputs at 90 and 130 decimal digits.
use bessel as b;

#[test]
fn high_precision_audit_references() {
    let mut count = 0;
    for line in include_str!("fixtures/accuracy_audit.csv").lines().skip(1) {
        let fields: Vec<_> = line.split(',').collect();
        let api = fields[0];
        let v: f64 = fields[1].parse().unwrap();
        let x: f64 = fields[2].parse().unwrap();
        let expected: f64 = fields[3].parse().unwrap();
        let actual = match api {
            "jv" => b::jv(v, x),
            "yv" => b::yv(v, x),
            "struve" => b::struve(v, x),
            "si" | "sid" => b::spherical_in(v as i64, x, api.ends_with('d')),
            "sjd" => b::spherical_jn(v as i64, x, true),
            _ => panic!("unknown audit API: {api}"),
        };
        let tolerance = 2e-11 * expected.abs() + 4.0 * f64::from_bits(1);
        assert!(
            actual.is_finite() && (actual - expected).abs() <= tolerance,
            "{api}({v}, {x}): actual={actual:e}, expected={expected:e}, bound={tolerance:e}"
        );
        count += 1;
    }
    assert_eq!(count, 178);
}

#[test]
fn reflected_large_orders_preserve_overflow_and_half_integer_limits() {
    // Exercise overflow without retrying the mutually dependent J/Y paths.
    assert_eq!(b::yv(300.25, 0.1), f64::NEG_INFINITY);
    assert_eq!(b::jv(-300.25, 0.1), f64::INFINITY);
    assert_eq!(b::yv(-300.25, 0.1), f64::NEG_INFINITY);
    // cos(pi*v) is exactly zero: the overflowing Y contribution vanishes.
    assert_eq!(b::yv(-300.5, 0.1), 0.0);
    assert_eq!(b::jv(-300.5, 0.1), f64::INFINITY);
}
