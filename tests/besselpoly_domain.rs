use bessel::besselpoly;

#[test]
fn negative_integer_orders_use_bessel_j_reflection_before_domain_validation() {
    for n in 1..=5 {
        for a in [0.0, 0.01, 0.5, 2.0] {
            let positive = besselpoly(a, 0.0, n as f64);
            let negative = besselpoly(a, 0.0, -(n as f64));
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert!(positive.is_finite(), "n={n}, a={a}");
            if a == 0.0 {
                assert_eq!(negative.to_bits(), (sign * 0.0_f64).to_bits());
            }
            assert!(
                (negative - sign * positive).abs() <= 3e-14 * positive.abs().max(1.0),
                "n={n}, a={a}, positive={positive}, negative={negative}"
            );
        }
    }
}

#[test]
fn singular_noninteger_orders_at_zero_have_the_continuous_limit() {
    assert_eq!(besselpoly(0.0, 0.0, -0.5), f64::INFINITY);
    assert_eq!(besselpoly(0.0, 3.0, -1.5), f64::NEG_INFINITY);
    assert!(besselpoly(0.0, 0.0, -1.5).is_nan());
}
