//! Verify the positive-real 3/2 Hankel closed forms against 80-digit mpmath.
use bessel::{hankel1, hankel1e, hankel2, hankel2e};
use num_complex::Complex64;

#[test]
fn real_half_integer_hankel_matches_independent_reference() {
    for line in include_str!("fixtures/hankel_three_halves.csv")
        .lines()
        .skip(1)
    {
        let columns: Vec<_> = line.split(',').collect();
        let x: f64 = columns[1].parse().unwrap();
        let expected = Complex64::new(columns[2].parse().unwrap(), columns[3].parse().unwrap());
        let actual = match columns[0] {
            "hankel1" => hankel1(1.5, x),
            "hankel2" => hankel2(1.5, x),
            "hankel1e" => hankel1e(1.5, x),
            "hankel2e" => hankel2e(1.5, x),
            _ => panic!("unexpected reference kind"),
        };
        let relative_error = (actual - expected).norm() / expected.norm();
        assert!(
            relative_error <= 4e-14,
            "{}({x}) got {actual:?}, expected {expected:?}, rel={relative_error:e}",
            columns[0]
        );
    }
}
