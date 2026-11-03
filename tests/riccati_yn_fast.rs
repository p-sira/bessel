//! 90-digit mpmath spherical Neumann values and derivatives, independent of Sci.
use bessel::riccati_yn;

#[test]
fn riccati_yn_sequence_matches_independent_reference() {
    for line in include_str!("fixtures/riccati_yn_fast.csv").lines().skip(1) {
        let fields: Vec<_> = line.split(',').collect();
        let n: usize = fields[0].parse().unwrap();
        let x: f64 = fields[1].parse().unwrap();
        let expected_value: f64 = fields[2].parse().unwrap();
        let expected_derivative: f64 = fields[3].parse().unwrap();
        let (values, derivatives) = riccati_yn(n, x).unwrap();
        for (actual, reference, name) in [
            (values[n], expected_value, "value"),
            (derivatives[n], expected_derivative, "derivative"),
        ] {
            let error = (actual / reference - 1.0).abs();
            assert!(
                error < 3e-11,
                "riccati_yn({n},{x}) {name}: got={actual:e}, ref={reference:e}, rel={error:e}"
            );
        }
    }
}

#[test]
fn riccati_yn_work_limit_and_endpoint_contract() {
    assert!(riccati_yn(4097, 1.0).is_err());
    assert!(riccati_yn(32, 0.0).is_err());
    assert!(riccati_yn(32, f64::NAN).is_err());
    let (values, derivatives) = riccati_yn(4096, 7.5).unwrap();
    assert_eq!(values.len(), 4097);
    assert_eq!(derivatives.len(), 4097);
    assert!(values[0].is_finite() && derivatives[0].is_finite());
}
