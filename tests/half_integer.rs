//! 100-digit mpmath references for low and moderate half-integer orders,
//! including the 32-point Sci derivative grid.

#[test]
fn half_integer_and_derivatives_match_independent_references() {
    let mut count = 0;
    for line in include_str!("fixtures/half_integer_reference.csv")
        .lines()
        .skip(1)
    {
        let mut fields = line.split(',');
        let api = fields.next().unwrap();
        let v: f64 = fields.next().unwrap().parse().unwrap();
        let x: f64 = fields.next().unwrap().parse().unwrap();
        let expected: f64 = fields.next().unwrap().parse().unwrap();
        let actual = match api {
            "jv" => bessel::jv(v, x),
            "yv" => bessel::yv(v, x),
            "jvp" => bessel::jvp(v, x, 2),
            "yvp" => bessel::yvp(v, x, 2),
            "jvp1" => bessel::jvp(v, x, 1),
            "yvp1" => bessel::yvp(v, x, 1),
            "jvp3" => bessel::jvp(v, x, 3),
            "yvp3" => bessel::yvp(v, x, 3),
            "jvp4" => bessel::jvp(v, x, 4),
            "yvp4" => bessel::yvp(v, x, 4),
            _ => panic!("unexpected API {api}"),
        };
        let tolerance = 5e-12 * expected.abs() + 1e-14;
        assert!(
            actual.is_finite() && (actual - expected).abs() <= tolerance,
            "{api}({v}, {x}): got={actual:e}, expected={expected:e}, bound={tolerance:e}"
        );
        count += 1;
    }
    assert_eq!(count, 496);
}
