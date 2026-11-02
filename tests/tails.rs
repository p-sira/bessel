//! Nonzero tiny-value and subnormal tests with pinned independent SciPy oracles.
//!
//! Regenerate: `uv run --offline --script tests/scipy/generate_tails.py`.
//! The input files and expected SciPy values carry matching immutable case IDs.
use bessel as b;

fn evaluate(api: &str, order: f64, x: f64) -> f64 {
    match api {
        "jv" => b::jv(order, x),
        "spherical_jn" => b::spherical_jn(order as i64, x, false),
        "kv" => b::kv(order, x),
        "k0" => b::k0(x),
        "ker" => b::ker(x),
        "kei" => b::kei(x),
        "iv" => b::iv(order, x),
        "ive" => b::ive(order, x),
        _ => panic!("unexpected API {api}"),
    }
}

#[test]
fn tiny_nonzero_values_match_scipy_by_relative_error() {
    let inputs = include_str!("fixtures/tails_inputs.csv");
    let expected = include_str!("fixtures/tails_expected.csv");
    let mut count = 0;
    for (input, reference) in inputs.lines().skip(1).zip(expected.lines().skip(1)) {
        let fields: Vec<_> = input.split(',').collect();
        let values: Vec<_> = reference.split(',').collect();
        assert_eq!(fields.len(), 4);
        assert_eq!(values.len(), 2);
        assert_eq!(fields[0], values[0], "fixture case IDs changed");
        let (api, order, x) = (
            fields[1],
            fields[2].parse::<f64>().unwrap(),
            fields[3].parse::<f64>().unwrap(),
        );
        let want = values[1].parse::<f64>().unwrap();
        assert!(
            want.is_finite() && want != 0.0,
            "{api} reference must be nonzero"
        );
        let got = evaluate(api, order, x);
        assert!(
            got.is_finite() && got != 0.0,
            "case {} {api}({order},{x}) returned zero/NaN",
            fields[0]
        );
        let relative = ((got - want) / want).abs();
        assert!(
            relative < 2e-10,
            "case {} {api}({order},{x}): actual={got:e}, SciPy={want:e}, rel={relative:e}",
            fields[0]
        );
        count += 1;
    }
    assert_eq!(count, 14);
    assert_eq!(inputs.lines().count(), expected.lines().count());
}

#[test]
fn near_origin_power_series_has_correct_scale_and_sign() {
    // Independent first two Taylor coefficients from DLMF 10.2.2 and 10.53.1.
    let x = 0.01_f64;
    let factorial_8 = (1..=8).product::<u64>() as f64;
    let j8 = (x / 2.0).powi(8) / factorial_8 * (1.0 - x * x / 36.0);
    let got_j8 = b::jv(8.0, x);
    assert!(got_j8 > 0.0 && ((got_j8 - j8) / j8).abs() < 2e-10);

    let x = 0.1_f64;
    let double_factorial_21: f64 = (1..=21).step_by(2).map(f64::from).product();
    let j10 = x.powi(10) / double_factorial_21 * (1.0 - x * x / 46.0);
    let got_j10 = b::spherical_jn(10, x, false);
    assert!(got_j10 > 0.0 && ((got_j10 - j10) / j10).abs() < 5e-8);
}

#[test]
fn nonzero_binary32_subnormals_are_checked_in_ulps() {
    // The Bessel API evaluates in binary64. We quantize its inputs to f32,
    // then compare its f32-rounded output to SciPy computed at those identical
    // quantized arguments; no absolute floor that swallows subnormal results.
    const ULP: f64 = f32::from_bits(1) as f64;
    let inputs = include_str!("fixtures/subnormal_inputs.csv");
    let expected = include_str!("fixtures/subnormal_expected.csv");
    let mut count = 0;
    let mut especially_tiny = 0;
    for (input, reference) in inputs.lines().skip(1).zip(expected.lines().skip(1)) {
        let fields: Vec<_> = input.split(',').collect();
        let values: Vec<_> = reference.split(',').collect();
        assert_eq!(fields[0], values[0], "fixture IDs changed");
        let api = fields[1];
        let v = fields[2].parse::<f32>().unwrap();
        let x = fields[3].parse::<f32>().unwrap();
        let want = values[1].parse::<f32>().unwrap();
        assert!(
            want.is_subnormal() && want > 0.0,
            "expected finite positive subnormal"
        );
        let raw = match api {
            "iv" => b::iv(f64::from(v), f64::from(x)),
            "ive" => b::ive(f64::from(v), f64::from(x)),
            _ => panic!("unexpected API {api}"),
        };
        let got = raw as f32;
        assert!(
            got.is_subnormal() && got > 0.0,
            "case {} {api}({v},{x}) returned {got:e}, wanted {want:e}",
            fields[0]
        );
        let ulp_distance = got.to_bits().abs_diff(want.to_bits());
        // At least the last few bits are comparable even in the lowest part
        // of the f32 subnormal interval. Zero always fails the guard above.
        assert!(
            ulp_distance <= 2 && ((f64::from(got) - f64::from(want)).abs() <= 2.0 * ULP),
            "case {} {api}({v},{x}): got={got:e}, want={want:e}, ulps={ulp_distance}",
            fields[0]
        );
        if want.to_bits() <= 64 {
            especially_tiny += 1;
        }
        count += 1;
    }
    assert_eq!(count, 34);
    assert!(especially_tiny >= 3);
    assert_eq!(inputs.lines().count(), expected.lines().count());
}
