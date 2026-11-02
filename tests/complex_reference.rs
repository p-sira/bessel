//! Independent SciPy principal-branch complex Bessel/Hankel references.
//!
//! Twelve complex APIs × twelve points include both sides of the negative-real
//! branch cut. A further 96 SciPy Hankel derivative references exercise the
//! real-argument `h1vp` and `h2vp` interfaces.
use bessel::{self as b, Complex64};

fn function(name: &str, v: f64, z: Complex64) -> Complex64 {
    let result = match name {
        "jv" => b::complex::jv(v, z),
        "yv" => b::complex::yv(v, z),
        "iv" => b::complex::iv(v, z),
        "kv" => b::complex::kv(v, z),
        "jve" => b::complex::jve(v, z),
        "yve" => b::complex::yve(v, z),
        "ive" => b::complex::ive(v, z),
        "kve" => b::complex::kve(v, z),
        "hankel1" => b::complex::hankel1(v, z),
        "hankel2" => b::complex::hankel2(v, z),
        "hankel1e" => b::complex::hankel1e(v, z),
        "hankel2e" => b::complex::hankel2e(v, z),
        _ => panic!("unknown complex API {name}"),
    };
    result.unwrap_or_else(|error| panic!("{name}({v},{z}) returned {error:?}"))
}

fn relative_complex(got: Complex64, expected: Complex64, tol: f64, label: &str) {
    assert!(expected.re.is_finite() && expected.im.is_finite());
    assert!(
        got.re.is_finite() && got.im.is_finite(),
        "{label}: nonfinite actual {got:?}"
    );
    let scale = expected.norm();
    assert!(scale > 0.0, "{label}: reference must be nonzero");
    let rel = (got - expected).norm() / scale;
    assert!(
        rel <= tol,
        "{label}: actual={got:?} expected={expected:?} rel={rel:e} tolerance={tol:e}"
    );
}

#[test]
fn twelve_complex_apis_match_independent_scipy_branch_references() {
    let inputs = include_str!("fixtures/complex_inputs.csv");
    let expected = include_str!("fixtures/complex_expected.csv");
    let mut count = 0;
    let mut checked_apis = std::collections::BTreeSet::new();
    for (input, reference) in inputs.lines().skip(1).zip(expected.lines().skip(1)) {
        let fields: Vec<_> = input.split(',').collect();
        let values: Vec<_> = reference.split(',').collect();
        assert_eq!(fields.len(), 5);
        assert_eq!(values.len(), 3);
        assert_eq!(fields[0], values[0], "changed fixture ID");
        let api = fields[1];
        let v = fields[2].parse::<f64>().unwrap();
        let z = Complex64::new(fields[3].parse().unwrap(), fields[4].parse().unwrap());
        let want = Complex64::new(values[1].parse().unwrap(), values[2].parse().unwrap());
        let actual = function(api, v, z);
        relative_complex(
            actual,
            want,
            4e-10,
            &format!("case {} {api}({v},{z})", fields[0]),
        );
        checked_apis.insert(api);
        count += 1;
    }
    assert_eq!(count, 144);
    assert_eq!(checked_apis.len(), 12);
    assert_eq!(inputs.lines().count(), expected.lines().count());
}

#[test]
fn principal_branches_obey_conjugation_relations_away_from_cuts() {
    // Ordinary Bessel J, Y, I and K with real order have real power-series
    // coefficients on their principal branches. The Hankel conjugation pairs
    // H1 with H2; their exponential scaling likewise forms a conjugate pair.
    for &(v, z) in &[
        (0.5, Complex64::new(0.4, 0.8)),
        (1.3, Complex64::new(-1.4, 0.6)),
        (2.5, Complex64::new(-2.0, 1e-8)),
        (0.0, Complex64::new(3.2, -4.1)),
    ] {
        for name in ["jv", "yv", "iv", "kv", "jve", "yve", "ive", "kve"] {
            let upper = function(name, v, z);
            let lower = function(name, v, z.conj());
            relative_complex(
                lower,
                upper.conj(),
                5e-11,
                &format!("{name} conjugacy {v} {z}"),
            );
        }
        for (upper, lower) in [("hankel1", "hankel2"), ("hankel1e", "hankel2e")] {
            let actual = function(lower, v, z.conj());
            let expected = function(upper, v, z).conj();
            relative_complex(
                actual,
                expected,
                5e-11,
                &format!("{upper}/{lower} conjugacy {z}"),
            );
        }
    }
}

#[test]
fn branch_cut_j_and_k_have_correct_rim_limits() {
    // DLMF 10.11.1: J_v(x exp(±iπ)) = exp(±iπv) J_v(x).
    // DLMF 10.34.2: K_0(-x±i0) = K_0(x) ∓ iπ I_0(x).
    let x = 2.0_f64;
    let eps = 1e-10_f64;
    let v = 0.5_f64;
    let j_positive = b::jv(v, x);
    for &side in &[-1.0, 1.0] {
        let z = Complex64::new(-x, side * eps);
        let phase = Complex64::new(0.0, std::f64::consts::PI * v * side).exp();
        relative_complex(
            function("jv", v, z),
            phase * j_positive,
            2e-9,
            &format!("J_v negative cut rim {side}"),
        );
        let expected_k = Complex64::new(b::k0(x), -side * std::f64::consts::PI * b::i0(x));
        relative_complex(
            function("kv", 0.0, z),
            expected_k,
            2e-9,
            &format!("K_0 negative cut rim {side}"),
        );
    }
}

#[test]
fn hankel_derivatives_match_scipy_at_four_orders() {
    let inputs = include_str!("fixtures/hankel_derivatives_inputs.csv");
    let expected = include_str!("fixtures/hankel_derivatives_expected.csv");
    let mut count = 0;
    for (input, reference) in inputs.lines().skip(1).zip(expected.lines().skip(1)) {
        let fields: Vec<_> = input.split(',').collect();
        let values: Vec<_> = reference.split(',').collect();
        assert_eq!(fields[0], values[0], "changed derivative fixture ID");
        let name = fields[1];
        let v = fields[2].parse::<f64>().unwrap();
        let x = fields[3].parse::<f64>().unwrap();
        let order = fields[4].parse::<u32>().unwrap();
        let want = Complex64::new(values[1].parse().unwrap(), values[2].parse().unwrap());
        let actual = match name {
            "h1vp" => b::h1vp(v, x, order),
            "h2vp" => b::h2vp(v, x, order),
            _ => panic!("unexpected derivative {name}"),
        };
        relative_complex(
            actual,
            want,
            4e-10,
            &format!("case {} {name}({v},{x},{order})", fields[0]),
        );
        count += 1;
    }
    assert_eq!(count, 96);
    assert_eq!(inputs.lines().count(), expected.lines().count());
}
