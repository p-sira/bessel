use bessel::{self as b, Complex64};

fn parse(s: &str) -> f64 {
    match s.trim().to_ascii_lowercase().as_str() {
        "nan" => f64::NAN,
        "inf" | "infinity" => f64::INFINITY,
        "-inf" | "-infinity" => f64::NEG_INFINITY,
        other => other.parse().unwrap(),
    }
}

fn check(name: &str, order: f64, x: f64, actual: f64, expected: f64, rtol: f64) {
    if expected.is_nan() {
        assert!(
            actual.is_nan(),
            "{name}({order},{x}) expected NaN, got {actual}"
        );
    } else if expected.is_infinite() {
        assert_eq!(actual, expected, "{name}({order},{x})");
    } else {
        let bound = if name.starts_with("hankel") && expected.abs() < 1e-12 {
            1e-12 // scaled Hankel half-integer phases cancel to machine zero
        } else if expected == 0.0 {
            2e-15
        } else {
            rtol * expected.abs() + 2.0 * f64::from_bits(1)
        };
        assert!(
            (actual - expected).abs() <= bound,
            "{name}({order},{x}) expected {expected:e}, got {actual:e}, error {:e}, bound {bound:e}",
            (actual - expected).abs()
        );
    }
}

fn evaluate(name: &str, v: f64, x: f64) -> Complex64 {
    match name {
        "itj0y0" => {
            let (a, c) = b::itj0y0(x);
            return Complex64::new(a, c);
        }
        "it2j0y0" => {
            let (a, c) = b::it2j0y0(x);
            return Complex64::new(a, c);
        }
        "iti0k0" => {
            let (a, c) = b::iti0k0(x);
            return Complex64::new(a, c);
        }
        "it2i0k0" => {
            let (a, c) = b::it2i0k0(x);
            return Complex64::new(a, c);
        }
        "itstruve0" => return Complex64::new(b::itstruve0(x), 0.0),
        "it2struve0" => return Complex64::new(b::it2struve0(x), 0.0),
        "itmodstruve0" => return Complex64::new(b::itmodstruve0(x), 0.0),
        _ => {}
    }
    let (name, _, derivative) = name.partition_colon();
    if name == "besselpoly" {
        let lambda: f64 = derivative.parse().unwrap();
        return Complex64::new(b::besselpoly(x, lambda, v), 0.0);
    }
    let n: u32 = derivative.parse().unwrap_or(0);
    let value = match name {
        "jv" => b::jv(v, x),
        "yv" => b::yv(v, x),
        "iv" => b::iv(v, x),
        "ive" => b::ive(v, x),
        "kv" => b::kv(v, x),
        "kve" => b::kve(v, x),
        "struve" => b::struve(v, x),
        "modstruve" => b::modstruve(v, x),
        "ber" => b::ber(x),
        "bei" => b::bei(x),
        "ker" => b::ker(x),
        "kei" => b::kei(x),
        "berp" => b::berp(x),
        "beip" => b::beip(x),
        "kerp" => b::kerp(x),
        "keip" => b::keip(x),
        "spherical_jn" => b::spherical_jn(v as i64, x, n != 0),
        "spherical_yn" => b::spherical_yn(v as i64, x, n != 0),
        "spherical_in" => b::spherical_in(v as i64, x, n != 0),
        "spherical_kn" => b::spherical_kn(v as i64, x, n != 0),
        "jvp" => b::jvp(v, x, n),
        "yvp" => b::yvp(v, x, n),
        "ivp" => b::ivp(v, x, n),
        "kvp" => b::kvp(v, x, n),
        "iv_ratio" => b::iv_ratio(v, x),
        "iv_ratio_c" => b::iv_ratio_c(v, x),
        "iv_ratioinv" => b::iv_ratioinv(v, x),
        "hankel1" => return b::hankel1(v, x),
        "hankel2" => return b::hankel2(v, x),
        "hankel1e" => return b::hankel1e(v, x),
        "hankel2e" => return b::hankel2e(v, x),
        _ => panic!("unexpected function {name}"),
    };
    Complex64::new(value, 0.0)
}

trait PartitionColon {
    fn partition_colon(&self) -> (&str, &str, &str);
}
impl PartitionColon for str {
    fn partition_colon(&self) -> (&str, &str, &str) {
        if let Some((a, b)) = self.split_once(':') {
            (a, ":", b)
        } else {
            (self, "", "")
        }
    }
}

#[test]
fn independent_scipy_families() {
    let inputs = include_str!("scipy/families_inputs.csv");
    let expected = include_str!("fixtures/families_expected.csv");
    let ins = inputs.lines().skip(1);
    let mut exs = expected.lines().skip(1);
    let mut n = 0;
    for input in ins {
        let reference = exs.next().expect("reference shorter than inputs");
        let args: Vec<_> = input.split(',').collect();
        let name = args[0];
        let v = parse(args[1]);
        let x = parse(args[2]);
        let expected: Vec<_> = reference.split(',').map(parse).collect();
        let actual = evaluate(name, v, x);
        let rtol = if name.starts_with("struve")
            || name == "modstruve"
            || name.contains("spherical")
            || name.starts_with("hankel")
            || name == "ber"
            || name == "bei"
            || name == "ker"
            || name == "kei"
            || name.ends_with('p')
        {
            2e-8
        } else {
            2e-10
        };
        check(name, v, x, actual.re, expected[0], rtol);
        check(name, v, x, actual.im, expected[1], rtol);
        n += 1;
    }
    assert!(exs.next().is_none(), "extra expected rows");
    assert_eq!(n, 1192);
}

#[test]
fn independently_generated_iv_extremes() {
    let fixture = include_str!("fixtures/bessel_iv.csv");
    let mut n = 0usize;
    for line in fixture.lines().skip(1) {
        let row: Vec<_> = line.split(',').collect();
        let v = parse(row[1]);
        let x = parse(row[2]);
        let want_i = parse(row[3]);
        let want_ie = parse(row[4]);
        if row[0] == "f32" {
            let vi = v as f32;
            let xi = x as f32;
            for (name, result, want) in [
                (
                    "iv",
                    b::iv(f64::from(vi), f64::from(xi)) as f32,
                    want_i as f32,
                ),
                (
                    "ive",
                    b::ive(f64::from(vi), f64::from(xi)) as f32,
                    want_ie as f32,
                ),
            ] {
                let tol = 3e-6f32 * want.abs() + 2.0 * f32::from_bits(1);
                assert!(
                    (result - want).abs() <= tol
                        || (result.is_nan() && want.is_nan())
                        || (result.is_infinite() && result == want),
                    "{name}({vi}, {xi}) got {result}, want {want}"
                );
            }
        } else {
            for (name, actual, expected) in
                [("iv", b::iv(v, x), want_i), ("ive", b::ive(v, x), want_ie)]
            {
                if expected.is_nan() {
                    assert!(
                        actual.is_nan(),
                        "{name}({v},{x}) expected NaN, got {actual}"
                    );
                } else if expected.is_infinite() {
                    assert_eq!(actual, expected, "{name}({v},{x})");
                } else if expected == 0.0 && actual == 0.0 {
                    assert_eq!(
                        actual.is_sign_negative(),
                        expected.is_sign_negative(),
                        "{name}({v},{x}) signed zero mismatch"
                    );
                } else {
                    let tolerance = 5e-12 * expected.abs() + 2.0 * f64::from_bits(1);
                    assert!(
                        (actual - expected).abs() <= tolerance,
                        "{name}({v},{x}) got {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
                    );
                }
            }
        }
        n += 1;
    }
    assert_eq!(n, 1259);
}

#[test]
fn zeros_and_recurrence() {
    let zeros = b::jn_zeros(0, 3).unwrap();
    for (got, reference) in zeros.iter().zip([
        2.404_825_557_695_773,
        5.520_078_110_286_311,
        8.653_727_912_911_013,
    ]) {
        assert!((got - reference).abs() < 1e-12);
    }
    let roots = b::ber_zeros(2).unwrap();
    assert!((roots[0] - 2.848_917_820_795_139).abs() < 2e-10);
    assert!((roots[1] - 7.238_829_447_632_408).abs() < 2e-10);
    for x in [0.25, 1.0, 4.0, 10.0] {
        let (values, derivatives) = b::riccati_jn(5, x).unwrap();
        for k in 0..6 {
            assert!((values[k] - x * b::spherical_jn(k as i64, x, false)).abs() <= 1e-14);
            assert!(
                (derivatives[k] - (values[k] / x + x * b::spherical_jn(k as i64, x, true))).abs()
                    <= 1e-14
            );
        }
        let wronskian = b::jv(0.0, x) * b::yv(1.0, x) - b::jv(1.0, x) * b::yv(0.0, x);
        assert!((wronskian + 2.0 / (std::f64::consts::PI * x)).abs() < 1e-12);
    }
}
