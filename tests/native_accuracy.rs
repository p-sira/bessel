//! Adversarial numerical references for the native complex Bessel/Hankel and
//! Airy kernels, plus real modified K and scaled K. All 2,299 fixtures are
//! static, so `cargo test` needs no Python.
//!
//! Regenerate with `uv run --script tests/fixtures/native_generate.py`.
//! Most complex values are independently computed with 85-digit mpmath,
//! cross-checked against SciPy AMOS (maximum relative difference ~5.5e-14).
//! Extreme scaled inputs use SciPy AMOS because J +/- iY loses hundreds of
//! digits to cancellation in the unscaled Hankel representation.
//! SciPy itself ignores -0 on the negative-real cut; lower-bank references
//! are generated via Schwarz reflection and the Hankel-kind exchange.
//! Real K references use 85-digit mpmath and 1,416 independent SciPy checks;
//! 139 finite mpmath values exceed SciPy's usable overflow/underflow range.
//! This includes K at subnormal x=2^-1074, where SciPy and intermediate
//! exponentiation in naive quadrature can return infinity.
//!
//! Per-case tolerance is a relative bound on the complex magnitude with a
//! tiny absolute floor for nonzero subnormal references. It never scales with
//! the (possibly enormous) input or with a result from the implementation.
use bessel::{self as b, Complex64};

const REFERENCES: &str = include_str!("fixtures/native_reference.csv");

fn complex_value(api: &str, order: f64, z: Complex64) -> b::Result<Complex64> {
    match api {
        "jv" => b::complex::jv(order, z),
        "yv" => b::complex::yv(order, z),
        "iv" => b::complex::iv(order, z),
        "kv" => b::complex::kv(order, z),
        "jve" => b::complex::jve(order, z),
        "yve" => b::complex::yve(order, z),
        "ive" => b::complex::ive(order, z),
        "kve" => b::complex::kve(order, z),
        "hankel1" => b::complex::hankel1(order, z),
        "hankel2" => b::complex::hankel2(order, z),
        "hankel1e" => b::complex::hankel1e(order, z),
        "hankel2e" => b::complex::hankel2e(order, z),
        _ => panic!("unexpected complex function {api}"),
    }
}

fn airy_value(family: &str, api: &str, x: f64) -> b::Result<Complex64> {
    let (ai, aip, bi, bip) = match family {
        "airy" => b::airy(x)?,
        "airye" => b::airye(x)?,
        _ => panic!("unexpected Airy family {family}"),
    };
    let value = match api {
        "ai" => ai,
        "aip" => aip,
        "bi" => bi,
        "bip" => bip,
        _ => panic!("unexpected Airy component {api}"),
    };
    Ok(Complex64::new(value, 0.0))
}

fn check_references(
    selected: &str,
    api_filter: Option<&str>,
    tag_prefix: Option<&str>,
    expected_count: usize,
) {
    let mut checked = 0;
    let mut failures: Vec<(f64, String)> = Vec::new();
    for line in REFERENCES.lines().skip(1) {
        let fields: Vec<_> = line.split(',').collect();
        assert_eq!(fields.len(), 10, "invalid fixture row: {line}");
        let [
            family,
            tag,
            api,
            order,
            re,
            im,
            expected_re,
            expected_im,
            rtol,
            atol,
        ] = fields.as_slice()
        else {
            unreachable!()
        };
        if *family != selected
            || api_filter.is_some_and(|filter| filter != *api)
            || tag_prefix.is_some_and(|prefix| !tag.starts_with(prefix))
        {
            continue;
        }

        let parse = |s: &&str| s.parse::<f64>().unwrap();
        let order = parse(order);
        let z = Complex64::new(parse(re), parse(im));
        let expected = Complex64::new(parse(expected_re), parse(expected_im));
        let rel_tol = parse(rtol);
        let abs_tol = parse(atol);
        assert!(expected.re.is_finite() && expected.im.is_finite());
        assert!(expected.norm().is_finite() && expected.norm() > 0.0);
        assert!(rel_tol > 0.0 && abs_tol >= 0.0);

        checked += 1;
        let actual = match selected {
            "complex" => complex_value(api, order, z),
            "airy" | "airye" => airy_value(family, api, z.re),
            "real_k" => {
                assert!(
                    order.fract() != 0.0,
                    "real K case should exercise noninteger order"
                );
                assert!(z.re > 0.0 && z.im == 0.0);
                assert!(expected.re > 0.0 && expected.im == 0.0);
                let value = match *api {
                    "kv" => b::kv(order, z.re),
                    "kve" => b::kve(order, z.re),
                    _ => panic!("unexpected real K api: {api}"),
                };
                Ok(Complex64::new(value, 0.0))
            }
            _ => panic!("unexpected numerical family {selected}"),
        };
        let actual = match actual {
            Ok(value) => value,
            Err(error) => {
                failures.push((
                    f64::INFINITY,
                    format!("{tag} {api}({order}, {z:?}): returned error {error:?}"),
                ));
                continue;
            }
        };
        let absolute_error = (actual - expected).norm();
        let bound = rel_tol * expected.norm() + abs_tol;
        if !actual.re.is_finite()
            || !actual.im.is_finite()
            || !absolute_error.is_finite()
            || absolute_error > bound
        {
            let ratio = absolute_error / bound;
            failures.push((
                if ratio.is_finite() {
                    ratio
                } else {
                    f64::INFINITY
                },
                format!(
                    "{tag} {api}({order}, {z:?}): got={actual:?}, expected={expected:?}, \
                     abs_error={absolute_error:.3e}, bound={bound:.3e}, ratio={ratio:.3e}"
                ),
            ));
        }
    }

    assert_eq!(checked, expected_count, "fixture row count changed");
    if !failures.is_empty() {
        failures.sort_by(|a, b| b.0.total_cmp(&a.0));
        let excerpt = failures
            .iter()
            .take(30)
            .map(|(_, message)| message.as_str())
            .collect::<Vec<_>>()
            .join("\n");
        panic!(
            "{}/{} {selected} reference cases failed (30 worst shown):\n{excerpt}",
            failures.len(),
            checked,
        );
    }
}

#[test]
fn complex_bessel_hankel_native_accuracy() {
    // Every public complex API across small arguments, all quadrants, signed
    // cuts, near-integer/half-integer orders, turning points, and |z|~750.
    check_references("complex", None, None, 600);
}

#[test]
fn airy_unscaled_native_accuracy() {
    // Negative oscillatory x, both sides of the real kernel boundaries at
    // x=10, and positive x through the exponentially unbalanced regime.
    check_references("airy", None, None, 88);
}

#[test]
fn airy_scaled_native_accuracy() {
    // Checks all four components through x=1000, where an unscaled Ai
    // underflows and unscaled Bi overflows in binary64.
    check_references("airye", None, None, 56);
}

#[test]
fn real_noninteger_kv_native_accuracy() {
    // Noninteger |v| from 1e-8 to 250.25 at 2^-1074 <= x <= ~740.
    // Finite nonzero subnormal K values and order-related transitions
    // remain checked instead of being silently replaced with zero.
    check_references("real_k", Some("kv"), None, 717);
}

#[test]
fn real_noninteger_kve_native_accuracy() {
    // Scaled K stays representable through x=1e6; include transition
    // points x~v, x~v^2 and the inverse-power asymptotic switch.
    check_references("real_k", Some("kve"), None, 838);
}

#[test]
fn real_k_ultratiny_argument_accuracy() {
    // mpmath, 85 digits: 6 positive noninteger orders × x=1e-308,1e-300
    // × both kv/kve; omit the two overflowing v=1.01,x=1e-308 results.
    // K_0.01(1e-308) ~= 60187.83686672893 is finite even though the
    // mode of the direct K integral lies beyond t=700.
    check_references("real_k", None, Some("ultratiny_"), 22);
}

#[test]
fn real_k_subnormal_argument_accuracy() {
    // 4 small noninteger orders x 3 subnormal arguments x both kv/kve.
    // The smallest input is the exact positive f64 with bits=1; mpmath
    // references were evaluated at that *binary* value, not decimal 5e-324.
    // K_0.5(f64::from_bits(1)) ~= 5.6385522612647099e161, finite.
    assert_eq!(5e-324_f64.to_bits(), 1);
    check_references("real_k", None, Some("subnormal_"), 24);
}
