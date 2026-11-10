//! Bounded Perron continued fraction for adjacent modified Bessel functions.
use super::double_double::Dd;

pub(crate) fn iv_ratio_pair(v: f64, x: f64) -> (f64, f64) {
    if v.is_nan() || x.is_nan() || v < 0.5 || x < 0.0 || (v.is_infinite() && x.is_infinite()) {
        return (f64::NAN, f64::NAN);
    }
    if x == 0.0 {
        return (x, 1.0);
    }
    if v.is_infinite() {
        return (0.0, 1.0);
    }
    if x.is_infinite() {
        return (1.0, 0.0);
    }
    if v == 0.5 {
        let e = (-x).exp();
        let e2 = e * e;
        return (x.tanh(), 2.0 * e2 / (1.0 + e2));
    }
    // Scale the coefficients before products; DD preserves the small complement
    // when the continued-fraction tail nearly cancels its leading 2v term.
    let c = Dd::from(1.0 / v.max(x));
    let vc = Dd::from(v) * c;
    let xc = Dd::from(x) * c;
    let two = Dd::from(2.0);
    let one = Dd::from(1.0);
    let coefficients = |n: usize| {
        (
            -(two * vc + Dd::from((2 * n - 1) as f64) * c) * xc,
            two * (vc + xc) + Dd::from(n as f64) * c,
        )
    };
    let (a, mut previous_b) = coefficients(1);
    let mut term = a.div_dd(previous_b);
    let mut ratio = Dd::from(0.0);
    let mut sum = two * vc;
    for n in 1..=1000 {
        sum = sum + term;
        if term.hi.abs() <= f64::EPSILON * sum.hi.abs() {
            let denominator = xc + sum;
            return (xc.div_dd(denominator).hi, sum.div_dd(denominator).hi);
        }
        let (a, b) = coefficients(n + 1);
        let p = a * (one + ratio);
        ratio = (-p).div_dd(p + b * previous_b);
        term = term * ratio;
        previous_b = b;
    }
    (f64::NAN, f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::{Dd, iv_ratio_pair};

    // The original full-DD quotient version, retained in tests to catch
    // changes to the cancellation-sensitive Perron iteration.
    fn original_pair(v: f64, x: f64) -> (f64, f64) {
        let c = Dd::from(1.0 / v.max(x));
        let vc = Dd::from(v) * c;
        let xc = Dd::from(x) * c;
        let two = Dd::from(2.0);
        let one = Dd::from(1.0);
        let coefficients = |n: usize| {
            (
                -(two * vc + Dd::from((2 * n - 1) as f64) * c) * xc,
                two * (vc + xc) + Dd::from(n as f64) * c,
            )
        };
        let (a, mut previous_b) = coefficients(1);
        let mut term = a / previous_b;
        let mut ratio = Dd::from(0.0);
        let mut sum = two * vc;
        for n in 1..=1000 {
            sum = sum + term;
            if term.hi.abs() <= f64::EPSILON * sum.hi.abs() {
                let denominator = xc + sum;
                return ((xc / denominator).hi, (sum / denominator).hi);
            }
            let (a, b) = coefficients(n + 1);
            let p = a * (one + ratio);
            ratio = -p / (p + b * previous_b);
            term = term * ratio;
            previous_b = b;
        }
        (f64::NAN, f64::NAN)
    }

    #[test]
    fn fast_dd_quotients_preserve_perron_outputs() {
        for v in [0.500_001, 0.51, 1.0, 1.5, 3.0, 20.0, 100.0] {
            for x in [0.01, 0.5, 1.0, 7.0, 20.0, 100.0, 1000.0] {
                let (ratio, complement) = iv_ratio_pair(v, x);
                let (old_ratio, old_complement) = original_pair(v, x);
                for (actual, expected) in [(ratio, old_ratio), (complement, old_complement)] {
                    assert!(
                        (actual.is_nan() && expected.is_nan())
                            || (actual - expected).abs() <= 3.0e-13 * expected.abs().max(1.0e-300),
                        "v={v} x={x}: actual={actual}, previous={expected}"
                    );
                }
            }
        }
    }

    #[test]
    fn half_integer_ratio_and_complement() {
        // I_(3/2)/I_(1/2) = coth(x) - 1/x.  Evaluating the complement
        // directly avoids subtracting nearly equal numbers for large x.
        for x in [0.5, 1.0, 2.0, 7.0, 20.0, 100.0] {
            let q = (-2.0_f64 * x).exp();
            let expected_complement = 1.0 / x - 2.0 * q / (1.0 - q);
            let (ratio, complement) = iv_ratio_pair(1.5, x);
            assert!(
                (ratio - (1.0 - expected_complement)).abs() <= 2.0e-13,
                "ratio at x={x}: {ratio}"
            );
            assert!(
                (complement - expected_complement).abs() <= 2.0e-13 * expected_complement,
                "complement at x={x}: {complement}"
            );
        }
    }
}
