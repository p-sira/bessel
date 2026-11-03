//! Independent positive-real modified Bessel K numerical kernel.
//!
//! For finite `x > 0`, K_v(x) = ∫ exp(-x cosh(t)) cosh(vt) dt is
//! positive and even in the real order v. The implementation combines
//! elementary and upward-recurrence cases, a tiny-argument expansion,
//! a descending large-x asymptotic series, a uniform large-parameter saddle
//! expansion, and logarithmically scaled quadrature. The `scaled` variant
//! evaluates exp(x) K_v(x).
use core::f64::consts::FRAC_PI_2;

#[inline]
fn log_cosh(y: f64) -> f64 {
    let a = y.abs();
    a + (-2.0 * a).exp().ln_1p() - core::f64::consts::LN_2
}

/// Positive real K_v(x), with the optional exp(x) scaling.
pub(crate) fn kv(v: f64, x: f64, scaled: bool) -> f64 {
    let v = v.abs();
    if !x.is_finite() || !v.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    if v == 0.0 {
        return if scaled {
            super::cephes64::k0e(x)
        } else {
            super::cephes64::k0(x)
        };
    }
    if v == 1.0 {
        return if scaled {
            super::cephes64::k1e(x)
        } else {
            super::cephes64::k1(x)
        };
    }

    // K_(1/2)(x) = sqrt(pi/(2x)) exp(-x). Divide by sqrt(x) rather than
    // taking sqrt(pi/(2x)), which would overflow prematurely at subnormal x.
    // Check this before the tiny-x approximation to retain the exact form.
    if v == 0.5 {
        let result = FRAC_PI_2.sqrt() / x.sqrt();
        return result * if scaled { 1.0 } else { (-x).exp() };
    }

    // Upward recursion for integer and higher half-integer orders is
    // stable at every positive x.
    if (v == v.floor() || v.fract() == 0.5) && v <= 4096.0 {
        let mut previous;
        let mut current;
        let first;
        let mut order;
        if v.fract() == 0.5 {
            first = 0.5;
            previous = FRAC_PI_2.sqrt() / x.sqrt();
            current = previous * (1.0 + 1.0 / x);
            order = 1.5;
        } else {
            first = 0.0;
            previous = super::cephes64::k0e(x);
            current = super::cephes64::k1e(x);
            order = 1.0;
        }
        if v == first {
            return previous * if scaled { 1.0 } else { (-x).exp() };
        }
        if v == first + 1.0 {
            return current * if scaled { 1.0 } else { (-x).exp() };
        }
        while order < v {
            let next = previous + (2.0 * order / x) * current;
            previous = current;
            current = next;
            order += 1.0;
            if current.is_infinite() {
                return current;
            }
        }
        return current * if scaled { 1.0 } else { (-x).exp() };
    }

    // The quadrature peak moves to t ≈ log(2v/x) as x approaches zero.
    // At subnormal x this is too close to the overflow threshold of cosh(t),
    // even when K_v(x) is easily representable. The leading two terms of
    // DLMF 10.30.2 converge to machine precision in this range:
    // K_v(x) ≈ Γ(v) (x/2)^(-v)/2 + Γ(-v) (x/2)^v/2.
    if x < 1e-250 {
        const EULER: f64 = 0.577_215_664_901_532_9;
        let log_two_over_x = core::f64::consts::LN_2 - x.ln();
        if v <= 1e-5 {
            // The two gamma terms nearly cancel as v -> 0. Expand their
            // logarithms instead, retaining the continuous K_0 limit.
            let a = log_two_over_x - EULER;
            let t = v * a;
            let sinhc = if t.abs() < 1e-4 {
                1.0 + t * t / 6.0 + t.powi(4) / 120.0
            } else {
                t.sinh() / t
            };
            return a * sinhc * (core::f64::consts::PI.powi(2) * v * v / 12.0).exp();
        }
        let gamma = super::cephes64::gamma(v);
        let lead = (gamma.ln() - core::f64::consts::LN_2 + v * log_two_over_x).exp();
        if v < 0.25 {
            return lead + 0.5 * super::cephes64::gamma(-v) * (-v * log_two_over_x).exp();
        }
        return lead;
    }

    // For x much larger than v², the descending asymptotic series has a
    // small first omitted term, including for the exponentially scaled K.
    if x > 100.0 * (v + 1.0).powi(2) {
        let mut sum = 1.0;
        let mut term = 1.0;
        for k in 1..=24 {
            let odd = (2 * k - 1) as f64;
            term *= (4.0 * v * v - odd * odd) / (8.0 * x * k as f64);
            sum += term;
            if term.abs() <= f64::EPSILON * sum.abs() {
                break;
            }
        }
        let result = (FRAC_PI_2 / x).sqrt() * sum;
        return if scaled { result } else { result * (-x).exp() };
    }

    // At large x or v the integral is localized to a saddle of width
    // O(1/sqrt(hypot(x,v))). The fixed [0,end] trapezoid grid below can
    // miss such a narrow peak even at its maximum 8192 subdivisions.
    //
    // Use the uniform Debye expansion with q=hypot(x,v), t=v/q. For
    // exp(x)*K_v(x) its leading log phase is
    //   v*asinh(v/x) - (q-x)
    // where q-x = v^2/(q+x) avoids subtracting nearly equal large values.
    // The correction factors are the first three Debye polynomials for K,
    // with alternating signs. At q>=1e5 the next term is negligible at
    // binary64 precision throughout 0<=t<=1.
    let q = x.hypot(v);
    if q >= 1e5 {
        let t2 = (v / q).powi(2);
        let inv_q = q.recip();
        let u1 = (5.0 * t2 - 3.0) * (inv_q / 24.0);
        let u2 = (81.0 + t2 * (-462.0 + 385.0 * t2)) * (inv_q * inv_q / 1152.0);
        let u3 = (-30375.0 + t2 * (369603.0 + t2 * (-765765.0 + 425425.0 * t2)))
            * (inv_q * inv_q * inv_q / 414720.0);
        let correction = 1.0 + u1 + u2 + u3;
        let log_phase = if scaled {
            v * (v / x).asinh() - v * (v / (q + x))
        } else {
            // Avoid exp(phase)*exp(-x), which can overflow its first
            // factor even when unscaled K_v(x) is representable.
            v * (v / x).asinh() - q
        };
        return (log_phase + 0.5 * (FRAC_PI_2.ln() - q.ln()) + correction.ln()).exp();
    }

    // Find the maximum of -x(cosh(t)-1)+log(cosh(vt)). Its derivative
    // is v*tanh(vt)-x*sinh(t), with an interior maximum when v²>x.
    // Subtract its value before exponentiation to keep the quadrature scaled.
    let mode = if v > x.sqrt() {
        let mut low = 0.0;
        let mut high = (v / x).asinh();
        for _ in 0..60 {
            let mid = low + (high - low) * 0.5;
            if v * (v * mid).tanh() > x * mid.sinh() {
                low = mid;
            } else {
                high = mid;
            }
        }
        low + (high - low) * 0.5
    } else {
        0.0
    };
    if !mode.is_finite() || mode >= 700.0 {
        return f64::INFINITY;
    }
    let log_integrand = |t: f64| -x * (t.cosh() - 1.0) + log_cosh(v * t);
    let peak = log_integrand(mode).max(0.0);
    let mut end = mode + 1.0;
    while end < 700.0 && log_integrand(end) > peak - 45.0 {
        end += 0.5;
    }

    // Periodic trapezoid refinements converge rapidly for this smooth even
    // integrand. The endpoint is placed deep in the decaying tail.
    let integrand = |t: f64| (log_integrand(t) - peak).exp();
    let mut steps = 32_usize;
    let mut previous = 0.0_f64;
    let mut area;
    loop {
        let h = end / steps as f64;
        let mut sum = 0.5 * (integrand(0.0) + integrand(end));
        for i in 1..steps {
            sum += integrand(i as f64 * h);
        }
        area = sum * h;
        if steps >= 128 && (area - previous).abs() <= 8.0 * f64::EPSILON * area {
            break;
        }
        if steps >= 8192 {
            break;
        }
        previous = area;
        steps *= 2;
    }
    (area.ln() + peak - if scaled { 0.0 } else { x }).exp()
}

#[cfg(test)]
mod tests {
    use super::kv;

    #[test]
    fn half_integer_closed_forms() {
        for x in [0.01, 0.5, 1.0, 10.0, 100.0] {
            let base = (core::f64::consts::FRAC_PI_2 / x).sqrt();
            assert!((kv(1.5, x, true) / (base * (1.0 + 1.0 / x)) - 1.0).abs() < 2e-15);
            assert!(
                (kv(2.5, x, true) / (base * (1.0 + 3.0 / x + 3.0 / x.powi(2))) - 1.0).abs() < 2e-15
            );
        }
    }

    #[test]
    fn public_half_integer_k_is_finite_at_subnormal_arguments() {
        for x in [f64::from_bits(1), 1e-320, 1e-308] {
            let expected = core::f64::consts::FRAC_PI_2.sqrt() / x.sqrt();
            for actual in [crate::kv(0.5, x), crate::kve(0.5, x)] {
                assert!(actual.is_finite());
                assert!(
                    (actual / expected - 1.0).abs() < 4.0 * f64::EPSILON,
                    "K_0.5({x}) = {actual}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn large_order_narrow_saddle_matches_mpmath() {
        // Independently evaluated with mpmath at 75 decimal digits. SciPy
        // AMOS cannot evaluate the two largest-argument cases (returns NaN).
        for (v, x, expected) in [
            (10_000.0, 1e10, 1.259_596_400_556_116e-5),
            (1e7, 1e15, 4.166_531_433_451_295e-8),
            (10_000.0, 1e8, 2.066_365_668_451_389_6e-4),
            (100_000.0, 1e8, 6.498_036_028_959_877e17),
        ] {
            let actual = kv(v, x, true);
            assert!(
                actual.is_finite() && (actual / expected - 1.0).abs() < 2e-10,
                "kve({v},{x})={actual:e}, mpmath={expected:e}"
            );
        }
    }
}
