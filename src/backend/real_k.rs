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
    // log(1 + exp(-2a)) is below half an ulp of a for a >= 20.
    // Quadrature at large order visits this region frequently.
    if a >= 20.0 {
        return a - core::f64::consts::LN_2;
    }
    a + (-2.0 * a).exp().ln_1p() - core::f64::consts::LN_2
}

// Real K for a reduced order |u| <= 1/2. Temme's power series and Steed's
// CF2 are adapted from the Sci-special Cephes I_v implementation, in turn
// based on SciPy's XSF cephes/scipy_iv.h. Both return exp(x) K_u and
// exp(x) K_(u+1); stable upward recurrence then reaches the requested order.
// If either iteration fails, the quadrature below remains the fallback.
fn temme_k(u: f64, x: f64) -> Option<(f64, f64)> {
    const EPS: f64 = f64::EPSILON;
    let gp = super::cephes64::gamma(1.0 + u) - 1.0;
    let gm = super::cephes64::gamma(1.0 - u) - 1.0;
    let a = (0.5 * x).ln();
    let b = (u * a).exp();
    let sigma = -a * u;
    let c = if u.abs() < EPS {
        1.0
    } else {
        (core::f64::consts::PI * u).sin() / (u * core::f64::consts::PI)
    };
    let d = if sigma.abs() < EPS {
        1.0
    } else {
        sigma.sinh() / sigma
    };
    let gamma1 = if u.abs() < EPS {
        -0.577_215_664_901_532_9
    } else {
        0.5 * (gp - gm) * c / u
    };
    let gamma2 = 0.5 * (2.0 + gp + gm) * c;
    let mut p = (1.0 + gp) / (2.0 * b);
    let mut q = (1.0 + gm) * b * 0.5;
    let mut f = (sigma.cosh() * gamma1 - d * a * gamma2) / c;
    let mut coefficient = 1.0;
    let mut sum = f;
    let mut sum1 = p;
    for k in 1..1000 {
        let k = k as f64;
        f = (k * f + p + q) / (k * k - u * u);
        p /= k - u;
        q /= k + u;
        let h = p - k * f;
        coefficient *= x * x / (4.0 * k);
        sum += coefficient * f;
        sum1 += coefficient * h;
        if (coefficient * f).abs() < sum.abs() * EPS {
            let ex = x.exp();
            return Some((sum * ex, 2.0 * sum1 / x * ex));
        }
    }
    None
}

fn cf2_k(u: f64, x: f64) -> Option<(f64, f64)> {
    let mut a = u * u - 0.25;
    let mut b = 2.0 * (x + 1.0);
    let mut d = b.recip();
    let mut f = d;
    let mut delta = d;
    let mut previous = 0.0;
    let mut current = 1.0;
    let mut c = -a;
    let mut q_sum = c;
    let mut sum = 1.0 + q_sum * delta;
    let mut k = 1.0;
    for _ in 2..1000 {
        a -= 2.0 * k;
        k += 1.0;
        b += 2.0;
        d = (b + a * d).recip();
        delta *= b * d - 1.0;
        f += delta;
        let q = (previous - (b - 2.0) * current) / a;
        previous = current;
        current = q;
        c *= -a / k;
        q_sum += c * q;
        sum += q_sum * delta;
        if (q_sum * delta).abs() < sum.abs() * f64::EPSILON {
            let value = FRAC_PI_2.sqrt() / x.sqrt() / sum;
            return Some((value, value * (0.5 + u + x + (u * u - 0.25) * f) / x));
        }
    }
    None
}

#[inline]
fn reduced_order_k(v: f64, x: f64, scaled: bool) -> Option<f64> {
    // Moderate noninteger orders are cheap via a reduced-order K and stable
    // recurrence. Restrict the fast path so near-integer cancellation, very
    // high-order overflow, and small-x special handling use their existing
    // vetted implementations.
    // The reduced-order algorithms also work for positive orders below one:
    // n=0 directly returns K_u; n=1 directly returns K_(u+1). Previously
    // those inputs fell through to expensive adaptive quadrature.
    if !(1e-4..=32.0).contains(&v) || (v - v.round()).abs() < 1e-4 {
        return None;
    }
    let n = v.round() as usize;
    let u = v - n as f64;
    let (mut previous, mut current) = if x <= 2.0 {
        temme_k(u, x)?
    } else {
        cf2_k(u, x)?
    };
    for k in 1..n {
        let next = previous + 2.0 * (u + k as f64) / x * current;
        previous = current;
        current = next;
    }
    let result = if n == 0 { previous } else { current };
    if !result.is_finite() || result <= 0.0 {
        return None;
    }
    Some(if scaled { result } else { result * (-x).exp() })
}

/// Positive saddle of the K integrand when v² > x. The derivative of its
/// logarithm is g(t) = v tanh(vt) - x sinh(t). For t > 0, g(t)/t decreases
/// strictly from v²-x to a negative value at asinh(v/x), giving one root.
fn saddle_mode(v: f64, x: f64) -> f64 {
    if v <= x.sqrt() {
        return 0.0;
    }
    let mut high = (v / x).asinh();
    // tanh(vt) rounds to one at this argument. The endpoint then agrees
    // with the positive saddle to floating-point precision.
    if v * high >= 32.0 {
        return high;
    }
    let delta = v.mul_add(v, -x);
    if delta <= 0.0 {
        return 0.0;
    }

    // Avoid subtracting nearly equal values when the interior saddle
    // emerges from t=0. This is g(t)/t and its derivative, using Taylor
    // expansions valid for |vt|, |t| < 0.02; the omitted terms are below
    // binary64 rounding in that region.
    let residual = |t: f64| {
        let vt = v * t;
        if vt.abs() < 0.02 && t.abs() < 0.02 {
            let y2 = vt * vt;
            let t2 = t * t;
            let tanh_correction = y2
                * (1.0 / 3.0
                    + y2 * (-2.0 / 15.0
                        + y2 * (17.0 / 315.0 + y2 * (-62.0 / 2835.0 + y2 * 1382.0 / 155925.0))));
            let sinh_correction =
                t2 * (1.0 / 6.0 + t2 * (1.0 / 120.0 + t2 * (1.0 / 5040.0 + t2 / 362880.0)));
            let tanh_derivative = v
                * vt
                * (2.0 / 3.0
                    + y2 * (-8.0 / 15.0
                        + y2 * (102.0 / 315.0 + y2 * (-496.0 / 2835.0 + y2 * 13820.0 / 155925.0))));
            let sinh_derivative =
                t * (1.0 / 3.0 + t2 * (1.0 / 30.0 + t2 * (1.0 / 840.0 + t2 / 45360.0)));
            let squared_order = v * v;
            (
                delta - squared_order * tanh_correction - x * sinh_correction,
                -squared_order * tanh_derivative - x * sinh_derivative,
            )
        } else {
            let tanh = vt.tanh();
            (
                v * tanh - x * t.sinh(),
                v * v * (1.0 - tanh * tanh) - x * t.cosh(),
            )
        }
    };

    // The known interval [0, asinh(v/x)] contains the positive saddle.
    // Newton steps are accepted only if they stay inside that bracket;
    // otherwise bisect it. Limit Newton to 16 iterations and retain a
    // 60-step bisection fail-safe for difficult or nearly flat roots.
    let mut low = 0.0;
    let mut current = high;
    if delta < 1e-3 * x {
        let a = (v * v) * (v * v / 3.0) + x / 6.0;
        let guess = (delta / a).sqrt();
        if guess > low && guess < high {
            current = guess;
        }
    }
    for _ in 0..16 {
        let (value, slope) = residual(current);
        if value == 0.0 {
            return current;
        }
        if !value.is_finite() || !slope.is_finite() {
            break;
        }
        if value > 0.0 {
            low = current;
        } else {
            high = current;
        }
        if high - low <= 4.0 * f64::EPSILON * high {
            return low + 0.5 * (high - low);
        }
        let correction = if slope < 0.0 { value / slope } else { f64::NAN };
        if correction.is_finite() && correction.abs() <= 4.0 * f64::EPSILON * current {
            return current;
        }
        let proposal = current - correction;
        let next = if proposal > low && proposal < high {
            proposal
        } else {
            low + 0.5 * (high - low)
        };
        if next == current {
            return current;
        }
        current = next;
    }
    for _ in 0..60 {
        let middle = low + 0.5 * (high - low);
        if middle == low || middle == high {
            return middle;
        }
        let (value, _) = residual(middle);
        if value == 0.0 {
            return middle;
        }
        if value > 0.0 {
            low = middle;
        } else {
            high = middle;
        }
        if high - low <= 4.0 * f64::EPSILON * high {
            return low + 0.5 * (high - low);
        }
    }
    low + 0.5 * (high - low)
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

    // Keep the vetted large-x and uniform Debye paths above in charge of
    // their original domains; use reduced-order Temme/CF2 only where the
    // remaining alternative is the much slower general quadrature.
    if let Some(result) = reduced_order_k(v, x, scaled) {
        return result;
    }

    // Find the maximum of -x(cosh(t)-1)+log(cosh(vt)). Its derivative
    // is v*tanh(vt)-x*sinh(t), with an interior maximum when v²>x.
    // Subtract its value before exponentiation to keep the quadrature scaled.
    let mode = saddle_mode(v, x);
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
    // integrand. The endpoint is placed deep in the decaying tail. Every
    // refined grid contains all nodes of the previous grid, so retain their
    // accumulated integral and evaluate only the new midpoint nodes.
    let integrand = |t: f64| (log_integrand(t) - peak).exp();
    let mut steps = 32_usize;
    let initial_h = end / steps as f64;
    let mut sum = 0.5 * (integrand(0.0) + integrand(end));
    for i in 1..steps {
        sum += integrand(i as f64 * initial_h);
    }
    let mut area = sum * initial_h;
    loop {
        let refined_steps = steps * 2;
        let h = end / refined_steps as f64;
        let mut midpoint_sum = 0.0;
        for i in 0..steps {
            midpoint_sum += integrand((2 * i + 1) as f64 * h);
        }
        let refined_area = 0.5 * area + midpoint_sum * h;
        steps = refined_steps;
        if steps >= 128 && (refined_area - area).abs() <= 8.0 * f64::EPSILON * refined_area {
            area = refined_area;
            break;
        }
        area = refined_area;
        if steps >= 8192 {
            break;
        }
    }
    (area.ln() + peak - if scaled { 0.0 } else { x }).exp()
}

#[cfg(test)]
mod tests {
    use super::{kv, log_cosh, saddle_mode};

    #[test]
    fn saddle_modes_match_high_precision_near_both_switches() {
        // Reference positive roots of v*tanh(v*t) = x*sinh(t), computed
        // independently with 105-digit Decimal exp, sinh, tanh and 280
        // bisection iterations using the exact binary64 input values.
        // Cases straddle v²=x and the v*asinh(v/x)=32 shortcut.
        for (v, x, expected) in [
            (20.0, 399.999_999_999_999_1, 4.126_952_385_024_223e-9),
            (20.0, 399.999_999_999_767_17, 6.603_123_816_041_052e-8),
            (20.0, 399.999_999_046_325_7, 4.225_999_248_301_342e-6),
            (20.0, 399.96, 0.000_865_536_502_905_355_6),
            (20.0, 400.0, 0.0),
            (20.0, 400.000_000_953_674_3, 0.0),
            (2.0, 3.999_999_99, 4.082_482_897_619_654e-5),
            (0.5, 0.249_999_9, 0.001_264_911_232_740_387_4),
            (35.25, 1_242.562_499_999_999, 1.329_092_202_813_063_4e-9),
            (35.25, 1_242.562_499_046_325_7, 1.360_990_416_307_023_4e-6),
            (35.25, 1242.52, 0.000_287_315_064_254_099_4),
            (35.25, 1_242.562_5, 0.0),
            (35.25, 400.0, 0.087_648_286_356_781_26),
            (35.25, 33.968_432_725_882_57, 0.907_801_418_439_716_3),
            (35.25, 33.968_432_725_9, 0.907_801_418_439_346_8),
            (35.25, 33.968_432_725_8, 0.907_801_418_441_466_7),
            (100.25, 705.0, 0.141_723_670_734_642_76),
            (150.25, 10_000.0, 0.014_662_161_810_281_675),
        ] {
            let mode = saddle_mode(v, x);
            assert!(
                mode.is_finite() && (mode - expected).abs() <= 5e-12 * expected.abs().max(1e-16),
                "saddle_mode({v},{x}) = {mode:.17e}, decimal reference {expected:.17e}"
            );
        }
    }

    #[test]
    fn large_log_cosh_matches_direct_expression() {
        for y in [20.0_f64, 32.0, 100.0, 1e6, 1e200] {
            let a = y.abs();
            let direct = a + (-2.0 * a).exp().ln_1p() - core::f64::consts::LN_2;
            assert_eq!(log_cosh(y).to_bits(), direct.to_bits());
            assert_eq!(log_cosh(-y).to_bits(), direct.to_bits());
        }
    }

    #[test]
    fn general_quadrature_matches_high_precision_with_and_without_saddle_shortcut() {
        // Independently computed with mpmath at 75 decimal digits. The
        // large-order cases hit the saddle shortcut, while (35.25,400)
        // retains the bisection path.
        for (v, x, expected) in [
            (40.25, 30.0, 5_838_319_055.265_112),
            (70.25, 1.0, 9.398_380_971_653_315e119),
            (80.25, 80.0, 2_497_200_778_473_931.5),
            (35.25, 400.0, 0.295_235_606_351_876_1),
        ] {
            let actual = kv(v, x, true);
            assert!(
                (actual / expected - 1.0).abs() < 2e-10,
                "kve({v},{x})={actual:e}, mpmath={expected:e}"
            );
        }
    }

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

    #[test]
    fn reduced_order_fast_path_matches_scipy() {
        // Independent SciPy 1.18.1 exp(x)*K_7.25(x) references across every
        // CodSpeed cylindrical benchmark argument. Also check K_-v = K_v
        // and the unscaled result where exp(-x) is representable.
        for (x, expected) in [
            (0.05, 250_208_156_526_753.7),
            (0.5, 21_849_353.312_655_04),
            (1.0, 229_694.590_696_991_5),
            (2.5, 1_093.115_229_492_91),
            (5.0, 44.296_120_801_812_414),
            (10.0, 4.447_085_082_482_836),
            (30.0, 0.537_524_873_380_052_6),
            (100.0, 0.162_571_610_063_528_16),
        ] {
            for v in [7.25, -7.25] {
                let scaled = kv(v, x, true);
                let unscaled = kv(v, x, false);
                assert!(
                    (scaled / expected - 1.0).abs() <= 2e-10,
                    "kve({v},{x})={scaled:e}, SciPy={expected:e}"
                );
                assert!(
                    (unscaled / (expected * (-x).exp()) - 1.0).abs() <= 2e-10,
                    "kv({v},{x})={unscaled:e}, SciPy scaled={expected:e}"
                );
            }
        }
    }
}
