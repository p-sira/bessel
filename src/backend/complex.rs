//! Native binary64 complex cylindrical Bessel and Hankel functions.
//!
//! Small and moderate arguments use the hypergeometric power series for J
//! and I, followed by the exact reflection identities for Y and K.  The
//! removable integer-order singularities in those identities are evaluated
//! by Cauchy's integral formula in the order plane.  For large arguments the
//! Hankel and modified-Bessel inverse-power expansions avoid cancellation
//! between exponentially growing solutions.  Continuation across the
//! negative-real axis follows the DLMF principal branches.
use crate::backend::{KernelError as Error, StrErr};
use num_complex::{Complex, Complex64};
use num_traits::Float;
use std::f64::consts::{FRAC_PI_2, PI};
pub(crate) fn promote<T: Float>(z: Complex<T>) -> Result<Complex64, StrErr> {
    if z.re.is_nan() || z.im.is_nan() {
        return Err("complex: NaN argument");
    }
    Ok(Complex64::new(
        z.re.to_f64().ok_or("complex: invalid argument")?,
        z.im.to_f64().ok_or("complex: invalid argument")?,
    ))
}
pub(crate) fn finish<T: Float>(z: Complex64) -> Result<Complex<T>, StrErr> {
    if z.re.is_nan() || z.im.is_nan() {
        return Err("complex: undefined result");
    }
    Ok(Complex::new(
        T::from(z.re).ok_or("complex: invalid result")?,
        T::from(z.im).ok_or("complex: invalid result")?,
    ))
}

type C = Complex64;

const ZERO: C = C::new(0.0, 0.0);
const ONE: C = C::new(1.0, 0.0);
const I: C = C::new(0.0, 1.0);
const MAX_SERIES: usize = 1500;

fn finite(z: C) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn checked(value: C) -> Result<C, Error> {
    if finite(value) {
        Ok(value)
    } else {
        Err(Error::Overflow)
    }
}

/// Apply an exponential scale in polar form, so an intermediate exp(±Re z)
/// need not overflow or underflow while the final Bessel value is finite.
fn restore_scale(value: C, log_scale: f64, phase: f64) -> Result<C, Error> {
    if value == ZERO {
        return Ok(ZERO);
    }
    let logarithm = value.norm().ln() + log_scale;
    if logarithm > f64::MAX.ln() {
        return Err(Error::Overflow);
    }
    checked(C::from_polar(logarithm.exp(), value.arg() + phase))
}

fn validate(v: f64, z: C) -> Result<(), Error> {
    if !v.is_finite() || !finite(z) {
        Err(Error::InvalidInput)
    } else if v.abs() >= 1.0e6 || z.norm() > 1.0e12 {
        Err(Error::TotalPrecisionLoss)
    } else {
        Ok(())
    }
}

/// Lanczos logarithm of Gamma for a complex argument.  The coefficients are
/// those of the standard g=7, nine-term double-precision approximation.
fn log_gamma(z: C) -> C {
    const P: [f64; 8] = [
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if z.re < 0.5 {
        return C::new(PI.ln(), 0.0) - (z * PI).sin().ln() - log_gamma(ONE - z);
    }
    let z = z - ONE;
    let mut x = C::new(0.999_999_999_999_809_9, 0.0);
    for (k, p) in P.iter().enumerate() {
        x += p / (z + (k + 1) as f64);
    }
    let t = z + 7.5;
    C::new(0.918_938_533_204_672_7, 0.0) + (z + 0.5) * t.ln() - t + x.ln()
}

fn inv_gamma(z: C) -> C {
    if z.im == 0.0 {
        if z.re <= 0.0 && z.re == z.re.floor() {
            return ZERO;
        }
        if z.re > 0.0 && z.re < 171.0 {
            return C::new(1.0 / crate::backend::cephes64::gamma(z.re), 0.0);
        }
    }
    (-log_gamma(z)).exp()
}

/// Sum (±z²/4)^k / (k! Gamma(k+nu+1)), with (z/2)^nu included.
/// A complex order is used for stable integer-order limits.
fn series(nu: C, z: C, modified: bool) -> Result<C, Error> {
    if z == ZERO {
        return if nu == ZERO {
            Ok(ONE)
        } else if nu.im == 0.0 && nu.re > 0.0 {
            Ok(ZERO)
        } else {
            Err(Error::Overflow)
        };
    }
    if nu.im == 0.0 && nu.re < 0.0 && nu.re == nu.re.floor() {
        let sign = if !modified && nu.re.abs() % 2.0 != 0.0 {
            -1.0
        } else {
            1.0
        };
        return Ok(sign * series(-nu, z, modified)?);
    }
    let u = z * 0.5;
    let t = if modified { u * u } else { -(u * u) };
    let mut term = (nu * u.ln()).exp() * inv_gamma(nu + ONE);
    let mut sum = term;
    let mut converged = false;
    // The threshold on k prevents a premature return when the early terms
    // are tiny but grow before reaching their eventual maximum.
    let minimum = (z.norm() * 0.6).ceil() as usize + 5;
    for k in 1..=MAX_SERIES {
        let divisor = C::new(k as f64, 0.0) * (nu + k as f64);
        if divisor == ZERO {
            return Err(Error::ConvergenceFailure);
        }
        term *= t / divisor;
        sum += term;
        if !finite(term) || !finite(sum) {
            return Err(Error::Overflow);
        }
        if k > minimum && term.norm() <= 5.0e-16 * sum.norm() {
            converged = true;
            break;
        }
        if term == ZERO && k > minimum {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(Error::ConvergenceFailure);
    }
    checked(sum)
}

/// Noninteger-order Y or K, allowing a complex order for order-plane
/// contour integration.  Never call this at a pole of sin(pi*nu).
fn reflection(nu: C, z: C, modified: bool) -> Result<C, Error> {
    let sin = (nu * PI).sin();
    if sin.norm() < 1.0e-9 {
        return Err(Error::TotalPrecisionLoss);
    }
    let positive = series(nu, z, modified)?;
    let negative = series(-nu, z, modified)?;
    let result = if modified {
        (negative - positive) * FRAC_PI_2 / sin
    } else {
        ((nu * PI).cos() * positive - negative) / sin
    };
    checked(result)
}

/// The quotient in the reflection identities has a removable singularity
/// at each integer order.  The trapezoidal Cauchy integral converges
/// geometrically for this analytic function of order, without subtracting
/// nearly identical J or I values or taking finite differences.
fn near_integer(v: f64, z: C, modified: bool) -> Result<C, Error> {
    const N: usize = 48;
    const R: f64 = 0.375;
    let n = v.round();
    let delta = v - n;
    let mut sum = ZERO;
    for k in 0..N {
        let theta = 2.0 * PI * (k as f64 + 0.5) / N as f64;
        let circle = C::from_polar(R, theta);
        let order = C::new(n, 0.0) + circle;
        let value = reflection(order, z, modified)?;
        sum += value * circle / (circle - delta);
    }
    checked(sum / N as f64)
}

/// Integer-order Y and K from their convergent harmonic-number expansions.
///
/// Differentiating the J/I power series with respect to the order at zero
/// gives Y0,Y1/K0,K1 without a removable-pole limit. The higher integer
/// orders follow stable dominant-solution recurrences. In particular, this
/// replaces 48 Cauchy contour nodes * two full complex-order series for
/// each common integer-order Hankel call.
fn integer_second(n: usize, z: C, modified: bool) -> Result<C, Error> {
    let u = z * 0.5;
    let t = if modified { u * u } else { -u * u };
    let logarithm = u.ln() + crate::backend::cephes64::EULER;
    let mut term = ONE;
    let mut sum = ONE;
    let mut next_sum = ONE;
    let mut harmonic = 0.0;
    let mut harmonic_sum = ZERO;
    let mut harmonic_derivative = ZERO;
    for k in 1..=160 {
        let kf = k as f64;
        term *= t / (kf * kf);
        harmonic += 1.0 / kf;
        sum += term;
        next_sum += term / (kf + 1.0);
        harmonic_sum += harmonic * term;
        harmonic_derivative += kf * harmonic * term;
        if k > 8 && term.norm() * (kf + 1.0) < 1e-17 * sum.norm().max(1e-250) {
            break;
        }
    }
    let first = if modified {
        -logarithm * sum + harmonic_sum
    } else {
        (2.0 / PI) * (logarithm * sum - harmonic_sum)
    };
    if n == 0 {
        return checked(first);
    }
    let next = if modified {
        sum / z + logarithm * u * next_sum - 2.0 * harmonic_derivative / z
    } else {
        (2.0 / PI) * (logarithm * u * next_sum - sum / z + 2.0 * harmonic_derivative / z)
    };
    if n == 1 {
        return checked(next);
    }
    let (mut a, mut b) = (first, next);
    for k in 1..n {
        let c = if modified {
            a + 2.0 * k as f64 / z * b
        } else {
            2.0 * k as f64 / z * b - a
        };
        a = b;
        b = c;
    }
    checked(b)
}

fn small_second(v: f64, z: C, modified: bool) -> Result<C, Error> {
    if z == ZERO {
        return Err(Error::Overflow);
    }
    if v == v.round() && v.abs() <= 7.0 && z.norm() <= 10.5 {
        let n = v.abs() as usize;
        let result = integer_second(n, z, modified)?;
        let sign = if !modified && v < 0.0 && n % 2 == 1 {
            -1.0
        } else {
            1.0
        };
        return checked(sign * result);
    }
    if (v - v.round()).abs() <= 0.065 {
        return near_integer(v, z, modified);
    }
    reflection(C::new(v, 0.0), z, modified)
}

/// The asymptotic coefficient a_k=(prod_j (4v²-(2j-1)²))/(k!8^k).
/// `phase` is +1 for K, +i for H1, and -i for H2.
fn asympt_sum(v: f64, z: C, phase: C) -> Result<C, Error> {
    let mut sum = ONE;
    let mut term = ONE;
    let mut previous = f64::INFINITY;
    for k in 1..=200 {
        let odd = (2 * k - 1) as f64;
        term *= phase * ((4.0 * v * v - odd * odd) / (8.0 * k as f64)) / z;
        let mag = term.norm();
        // An asymptotic expansion must stop at its smallest term.
        if mag > previous {
            break;
        }
        sum += term;
        if mag < 3.0e-16 * sum.norm() || term == ZERO {
            break;
        }
        previous = mag;
    }
    checked(sum)
}

fn hankel_asym(v: f64, z: C, first: bool, scale: f64) -> Result<C, Error> {
    let sign = if first { 1.0 } else { -1.0 };
    let polynomial = asympt_sum(v, z, I * sign)?;
    let phase = I * sign * (z - (FRAC_PI_2 * v + PI / 4.0)) + scale;
    checked((C::new(2.0 / PI, 0.0) / z).sqrt() * phase.exp() * polynomial)
}

fn k_asym(v: f64, z: C, scaled: bool) -> Result<C, Error> {
    let polynomial = asympt_sum(v, z, ONE)?;
    let exponent = if scaled { ZERO } else { -z };
    checked((C::new(FRAC_PI_2, 0.0) / z).sqrt() * exponent.exp() * polynomial)
}

/// Use the recessive K expansion near the imaginary axis only when its
/// smallest term is sufficiently small. For orders near eight and |z|≈15,
/// the first term can already be increasing, so the ordinary radius cutoff
/// alone is insufficient.
fn k_asym_precise(v: f64, z: C, scaled: bool) -> Option<C> {
    let mut sum = ONE;
    let mut term = ONE;
    let mut previous = f64::INFINITY;
    for k in 1..=200 {
        let odd = (2 * k - 1) as f64;
        term *= ((4.0 * v * v - odd * odd) / (8.0 * k as f64)) / z;
        let magnitude = term.norm();
        if magnitude > previous {
            if previous <= 2.0e-12 * sum.norm() {
                return k_asym(v, z, scaled).ok();
            }
            return None;
        }
        sum += term;
        if magnitude <= 2.0e-12 * sum.norm() {
            return checked(
                (C::new(FRAC_PI_2, 0.0) / z).sqrt() * if scaled { ONE } else { (-z).exp() } * sum,
            )
            .ok();
        }
        previous = magnitude;
    }
    None
}

fn use_asympt(v: f64, z: C) -> bool {
    let radius = z.norm();
    // Fixed-order asymptotics are inappropriate close to a large-order
    // turning point.  Revert to the convergent series there.
    radius >= 22.0_f64.max(1.6 * v * v + 8.0)
        || (v.abs() <= 2.0 && radius >= 12.0)
        || (v.abs().fract() == 0.5 && radius >= 7.0 && v.abs() <= 16.5)
}

/// Stable backward Miller recurrence for the minimal solution J or I at
/// orders above the argument.  A low-order reference fixes normalization.
/// This is essential near the turning point z ≈ nu, where direct series
/// terms can exceed their final sum by many orders of magnitude.
fn minimal_recurrence(v: f64, z: C, modified: bool, scaled: bool) -> Result<C, Error> {
    let n = v.floor() as usize;
    let fraction = v - n as f64;
    let base = if modified {
        i_right(fraction, z, scaled)?
    } else {
        first_second_right(fraction, z, false, scaled)?
    };
    if n == 0 {
        return Ok(base);
    }
    let m = n.max(z.norm().ceil() as usize) + (8.0 * z.norm().sqrt()).ceil() as usize + 45;
    let mut higher = ZERO;
    let mut current = ONE;
    let mut target = ZERO;
    let mut at_one = ZERO;
    // Keep the normalization of saved order values separately. Repeated
    // rescaling can otherwise underflow the saved high-order value before
    // the recurrence reaches order zero (e.g. I_1200(1200)).
    let mut log_shift = 0.0;
    let mut target_shift = 0.0;
    let mut at_one_shift = 0.0;
    for k in (1..=m).rev() {
        if k == n {
            target = current;
            target_shift = log_shift;
        }
        if k == 1 {
            at_one = current;
            at_one_shift = log_shift;
        }
        let ratio = 2.0 * (fraction + k as f64) / z;
        let previous = if modified {
            ratio * current + higher
        } else {
            ratio * current - higher
        };
        higher = current;
        current = previous;
        if current.norm() > 1e190 {
            current *= 1e-190;
            higher *= 1e-190;
            log_shift += 190.0 * std::f64::consts::LN_10;
        }
    }
    let one = if modified {
        i_right(fraction + 1.0, z, scaled)?
    } else {
        first_second_right(fraction + 1.0, z, false, scaled)?
    };
    let (denominator, denominator_shift, reference) = if base.norm() >= one.norm() {
        (current, log_shift, base)
    } else {
        (at_one, at_one_shift, one)
    };
    if target == ZERO || reference == ZERO {
        return Ok(ZERO);
    }
    let log_size = target.norm().ln() - denominator.norm().ln() + target_shift - denominator_shift
        + reference.norm().ln();
    if log_size > f64::MAX.ln() {
        return Err(Error::Overflow);
    }
    checked(C::from_polar(
        log_size.exp(),
        target.arg() - denominator.arg() + reference.arg(),
    ))
}

/// Forward recurrence for the dominant solutions Y and K.  Low fractional
/// order evaluations have no large-order cancellation.
fn dominant_recurrence(v: f64, z: C, modified: bool, scaled: bool) -> Result<C, Error> {
    let n = v.floor() as usize;
    let f = v - n as f64;
    let mut a = if modified {
        k_right(f, z, scaled)?
    } else {
        first_second_right(f, z, true, scaled)?
    };
    if n == 0 {
        return Ok(a);
    }
    let mut b = if modified {
        k_right(f + 1.0, z, scaled)?
    } else {
        first_second_right(f + 1.0, z, true, scaled)?
    };
    for k in 1..n {
        let ratio = 2.0 * (f + k as f64) / z;
        let next = if modified {
            a + ratio * b
        } else {
            ratio * b - a
        };
        a = b;
        b = checked(next)?;
    }
    Ok(b)
}

// Nodes and weights of the 16-point Gauss-Legendre rule, positive half.
const GL16_X: [f64; 8] = [
    0.09501250983763744,
    0.2816035507792589,
    0.4580167776572274,
    0.6178762444026437,
    0.755404408355003,
    0.8656312023878317,
    0.9445750230732326,
    0.9894009349916499,
];
const GL16_W: [f64; 8] = [
    0.1894506104550685,
    0.1826034150449236,
    0.1691565193950025,
    0.1495959888165767,
    0.1246289712555339,
    0.0951585116824928,
    0.06225352393864789,
    0.02715245941175409,
];

fn k_integrand(t: f64, v: f64, z: C) -> C {
    // exp(z) K_v(z) = 1/2 integral_0^inf
    //    [exp(-z(cosh(t)-1)+vt)+exp(-z(cosh(t)-1)-vt)] dt.
    let exponent = -z * (t.cosh() - 1.0);
    ((exponent + v * t).exp() + (exponent - v * t).exp()) * 0.5
}

fn k_integral(v: f64, z: C) -> Result<C, Error> {
    let v = v.abs();
    let x = z.re;
    let t_peak = (v / x).asinh();
    let peak = v * t_peak - x * (t_peak.cosh() - 1.0);
    let mut end = t_peak + 0.5;
    while end < 30.0 && v * end - x * (end.cosh() - 1.0) > peak - 65.0 {
        end += 0.25;
    }
    if end >= 30.0 {
        return Err(Error::ConvergenceFailure);
    }
    let mut last = ZERO;
    // Start with one 16-node Gauss panel rather than 16 panels. For the
    // typical nonoscillatory right-half-plane K integral the first 2–8
    // panels already reach binary64 accuracy. Starting with 16 and 32
    // wasted hundreds of complex exponentials on every invocation.
    for iter in 0..14 {
        let intervals = 1 << iter;
        let half = end / (2 * intervals) as f64;
        let mut sum = ZERO;
        for k in 0..intervals {
            let center = (k as f64 + 0.5) * end / intervals as f64;
            for j in 0..8 {
                let delta = half * GL16_X[j];
                sum += GL16_W[j]
                    * (k_integrand(center + delta, v, z) + k_integrand(center - delta, v, z))
                    * half;
            }
        }
        checked(sum)?;
        if iter >= 1 && (sum - last).norm() < 2.0e-13 * sum.norm() {
            return Ok(sum);
        }
        last = sum;
    }
    Err(Error::ConvergenceFailure)
}

/// The direct real-axis K integral oscillates too rapidly for nearly
/// imaginary z. Shift the t contour by -sign(Im z)*i*pi/2 instead: the
/// horizontal part then decays as exp(-|Im z| sinh(t)), while the short
/// vertical part is a smooth integral over [0, pi/2].
fn k_rotated_integral(v: f64, z: C) -> Result<C, Error> {
    let sign = z.im.signum();
    let decay = z.im.abs();
    let angle = v * FRAC_PI_2;
    let mut end: f64 = 1.0;
    while decay * end.sinh() - v * end < 65.0 && end < 12.0 {
        end += 0.25;
    }
    if end >= 12.0 {
        return Err(Error::ConvergenceFailure);
    }
    let mut last = ZERO;
    for iter in 0..9 {
        let panels = 1 << iter;
        let mut sum = ZERO;
        for panel in 0..panels {
            for (span, vertical) in [(FRAC_PI_2, true), (end, false)] {
                let half = span / (2 * panels) as f64;
                let center = (panel as f64 + 0.5) * span / panels as f64;
                for j in 0..8 {
                    for t in [center - half * GL16_X[j], center + half * GL16_X[j]] {
                        let integrand = if vertical {
                            -sign * I * (-z * t.cos()).exp() * (v * t).cos()
                        } else {
                            let hyperbolic = C::new(
                                (v * t).cosh() * angle.cos(),
                                -sign * (v * t).sinh() * angle.sin(),
                            );
                            (sign * I * z * t.sinh()).exp() * hyperbolic
                        };
                        sum += half * GL16_W[j] * integrand;
                    }
                }
            }
        }
        checked(sum)?;
        if iter >= 1 && (sum - last).norm() <= 1.0e-13 * sum.norm() {
            return Ok(sum);
        }
        last = sum;
    }
    Err(Error::ConvergenceFailure)
}

/// Complex continuation of Steed's CF2 for exp(z)*K_u(z), |u| <= 1/2.
/// The real version is proven in sci-special's Cephes/XSF I_v kernel.
/// In the moderately right-half-plane sector, it avoids repeated complex
/// exponentials in the K-integral and cancellation of I_{-u}-I_u.
fn k_cf2_base(u: f64, z: C) -> Option<(C, C)> {
    let mut a = u * u - 0.25;
    let mut b = 2.0 * (z + ONE);
    let mut d = ONE / b;
    let mut f = d;
    let mut delta = d;
    let mut previous = ZERO;
    let mut current = ONE;
    let mut c = -a;
    let mut q_sum = C::new(c, 0.0);
    let mut sum = ONE + q_sum * delta;
    let mut k = 1.0;
    for _ in 2..400 {
        a -= 2.0 * k;
        k += 1.0;
        b += 2.0;
        d = ONE / (b + a * d);
        delta *= b * d - ONE;
        f += delta;
        let q = (previous - (b - 2.0) * current) / a;
        previous = current;
        current = q;
        c *= -a / k;
        q_sum += c * q;
        let correction = q_sum * delta;
        sum += correction;
        if !finite(sum) || !finite(f) {
            return None;
        }
        if correction.norm() < f64::EPSILON * sum.norm() {
            let value = C::new(FRAC_PI_2.sqrt(), 0.0) / z.sqrt() / sum;
            let next = value * (0.5 + u + z + (u * u - 0.25) * f) / z;
            return (finite(value) && finite(next)).then_some((value, next));
        }
    }
    None
}

/// Reduced-order scaled K with stable recurrence to the requested order.
fn k_cf2(v: f64, z: C, scaled: bool) -> Option<C> {
    // Conservative sector: the continued fraction converges quickly and
    // accurately here, whereas it is not reliable near imaginary axes.
    if !(z.re > 2.0 && z.im.abs() <= z.re * 0.5 && v <= 8.0 && z.norm() < 20.0) || v.fract() == 0.5
    {
        return None;
    }
    let n = v.round() as usize;
    let u = v - n as f64;
    let (mut previous, mut current) = k_cf2_base(u, z)?;
    for j in 1..n {
        let next = previous + (2.0 * (u + j as f64) / z) * current;
        previous = current;
        current = next;
    }
    let value = if n == 0 { previous } else { current };
    let result = if scaled { value } else { value * (-z).exp() };
    finite(result).then_some(result)
}

/// Unscaled/scaled J and Y in the closed right half-plane.
fn first_second_right(v: f64, z: C, second: bool, scaled: bool) -> Result<C, Error> {
    let factor = if scaled { -z.im.abs() } else { 0.0 };
    // These pure-real Cephes kernels are also the independently validated
    // fast paths used by sci-special. Do not use them on the negative-real
    // cut: first_second() applies signed-rim continuation there.
    if z.im == 0.0 && z.re > 0.0 && v.abs() <= 8.0 {
        let value = if second {
            crate::backend::cephes64::yv(v, z.re)
        } else {
            crate::backend::cephes64::jv(v, z.re)
        };
        return checked(C::new(value * factor.exp(), 0.0));
    }
    if v < -8.0 || (v < 0.0 && (v - v.round()).abs() < 0.1) {
        let positive = -v;
        let j = first_second_right(positive, z, false, scaled)?;
        let y = first_second_right(positive, z, true, scaled)?;
        return checked(if second {
            crate::backend::cephes64::cospi(positive) * y
                + crate::backend::cephes64::sinpi(positive) * j
        } else {
            crate::backend::cephes64::cospi(positive) * j
                - crate::backend::cephes64::sinpi(positive) * y
        });
    }
    // Use order recurrences through large-order turning points, even when
    // |z| exceeds the previous 800 cutoff.  The fixed-order expansion
    // cannot describe nu ~= |z|, and the raw series overflows there.
    if v >= 8.0 && z != ZERO && (z.norm() < 800.0 || (z.norm() <= 50_000.0 && !use_asympt(v, z))) {
        return if second {
            if z.im.abs() > 7.0 && z.norm() < 80.0 {
                // At moderate imaginary arguments, upward Y recurrence
                // magnifies the roundoff of its two low-order starting
                // values.  Recover Y from the accurately computed minimal
                // J and the decaying Hankel solution, represented by K on
                // a rotated argument.
                let sign = if z.im > 0.0 { -1.0 } else { 1.0 };
                let k = k_right(v, sign * I * z, false)?;
                let h = sign * (2.0 / PI) * I * C::from_polar(1.0, sign * FRAC_PI_2 * v) * k;
                let j = minimal_recurrence(v, z, false, false)?;
                let y = sign * I * (h - j);
                checked(y * factor.exp())
            } else if z.im.abs() >= 1.0 && z.norm() <= 50_000.0 {
                // Upward Y recurrence is poorly conditioned in complex
                // turning regions: cancellation can amplify errors in
                // low-order seeds by many orders of magnitude. Construct
                // the recessive Hankel solution using scaled K instead.
                let sign = if z.im > 0.0 { -1.0 } else { 1.0 };
                let w = sign * I * z;
                let k = k_right(v, w, true)?;
                let h_factor =
                    C::from_polar((-w.re - if scaled { z.im.abs() } else { 0.0 }).exp(), -w.im);
                let h =
                    sign * (2.0 / PI) * I * C::from_polar(1.0, sign * FRAC_PI_2 * v) * k * h_factor;
                let j = minimal_recurrence(v, z, false, scaled)?;
                checked(sign * I * (h - j))
            } else {
                dominant_recurrence(v, z, false, scaled)
            }
        } else {
            minimal_recurrence(v, z, false, scaled)
        };
    }
    if use_asympt(v, z) {
        let a = hankel_asym(v, z, true, factor)?;
        let b = hankel_asym(v, z, false, factor)?;
        return checked(if second {
            (a - b) / (2.0 * I)
        } else {
            (a + b) * 0.5
        });
    }
    let result = if second {
        small_second(v, z, false)?
    } else {
        series(C::new(v, 0.0), z, false)?
    };
    checked(result * factor.exp())
}

fn first_second(v: f64, z: C, second: bool, scaled: bool) -> Result<C, Error> {
    if z.re >= 0.0 {
        return first_second_right(v, z, second, scaled);
    }
    // For z=-w, arg(w) is opposite the cut rim; the rotation from w to
    // z is +pi on the upper rim and -pi on the lower rim.
    let side = if z.im.is_sign_negative() { -1.0 } else { 1.0 };
    let w = -z;
    let j = first_second_right(v, w, false, scaled)?;
    if !second {
        return checked(C::from_polar(1.0, side * PI * v) * j);
    }
    let y = first_second_right(v, w, true, scaled)?;
    checked(
        C::from_polar(1.0, -side * PI * v) * y
            + 2.0 * side * I * crate::backend::cephes64::cospi(v) * j,
    )
}

/// I on the closed right half-plane, scaled by exp(-Re z) as needed.
fn i_right(v: f64, z: C, scaled: bool) -> Result<C, Error> {
    let real_scale = if scaled { -z.re } else { 0.0 };
    if v < -8.0 || (v < 0.0 && (v - v.round()).abs() < 0.1) {
        let positive = -v;
        // I_-v = I_v + (2/pi) sin(pi*v) K_v, real v.
        return checked(
            i_right(positive, z, scaled)?
                + (2.0 / PI)
                    * crate::backend::cephes64::sinpi(positive)
                    * k_right(positive, z, false)?
                    * real_scale.exp(),
        );
    }
    if v >= 8.0 && z != ZERO && (z.norm() < 800.0 || (z.norm() <= 50_000.0 && !use_asympt(v, z))) {
        if scaled || z.re <= 700.0 {
            return minimal_recurrence(v, z, true, scaled);
        }
        // I_0(z) grows like exp(Re z), even where I_nu(z) remains
        // representable. Normalize Miller with I_0 exp(-Re z), then
        // restore the real exponential only after the order ratio.
        return restore_scale(minimal_recurrence(v, z, true, true)?, z.re, 0.0);
    }
    if use_asympt(v, z) {
        // Rotate into a right-half-plane J argument.  Its two Hankel
        // exponentials are explicitly rescaled before they are summed.
        let side = if z.im < 0.0 { 1.0 } else { -1.0 };
        let w = side * I * z;
        let prefactor = C::from_polar(1.0, -side * FRAC_PI_2 * v);
        let j = first_second_right(v, w, false, scaled)?;
        return checked(prefactor * j);
    }
    checked(series(C::new(v, 0.0), z, true)? * real_scale.exp())
}

fn i_general(v: f64, z: C, scaled: bool) -> Result<C, Error> {
    if z.re >= 0.0 {
        return i_right(v, z, scaled);
    }
    let side = if z.im.is_sign_negative() { -1.0 } else { 1.0 };
    let w = -z;
    checked(C::from_polar(1.0, side * PI * v) * i_right(v, w, scaled)?)
}

/// K on the right half-plane; its analytic continuation into the left half
/// plane is handled by k_general.
fn k_right(v: f64, z: C, scaled: bool) -> Result<C, Error> {
    let v = v.abs();
    if z.im == 0.0 && z.re > 0.0 {
        if !scaled && z.re > 700.0 {
            // At order near the argument, K_v is representable even when
            // exp(-Re z) underflows before multiplication by Kve.
            let scaled_value = crate::backend::real_k::kv(v, z.re, true);
            return restore_scale(C::new(scaled_value, 0.0), -z.re, 0.0);
        }
        return checked(C::new(crate::backend::real_k::kv(v, z.re, scaled), 0.0));
    }
    if v == 0.5 {
        // Avoid overflow in pi/(2*z) when the final half-order K is finite.
        let prefactor = FRAC_PI_2.sqrt() / z.sqrt();
        return checked(prefactor * if scaled { ONE } else { (-z).exp() });
    }
    if let Some(value) = k_cf2(v, z, scaled) {
        return Ok(value);
    }
    if v >= 8.0 && z.norm() >= 11.5 {
        if !scaled && z.re > 700.0 {
            // K_0(z) underflows before the upward recurrence can reach a
            // large order for which K_nu(z) is still representable.
            return restore_scale(dominant_recurrence(v, z, true, true)?, -z.re, -z.im);
        }
        return dominant_recurrence(v, z, true, scaled);
    }
    // The inverse-power expansion is extremely accurate once |z| >= 20
    // for orders below eight, including arguments close to the imaginary
    // axis where the I_-v - I_v reflection loses many digits. For example
    // K_7.25(0.1+25i) is accurate to binary64 here, not ~5e-7.
    if use_asympt(v, z) || (v < 8.0 && z.norm() >= 20.0) {
        return k_asym(v, z, scaled);
    }
    if v < 9.0 && z.im.abs() >= 8.0 && z.re <= 0.5 * z.im.abs() {
        if let Some(value) = k_asym_precise(v, z, scaled) {
            return Ok(value);
        }
        let value = k_rotated_integral(v, z)?;
        return checked(if scaled { value * z.exp() } else { value });
    }
    if z.re > 0.1 && (z.norm() > 6.0 || v > 8.0) {
        let value = k_integral(v, z)?;
        return checked(if scaled { value } else { value * (-z).exp() });
    }
    if z.im == 0.0 && z.re > 0.0 && (v == 0.0 || v.abs() == 1.0) {
        let value = if v == 0.0 {
            if scaled {
                crate::backend::cephes64::k0e(z.re)
            } else {
                crate::backend::cephes64::k0(z.re)
            }
        } else if scaled {
            crate::backend::cephes64::k1e(z.re)
        } else {
            crate::backend::cephes64::k1(z.re)
        };
        return checked(C::new(value, 0.0));
    }
    checked(small_second(v, z, true)? * if scaled { z.exp() } else { ONE })
}

fn k_general(v: f64, z: C, scaled: bool) -> Result<C, Error> {
    if z.re >= 0.0 {
        return k_right(v, z, scaled);
    }
    let side = if z.im.is_sign_negative() { -1.0 } else { 1.0 };
    let w = -z;
    let rotation = C::from_polar(1.0, -side * PI * v);
    if scaled {
        // exp(z)K(z) = exp(-2w-i*pi*v*side) K_scaled(w)
        //              - side*i*pi*exp(-i*Im(w)) I_scaled(w).
        let k = k_right(v, w, true)?;
        let i = i_right(v, w, true)?;
        return checked(
            (-2.0 * w).exp() * rotation * k - side * PI * I * C::from_polar(1.0, -w.im) * i,
        );
    }
    checked(rotation * k_right(v, w, false)? - side * PI * I * i_right(v, w, false)?)
}

fn hankel(v: f64, z: C, first: bool, scaled: bool) -> Result<C, Error> {
    // The decaying Hankel solution is evaluated through K in its decaying
    // sector.  This also works when exp(±Im z) is far beyond f64 range.
    if ((first && z.im > 0.0) || (!first && z.im < 0.0))
        // At moderate arguments, summing J +/- iY is both faster and
        // stabler than the rotated K quadrature. In particular the latter
        // fails to converge for H1(7.25, 10+2i).
        && !(z.norm() <= 14.0 && v.abs() <= 8.0 && z.re.abs() >= z.im.abs())
        && !(v.abs() >= 8.0 && z.re.abs() >= z.im.abs() && z.im.abs() < 12.0)
    {
        let sign = if first { -1.0 } else { 1.0 };
        let w = sign * I * z;
        let k = k_right(v, w, scaled)?;
        let prefactor = sign * (2.0 / PI) * I * C::from_polar(1.0, sign * FRAC_PI_2 * v);
        return checked(prefactor * k);
    }
    let j = first_second(v, z, false, scaled)?;
    let y = first_second(v, z, true, scaled)?;
    let h = j + if first { I * y } else { -I * y };
    if scaled {
        let sign = if first { -1.0 } else { 1.0 };
        let phase = C::from_polar(1.0, sign * z.re);
        // j and y already contain exp(-|Im z|).
        checked(h * phase * (z.im.abs() - sign * z.im).exp())
    } else {
        checked(h)
    }
}

/// The cylindrical solution to evaluate on the principal complex branch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Kind {
    J,
    Y,
    I,
    K,
    Hankel1,
    Hankel2,
}

/// Evaluate a complex cylindrical Bessel or Hankel function.
///
/// Scaling matches the SciPy/AMOS conventions: exp(-|Im(z)|) for J/Y,
/// exp(-|Re(z)|) for I, exp(z) for K, and exp(∓i*z) for Hankel 1/2.
/// Negative-real branch cut rims are distinguished by the sign of Im(z),
/// including signed zero.
pub(crate) fn cylindrical(
    kind: Kind,
    v: f64,
    z: Complex64,
    scaled: bool,
) -> Result<Complex64, Error> {
    validate(v, z)?;
    if z == ZERO && v == v.round() && matches!(kind, Kind::J | Kind::I) {
        return Ok(if v == 0.0 { ONE } else { ZERO });
    }
    match kind {
        Kind::J => first_second(v, z, false, scaled),
        Kind::Y => first_second(v, z, true, scaled),
        Kind::I => i_general(v, z, scaled),
        Kind::K => k_general(v, z, scaled),
        Kind::Hankel1 => hankel(v, z, true, scaled),
        Kind::Hankel2 => hankel(v, z, false, scaled),
    }
}
