//! Complex Airy functions, including AMOS-compatible exponential scaling.
//!
//! The power series is evaluated with double-double intermediates because Ai
//! on the positive real axis is the difference of two exponentially growing
//! series. For large arguments, sector-rotated recessive asymptotic solutions
//! account for both traveling waves across the negative-real Stokes line.

use core::f64::consts::PI;
use num_complex::Complex64;
use num_traits::Float;

use super::KernelError as Error;
use super::StrErr;
use super::double_double::Dd;

/// Tuple (Ai, Ai′, Bi, Bi′), preserving XSF's output ordering.
pub type Airy<T> = (T, T, T, T);

/// Real Ai, Ai′, Bi and Bi′, with optional exponential scaling.
pub(crate) fn real<T: Float>(x: T, scaled: bool) -> Result<Airy<T>, StrErr> {
    if !x.is_finite() {
        return Err("airy: Argument must be finite.");
    }
    let x = x.to_f64().ok_or("airy: Invalid argument.")?;
    if scaled && x < 0.0 {
        return Err("airye: Real scaled Ai requires x >= 0; use complex::airye for negative x.");
    }
    if !scaled && (-10.0..=10.0).contains(&x) {
        let (a, b, c, d) = super::cephes64::airy(x);
        return Ok((
            super::utils::finish(a, "airy: Undefined result.")?,
            super::utils::finish(b, "airy: Undefined result.")?,
            super::utils::finish(c, "airy: Undefined result.")?,
            super::utils::finish(d, "airy: Undefined result.")?,
        ));
    }
    // The Cephes real Airy kernels are accurate through this moderate positive
    // interval. Apply the real exponential factors once to all four outputs;
    // this avoids four separate complex-Bessel evaluations for airye(x).
    // Beyond it, keep the natively scaled kernels to avoid cancellation and
    // overflow in the unscaled Ai/Bi values.
    if scaled && x <= 20.0 {
        let (ai, aip, bi, bip) = super::cephes64::airy(x);
        let zeta = (2.0 / 3.0) * x * x.sqrt();
        let exp = zeta.exp();
        let inv_exp = (-zeta).exp();
        return Ok((
            super::utils::finish(ai * exp, "airye: Undefined result.")?,
            super::utils::finish(aip * exp, "airye: Undefined result.")?,
            super::utils::finish(bi * inv_exp, "airye: Undefined result.")?,
            super::utils::finish(bip * inv_exp, "airye: Undefined result.")?,
        ));
    }
    let (a, b, c, d) = values(Complex64::new(x, 0.0), scaled).map_err(Error::message)?;
    Ok((
        super::utils::finish(a.re, "airy: Undefined result.")?,
        super::utils::finish(b.re, "airy: Undefined result.")?,
        super::utils::finish(c.re, "airy: Undefined result.")?,
        super::utils::finish(d.re, "airy: Undefined result.")?,
    ))
}

const AI0: Dd = Dd {
    hi: 0.35502805388781722,
    lo: 2.0523363243621199e-17,
};
const AIP0: Dd = Dd {
    hi: -0.25881940379280682,
    lo: 2.5222431116108321e-17,
};
const BI0: Dd = Dd {
    hi: 0.61492662744600068,
    lo: 5.0899207794891416e-17,
};
const BIP0: Dd = Dd {
    hi: 0.44828835735382638,
    lo: -2.5363237774417305e-17,
};

#[derive(Clone, Copy)]
struct ComplexDd {
    re: Dd,
    im: Dd,
}

impl ComplexDd {
    #[inline]
    fn real(re: f64) -> Self {
        Self {
            re: Dd::from(re),
            im: Dd::from(0.0),
        }
    }

    #[inline]
    fn from_complex(value: Complex64) -> Self {
        Self {
            re: Dd::from(value.re),
            im: Dd::from(value.im),
        }
    }

    #[inline]
    fn scale_dd(self, factor: Dd) -> Self {
        Self {
            re: self.re * factor,
            im: self.im * factor,
        }
    }

    #[inline]
    fn divided(self, divisor: f64) -> Self {
        Self {
            re: self.re.div_small(divisor),
            im: self.im.div_small(divisor),
        }
    }

    #[inline]
    fn value(self) -> Complex64 {
        Complex64::new(self.re.hi + self.re.lo, self.im.hi + self.im.lo)
    }

    #[inline]
    fn magnitude(self) -> f64 {
        self.re.hi.hypot(self.im.hi)
    }
}

impl core::ops::Add for ComplexDd {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl core::ops::Mul for ComplexDd {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

/// Two entire solutions f(0)=1, f'(0)=0 and g(0)=0, g'(0)=1.
/// They contain only z^(3k) and z^(3k+1), respectively.
fn series(z: Complex64) -> Result<(Complex64, Complex64, Complex64, Complex64), Error> {
    let z = ComplexDd::from_complex(z);
    let z2 = z * z;
    let z3 = z2 * z;
    let mut f = ComplexDd::real(1.0);
    let mut g = z;
    let mut fp = ComplexDd::real(0.0);
    let mut gp = ComplexDd::real(1.0);
    let mut ft = f;
    let mut gt = g;
    for k in 1..=180 {
        let k = k as f64;
        fp = fp + (ft * z2).divided(3.0 * k - 1.0);
        gp = gp + (gt * z2).divided(3.0 * k);
        ft = (ft * z3).divided((3.0 * k) * (3.0 * k - 1.0));
        gt = (gt * z3).divided((3.0 * k + 1.0) * (3.0 * k));
        f = f + ft;
        g = g + gt;
        if ft.magnitude().max(gt.magnitude()) < 5e-32 * f.magnitude().max(g.magnitude()).max(1.0) {
            break;
        }
        if k == 180.0 {
            return Err(Error::ConvergenceFailure);
        }
    }
    let ai = f.scale_dd(AI0) + g.scale_dd(AIP0);
    let aip = fp.scale_dd(AI0) + gp.scale_dd(AIP0);
    let bi = f.scale_dd(BI0) + g.scale_dd(BIP0);
    let bip = fp.scale_dd(BI0) + gp.scale_dd(BIP0);
    Ok((ai.value(), aip.value(), bi.value(), bip.value()))
}

const OMEGA: Complex64 = Complex64 {
    re: -0.5,
    im: 0.8660254037844386,
};
const OMEGA2: Complex64 = Complex64 {
    re: -0.5,
    im: -0.8660254037844386,
};
const ROOT_PI_FACTOR: f64 = 0.28209479177387814; // 1/(2 sqrt(pi))
const PHASE_PLUS: Complex64 = Complex64 {
    re: 0.8660254037844386,
    im: 0.5,
};
const PHASE_MINUS: Complex64 = Complex64 {
    re: 0.8660254037844386,
    im: -0.5,
};

/// Principal-branch single recessive Airy wave, with the scaling already
/// absorbed into its exponent. Only invoked in the sector |arg(z)| <= 2pi/3.
fn wave(
    z: Complex64,
    zeta: Complex64,
    adjustment: Complex64,
) -> Result<(Complex64, Complex64), Error> {
    let inv = Complex64::new(1.0, 0.0) / zeta;
    let mut term = Complex64::new(1.0, 0.0);
    let mut sum = term;
    let mut derivative_sum = term;
    let mut old_size = f64::INFINITY;
    for k in 1..=180 {
        let kf = k as f64;
        let next = term * inv * (-((6.0 * kf - 5.0) * (6.0 * kf - 1.0) / (72.0 * kf)));
        let size = next.norm();
        if size >= old_size {
            break;
        }
        term = next;
        sum += term;
        derivative_sum -= term * ((6.0 * kf + 1.0) / (6.0 * kf - 1.0));
        if size < 3e-17 * sum.norm().max(derivative_sum.norm()) {
            break;
        }
        old_size = size;
    }
    let exponent = -zeta + adjustment;
    if exponent.re > 709.0 {
        return Err(Error::Overflow);
    }
    let exponential = exponent.exp() * ROOT_PI_FACTOR;
    let quartic = z.sqrt().sqrt();
    let ai = exponential * sum / quartic;
    let aip = -exponential * derivative_sum * quartic;
    if !ai.re.is_finite() || !ai.im.is_finite() || !aip.re.is_finite() || !aip.im.is_finite() {
        return Err(Error::Overflow);
    }
    Ok((ai, aip))
}

/// Ai and Ai' in all sectors. On the negative axis the two rotated solutions
/// represent both oscillatory waves, including their Stokes multipliers.
fn rotated_zeta(z: Complex64, zeta: Complex64, rotated: Complex64) -> Complex64 {
    // Rotating z by 2pi/3 changes z^(3/2) by a sign. Determine that sign
    // from the principal arguments instead of recomputing a large zeta and
    // losing digits when the scaled exponential factors cancel.
    if (rotated.arg() - z.arg()).abs() > PI {
        zeta
    } else {
        -zeta
    }
}

fn asympt_ai(
    z: Complex64,
    zeta: Complex64,
    adjustment: Complex64,
) -> Result<(Complex64, Complex64), Error> {
    if z.arg().abs() <= 2.0 * PI / 3.0 {
        wave(z, zeta, adjustment)
    } else {
        // Ai(z) + omega Ai(omega*z) + omega² Ai(omega²*z) = 0.
        let wa = OMEGA * z;
        let wb = OMEGA2 * z;
        let (a, ap) = wave(wa, rotated_zeta(z, zeta, wa), adjustment)?;
        let (b, bp) = wave(wb, rotated_zeta(z, zeta, wb), adjustment)?;
        Ok((-OMEGA * a - OMEGA2 * b, -OMEGA2 * ap - OMEGA * bp))
    }
}

/// Bi(z) = exp(-i*pi/6) Ai(omega²*z) + exp(i*pi/6) Ai(omega*z).
fn asympt_bi(
    z: Complex64,
    zeta: Complex64,
    adjustment: Complex64,
) -> Result<(Complex64, Complex64), Error> {
    let wa = OMEGA2 * z;
    let wb = OMEGA * z;
    let (a, ap) = asympt_ai(wa, rotated_zeta(z, zeta, wa), adjustment)?;
    let (b, bp) = asympt_ai(wb, rotated_zeta(z, zeta, wb), adjustment)?;
    Ok((
        PHASE_MINUS * a + PHASE_PLUS * b,
        PHASE_MINUS * OMEGA2 * ap + PHASE_PLUS * OMEGA * bp,
    ))
}

/// Ai, Ai', Bi, Bi'. For `scaled`, Ai and Ai' carry `exp(zeta)` while Bi
/// and Bi' carry `exp(-abs(Re(zeta)))`, zeta = 2*z*sqrt(z)/3.
pub(crate) fn values(
    z: Complex64,
    scaled: bool,
) -> Result<(Complex64, Complex64, Complex64, Complex64), Error> {
    if !z.re.is_finite() || !z.im.is_finite() {
        return Err(Error::InvalidInput);
    }
    let radius = z.norm();
    if radius > 1e10 {
        return Err(Error::TotalPrecisionLoss);
    }
    let zeta = (2.0 / 3.0) * z * z.sqrt();
    if radius < 8.0 {
        let (mut ai, mut aip, mut bi, mut bip) = series(z)?;
        if scaled {
            let ai_factor = zeta.exp();
            let bi_factor = (-zeta.re.abs()).exp();
            ai *= ai_factor;
            aip *= ai_factor;
            bi *= bi_factor;
            bip *= bi_factor;
        }
        return Ok((ai, aip, bi, bip));
    }
    let (ai, aip) = asympt_ai(
        z,
        zeta,
        if scaled {
            zeta
        } else {
            Complex64::new(0.0, 0.0)
        },
    )?;
    let (bi, bip) = asympt_bi(
        z,
        zeta,
        if scaled {
            Complex64::new(-zeta.re.abs(), 0.0)
        } else {
            Complex64::new(0.0, 0.0)
        },
    )?;
    Ok((ai, aip, bi, bip))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn near_origin_and_negative_real_wronskian() {
        let (a, ap, b, bp) = values(Complex64::new(0.0, 0.0), false).unwrap();
        assert!((a.re - AI0.hi).abs() < 1e-16);
        assert!((ap.re - AIP0.hi).abs() < 1e-16);
        assert!((b.re - BI0.hi).abs() < 1e-16);
        assert!((bp.re - BIP0.hi).abs() < 1e-16);
        for z in [
            Complex64::new(8.0, 0.0),
            Complex64::new(-8.0, 0.0),
            Complex64::new(-20.0, 0.0),
            Complex64::new(11.0, 13.0),
            Complex64::new(-15.0, 0.5),
        ] {
            let (a, ap, b, bp) = values(z, false).unwrap();
            let det = a * bp - ap * b;
            assert!(
                (det.re - 1.0 / PI).abs() < 2e-9 && det.im.abs() < 2e-9,
                "Wronskian at {z}: {det}"
            );
        }
    }

    #[test]
    fn scaled_large_real_reference() {
        for (x, ai_ref, aip_ref, bi_ref, bip_ref) in [
            (
                12.0,
                0.15119256068463708,
                -0.5268505009124518,
                0.3039054138807329,
                1.0463290385080119,
            ),
            (
                30.0,
                0.12045939663973669,
                -0.6607833325212036,
                0.24122445526882688,
                1.3192228350679096,
            ),
            (
                100.0,
                0.08919692093633041,
                -0.8921920625040315,
                0.1784310111708354,
                1.7838637549628087,
            ),
        ] {
            let (ai, aip, bi, bip) = values(Complex64::new(x, 0.0), true).unwrap();
            for (value, reference) in [
                (ai.re, ai_ref),
                (aip.re, aip_ref),
                (bi.re, bi_ref),
                (bip.re, bip_ref),
            ] {
                assert!(
                    (value / reference - 1.0).abs() < 3e-11,
                    "airye({x}): {value} vs {reference}"
                );
            }
        }
    }

    #[test]
    fn complex_references_across_asymptotic_sectors() {
        // Independent AMOS values (the former complex-bessel implementation).
        // Include both sides of the negative-real Stokes line and the series/
        // asymptotic transition, with all four scaled quantities checked.
        let cases = [
            (
                Complex64::new(3.0, -4.0),
                [
                    Complex64::new(0.18357553003565633, 0.0416176487393803),
                    Complex64::new(-0.4129909058014225, 0.09208370774312627),
                    Complex64::new(0.2731914926203999, -0.2771397791581104),
                    Complex64::new(0.20768534826166007, -0.7907563262094387),
                ],
            ),
            (
                Complex64::new(8.1, 1.0),
                [
                    Complex64::new(0.16611374675853574, -0.004974539243459435),
                    Complex64::new(-0.4794463188231269, -0.014209897925746479),
                    Complex64::new(-0.3177260152843692, 0.10717216239989037),
                    Complex64::new(-0.9152792044240922, 0.24525456963640316),
                ],
            ),
            (
                Complex64::new(-10.0, 1e-8),
                [
                    Complex64::new(-0.02472019295539314, -0.03175325546701553),
                    Complex64::new(-0.6120058182017353, -0.7861252151369666),
                    Complex64::new(-0.31467981969278364, 1.1941409870705769e-9),
                    Complex64::new(0.11941410962372998, 3.146798310815391e-8),
                ],
            ),
            (
                Complex64::new(-10.0, -1e-8),
                [
                    Complex64::new(-0.02472019295539314, 0.03175325546701553),
                    Complex64::new(-0.6120058182017353, 0.7861252151369666),
                    Complex64::new(-0.31467981969278364, -1.1941409870705769e-9),
                    Complex64::new(0.11941410962372998, -3.146798310815391e-8),
                ],
            ),
            (
                Complex64::new(-12.0, 5.0),
                [
                    Complex64::new(0.11485336302626185, -0.0945182710326984),
                    Complex64::new(-0.4126378427959292, -0.34009917718435384),
                    Complex64::new(-0.01572707360216159, 0.14791097901563985),
                    Complex64::new(0.5327105987123407, -0.04644198152849127),
                ],
            ),
        ];
        for (z, expected) in cases {
            let (a, ap, b, bp) = values(z, true).unwrap();
            for (actual, reference) in [a, ap, b, bp].into_iter().zip(expected) {
                let relative_error = (actual - reference).norm() / reference.norm();
                assert!(
                    relative_error < 2e-11,
                    "airye({z}): actual={actual}, ref={reference}, relative={relative_error:e}"
                );
            }
        }
    }

    #[test]
    fn scaled_values_remain_finite_past_unscaled_overflow() {
        for x in [100.0, 1_000.0, 10_000.0] {
            let (ai, aip, bi, bip) = values(Complex64::new(x, 0.0), true).unwrap();
            for result in [ai, aip, bi, bip] {
                assert!(
                    result.re.is_finite() && result.im.is_finite(),
                    "x={x}, value={result}"
                );
            }
            let wronskian = ai * bip - aip * bi;
            assert!(
                (wronskian.re - 1.0 / PI).abs() < 2e-13,
                "x={x}, W={wronskian}"
            );
        }
    }
}

#[cfg(test)]
mod real_tests {
    use super::real as airy_real;

    #[test]
    fn scaled_real_airy_matches_high_precision_at_kernel_transitions() {
        // mpmath 1.4, 70-digit Ai, Ai', Bi, Bi' with exp(±2x^(3/2)/3).
        let references: [(f64, [f64; 4]); 14] = [
            (
                0.0,
                [
                    0.3550280538878172,
                    -0.2588194037928068,
                    0.6149266274460007,
                    0.44828835735382636,
                ],
            ),
            (
                0.5,
                [
                    0.29327715912994736,
                    -0.28469116209194257,
                    0.6748924111156302,
                    0.4302209614637694,
                ],
            ),
            (
                1e-12,
                [
                    0.3550280538875584,
                    -0.2588194037928068,
                    0.614926627446449,
                    0.4482883573538264,
                ],
            ),
            (
                2.09,
                [
                    0.22801909727853855,
                    -0.35295656709963385,
                    0.49299852416324267,
                    0.6328540952588387,
                ],
            ),
            (
                2.1,
                [
                    0.22778653656129074,
                    -0.3532931777439597,
                    0.492202880303034,
                    0.6340057173130742,
                ],
            ),
            (
                4.0,
                [
                    0.1970948026430665,
                    -0.40568791157876305,
                    0.4048094678892981,
                    0.7817739306272438,
                ],
            ),
            (
                8.3203353,
                [
                    0.16539804113168443,
                    -0.48193818011893056,
                    0.33368593640111104,
                    0.9522113574366132,
                ],
            ),
            (
                8.33,
                [
                    0.16535122268088864,
                    -0.482073267100383,
                    0.33358641257956495,
                    0.9524985171426074,
                ],
            ),
            (
                12.0,
                [
                    0.15119256068463708,
                    -0.5268505009124518,
                    0.3039054138807329,
                    1.0463290385080119,
                ],
            ),
            (
                19.999,
                [
                    0.1332420614603509,
                    -0.5975158764105255,
                    0.2671056874514886,
                    1.1911403614781674,
                ],
            ),
            (
                20.0,
                [
                    0.13324040735181362,
                    -0.5975232737163193,
                    0.2671023248341684,
                    1.1911553993720082,
                ],
            ),
            (
                20.0001,
                [
                    0.13324024194659163,
                    -0.5975240134318844,
                    0.26710198858411147,
                    1.1911569031298575,
                ],
            ),
            (
                30.0,
                [
                    0.12045939663973669,
                    -0.6607833325212036,
                    0.24122445526882688,
                    1.3192228350679096,
                ],
            ),
            (
                100.0,
                [
                    0.08919692093633041,
                    -0.8921920625040315,
                    0.1784310111708354,
                    1.7838637549628087,
                ],
            ),
        ];
        for (x, expected) in references {
            let (ai, aip, bi, bip) = airy_real(x, true).unwrap();
            for (actual, reference) in [ai, aip, bi, bip].into_iter().zip(expected) {
                assert!(
                    (actual - reference).abs() <= 1e-12 * reference.abs(),
                    "airye({x}) = {actual}, expected {reference}"
                );
            }
            // W[Ai,Bi] = 1/pi; opposing exponential scalings cancel.
            assert!(
                (ai * bip - aip * bi - core::f64::consts::FRAC_1_PI).abs() < 1e-12,
                "Wronskian at x={x}"
            );
        }
    }
}
