#[allow(unused_imports)]
use num_traits::Float;
// Rust translation adapted from spec_math 0.1.6 (MIT); see LICENSES.
/*
* Cephes Math Library Release 2.8: June, 2000
* Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

use crate::backend::cephes64::jv::jv;
use crate::backend::cephes64::yn::yn;
use num_complex::Complex64;

/// Leading small-x Neumann terms for |v| < 1, evaluated with logarithms so
/// even the least subnormal x does not round x/2 to zero.
fn small_fractional(v: f64, x: f64) -> f64 {
    let log_half = x.ln() - core::f64::consts::LN_2;
    let positive = (v * log_half - super::gamma(1.0 + v).ln()).exp();
    let negative = (-v * log_half - super::gamma(1.0 - v).ln()).exp();
    (super::cospi(v) * positive - negative) / super::sinpi(v)
}

/// Y_v for near-integer orders and x <= 1e-6. A fixed-radius Cauchy
/// contour loses precision as x -> 0: its order-plane growth eventually
/// overwhelms the finite number of quadrature points. Evaluate the two
/// leading small-x solutions instead and use stable upward Y recurrence.
fn small_near_integer(v: f64, x: f64) -> f64 {
    let a = v.abs();
    let n = a.round() as usize;
    let u = a - n as f64;
    let log_half = x.ln() - core::f64::consts::LN_2;
    let pi = core::f64::consts::PI;
    let euler = super::EULER;
    // log(Gamma(1+u)) = -gamma*u + zeta(2)*u²/2 - zeta(3)*u³/3 + ...
    // Write the leading Y_u reflection identity as a hyperbolic sine;
    // division by sin(pi*u) then remains stable even for subnormal u.
    let t_over_u = log_half + euler + 1.202_056_903_159_594_2 * u * u / 3.0;
    let t = u * t_over_u;
    let sinhc = if t.abs() < 1e-4 {
        1.0 + t * t / 6.0 + t.powi(4) / 120.0
    } else {
        t.sinh() / t
    };
    let pi_u = pi * u;
    let sinc = 1.0 - pi_u * pi_u / 6.0;
    let common = (-(pi_u * pi_u) / 12.0).exp();
    let mut previous =
        common * (2.0 * t_over_u * sinhc - 0.5 * pi * pi * u * t.exp()) / (pi * sinc);
    if n == 0 {
        return if v >= 0.0 {
            previous
        } else {
            // Y_-a = cos(pi*a)*Y_a + sin(pi*a)*J_a.
            super::cospi(a) * previous + super::sinpi(a) * jv(a, x)
        };
    }
    // For Y_(u+1), the (x/2)^(-u-1) term dominates by at least 1e12.
    let mut current = -(super::gamma(1.0 + u).ln() - pi.ln() - (1.0 + u) * log_half).exp();
    for k in 1..n {
        let next = (2.0 * (k as f64 + u) / x) * current - previous;
        previous = current;
        current = next;
        if !current.is_finite() {
            break;
        }
    }
    if v >= 0.0 {
        current
    } else {
        super::cospi(a) * current + super::sinpi(a) * jv(a, x)
    }
}

pub fn yv(v: f64, x: f64) -> f64 {
    let n = v as isize;
    if n as f64 == v {
        yn(n, x)
    } else if v == v.floor() {
        /* Zero in denominator. */
        //sf_error("yv", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        if x > 0.0 && x <= 1e-6 && v.is_finite() {
            if (v - v.round()).abs() < 1e-5 {
                return small_near_integer(v, x);
            }
            if v.abs() < 1.0 {
                return small_fractional(v, x);
            }
        }
        // The reflection formula below loses all significant digits when
        // sin(pi*v) is tiny: J_v and J_-v nearly cancel. The native complex
        // kernel evaluates the removable integer-order limit by a stable
        // contour formula. Supply an infinitesimal imaginary component to
        // bypass its real Cephes fast path (which calls this function).
        //
        // Cephes J_v also rejects negative noninteger orders of magnitude
        // >= 500, making the reflection formula NaN for otherwise ordinary
        // Y_v values. The native kernel uses stable order recurrence there.
        if x > 0.0
            && x.is_finite()
            && v.is_finite()
            && ((v - v.round()).abs() < 1e-5 || v.abs() >= 500.0)
        {
            let imaginary = (x * f64::EPSILON).max(f64::from_bits(1));
            if let Ok(value) = crate::backend::complex::cylindrical(
                crate::backend::complex::Kind::Y,
                v,
                Complex64::new(x, imaginary),
                false,
            ) {
                return value.re;
            }
        }
        let c = crate::backend::cephes64::cospi(v);
        let s = crate::backend::cephes64::sinpi(v);
        let first = if c == 0.0 { 0.0 } else { c * jv(v, x) };
        let y = (first - jv(-v, x)) / s;

        if y.is_infinite() {
            if v > 0.0 {
                //sf_error("yv", SF_ERROR_OVERFLOW, NULL);
                -f64::INFINITY
            } else if v < -1e10 {
                /* Whether it's +inf or -inf is numerically ill-defined. */
                //sf_error("yv", SF_ERROR_DOMAIN, NULL);
                f64::NAN
            } else {
                y
            }
        } else {
            y
        }
    }
}
