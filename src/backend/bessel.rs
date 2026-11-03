//! Cylindrical real Bessel functions and checked complex Hankel adapters.
//!
//! Real J, Y, and I use the Cephes/Temme kernels; positive-real K has a
//! dedicated stable kernel. Complex families share the native complex backend.
use crate::backend::{
    StrErr,
    complex::{self, Kind, finish, promote},
    utils::real_kernel,
};
use num_complex::{Complex, Complex64};
#[allow(unused_imports)]
use num_traits::Float;

fn checked_hankel<T: Float>(
    v: T,
    z: Complex<T>,
    kind: Kind,
    scaled: bool,
) -> Result<Complex<T>, StrErr> {
    let z = promote(z)?;
    let v = v.to_f64().ok_or("Bessel: Invalid order.")?;
    if !v.is_finite() {
        return Err("Bessel: Order must be finite.");
    }
    let value =
        complex::cylindrical(kind, v, z, scaled).map_err(crate::backend::KernelError::message)?;
    finish(value)
}

macro_rules! checked_hankel {
    ($($name:ident, $kind:ident, $scaled:expr, $doc:literal;)+) => {$(
        #[doc = $doc]
        pub fn $name<T: Float>(v: T, z: Complex<T>) -> Result<Complex<T>, StrErr> {
            checked_hankel(v, z, Kind::$kind, $scaled)
        }
    )+};
}
fn modified_i_real(v: f64, x: f64, scaled: bool) -> f64 {
    if !v.is_finite() || x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 && v != v.floor() {
        return f64::NAN;
    }
    // `frem` is speculatable in LLVM, so keep the libm `fmod` call behind an
    // explicit branch rather than a short-circuiting condition.
    let sign = if x < 0.0 {
        if v % 2.0 != 0.0 { -1.0 } else { 1.0 }
    } else {
        1.0
    };
    let x = x.abs();
    if x == 0.0 {
        if v < 0.0 && v != v.floor() {
            return f64::INFINITY;
        }
        return if v == 0.0 { 1.0 } else { 0.0 };
    }
    if x.is_infinite() {
        return if scaled { 0.0 } else { sign * f64::INFINITY };
    }
    // DLMF 10.39.1: elementary half-integer I values are more accurate near zero.
    if v.abs() == 0.5 {
        let prefactor = (2.0 * core::f64::consts::PI).sqrt().recip() / x.sqrt();
        let difference = if v > 0.0 {
            -(-2.0 * x).exp_m1()
        } else {
            1.0 + (-2.0 * x).exp()
        };
        let value = prefactor * difference;
        return if scaled {
            value
        } else {
            let half = (0.5 * x).exp();
            (value * half) * half
        };
    }
    sign * crate::backend::cephes64::iv::iv(v, x, scaled)
}

fn modified_k_real(v: f64, x: f64, scaled: bool) -> f64 {
    if !v.is_finite() || x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::INFINITY;
    }
    if x.is_infinite() {
        return 0.0;
    }
    crate::backend::real_k::kv(v, x, scaled)
}
real_kernel! {
    /// Real modified Bessel I_v(x). Negative x requires integer v.
    iv(v,x) { modified_i_real(v,x,false) }
}
real_kernel! {
    /// Scaled modified Bessel exp(-|x|)*I_v(x). Negative x requires integer v.
    ive(v,x) { modified_i_real(v,x,true) }
}
real_kernel! {
    /// Real modified Bessel K_v(x), for x >= 0. K_-v = K_v.
    kv(v,x) { modified_k_real(v,x,false) }
}
real_kernel! {
    /// Scaled modified Bessel exp(x)*K_v(x), for x >= 0.
    kve(v,x) { modified_k_real(v,x,true) }
}

pub(crate) fn real_y(v: f64, x: f64) -> f64 {
    if !v.is_finite() || !x.is_finite() || x < 0.0 {
        return f64::NAN;
    }
    if v.abs() < 10000.0 {
        return crate::backend::cephes64::yv(v, x);
    }
    // Avoid converting huge orders into machine integers or performing a
    // recurrence whose running time is proportional to the order.
    match complex::cylindrical(Kind::Y, v, Complex64::new(x, 0.0), false) {
        Ok(z) => z.re,
        Err(_) => f64::NAN,
    }
}

checked_hankel! {
    cyl_hankel_1, Hankel1, false, "Hankel function H_v^(1)(z).";
    cyl_hankel_1e, Hankel1, true, "Scaled Hankel function exp(-i*z)*H_v^(1)(z).";
    cyl_hankel_2, Hankel2, false, "Hankel function H_v^(2)(z).";
    cyl_hankel_2e, Hankel2, true, "Scaled Hankel function exp(i*z)*H_v^(2)(z).";
}
