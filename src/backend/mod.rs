#![allow(clippy::excessive_precision, unused_imports)]
pub type StrErr = &'static str;
pub(crate) mod bessel;
pub(crate) mod cephes64;
pub(crate) mod complex;
pub(crate) mod double_double;
pub(crate) mod integrals;
pub(crate) mod kelvin;
pub(crate) mod ratios;
pub(crate) mod spherical;
pub(crate) mod struve;
pub(crate) mod utils;
pub use bessel::{cyl_hankel_1, cyl_hankel_1e, cyl_hankel_2, cyl_hankel_2e};
pub use bessel::{
    iv as cyl_bessel_i, ive as cyl_bessel_ie, kv as cyl_bessel_k, kve as cyl_bessel_ke,
};
pub use kelvin::*;
use num_traits::Float;
pub(crate) use ratios::iv_ratio_pair;
pub use spherical::*;
pub use struve::*;
pub fn gammasgn<T: Float>(x: T) -> Result<T, StrErr> {
    let x = x.to_f64().ok_or("gamma: invalid argument")?;
    if x.is_nan() || (x <= 0.0 && x.fract() == 0.0) {
        return Err("gamma: pole");
    }
    let sign = if x > 0.0 {
        1.0
    } else {
        cephes64::sinpi(x).signum()
    };
    utils::finish(sign, "gamma: invalid result")
}
utils::real_kernel! {
    /// Cylindrical first-kind Bessel J_v(x).
    cyl_bessel_j(v,x) {if v.is_finite() && x.is_finite() { cephes64::jv(v,x) } else {f64::NAN}}
}
utils::real_kernel! {
    /// Cylindrical second-kind Bessel Y_v(x), x >= 0.
    cyl_bessel_y(v,x) {bessel::real_y(v,x)}
}
