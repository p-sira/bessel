use crate::backend::StrErr;
use num_traits::Float;

pub(crate) fn frexp(x: f64) -> (f64, i32) {
    if x == 0.0 || !x.is_finite() {
        return (x, 0);
    }
    let (mantissa, exponent, sign) = Float::integer_decode(x);
    let shift = 64 - mantissa.leading_zeros() as i32;
    (
        sign as f64 * mantissa as f64 / (2.0_f64).powi(shift),
        exponent as i32 + shift,
    )
}

pub(crate) fn finish<T: Float>(value: f64, error: StrErr) -> Result<T, StrErr> {
    if value.is_nan() {
        Err(error)
    } else {
        T::from(value).ok_or("Unsupported floating-point representation.")
    }
}

/// Define checked real APIs over double-precision kernels, matching xsf's float promotion.
macro_rules! real_kernel {
    ($(#[$meta:meta])* $name:ident($($arg:ident),+) $body:block) => {
        $(#[$meta])*
        pub fn $name<T: num_traits::Float>($($arg: T),+) -> Result<T, crate::backend::StrErr> {
            $(if $arg.is_nan() { return Err(concat!(stringify!($name), ": Arguments cannot be NAN.")); }
            let $arg = $arg.to_f64().ok_or(concat!(stringify!($name), ": Cannot represent argument as f64."))?;)+
            let value: f64 = $body;
            crate::backend::utils::finish(value, concat!(stringify!($name), ": Result is undefined in the real domain."))
        }
    };
}
pub(crate) use real_kernel;
