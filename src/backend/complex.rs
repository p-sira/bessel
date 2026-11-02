use crate::backend::StrErr;
use num_complex::{Complex, Complex64};
use num_traits::Float;
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
