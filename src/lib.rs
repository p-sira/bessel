//! Cylindrical Bessel, Hankel, modified Bessel, spherical Bessel, Struve, and Kelvin functions.
//!
//! # Function families
//!
//! - **Cylindrical Bessel:** [`jv`], [`yv`], [`iv`], and [`kv`] evaluate J, Y, I,
//!   and K for real order and argument. Integer-order and order-zero/one
//!   shortcuts, scaled variants (`ive`, `kve`, etc.), and derivatives
//!   ([`jvp`], [`yvp`], [`ivp`], [`kvp`]) are also available.
//! - **Complex Bessel and Hankel:** [`complex`] evaluates J, Y, I, K, and both
//!   Hankel functions with a complex argument, including scaled variants.
//!   The root-level [`hankel1`] and [`hankel2`] accept real arguments, while
//!   [`checked_complex`] exposes the generic complex Hankel kernels.
//! - **Spherical Bessel:** [`spherical_jn`], [`spherical_yn`],
//!   [`spherical_in`], and [`spherical_kn`] optionally calculate derivatives.
//!   [`riccati_jn`] and [`riccati_yn`] return values and derivatives through
//!   a requested order.
//! - **Kelvin and Struve:** [`ber`], [`bei`], [`ker`], [`kei`], their derivatives
//!   and zeros, [`struve`] and [`modstruve`], with fallible equivalents for
//!   selected kernels in [`checked`].
//! - **Additional functions:** Cylindrical Bessel zeros ([`jn_zeros`],
//!   [`jnyn_zeros`], [`jnjnp_zeros`]), modified Bessel ratios ([`iv_ratio`],
//!   [`iv_ratio_c`], [`iv_ratioinv`]), integral families ([`itj0y0`],
//!   [`it2i0k0`], [`itstruve0`]), [`besselpoly`], and the Airy quartet ([`airy`],
//!   [`airye`]).
//!
//! # Domains and numerical behavior
//!
//! Arguments, orders, and normalizations are dimensionless. Real J and I
//! require integer order for negative real arguments; real Y and K generally
//! require positive arguments, with singular limits at zero. For invalid
//! domains, infallible real functions return IEEE NaN; some singular limits
//! return infinity. The [`checked`], [`checked_complex`], and [`complex`]
//! interfaces return errors for their respective fallible operations, including
//! invalid inputs, failure to converge, and overflow where supported.
//!
//! Numerical kernels use binary64 arithmetic. Scalar Bessel, Struve, and Kelvin
//! functions accept `f32` or `f64` and evaluate internally in `f64`, converting
//! the result back to the input type. Complex Bessel/Hankel functions also
//! support generic `Complex<T>` inputs. All iterations have finite limits, but
//! relative accuracy can deteriorate near zeros, branch cuts, turning points,
//! and derivatives with severe cancellation. The real derivatives support
//! orders up to 1024; spherical functions support integer orders through 4096.
//! Individual functions document other domain and count restrictions.
//!
//! This standalone library was extracted from the Bessel namespace in
//! [Sci](https://github.com/p-sira/sci). It covers all 79 original Bessel names.
//! [`kelvin_zeros`] returns all eight zero sequences, with
//! [`kelvin_zeros_for`] available to select one family. Fallible public
//! functions use crate-level [`Error`] and [`Result`] where applicable.
//!
//! # Mathematical conventions
//!
//! - The Bessel differential equation is `x² y″ + x y′ + (x² - ν²)y = 0`.
//! - `I_ν(z) = i^{-ν}J_ν(iz)` and `K_ν` is the decaying modified solution.
//! - `H_ν^{(1,2)} = J_ν ± iY_ν`.
//! - `j_n(x) = √(π/(2x)) J_{n+1/2}(x)`; `k_n` uses the
//!   SciPy convention `√(π/(2x)) K_{n+1/2}(x)`.
//! - Kelvin functions are the real and imaginary components of
//!   `J_0(x exp(3πi/4))` and `K_0(x exp(πi/4))`.
//!
//! Complex functions follow the principal branches used by AMOS. For defining formulas,
//! principal branches, and analytic continuation, see
//! [NIST DLMF Chapter 10](https://dlmf.nist.gov/10) and
//! [NIST DLMF Chapter 11](https://dlmf.nist.gov/11). For the modified Bessel
//! ratios, see Gautschi and Slavik (1978). Numerical validation uses independent
//! SciPy references and high-precision checks at extreme arguments.
//!
//! # Example
//!
//! ```
//! use bessel::{j0, jv, iv, spherical_jn};
//! assert!((j0(0.0_f64) - 1.0).abs() < 1e-14);
//! assert!((jv(0.0_f64, 1.0) - j0(1.0_f64)).abs() < 1e-14);
//! assert!((iv(0.0_f64, 0.0) - 1.0).abs() < 1e-14);
//! assert!((spherical_jn(0,1.0_f64,false) - 1.0_f64.sin()).abs() < 1e-14);
//! ```
#![forbid(unsafe_code)]
mod backend;
pub use backend::airy::Airy;
pub use backend::kelvin::{Kelvin, KelvinFunction};
pub use num_complex::{Complex, Complex64};

/// A validation or numerical evaluation failure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// A parameter or sequence size violates its documented domain.
    InvalidInput(&'static str),
    /// The numerical algorithm failed to converge or evaluate.
    NumericalFailure(&'static str),
    /// The result cannot be represented as a finite number.
    NonFiniteEvaluation,
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(message) | Self::NumericalFailure(message) => f.write_str(message),
            Self::NonFiniteEvaluation => f.write_str("evaluation returned a nonfinite value"),
        }
    }
}
impl std::error::Error for Error {}
/// A checked numerical result.
pub type Result<T> = std::result::Result<T, Error>;
/// Checked complex Hankel functions with generic real order and complex argument.
pub mod checked_complex {
    pub use crate::backend::bessel::{cyl_hankel_1, cyl_hankel_1e, cyl_hankel_2, cyl_hankel_2e};
}
/// The Airy quartet (Ai, Ai′, Bi, Bi′), using real Cephes and native complex kernels.
pub fn airy<T: num_traits::Float>(x: T) -> Result<Airy<T>> {
    if !x.is_finite() {
        return Err(Error::InvalidInput("airy requires a finite argument"));
    }
    backend::airy::real(x, false).map_err(Error::NumericalFailure)
}
/// Scaled Airy quartet (Ai, Ai′, Bi, Bi′). For real x, requires x>=0.
pub fn airye<T: num_traits::Float>(x: T) -> Result<Airy<T>> {
    if !x.is_finite() || x < T::zero() {
        return Err(Error::InvalidInput(
            "airye requires a finite nonnegative argument",
        ));
    }
    backend::airy::real(x, true).map_err(Error::NumericalFailure)
}
/// Checked Kelvin component and zero-finder functions.
pub mod checked {
    pub use crate::backend::kelvin::{
        bei, beip, ber, berp, kei, keip, kelvin, kelvin_zeros, ker, kerp,
    };
    pub use crate::backend::struve::{struve_h, struve_l};
}
macro_rules! scalar {
    ($($(#[$meta:meta])* $name:ident($($arg:ident: $ty:ty),*) => $body:expr;)+) => {$(
        $(#[$meta])*
        #[inline]
        pub fn $name($($arg: $ty),*) -> f64 { $body }
    )+};
}
macro_rules! generic_scalar {
    ($($(#[$meta:meta])* $name:ident($($arg:ident: f64),*) => $body:expr;)+) => {$(
        $(#[$meta])*
        #[inline]
        pub fn $name<T: num_traits::Float>($($arg: T),*) -> T {
            $(let $arg = $arg.to_f64().unwrap_or(f64::NAN);)*
            T::from($body).unwrap_or(T::nan())
        }
    )+};
}
generic_scalar! {
    /// First-kind Bessel J₀(x), for all finite real x.
    j0(x:f64) => backend::cephes64::j0(x);
    /// First-kind Bessel J₁(x), odd in x.
    j1(x:f64) => backend::cephes64::j1(x);
    /// Second-kind Bessel Y₀(x), for x>0; singular at zero.
    y0(x:f64) => backend::cephes64::y0(x);
    /// Second-kind Bessel Y₁(x), for x>0; singular at zero.
    y1(x:f64) => backend::cephes64::y1(x);
    /// First-kind Bessel Jν(x), where negative x requires integer ν.
    jv(v:f64,x:f64) => backend::cyl_bessel_j(v,x).unwrap_or(f64::NAN);
    /// Second-kind Bessel Yν(x), for x≥0.
    yv(v:f64,x:f64) => backend::cyl_bessel_y(v,x).unwrap_or(f64::NAN);
    /// Modified first-kind Bessel Iν(x), where negative x requires integer ν.
    iv(v:f64,x:f64) => backend::cyl_bessel_i(v,x).unwrap_or(f64::NAN);
    /// Exponentially scaled Iν(x): exp(-|x|) Iν(x).
    ive(v:f64,x:f64) => backend::cyl_bessel_ie(v,x).unwrap_or(f64::NAN);
    /// Modified second-kind Bessel Kν(x), for x≥0.
    kv(v:f64,x:f64) => backend::cyl_bessel_k(v,x).unwrap_or(f64::NAN);
    /// Exponentially scaled Kν(x): exp(x) Kν(x), for x≥0.
    kve(v:f64,x:f64) => backend::cyl_bessel_ke(v,x).unwrap_or(f64::NAN);
    /// Modified first-kind Bessel I₀(x).
    i0(x:f64) => backend::cephes64::i0(x);
    /// Modified first-kind Bessel I₁(x).
    i1(x:f64) => backend::cephes64::i1(x);
    /// Scaled I₀(x): exp(-|x|) I₀(x).
    i0e(x:f64) => backend::cephes64::i0e(x);
    /// Scaled I₁(x): exp(-|x|) I₁(x).
    i1e(x:f64) => backend::cephes64::i1e(x);
    /// Modified second-kind Bessel K₀(x), x≥0.
    k0(x:f64) => backend::cephes64::k0(x);
    /// Modified second-kind Bessel K₁(x), x≥0.
    k1(x:f64) => backend::cephes64::k1(x);
    /// Scaled K₀(x): exp(x) K₀(x), x≥0.
    k0e(x:f64) => backend::cephes64::k0e(x);
    /// Scaled K₁(x): exp(x) K₁(x), x≥0.
    k1e(x:f64) => backend::cephes64::k1e(x);
    /// Scaled first-kind cylindrical Bessel Jν(x), identical on the real axis.
    jve(v:f64,x:f64) => jv(v,x);
    /// Scaled second-kind cylindrical Bessel Yν(x), identical on the real axis.
    yve(v:f64,x:f64) => yv(v,x);
    /// Real Struve Hν(x), defining series in DLMF §11.2.
    struve(v:f64,x:f64) => backend::struve_h(v,x).unwrap_or(f64::NAN);
    /// Modified Struve Lν(x), defining series in DLMF §11.2.
    modstruve(v:f64,x:f64) => backend::struve_l(v,x).unwrap_or(f64::NAN);
    /// Kelvin ber(x), Re[J₀(x e^(3πi/4))].
    ber(x:f64) => backend::kelvin::ber(x).unwrap_or(f64::NAN);
    /// Kelvin bei(x), Im[J₀(x e^(3πi/4))].
    bei(x:f64) => backend::kelvin::bei(x).unwrap_or(f64::NAN);
    /// Kelvin ker(x), Re[K₀(x e^(πi/4))].
    ker(x:f64) => backend::kelvin::ker(x).unwrap_or(f64::NAN);
    /// Kelvin kei(x), Im[K₀(x e^(πi/4))].
    kei(x:f64) => backend::kelvin::kei(x).unwrap_or(f64::NAN);
    /// Derivative ber′(x).
    berp(x:f64) => backend::kelvin::berp(x).unwrap_or(f64::NAN);
    /// Derivative bei′(x).
    beip(x:f64) => backend::kelvin::beip(x).unwrap_or(f64::NAN);
    /// Derivative ker′(x).
    kerp(x:f64) => backend::kelvin::kerp(x).unwrap_or(f64::NAN);
    /// Derivative kei′(x).
    keip(x:f64) => backend::kelvin::keip(x).unwrap_or(f64::NAN);
}
scalar! {
    /// Bessel J of integer order n.
    jn(n:i32,x:f64) => jv(n as f64,x);
    /// Bessel Y of integer order n.
    yn(n:i32,x:f64) => yv(n as f64,x);
    /// Modified Bessel K of integer order n.
    kn(n:i32,x:f64) => kv(n as f64,x);
}
/// All Kelvin values and derivatives, in the order
/// `ber+i bei, ker+i kei, berp+i beip, kerp+i keip`.
/// Returns NaN components on a backend evaluation failure.
pub fn kelvin(x: f64) -> [Complex64; 4] {
    let nan = Complex64::new(f64::NAN, f64::NAN);
    backend::kelvin::kelvin(x)
        .map(|(a, b, c, d)| [a, b, c, d])
        .unwrap_or([nan; 4])
}
/// Return the first `nt` positive zeros of one Kelvin family.
///
/// `function` selects one of the eight real Kelvin functions and derivatives.
/// The supported domain is `1..=4096` zeros. Evaluations use scaled kernels
/// with bounded bracketing and safeguarded Halley refinement.
///
/// # Errors
///
/// [`Error::InvalidInput`] for an empty or excessive count;
/// [`Error::NumericalFailure`] when the solver cannot find a requested root.
///
/// # Example
///
/// ```
/// let roots = bessel::kelvin_zeros_for(bessel::KelvinFunction::Ber, 2)?;
/// assert!(roots[0] > 0.0 && roots[1] > roots[0]);
/// # Ok::<(), bessel::Error>(())
/// ```
pub fn kelvin_zeros_for(function: KelvinFunction, nt: usize) -> Result<Vec<f64>> {
    if !(1..=4096).contains(&nt) {
        return Err(Error::InvalidInput("Kelvin zero count must be in 1..=4096"));
    }
    let mut zeros = vec![0.0; nt];
    backend::kelvin::kelvin_zeros(function, &mut zeros).map_err(Error::NumericalFailure)?;
    Ok(zeros)
}

/// First `nt` positive zeros of ber, bei, ker, kei, ber′, bei′, ker′, kei′.
///
/// Returns eight vectors of equal length in the order listed above.
/// `nt` must lie in `1..=4096`; an invalid count or a numerical failure
/// returns [`Error`]. See [`kelvin_zeros_for`] for the per-family solver.
///
/// # Example
///
/// ```
/// let roots = bessel::kelvin_zeros(2)?;
/// assert_eq!(roots.len(), 8);
/// assert!(roots.iter().all(|family| family.len() == 2));
/// # Ok::<(), bessel::Error>(())
/// ```
pub fn kelvin_zeros(nt: usize) -> Result<[Vec<f64>; 8]> {
    Ok([
        ber_zeros(nt)?,
        bei_zeros(nt)?,
        ker_zeros(nt)?,
        kei_zeros(nt)?,
        berp_zeros(nt)?,
        beip_zeros(nt)?,
        kerp_zeros(nt)?,
        keip_zeros(nt)?,
    ])
}
/// First positive roots of ber(x).
pub fn ber_zeros(nt: usize) -> Result<Vec<f64>> {
    kelvin_zeros_for(KelvinFunction::Ber, nt)
}
/// First positive roots of bei(x).
pub fn bei_zeros(nt: usize) -> Result<Vec<f64>> {
    kelvin_zeros_for(KelvinFunction::Bei, nt)
}
/// First positive roots of ker(x).
pub fn ker_zeros(nt: usize) -> Result<Vec<f64>> {
    kelvin_zeros_for(KelvinFunction::Ker, nt)
}
/// First positive roots of kei(x).
pub fn kei_zeros(nt: usize) -> Result<Vec<f64>> {
    kelvin_zeros_for(KelvinFunction::Kei, nt)
}
/// First positive roots of ber′(x).
pub fn berp_zeros(nt: usize) -> Result<Vec<f64>> {
    kelvin_zeros_for(KelvinFunction::BerPrime, nt)
}
/// First positive roots of bei′(x).
pub fn beip_zeros(nt: usize) -> Result<Vec<f64>> {
    kelvin_zeros_for(KelvinFunction::BeiPrime, nt)
}
/// First positive roots of ker′(x).
pub fn kerp_zeros(nt: usize) -> Result<Vec<f64>> {
    kelvin_zeros_for(KelvinFunction::KerPrime, nt)
}
/// First positive roots of kei′(x).
pub fn keip_zeros(nt: usize) -> Result<Vec<f64>> {
    kelvin_zeros_for(KelvinFunction::KeiPrime, nt)
}

/// Checked complex cylindrical Bessel and Hankel functions with real order.
///
/// Signed imaginary zero selects the upper or lower rim of the negative-real
/// branch cut. All kernels evaluate internally in binary64 and convert to the
/// caller's floating-point type. Invalid, overflowing, and unconverged
/// evaluations return the shared [`Error`] type.
pub mod complex {
    use crate::backend::complex::{self as kernel, Kind};
    use crate::{Complex, Complex64, Error, Result};
    use num_traits::Float;

    fn evaluate<T: Float>(v: T, z: Complex<T>, family: Kind, scaled: bool) -> Result<Complex<T>> {
        if !v.is_finite() || !z.re.is_finite() || !z.im.is_finite() {
            return Err(Error::InvalidInput(
                "Bessel order and argument must be finite",
            ));
        }
        let v = v.to_f64().ok_or(Error::NonFiniteEvaluation)?;
        let z = Complex64::new(
            z.re.to_f64().ok_or(Error::NonFiniteEvaluation)?,
            z.im.to_f64().ok_or(Error::NonFiniteEvaluation)?,
        );
        let result = kernel::cylindrical(family, v, z, scaled)
            .map_err(|error| Error::NumericalFailure(error.message()))?;
        if !result.re.is_finite() || !result.im.is_finite() {
            return Err(Error::NonFiniteEvaluation);
        }
        let re = T::from(result.re).ok_or(Error::NonFiniteEvaluation)?;
        let im = T::from(result.im).ok_or(Error::NonFiniteEvaluation)?;
        if !re.is_finite() || !im.is_finite() {
            return Err(Error::NonFiniteEvaluation);
        }
        Ok(Complex::new(re, im))
    }

    macro_rules! bessel_complex {
        ($($(#[$meta:meta])* $name:ident, $kind:ident, $scaled:expr;)+) => {$(
            $(#[$meta])*
            pub fn $name<T: Float>(v: T, z: Complex<T>) -> Result<Complex<T>> {
                evaluate(v, z, Kind::$kind, $scaled)
            }
        )+};
    }

    bessel_complex! {
        /// Jν(z) on the principal branch.
        jv, J, false;
        /// Yν(z) on the principal branch.
        yv, Y, false;
        /// Iν(z) on the principal branch.
        iv, I, false;
        /// Kν(z) on the principal branch.
        kv, K, false;
        /// Scaled Jν(z), exp(-|Im z|) Jν(z).
        jve, J, true;
        /// Scaled Yν(z), exp(-|Im z|) Yν(z).
        yve, Y, true;
        /// Scaled Iν(z), exp(-|Re z|) Iν(z).
        ive, I, true;
        /// Scaled Kν(z), exp(z) Kν(z).
        kve, K, true;
        /// Hankel Hν¹(z).
        hankel1, Hankel1, false;
        /// Hankel Hν²(z).
        hankel2, Hankel2, false;
        /// Scaled Hν¹(z) = exp(-iz) Hν¹(z).
        hankel1e, Hankel1, true;
        /// Scaled Hν²(z) = exp(iz) Hν²(z).
        hankel2e, Hankel2, true;
    }
}
/// First-kind Hankel Hν¹(x), on the upper rim of the negative real cut.
pub fn hankel1(v: f64, x: f64) -> Complex64 {
    complex::hankel1(v, Complex64::new(x, 0.0)).unwrap_or(Complex64::new(f64::NAN, f64::NAN))
}
/// Second-kind Hankel Hν²(x), on the upper rim of the negative real cut.
pub fn hankel2(v: f64, x: f64) -> Complex64 {
    complex::hankel2(v, Complex64::new(x, 0.0)).unwrap_or(Complex64::new(f64::NAN, f64::NAN))
}
/// Scaled Hankel exp(-ix) Hν¹(x).
pub fn hankel1e(v: f64, x: f64) -> Complex64 {
    complex::hankel1e(v, Complex64::new(x, 0.0)).unwrap_or(Complex64::new(f64::NAN, f64::NAN))
}
/// Scaled Hankel exp(ix) Hν²(x).
pub fn hankel2e(v: f64, x: f64) -> Complex64 {
    complex::hankel2e(v, Complex64::new(x, 0.0)).unwrap_or(Complex64::new(f64::NAN, f64::NAN))
}
/// Derivative dⁿ/dxⁿ fν(x) from adjacent orders:
/// `2^{-n}Σ(k=0..n) C(n,k) s^k f_{ν-n+2k}(x)`,
/// with `s=-1` for J,Y and `+1` for I,K.
/// Orders above 1024 return NaN; cancellation may limit accuracy.
fn derivative(v: f64, x: f64, n: u32, f: fn(f64, f64) -> f64, sign: f64) -> f64 {
    if n > 1024 {
        return f64::NAN;
    }
    let mut weight = 2.0_f64.powi(-(n as i32));
    let mut sum = weight * f(v - f64::from(n), x);
    for i in 1..=n {
        weight *= sign * f64::from(n - i + 1) / f64::from(i);
        sum += weight * f(v - f64::from(n) + 2.0 * f64::from(i), x);
    }
    sum
}
/// n-th derivative of Jν(x) with respect to x (n≤1024).
pub fn jvp(v: f64, x: f64, n: u32) -> f64 {
    derivative(v, x, n, jv, -1.0)
}
/// n-th derivative of Yν(x) with respect to x (n≤1024).
pub fn yvp(v: f64, x: f64, n: u32) -> f64 {
    derivative(v, x, n, yv, -1.0)
}
/// n-th derivative of Iν(x) with respect to x (n≤1024).
pub fn ivp(v: f64, x: f64, n: u32) -> f64 {
    derivative(v, x, n, iv, 1.0)
}
/// n-th derivative of Kν(x) with respect to x (n≤1024).
pub fn kvp(v: f64, x: f64, n: u32) -> f64 {
    derivative(v, x, n, kv, 1.0) * if n.is_multiple_of(2) { 1.0 } else { -1.0 }
}
/// n-th derivative of Hankel Hν¹(x) with respect to x.
pub fn h1vp(v: f64, x: f64, n: u32) -> Complex64 {
    Complex64::new(jvp(v, x, n), yvp(v, x, n))
}
/// n-th derivative of Hankel Hν²(x) with respect to x.
pub fn h2vp(v: f64, x: f64, n: u32) -> Complex64 {
    Complex64::new(jvp(v, x, n), -yvp(v, x, n))
}

macro_rules! spherical {
    ($($(#[$meta:meta])* $name:ident($fun:ident,$jac:ident);)+)=>{$(
       $(#[$meta])*
       pub fn $name<T:num_traits::Float>(n:i64,x:T,derivative:bool)->T {
           if !(0..=4096).contains(&n) {return T::nan();}
           if derivative {backend::$jac(n as i32,x).unwrap_or(T::nan())}
           else {backend::$fun(n as i32,x).unwrap_or(T::nan())}
       }
    )+}
}
spherical! {
    /// Spherical j_n(x)=sqrt(π/(2x)) J_(n+1/2)(x).
    spherical_jn(sph_bessel_j,sph_bessel_j_jac);
    /// Spherical y_n(x)=sqrt(π/(2x)) Y_(n+1/2)(x).
    spherical_yn(sph_bessel_y,sph_bessel_y_jac);
    /// Spherical i_n(x)=sqrt(π/(2x)) I_(n+1/2)(x).
    spherical_in(sph_bessel_i,sph_bessel_i_jac);
    /// Spherical k_n(x)=sqrt(π/(2x)) K_(n+1/2)(x).
    spherical_kn(sph_bessel_k,sph_bessel_k_jac);
}
/// Riccati-Bessel functions x*j_n(x), and their first derivatives, through n.
pub fn riccati_jn(n: usize, x: f64) -> Result<(Vec<f64>, Vec<f64>)> {
    riccati(n, x, true)
}
/// Riccati-Bessel functions x*y_n(x), and their first derivatives, through n.
pub fn riccati_yn(n: usize, x: f64) -> Result<(Vec<f64>, Vec<f64>)> {
    riccati(n, x, false)
}
fn riccati(n: usize, x: f64, first: bool) -> Result<(Vec<f64>, Vec<f64>)> {
    if n > 4096 || !x.is_finite() || x <= 0.0 {
        return Err(Error::InvalidInput(
            "Riccati requires positive finite x and n<=4096",
        ));
    }
    Ok((0..=n)
        .map(|i| {
            let f = if first { spherical_jn } else { spherical_yn };
            let y = f(i as i64, x, false);
            (x * y, y + x * f(i as i64, x, true))
        })
        .unzip())
}
/// Bessel zero search kind (Jₙ, Jₙ′, Yₙ, Yₙ′).
#[derive(Clone, Copy)]
enum RootKind {
    J,
    JPrime,
    Y,
    YPrime,
}
fn bessel_roots(n: u32, nt: usize, kind: RootKind) -> Result<Vec<f64>> {
    if n > 256 || !(1..=256).contains(&nt) {
        return Err(Error::InvalidInput(
            "Bessel roots require n<=256 and count in 1..=256",
        ));
    }
    let v = f64::from(n);
    let f = |x: f64| match kind {
        RootKind::J => jv(v, x),
        RootKind::JPrime => jvp(v, x, 1),
        RootKind::Y => yv(v, x),
        RootKind::YPrime => yvp(v, x, 1),
    };
    let mut a = v.max(0.01);
    let mut fa = f(a);
    let mut zeros = Vec::with_capacity(nt);
    for _ in 0..(nt + 256) * 32 {
        let b = a + std::f64::consts::FRAC_PI_4;
        let fb = f(b);
        if fa.is_finite()
            && fb.is_finite()
            && fa != 0.0
            && fa.is_sign_positive() != fb.is_sign_positive()
        {
            let (mut low, mut high) = (a, b);
            let mut fl = fa;
            for _ in 0..100 {
                let mid = low + (high - low) * 0.5;
                let fm = f(mid);
                if !fm.is_finite() {
                    return Err(Error::NonFiniteEvaluation);
                }
                if fm == 0.0 || high - low <= 4.0 * f64::EPSILON * mid.max(1.0) {
                    low = mid;
                    high = mid;
                    break;
                }
                if fm.is_sign_positive() == fl.is_sign_positive() {
                    low = mid;
                    fl = fm;
                } else {
                    high = mid;
                }
            }
            zeros.push(low + (high - low) * 0.5);
            if zeros.len() == nt {
                return Ok(zeros);
            }
        }
        a = b;
        fa = fb;
    }
    Err(Error::NumericalFailure(
        "unable to bracket requested Bessel roots",
    ))
}
/// First positive zeros of cylindrical Jₙ(x).
pub fn jn_zeros(n: u32, nt: usize) -> Result<Vec<f64>> {
    bessel_roots(n, nt, RootKind::J)
}
/// First positive zeros of Jₙ′(x).
pub fn jnp_zeros(n: u32, nt: usize) -> Result<Vec<f64>> {
    bessel_roots(n, nt, RootKind::JPrime)
}
/// First positive zeros of cylindrical Yₙ(x).
pub fn yn_zeros(n: u32, nt: usize) -> Result<Vec<f64>> {
    bessel_roots(n, nt, RootKind::Y)
}
/// First positive zeros of Yₙ′(x).
pub fn ynp_zeros(n: u32, nt: usize) -> Result<Vec<f64>> {
    bessel_roots(n, nt, RootKind::YPrime)
}
/// All four real zero sequences in the order Jₙ, Jₙ′, Yₙ, Yₙ′.
pub fn jnyn_zeros(n: u32, nt: usize) -> Result<[Vec<f64>; 4]> {
    Ok([
        jn_zeros(n, nt)?,
        jnp_zeros(n, nt)?,
        yn_zeros(n, nt)?,
        ynp_zeros(n, nt)?,
    ])
}
/// Positive real Y₀ zeros and corresponding Y₁ values.
pub fn y0_zeros(nt: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let roots = yn_zeros(0, nt)?;
    let values = roots.iter().map(|&x| y1(x)).collect();
    Ok((roots, values))
}
/// Positive real Y₁ zeros and corresponding Y₁ derivatives.
pub fn y1_zeros(nt: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let roots = yn_zeros(1, nt)?;
    let values = roots.iter().map(|&x| yvp(1.0, x, 1)).collect();
    Ok((roots, values))
}
/// Positive real Y₁′ zeros and corresponding Y₁ values.
pub fn y1p_zeros(nt: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let roots = ynp_zeros(1, nt)?;
    let values = roots.iter().map(|&x| y1(x)).collect();
    Ok((roots, values))
}
/// First merged real zeros of Jₙ and Jₙ′ across increasing n.
///
/// Returns parallel vectors `[root, order, root_index, derivative_flag]`;
/// the J₀′ root at zero is included.
pub fn jnjnp_zeros(nt: usize) -> Result<[Vec<f64>; 4]> {
    if !(1..=128).contains(&nt) {
        return Err(Error::InvalidInput("merged root count must be in 1..=128"));
    }
    let mut entries: Vec<(f64, f64, f64, f64)> = vec![(0.0, 0.0, 0.0, 1.0)];
    let cutoff = ((nt as f64) * 8.0).sqrt() + 4.0;
    for n in 0..=cutoff.ceil() as u32 {
        for (kind, flag) in [(RootKind::J, 0.0), (RootKind::JPrime, 1.0)] {
            for (index, x) in bessel_roots(n, nt, kind)?.into_iter().enumerate() {
                entries.push((x, f64::from(n), (index + 1) as f64, flag));
            }
        }
    }
    entries.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| b.3.total_cmp(&a.3)));
    let mut start = 0;
    while start < entries.len() {
        let mut end = start + 1;
        while end < entries.len() && (entries[end].0 - entries[start].0).abs() <= 1e-12 {
            end += 1;
        }
        entries[start..end].sort_by(|a, b| a.3.total_cmp(&b.3));
        start = end;
    }
    let mut out: [Vec<f64>; 4] = std::array::from_fn(|_| Vec::with_capacity(nt));
    for (x, n, index, flag) in entries.into_iter().take(nt) {
        out[0].push(x);
        out[1].push(n);
        out[2].push(index);
        out[3].push(flag);
    }
    Ok(out)
}

/// Ratio `I_v(x)/I_(v-1)(x)`, using a scaled Perron continued fraction.
///
/// Requires v≥1/2 and x≥0. Returns the continuous limit at x=0,
/// 1 at x=+infinity, and NaN on invalid input or failed convergence.
pub fn iv_ratio<T: num_traits::Float>(v: T, x: T) -> T {
    let (ratio, _) = backend::iv_ratio_pair(
        v.to_f64().unwrap_or(f64::NAN),
        x.to_f64().unwrap_or(f64::NAN),
    );
    T::from(ratio).unwrap_or(T::nan())
}
/// Complement `1 - I_v(x)/I_(v-1)(x)`, calculated directly to
/// preserve small values where subtracting a rounded ratio would lose them.
pub fn iv_ratio_c<T: num_traits::Float>(v: T, x: T) -> T {
    let (_, complement) = backend::iv_ratio_pair(
        v.to_f64().unwrap_or(f64::NAN),
        x.to_f64().unwrap_or(f64::NAN),
    );
    T::from(complement).unwrap_or(T::nan())
}
/// Inverse of `I_v(x)/I_(v-1)(x)` for v≥1/2 and ratio r in `\[0,1]`.
///
/// Selects nonnegative x. Bisection is bounded to 128 evaluations; r=1
/// returns positive infinity, r=0 returns zero, invalid parameters return NaN.
pub fn iv_ratioinv<T: num_traits::Float>(v: T, r: T) -> T {
    let (v, r) = (
        v.to_f64().unwrap_or(f64::NAN),
        r.to_f64().unwrap_or(f64::NAN),
    );
    if !v.is_finite() || v < 0.5 || !(0.0..=1.0).contains(&r) {
        return T::nan();
    }
    if r == 0.0 {
        return T::zero();
    }
    if r == 1.0 {
        return T::infinity();
    }
    if v == 0.5 {
        return T::from(r.atanh()).unwrap_or(T::nan());
    }
    let mut high = (v * r) / (1.0 - r) * (2.0 / (1.0 + r));
    if high.is_infinite() {
        return T::infinity();
    }
    let mut low = 0.0;
    for _ in 0..128 {
        let mid = low + (high - low) * 0.5;
        let (value, complement) = backend::iv_ratio_pair(v, mid);
        if value.is_nan() {
            return T::nan();
        }
        let residual = if r > 0.5 {
            (1.0 - r) - complement
        } else {
            value - r
        };
        if residual == 0.0 || high - low <= 4.0 * f64::EPSILON * mid {
            return T::from(mid).unwrap_or(T::nan());
        }
        if residual < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    T::nan()
}
/// Weighted Bessel integral `∫₀¹ t^lambda J_v(2 a t) dt`.
///
/// The defining series converges for a≥0 and lambda+v>−1.
/// Returns NaN for invalid nonfinite inputs. Uses a bounded series with
/// up to 1000 terms; loss of significance is possible for large a.
pub fn besselpoly(a: f64, lambda: f64, mut v: f64) -> f64 {
    if !a.is_finite() || !lambda.is_finite() || !v.is_finite() || a < 0.0 || lambda + v <= -1.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return if v == 0.0 { 1.0 / (lambda + 1.0) } else { 0.0 };
    }
    let mut sign = 1.0;
    if v < 0.0 && v == v.floor() {
        v = -v;
        if v % 2.0 != 0.0 {
            sign = -1.0;
        }
    }
    let mut term = (v * a.ln()).exp() / (backend::cephes64::gamma(v + 1.0) * (lambda + v + 1.0));
    let mut sum = 0.0;
    for m in 0..1000 {
        sum += term;
        let m = m as f64;
        term *= -a * a * (lambda + v + 1.0 + 2.0 * m)
            / ((v + m + 1.0) * (m + 1.0) * (lambda + v + 3.0 + 2.0 * m));
        if term == 0.0 {
            break;
        }
    }
    sign * sum
}

/// Integrals `(∫₀ˣ J₀(t) dt, ∫₀ˣ Y₀(t) dt)` for x ≥ 0.
/// Y₀'s logarithmic endpoint singularity is integrable.
pub fn itj0y0(x: f64) -> (f64, f64) {
    paired_integral(x, backend::integrals::itjya, true)
}
/// Integrals `(∫₀ˣ (1-J₀(t))/t dt, ∫ₓ^∞ Y₀(t)/t dt)` for x > 0.
pub fn it2j0y0(x: f64) -> (f64, f64) {
    paired_integral(x, backend::integrals::ittjya, false)
}
/// Integrals `(∫₀ˣ I₀(t) dt, ∫₀ˣ K₀(t) dt)` for x ≥ 0.
/// K₀'s logarithmic endpoint singularity is integrable.
pub fn iti0k0(x: f64) -> (f64, f64) {
    paired_integral(x, backend::integrals::itika, true)
}
/// Integrals `(∫₀ˣ (I₀(t)-1)/t dt, ∫ₓ^∞ K₀(t)/t dt)` for x > 0.
///
/// For x ≥ 4, evaluates the exponentially decaying K₀ tail with a scaled
/// kernel to avoid loss of accuracy in the original specfun approximation.
pub fn it2i0k0(x: f64) -> (f64, f64) {
    if x.is_finite() && x >= 4.0 {
        backend::integrals::it2i0k0_corrected(x)
    } else {
        paired_integral(x, backend::integrals::ittika, false)
    }
}
/// Integral of the Struve function H₀ from zero to `x`.
///
/// # Definition and conventions
///
/// Returns `∫₀ˣ H₀(t) dt`. Since H₀ is odd, its integral is even in `x`.
/// The argument and the result are dimensionless, with real Struve H₀.
///
/// # Domain and special values
///
/// All finite real arguments are supported; `itstruve0(±0.0) = +0.0`.
/// Either infinity returns positive infinity (logarithmic growth), and NaN
/// propagates. This infallible function has no separate error result.
///
/// # Numerical behavior
///
/// Uses a convergent series for `|x| <= 10`, quadrature for `10 < |x| <= 512`,
/// then an asymptotic expansion including the oscillatory Y₀ tail. The large
/// argument branch avoids intermediate overflow and phase loss in `x - π/4`.
///
/// # Example
///
/// ```
/// assert_eq!(bessel::itstruve0(0.0), 0.0);
/// assert_eq!(bessel::itstruve0(-600.0), bessel::itstruve0(600.0));
/// assert!(bessel::itstruve0(1.0e308).is_finite());
/// ```
///
/// # References
///
/// [DLMF 11.6.3](https://dlmf.nist.gov/11.6.E3) and
/// [DLMF 10.17.4](https://dlmf.nist.gov/10.17.E4).
pub fn itstruve0(x: f64) -> f64 {
    backend::integrals::itstruve0_corrected(x)
}
/// Integral `∫ₓ^∞ H₀(t)/t dt` for x ≥ 0.
pub fn it2struve0(x: f64) -> f64 {
    if !x.is_finite() {
        f64::NAN
    } else if x < 0.0 {
        std::f64::consts::PI - backend::integrals::itth0(-x)
    } else {
        backend::integrals::itth0(x)
    }
}
/// Integral `∫₀ˣ L₀(t) dt`, even in x; L₀ is the modified Struve function.
pub fn itmodstruve0(x: f64) -> f64 {
    if !x.is_finite() {
        f64::NAN
    } else {
        backend::integrals::itsl0(x.abs())
    }
}

/// Invoke a bundled specfun integral pair without backend aliases.
fn paired_integral(x: f64, evaluate: fn(f64, &mut [f64], &mut [f64]), odd: bool) -> (f64, f64) {
    if !x.is_finite() {
        return (f64::NAN, f64::NAN);
    }
    let (mut a, mut b) = ([0.0], [0.0]);
    evaluate(x.abs(), &mut a, &mut b);
    if x < 0.0 {
        (if odd { -a[0] } else { a[0] }, f64::NAN)
    } else {
        (a[0], b[0])
    }
}
