use bessel::{self as b, Complex, Complex64, Error, KelvinFunction};

#[test]
fn kelvin_all_families_match_the_individual_selectors() {
    let all = b::kelvin_zeros(3).unwrap();
    let families = [
        KelvinFunction::Ber,
        KelvinFunction::Bei,
        KelvinFunction::Ker,
        KelvinFunction::Kei,
        KelvinFunction::BerPrime,
        KelvinFunction::BeiPrime,
        KelvinFunction::KerPrime,
        KelvinFunction::KeiPrime,
    ];
    for (family, (kind, roots)) in families.into_iter().zip(all.iter()).enumerate() {
        assert_eq!(
            roots,
            &b::kelvin_zeros_for(kind, 3).unwrap(),
            "family {family}"
        );
        assert_eq!(roots.len(), 3);
        assert!(roots.iter().all(|x| *x > 0.0 && x.is_finite()));
        assert!(roots.windows(2).all(|w| w[0] < w[1]));
    }
    assert!((all[0][0] - 2.848_917_820_795_139).abs() < 5e-12);
    assert!(matches!(b::kelvin_zeros(0), Err(Error::InvalidInput(_))));
    assert!(matches!(b::kelvin_zeros(4097), Err(Error::InvalidInput(_))));
    assert!(matches!(
        b::kelvin_zeros_for(KelvinFunction::Ker, 0),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn max_kelvin_count_uses_scaled_bracketed_solver() {
    // Unscaled Kelvin functions overflow or underflow before the 4096th
    // root. Every family must keep the full sequence via scaled kernels.
    let all = b::kelvin_zeros(4096).unwrap();
    for (family, roots) in all.into_iter().enumerate() {
        assert_eq!(roots.len(), 4096, "family {family}");
        assert!(roots.iter().all(|x| x.is_finite()), "family {family}");
        assert!(roots.windows(2).all(|w| w[1] > w[0]), "family {family}");
        assert!(roots[4095] > 18_000.0, "family {family}");
        let spacing = roots[4095] - roots[4094];
        assert!(
            (spacing - std::f64::consts::PI * 2.0_f64.sqrt()).abs() < 1e-6,
            "family {family}"
        );
    }
}

#[test]
fn generic_scalar_functions_accept_f32_and_f64() {
    let x = 1.25_f32;
    let v = 0.75_f32;
    for (got, expected) in [
        (b::j0(x), b::j0(1.25_f64) as f32),
        (b::jv(v, x), b::jv(0.75_f64, 1.25) as f32),
        (b::yv(v, x), b::yv(0.75_f64, 1.25) as f32),
        (b::ive(v, x), b::ive(0.75_f64, 1.25) as f32),
        (b::kve(v, x), b::kve(0.75_f64, 1.25) as f32),
        (b::struve(v, x), b::struve(0.75_f64, 1.25) as f32),
        (b::ber(x), b::ber(1.25_f64) as f32),
        (b::keip(x), b::keip(1.25_f64) as f32),
    ] {
        assert_eq!(got, expected);
    }
    let zero = b::jv(8.0_f32, 0.01_f32);
    assert!(zero > 0.0 && zero.is_finite());
    assert!(matches!(b::airye(-1.0_f32), Err(Error::InvalidInput(_))));
    assert!(matches!(b::jn_zeros(0, 0), Err(Error::InvalidInput(_))));
    assert!(matches!(b::riccati_jn(5, 0.0), Err(Error::InvalidInput(_))));
}

#[test]
fn complex_branch_cut_uses_signed_imaginary_zero() {
    let order = 0.75;
    let upper = Complex64::new(-3.0, 0.0);
    let lower = Complex64::new(-3.0, -0.0);
    let check = |got: Complex64, want: Complex64| {
        assert!(
            (got - want).norm() <= 5e-13 * want.norm().max(1.0),
            "got {got:?}, want {want:?}"
        );
    };
    for kernel in [
        b::complex::jv,
        b::complex::yv,
        b::complex::iv,
        b::complex::kv,
        b::complex::jve,
        b::complex::yve,
        b::complex::ive,
        b::complex::kve,
    ] {
        check(
            kernel(order, lower).unwrap(),
            kernel(order, upper).unwrap().conj(),
        );
    }
    type BesselComplex = fn(f64, Complex64) -> b::Result<Complex64>;
    for (first, second) in [
        (
            b::complex::hankel1 as BesselComplex,
            b::complex::hankel2 as BesselComplex,
        ),
        (
            b::complex::hankel1e as BesselComplex,
            b::complex::hankel2e as BesselComplex,
        ),
    ] {
        check(
            first(order, lower).unwrap(),
            second(order, upper).unwrap().conj(),
        );
        check(
            second(order, lower).unwrap(),
            first(order, upper).unwrap().conj(),
        );
    }
    let z = Complex::new(1.0_f32, 0.5_f32);
    let y: Complex<f32> = b::complex::hankel1(0.5_f32, z).unwrap();
    assert!(y.re.is_finite() && y.im.is_finite());
    assert!(matches!(
        b::complex::jv(f64::NAN, upper),
        Err(Error::InvalidInput(_))
    ));
}
