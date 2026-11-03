//! Independent high-precision references for the moderate complex-argument
//! Steed CF2 K fast path (mpmath 1.4.1, 75-digit working precision).
use bessel::{self as b, Complex64};

#[test]
fn complex_k_reduced_order_cf2_matches_mpmath() {
    let z = Complex64::new(10.0, 2.0);
    for (v, expected, expected_scaled) in [
        (
            0.0,
            Complex64::new(-8.839_035_332_607_715e-6, -1.523_563_981_571_953e-5),
            Complex64::new(0.386_169_422_816_089_8, -0.037_379_986_717_081_34),
        ),
        (
            0.5,
            Complex64::new(-8.975_266_415_492_014e-6, -1.539_228_846_578_46e-5),
            Complex64::new(0.390_555_604_681_202, -0.038_672_626_006_468_58),
        ),
        (
            2.0,
            Complex64::new(-1.125_624_230_042_379e-5, -1.792_642_790_774_655e-5),
            Complex64::new(0.462_218_873_535_785_1, -0.061_128_858_780_655_394),
        ),
        (
            7.25,
            Complex64::new(-1.509_927_209_634_373_2e-4, -1.076_563_578_741_732_4e-4),
            Complex64::new(3.540_243_194_380_512_4, -2.037_369_669_112_786),
        ),
    ] {
        for (name, result, reference) in [
            ("kv", b::complex::kv(v, z), expected),
            ("kve", b::complex::kve(v, z), expected_scaled),
        ] {
            let actual = result.unwrap();
            let relative = (actual - reference).norm() / reference.norm();
            assert!(
                relative < 3e-12,
                "{name}({v}, {z}): got {actual}, reference {reference}, rel={relative:e}"
            );
        }
    }
}
