//! Independent 100-digit mpmath references for the Chebyshev fast path.
//! Values use x²/π · ₂F₃(1,1;2,3/2,3/2;-x²/4), not Sci outputs.
#![allow(clippy::excessive_precision)]

use bessel::itstruve0;

#[test]
fn moderate_struve_integrals_match_high_precision_reference() {
    for (x, expected) in [
        (10.0_f64, 2.5189576695087033443_f64),
        (10.00001_f64, 2.518958856932779578_f64),
        (10.125_f64, 2.5317996937640715336_f64),
        (12.75_f64, 2.2541237096425421946_f64),
        (15.5_f64, 2.6583579062487770507_f64),
        (19.9999_f64, 2.5484422537438187876_f64),
        (20.0_f64, 2.548451692293957187_f64),
        (23.5_f64, 2.9451387822967178826_f64),
        (23.9999_f64, 2.8918158419558812925_f64),
        (24.0_f64, 2.8918032068700988167_f64),
        (24.0001_f64, 2.8917905712427223017_f64),
    ] {
        let actual = itstruve0(x);
        assert!(
            (actual / expected - 1.0).abs() <= 1.5e-14,
            "itstruve0({x}) = {actual:.17e}; expected {expected:.17e}"
        );
        assert_eq!(itstruve0(-x), actual, "integral must be even at {x}");
    }
}
