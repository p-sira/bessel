//! Independent high-precision references across the integral's asymptotic
//! crossover. Values come from the exact identity
//! x²/pi * 2F3(1, 1; 2, 3/2, 3/2; -x²/4) evaluated at 80 decimal digits.
#![allow(clippy::excessive_precision)]

use bessel::{itstruve0, struve};

#[test]
fn integrated_struve_remains_accurate_across_asymptotic_crossover() {
    let reference = [
        (150.0, 3.99961290534207903e+00),
        (175.0, 4.10560902108119219e+00),
        (180.0, 4.17357150445117409e+00),
        (190.0, 4.10542748807864477e+00),
        (195.0, 4.11745738319333210e+00),
        (199.0, 4.23274811492877756e+00),
        (199.75, 4.20948578511452265e+00),
        (199.999, 4.19738348392023841e+00),
        (199.999999, 4.19733245991496418e+00),
        (200.0, 4.19733240883221637e+00),
        (200.000001, 4.19733235774945257e+00),
        (200.001, 4.19728131842645702e+00),
        (200.5, 4.17103318247960786e+00),
        (201.0, 4.14781116133519490e+00),
        (203.0, 4.16834197312979082e+00),
        (207.5, 4.16055514369328616e+00),
        (211.0, 4.26870669052992735e+00),
        (220.0, 4.20128233601793344e+00),
        (235.0, 4.29288505785656405e+00),
        (250.0, 4.34995382661899388e+00),
        (275.0, 4.41463210890577518e+00),
        (300.0, 4.47323362427752969e+00),
        (330.0, 4.53540296467747606e+00),
        (360.0, 4.53580870420704674e+00),
        (400.0, 4.66186132678924992e+00),
        (450.0, 4.73559754155757506e+00),
        (500.0, 4.79917170271420268e+00),
        (511.999, 4.80299140727448481e+00),
        (512.0, 4.80301950532947597e+00),
        (512.001, 4.80304758049014247e+00),
        (513.0, 4.81633488990365422e+00),
    ];
    for (x, expected) in reference {
        let got = itstruve0(x);
        let tolerance = if x < 200.0 { 4e-14 } else { 4e-15 };
        assert!(
            (got - expected).abs() <= tolerance * expected,
            "itstruve0({x})={got:.17e}, reference={expected:.17e}"
        );
        assert_eq!(itstruve0(-x), got, "evenness at {x}");
    }
}

#[test]
fn integrated_struve_derivative_is_continuous_at_crossover() {
    let x = 200.0;
    let h = 1e-3;
    let slope = (itstruve0(x + h) - itstruve0(x - h)) / (2.0 * h);
    let expected = struve(0.0, x);
    assert!(
        (slope - expected).abs() <= 5e-8,
        "slope across x=200 is {slope}, H0(200) is {expected}"
    );
}
