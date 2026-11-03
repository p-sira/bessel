//! Translated legacy specfun integral kernels; source attribution in LICENSES/scipy-specfun.txt.
#![allow(
    unused_variables,
    unused_assignments,
    unused_mut,
    unused_parens,
    clippy::all
)]
pub(crate) fn itika(mut x: f64, ti: &mut [f64], tk: &mut [f64]) {
    let mut k: i32 = 0;
    let mut rc1: f64 = 0.0;
    let mut rc2: f64 = 0.0;
    let mut e0: f64 = 0.0;
    let mut b1: f64 = 0.0;
    let mut b2: f64 = 0.0;
    let mut rs: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut tw: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut a = [
        0.625f64,
        1.0078125f64,
        2.5927734375f64,
        9.1868591308594f64,
        4.1567974090576e+1f64,
        2.2919635891914e+2f64,
        1.491504060477e+3f64,
        1.1192354495579e+4f64,
        9.515939374212e+4f64,
        9.0412425769041e+5f64,
    ];
    let mut pi: f64 = 3.141592653589793f64;
    let mut el: f64 = 0.5772156649015329f64;
    if (x == 0.0f64) {
        {
            ti[0] = 0.0f64;
            tk[0] = 0.0f64;
            return;
        }
    } else {
        if (x < 20.0f64) {
            {
                x2 = (x * x);
                ti[0] = 1.0f64;
                r = 1.0f64;
                k = 1;
                while (k <= 50) {
                    {
                        r = (((((0.25f64 * r) * (((2 * k) as f64) - 1.0f64))
                            / (((2 * k) as f64) + 1.0f64))
                            / ((k * k) as f64))
                            * x2);
                        ti[0] += r;
                        if ((r / ti[0]).abs() < 1.0e-12f64) {
                            {
                                break;
                            }
                        }
                    }
                    {
                        k += 1;
                        k
                    };
                }
                ti[0] *= x;
            }
        } else {
            {
                x2 = 0.0f64;
                ti[0] = 1.0f64;
                r = 1.0f64;
                k = 1;
                while (k <= 10) {
                    {
                        r = (r / x);
                        ti[0] += (a[(k - 1) as usize] * r);
                    }
                    {
                        k += 1;
                        k
                    };
                }
                rc1 = (1.0f64 / ((2.0f64 * pi) * x).sqrt());
                ti[0] = ((rc1 * (x).exp()) * ti[0]);
            }
        }
    }
    if (x < 12.0f64) {
        {
            e0 = (el + (x / 2.0f64).ln());
            b1 = (1.0f64 - e0);
            b2 = 0.0f64;
            rs = 0.0f64;
            r = 1.0f64;
            tw = 0.0f64;
            k = 1;
            while (k <= 50) {
                {
                    r = (((((0.25f64 * r) * (((2 * k) as f64) - 1.0f64))
                        / (((2 * k) as f64) + 1.0f64))
                        / ((k * k) as f64))
                        * x2);
                    b1 += (r * ((1.0f64 / (((2 * k) + 1) as f64)) - e0));
                    rs += (1.0f64 / ((k) as f64));
                    b2 += (r * rs);
                    tk[0] = (b1 + b2);
                    if (((tk[0] - tw) / tk[0]).abs() < 1.0e-12f64) {
                        {
                            break;
                        }
                    }
                    tw = tk[0];
                }
                {
                    k += 1;
                    k
                };
            }
            tk[0] *= x;
        }
    } else {
        {
            tk[0] = 1.0f64;
            r = 1.0f64;
            k = 1;
            while (k <= 10) {
                {
                    r = ((-r) / x);
                    tk[0] = (tk[0] + (a[(k - 1) as usize] * r));
                }
                {
                    k += 1;
                    k
                };
            }
            rc2 = (pi / (2.0f64 * x)).sqrt();
            tk[0] = ((pi / 2.0f64) - ((rc2 * tk[0]) * (-x).exp()));
        }
    }
    return;
}

pub(crate) fn itjya(mut x: f64, tj: &mut [f64], ty: &mut [f64]) {
    let mut k: i32 = 0;
    let mut a = [0.0; (18) as usize];
    let mut a0: f64 = 0.0;
    let mut a1: f64 = 0.0;
    let mut af: f64 = 0.0;
    let mut bf: f64 = 0.0;
    let mut bg: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut r2: f64 = 0.0;
    let mut rc: f64 = 0.0;
    let mut rs: f64 = 0.0;
    let mut ty1: f64 = 0.0;
    let mut ty2: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut xp: f64 = 0.0;
    let mut pi: f64 = 3.141592653589793f64;
    let mut el: f64 = 0.5772156649015329f64;
    let mut eps: f64 = 1.0e-12f64;
    if (x == 0.0f64) {
        {
            tj[0] = 0.0f64;
            ty[0] = 0.0f64;
        }
    } else {
        if (x <= 20.0f64) {
            {
                x2 = (x * x);
                tj[0] = x;
                r = x;
                k = 1;
                while (k < 61) {
                    {
                        r = ((((((-0.25f64) * r) * ((2.0f64 * ((k) as f64)) - 1.0f64))
                            / ((2.0f64 * ((k) as f64)) + 1.0f64))
                            / ((k * k) as f64))
                            * x2);
                        tj[0] += r;
                        if ((r).abs() < ((tj[0]).abs() * eps)) {
                            {
                                break;
                            }
                        }
                    }
                    {
                        k += 1;
                        k
                    };
                }
                ty1 = ((el + (x / 2.0f64).ln()) * tj[0]);
                rs = 0.0f64;
                ty2 = 1.0f64;
                r = 1.0f64;
                k = 1;
                while (k < 61) {
                    {
                        r = ((((((-0.25f64) * r) * ((2.0f64 * ((k) as f64)) - 1.0f64))
                            / ((2.0f64 * ((k) as f64)) + 1.0f64))
                            / ((k * k) as f64))
                            * x2);
                        rs += (1.0f64 / ((k) as f64));
                        r2 = (r * (rs + (1.0f64 / ((2.0f64 * ((k) as f64)) + 1.0f64))));
                        ty2 += r2;
                        if ((r2).abs() < ((ty2).abs() * eps)) {
                            {
                                break;
                            }
                        }
                    }
                    {
                        k += 1;
                        k
                    };
                }
                ty[0] = (((ty1 - (x * ty2)) * 2.0f64) / pi);
            }
        } else {
            {
                a0 = 1.0f64;
                a1 = (5.0f64 / 8.0f64);
                a[(0) as usize] = a1;
                let mut k: i32 = 1;
                while (k <= 16) {
                    {
                        af = (((((1.5f64 * (((k) as f64) + 0.5f64))
                            * (((k) as f64) + (5.0f64 / 6.0f64)))
                            * a1)
                            - ((((0.5f64 * (((k) as f64) + 0.5f64)) * (((k) as f64) + 0.5f64))
                                * (((k) as f64) - 0.5f64))
                                * a0))
                            / (((k) as f64) + 1.0f64));
                        a[(k) as usize] = af;
                        a0 = a1;
                        a1 = af;
                    }
                    {
                        k += 1;
                        k
                    };
                }
                bf = 1.0f64;
                r = 1.0f64;
                let mut k: i32 = 1;
                while (k <= 8) {
                    {
                        r = ((-r) / (x * x));
                        bf += (a[((2 * k) - 1) as usize] * r);
                    }
                    {
                        k += 1;
                        k
                    };
                }
                bg = (a[(0) as usize] / x);
                r = (1.0f64 / x);
                let mut k: i32 = 1;
                while (k <= 8) {
                    {
                        r = ((-r) / (x * x));
                        bg += (a[(2 * k) as usize] * r);
                    }
                    {
                        k += 1;
                        k
                    };
                }
                xp = (x + (0.25f64 * pi));
                rc = (2.0f64 / (pi * x)).sqrt();
                tj[0] = (1.0f64 - (rc * ((bf * (xp).cos()) + (bg * (xp).sin()))));
                ty[0] = (rc * ((bg * (xp).cos()) - (bf * (xp).sin())));
            }
        }
    }
    return;
}

pub(crate) fn ittika(mut x: f64, tti: &mut [f64], ttk: &mut [f64]) {
    let mut k: i32 = 0;
    let mut b1: f64 = 0.0;
    let mut e0: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut r2: f64 = 0.0;
    let mut rc: f64 = 0.0;
    let mut rs: f64 = 0.0;
    let mut pi: f64 = 3.141592653589793f64;
    let mut el: f64 = 0.5772156649015329f64;
    let mut c = [
        1.625f64,
        4.1328125f64,
        1.45380859375e1f64,
        6.553353881835e1f64,
        3.6066157150269e2f64,
        2.3448727161884e3f64,
        1.7588273098916e4f64,
        1.4950639538279e5f64,
    ];
    if (x == 0.0f64) {
        {
            tti[0] = 0.0f64;
            ttk[0] = 1.0e+300f64;
            return;
        }
    }
    if (x < 700.0f64) {
        {
            tti[0] = 1.0f64;
            r = 1.0f64;
            let mut k: i32 = 2;
            while (k <= 1000) {
                {
                    r = (((((0.25f64 * r) * (((k) as f64) - 1.0f64)) / (((k * k) * k) as f64))
                        * x)
                        * x);
                    tti[0] += r;
                    if ((r / tti[0]).abs() < f64::EPSILON) {
                        {
                            break;
                        }
                    }
                }
                {
                    k += 1;
                    k
                };
            }
            tti[0] *= ((0.125f64 * x) * x);
        }
    } else {
        {
            tti[0] = 1.0f64;
            r = 1.0f64;
            let mut k: i32 = 1;
            while (k <= 8) {
                {
                    r = (r / x);
                    tti[0] += (c[(k - 1) as usize] * r);
                }
                {
                    k += 1;
                    k
                };
            }
            rc = (x * ((2.0f64 * pi) * x).sqrt());
            tti[0] = ((tti[0] * (x).exp()) / rc);
        }
    }
    if (x <= 12.0f64) {
        {
            e0 = (((((0.5f64 * (x / 2.0f64).ln()) + el) * (x / 2.0f64).ln())
                + ((pi * pi) / 24.0f64))
                + ((0.5f64 * el) * el));
            b1 = (1.5f64 - (el + (x / 2.0f64).ln()));
            rs = 1.0f64;
            r = 1.0f64;
            k = 2;
            while (k <= 50) {
                {
                    r = (((((0.25f64 * r) * (((k) as f64) - 1.0f64)) / (((k * k) * k) as f64))
                        * x)
                        * x);
                    rs += (1.0f64 / ((k) as f64));
                    r2 = (r
                        * ((rs + (1.0f64 / (2.0f64 * ((k) as f64)))) - (el + (x / 2.0f64).ln())));
                    b1 += r2;
                    if ((r2 / b1).abs() < 1.0e-12f64) {
                        {
                            break;
                        }
                    }
                }
                {
                    k += 1;
                    k
                };
            }
            ttk[0] = (e0 - (((0.125f64 * x) * x) * b1));
        }
    } else {
        {
            ttk[0] = 1.0f64;
            r = 1.0f64;
            k = 1;
            while (k <= 8) {
                {
                    r = ((-r) / x);
                    ttk[0] += (c[(k - 1) as usize] * r);
                }
                {
                    k += 1;
                    k
                };
            }
            rc = (x * ((2.0f64 / pi) * x).sqrt());
            ttk[0] = ((ttk[0] * (-x).exp()) / rc);
        }
    }
    return;
}

pub(crate) fn ittjya(mut x: f64, ttj: &mut [f64], tty: &mut [f64]) {
    let mut k: i32 = 0;
    let mut l: i32 = 0;
    let mut a0: f64 = 0.0;
    let mut bj0: f64 = 0.0;
    let mut bj1: f64 = 0.0;
    let mut by0: f64 = 0.0;
    let mut by1: f64 = 0.0;
    let mut e0: f64 = 0.0;
    let mut b1: f64 = 0.0;
    let mut g0: f64 = 0.0;
    let mut g1: f64 = 0.0;
    let mut px: f64 = 0.0;
    let mut qx: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut r0: f64 = 0.0;
    let mut r1: f64 = 0.0;
    let mut r2: f64 = 0.0;
    let mut rs: f64 = 0.0;
    let mut t: f64 = 0.0;
    let mut vt: f64 = 0.0;
    let mut xk: f64 = 0.0;
    let mut pi: f64 = 3.141592653589793f64;
    let mut el: f64 = 0.5772156649015329f64;
    if (x == 0.0f64) {
        {
            ttj[0] = 0.0f64;
            tty[0] = (-1.0e+300f64);
        }
    } else {
        if (x <= 20.0f64) {
            {
                ttj[0] = 1.0f64;
                r = 1.0f64;
                k = 2;
                while (k <= 100) {
                    {
                        r = ((((((-0.25f64) * r) * (((k) as f64) - 1.0f64))
                            / (((k * k) * k) as f64))
                            * x)
                            * x);
                        ttj[0] += r;
                        if ((r).abs() < ((ttj[0]).abs() * 1.0e-12f64)) {
                            {
                                break;
                            }
                        }
                    }
                    {
                        k += 1;
                        k
                    };
                }
                ttj[0] *= ((0.125f64 * x) * x);
                e0 = ((0.5f64 * (((pi * pi) / 6.0f64) - (el * el)))
                    - (((0.5f64 * (x / 2.0f64).ln()) + el) * (x / 2.0f64).ln()));
                b1 = ((el + (x / 2.0f64).ln()) - 1.5f64);
                rs = 1.0f64;
                r = (-1.0f64);
                k = 2;
                while (k <= 100) {
                    {
                        r = ((((((-0.25f64) * r) * (((k) as f64) - 1.0f64))
                            / (((k * k) * k) as f64))
                            * x)
                            * x);
                        rs = (rs + (1.0f64 / ((k) as f64)));
                        r2 = (r
                            * ((rs + (1.0f64 / (2.0f64 * ((k) as f64))))
                                - (el + (x / 2.0f64).ln())));
                        b1 = (b1 + r2);
                        if ((r2).abs() < ((b1).abs() * 1.0e-12f64)) {
                            {
                                break;
                            }
                        }
                    }
                    {
                        k += 1;
                        k
                    };
                }
                tty[0] = ((2.0f64 / pi) * (e0 + (((0.125f64 * x) * x) * b1)));
            }
        } else {
            {
                a0 = (2.0f64 / (pi * x)).sqrt();
                bj0 = 0.0f64;
                by0 = 0.0f64;
                bj1 = 0.0f64;
                l = 0;
                while (l <= 1) {
                    {
                        vt = ((4.0f64 * ((l) as f64)) * ((l) as f64));
                        px = 1.0f64;
                        r = 1.0f64;
                        k = 1;
                        while (k <= 14) {
                            {
                                r = ((((((-0.0078125f64) * r)
                                    * (vt
                                        - ((4.0f64 * ((k) as f64)) - 3.0f64)
                                            .powf(((2) as f64))))
                                    / (x * ((k) as f64)))
                                    * (vt
                                        - ((4.0f64 * ((k) as f64)) - 1.0f64).powf(((2) as f64))))
                                    / (((2.0f64 * ((k) as f64)) - 1.0f64) * x));
                                px += r;
                                if ((r).abs() < ((px).abs() * 1.0e-12f64)) {
                                    {
                                        break;
                                    }
                                }
                            }
                            {
                                k += 1;
                                k
                            };
                        }
                        qx = 1.0f64;
                        r = 1.0f64;
                        k = 1;
                        while (k <= 14) {
                            {
                                r = ((((((-0.0078125f64) * r)
                                    * (vt
                                        - ((4.0f64 * ((k) as f64)) - 1.0f64)
                                            .powf(((2) as f64))))
                                    / (x * ((k) as f64)))
                                    * (vt
                                        - ((4.0f64 * ((k) as f64)) + 1.0f64).powf(((2) as f64))))
                                    / (((2.0f64 * ((k) as f64)) + 1.0f64) * x));
                                qx += r;
                                if ((r).abs() < ((qx).abs() * 1.0e-12f64)) {
                                    {
                                        break;
                                    }
                                }
                            }
                            {
                                k += 1;
                                k
                            };
                        }
                        qx = (((0.125f64 * (vt - 1.0f64)) / x) * qx);
                        xk = (x - ((0.25f64 + (0.5f64 * ((l) as f64))) * pi));
                        bj1 = (a0 * ((px * (xk).cos()) - (qx * (xk).sin())));
                        by1 = (a0 * ((px * (xk).sin()) + (qx * (xk).cos())));
                        if (l == 0) {
                            {
                                bj0 = bj1;
                                by0 = by1;
                            }
                        }
                    }
                    {
                        l += 1;
                        l
                    };
                }
                t = (2.0f64 / x);
                g0 = 1.0f64;
                r0 = 1.0f64;
                k = 1;
                while (k <= 10) {
                    {
                        r0 = (((((-r0) * ((k) as f64)) * ((k) as f64)) * t) * t);
                        g0 += r0;
                    }
                    {
                        k += 1;
                        k
                    };
                }
                g1 = 1.0f64;
                r1 = 1.0f64;
                k = 1;
                while (k <= 10) {
                    {
                        r1 = (((((-r1) * ((k) as f64)) * (((k) as f64) + 1.0f64)) * t) * t);
                        g1 += r1;
                    }
                    {
                        k += 1;
                        k
                    };
                }
                ttj[0] = ((((((2.0f64 * g1) * bj0) / (x * x)) - ((g0 * bj1) / x)) + el)
                    + (x / 2.0f64).ln());
                tty[0] = ((((2.0f64 * g1) * by0) / (x * x)) - ((g0 * by1) / x));
            }
        }
    }
    return;
}

pub(crate) fn itth0(mut x: f64) -> f64 {
    let mut k: i32 = 0;
    let mut f0: f64 = 0.0;
    let mut g0: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut t: f64 = 0.0;
    let mut tth: f64 = 0.0;
    let mut tty: f64 = 0.0;
    let mut xt: f64 = 0.0;
    let mut pi: f64 = 3.141592653589793f64;
    s = 1.0f64;
    r = 1.0f64;
    if (x < 24.5f64) {
        {
            k = 1;
            while (k < 61) {
                {
                    r = (((((-r) * x) * x) * ((2.0f64 * ((k) as f64)) - 1.0f64))
                        / ((2.0f64 * ((k) as f64)) + 1.0f64).powf(((3) as f64)));
                    s += r;
                    if ((r).abs() < ((s).abs() * 1.0e-12f64)) {
                        {
                            break;
                        }
                    }
                }
                {
                    k += 1;
                    k
                };
            }
            tth = ((pi / 2.0f64) - (((2.0f64 / pi) * x) * s));
        }
    } else {
        {
            k = 1;
            while (k < 11) {
                {
                    r = (((-r) * ((2.0f64 * ((k) as f64)) - 1.0f64).powf(((3) as f64)))
                        / ((((2.0f64 * ((k) as f64)) + 1.0f64) * x) * x));
                    s += r;
                    if ((r).abs() < ((s).abs() * 1.0e-12f64)) {
                        {
                            break;
                        }
                    }
                }
                {
                    k += 1;
                    k
                };
            }
            tth = ((2.0f64 / (pi * x)) * s);
            t = (8.0f64 / x);
            xt = (x + (0.25f64 * pi));
            f0 = ((((((((((((0.18118e-2f64 * t) - 0.91909e-2f64) * t) + 0.017033f64) * t)
                - 0.9394e-3f64)
                * t)
                - 0.051445f64)
                * t)
                - 0.11e-5f64)
                * t)
                + 0.7978846f64);
            g0 = ((((((((((((-0.23731e-2f64) * t) + 0.59842e-2f64) * t) + 0.24437e-2f64) * t)
                - 0.0233178f64)
                * t)
                + 0.595e-4f64)
                * t)
                + 0.1620695f64)
                * t);
            tty = (((f0 * (xt).sin()) - (g0 * (xt).cos())) / ((x).sqrt() * x));
            tth = (tth + tty);
        }
    }
    return tth;
}

pub(crate) fn itsl0(mut x: f64) -> f64 {
    let mut k: i32 = 0;
    let mut a = [0.0; (18) as usize];
    let mut a0: f64 = 0.0;
    let mut a1: f64 = 0.0;
    let mut af: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut rd: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut s0: f64 = 0.0;
    let mut ti: f64 = 0.0;
    let mut tl0: f64 = 0.0;
    let mut pi: f64 = 3.141592653589793f64;
    let mut el: f64 = 0.57721566490153f64;
    r = 1.0f64;
    if (x <= 20.0f64) {
        {
            s = 0.5f64;
            k = 1;
            while (k < 101) {
                {
                    rd = 1.0f64;
                    if (k == 1) {
                        {
                            rd = 0.5f64;
                        }
                    }
                    r = ((((r * rd) * ((k) as f64)) / (((k) as f64) + 1.0f64))
                        * (x / ((2.0f64 * ((k) as f64)) + 1.0f64)).powf(((2) as f64)));
                    s += r;
                    if ((r / s).abs() < 1.0e-12f64) {
                        {
                            break;
                        }
                    }
                }
                {
                    k += 1;
                    k
                };
            }
            tl0 = ((((2.0f64 / pi) * x) * x) * s);
        }
    } else {
        {
            s = 1.0f64;
            k = 1;
            while (k < 11) {
                {
                    r = (((r * ((k) as f64)) / (((k) as f64) + 1.0f64))
                        * (((2.0f64 * ((k) as f64)) + 1.0f64) / x).powf(((2) as f64)));
                    s += r;
                    if ((r / s).abs() < 1.0e-12f64) {
                        {
                            break;
                        }
                    }
                }
                {
                    k += 1;
                    k
                };
            }
            s0 = (((-s) / ((pi * x) * x)) + ((2.0f64 / pi) * ((2.0f64 * x).ln() + el)));
            a0 = 1.0f64;
            a1 = (5.0f64 / 8.0f64);
            a[(0) as usize] = a1;
            k = 1;
            while (k < 11) {
                {
                    af = (((((1.5f64 * (((k) as f64) + 0.50f64))
                        * (((k) as f64) + (5.0f64 / 6.0f64)))
                        * a1)
                        - (((0.5f64 * (((k) as f64) + 0.5f64).powf(((2) as f64)))
                            * (((k) as f64) - 0.5f64))
                            * a0))
                        / (((k) as f64) + 1.0f64));
                    a[(k) as usize] = af;
                    a0 = a1;
                    a1 = af;
                }
                {
                    k += 1;
                    k
                };
            }
            ti = 1.0f64;
            r = 1.0f64;
            k = 1;
            while (k < 11) {
                {
                    r = (r / x);
                    ti += (a[(k - 1) as usize] * r);
                }
                {
                    k += 1;
                    k
                };
            }
            tl0 = (((ti / ((((2) as f64) * pi) * x).sqrt()) * (x).exp()) + s0);
        }
    }
    return tl0;
}

/// Corrected exponentially scaled K0-tail integral (from Sci backend).
pub(crate) fn it2i0k0_corrected(x: f64) -> (f64, f64) {
    let first = if x >= 700.0 {
        // The positive series is no longer appropriate near overflow.
        {
            let mut a = [0.0];
            let mut b = [0.0];
            ittika(x, &mut a, &mut b);
            a[0]
        }
    } else {
        let square = x * x;
        let mut term = 1.0;
        let mut sum = 1.0;
        for k in 2..=1000 {
            let k = k as f64;
            term *= 0.25 * square * (k - 1.0) / (k * k * k);
            sum += term;
            if term <= f64::EPSILON * sum {
                break;
            }
        }
        0.125 * square * sum
    };
    // The scaled K tail, H(x) = exp(x) * x^(3/2) * integral_x^inf K0(t)/t dt,
    // is smooth in these variables. Chebyshev interpolation avoids 8-26
    // independent 25-term K0e evaluations per call. The coefficients were
    // formed by integrating the exponentially shifted K0e kernel with a
    // 128-point Gauss-Legendre rule on [0, 48] and cosine-node interpolation.
    // Against the defining integral, the relative discrepancy is below
    // 5e-13 throughout the two intervals, including their shared endpoint;
    // the quadrature tail beyond 48 is exponentially negligible.
    const TAIL_4_TO_20: [f64; 30] = [
        1.08149453318463951e+00,
        1.06404120700074878e-01,
        -3.33484655396161719e-02,
        1.05675055860129964e-02,
        -3.38232887722809137e-03,
        1.09249365161849066e-03,
        -3.55824511700227519e-04,
        1.16777067548016575e-04,
        -3.85929280831032200e-05,
        1.28362792801558122e-05,
        -4.29467611678576357e-06,
        1.44471726971998984e-06,
        -4.88447612504809061e-07,
        1.65910964078364032e-07,
        -5.65990551389579822e-08,
        1.93860449619146857e-08,
        -6.66494487748063095e-09,
        2.29944597437044722e-09,
        -7.95925244230962404e-10,
        2.76345442825250619e-10,
        -9.62261270665010241e-11,
        3.35964071377740454e-11,
        -1.17622948418253746e-11,
        4.12563316842806687e-12,
        -1.45310800429380521e-12,
        5.11600021605810643e-13,
        -1.83782618753032985e-13,
        6.24333917897956806e-14,
        -2.29538610341251115e-14,
        5.52012139702166354e-15,
    ];
    const TAIL_20_TO_INF: [f64; 16] = [
        1.20665583860878023e+00,
        -4.53290024763691179e-02,
        1.27891347997383109e-03,
        -4.80235938443152044e-05,
        2.22885246944914961e-06,
        -1.22074026914775757e-07,
        7.64871926905996702e-09,
        -5.36277064422385763e-10,
        4.13909462260164673e-11,
        -3.47147310897355510e-12,
        3.13582493305375465e-13,
        -2.96013213940682363e-14,
        1.37390099297363122e-15,
        2.06779038336435406e-15,
        8.84708972748171618e-16,
        4.03670152859803011e-15,
    ];
    let (t, coefficients): (f64, &[f64]) = if x < 20.0 {
        ((x - 12.0) * 0.125, &TAIL_4_TO_20)
    } else {
        (40.0 / x - 1.0, &TAIL_20_TO_INF)
    };
    let two_t = 2.0 * t;
    let (mut b0, mut b1) = (0.0, 0.0);
    for &coefficient in coefficients[1..].iter().rev() {
        let next = two_t * b0 + (coefficient - b1);
        b1 = b0;
        b0 = next;
    }
    let scaled_tail = t * b0 + (coefficients[0] - b1);
    let second = scaled_tail * (-x).exp() / (x * x.sqrt());
    (first, second)
}

use std::f64::consts::PI;
const GAUSS_32: [(f64, f64); 16] = [
    (0.9972638618494816, 0.007018610009469521),
    (0.9856115115452684, 0.01627439473090571),
    (0.9647622555875064, 0.02539206530926214),
    (0.9349060759377397, 0.03427386291302141),
    (0.8963211557660521, 0.042835898022226704),
    (0.84936761373257, 0.050998059262376154),
    (0.7944837959679424, 0.05868409347853558),
    (0.7321821187402897, 0.06582222277636181),
    (0.6630442669302152, 0.07234579410884862),
    (0.5877157572407623, 0.07819389578707044),
    (0.5068999089322294, 0.08331192422694672),
    (0.42135127613063533, 0.08765209300440374),
    (0.33186860228212767, 0.0911738786957639),
    (0.23928736225213706, 0.09384439908080441),
    (0.1444719615827965, 0.09563872007927485),
    (0.04830766568773831, 0.09654008851472785),
];
// Composite 32-point Gaussian quadrature. Intervals are at most four units;
// callers bound their integration domains before entering this kernel.
fn quadrature(left: f64, right: f64, evaluate: impl Fn(f64) -> f64) -> f64 {
    let pieces = ((right - left).abs() / 4.0).ceil().max(1.0) as usize;
    let width = (right - left) / pieces as f64;
    let mut total = 0.0;
    for piece in 0..pieces {
        let midpoint = left + (piece as f64 + 0.5) * width;
        for &(z, weight) in &GAUSS_32 {
            total += weight
                * width
                * 0.5
                * (evaluate(midpoint - z * width * 0.5) + evaluate(midpoint + z * width * 0.5));
        }
    }
    total
}
/// Large-x expansion of the integrated Struve H0 function.
///
/// H0 = Y0 + K0 (Struve K0 in DLMF notation). DLMF 11.6.3 gives the
/// logarithm and inverse even powers from the integral of K0; the remaining
/// oscillatory part is -integral_x^inf Y0(t) dt. Termwise integration of the
/// Hankel expansion of Y0 (DLMF 10.17.4) gives the cosine and sine polynomials.
/// For x >= 200, this expansion agrees with independent high-precision
/// hypergeometric references to within about 2e-15 relative error.
fn itstruve0_asymptotic(x: f64) -> f64 {
    const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;
    let inverse = 1.0 / x;
    let z = inverse * inverse;
    let algebraic = z * (0.5 + z * (-2.25 + z * 37.5));

    // Construct x - pi/4 from sin(x), cos(x): subtracting pi/4 from a huge
    // floating-point x would discard the phase and corrupt the oscillation.
    let (sin_x, cos_x) = x.sin_cos();
    let phase_cos = (cos_x + sin_x) * std::f64::consts::FRAC_1_SQRT_2;
    let phase_sin = (sin_x - cos_x) * std::f64::consts::FRAC_1_SQRT_2;
    let cosine_polynomial =
        1.0 + z * (-129.0 / 128.0 + z * (301_035.0 / 32_768.0 - z * 961_319_205.0 / 4_194_304.0));
    let sine_polynomial =
        inverse * (5.0 / 8.0 + z * (-2_655.0 / 1_024.0 + z * 10_896_795.0 / 262_144.0));
    let y0_tail =
        ((2.0 / PI) / x).sqrt() * (phase_cos * cosine_polynomial + phase_sin * sine_polynomial);

    // ln(x)+ln(2) avoids overflow in ln(2*x) near f64::MAX.
    (2.0 / PI) * (x.ln() + std::f64::consts::LN_2 + EULER_GAMMA + algebraic) - y0_tail
}

/// Integrated H0 series for small x, Gaussian quadrature for moderate x,
/// and a DLMF asymptotic from x = 200 onward. Earlier use of quadrature
/// through x = 512 required 4096 Struve evaluations at its upper endpoint.
pub(crate) fn itstruve0_corrected(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    let x = x.abs();
    if x >= 200.0 {
        return itstruve0_asymptotic(x);
    }
    if x <= 10.0 {
        let square = x * x;
        let (mut term, mut sum) = (1.0, 1.0);
        for n in 1..128 {
            let k = n as f64;
            term *= -square * k / ((k + 1.0) * (2.0 * k + 1.0).powi(2));
            sum += term;
            if square < (k + 2.0) * (2.0 * k + 3.0) && term.abs() <= 1e-16 * sum.abs() {
                break;
            }
        }
        square * sum / PI
    } else if x <= 24.0 {
        // Chebyshev approximation on [10,24], with coefficients from the
        // independently evaluated 100-digit integrated defining series.
        // Avoid repeated Struve evaluations inside Gaussian quadrature for
        // this well-conditioned moderate-argument range. 30-term Clenshaw
        // evaluation agrees with the defining series to binary64 precision.
        const COEFFS: [f64; 30] = [
            2.64043445892250972e+00,
            2.49257779660430639e-01,
            8.40097285831148816e-02,
            -4.47291815478741001e-02,
            6.37105452831351754e-02,
            -5.53446457174878309e-02,
            -1.17171559111163739e-01,
            4.91628886389838052e-02,
            4.06014366605065569e-02,
            -1.38077229284536557e-02,
            -6.86794559981943921e-03,
            2.06806073252160764e-03,
            7.11813275090017006e-04,
            -1.96884902295914126e-04,
            -5.05757532668807042e-05,
            1.30976606909748167e-05,
            2.63852449695768356e-06,
            -6.46932392889342221e-07,
            -1.05890399083151552e-07,
            2.47553628306170966e-08,
            3.38147571122969529e-09,
            -7.57380213532717924e-10,
            -8.81591380397010012e-11,
            1.89829318710235619e-11,
            1.91480515312742977e-12,
            -3.97407039435626742e-13,
            -3.52198413143367208e-14,
            7.05995285185939658e-15,
            5.55988095185977544e-16,
            -1.09271641296022452e-16,
        ];
        let y = (x - 17.0) / 7.0;
        let two_y = 2.0 * y;
        let (mut b0, mut b1) = (0.0, 0.0);
        for &coefficient in COEFFS[1..].iter().rev() {
            let next = two_y * b0 + (coefficient - b1);
            b1 = b0;
            b0 = next;
        }
        y * b0 + (COEFFS[0] - b1)
    } else {
        quadrature(0.0, x, |t| super::struve_h(0.0, t).unwrap_or(f64::NAN))
    }
}
