#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["scipy==1.18.1", "mpmath==1.4.1"]
# ///
"""Rebuild native_reference.csv from independent special-function libraries.

Most Bessel references use 85-digit mpmath evaluations, cross-checked against
SciPy AMOS. Extreme scaled arguments use SciPy AMOS directly because the
unscaled Hankel definition J +/- iY catastrophically cancels in mpmath at
hundreds of digits of imaginary argument. Signed -0 cut values use the exact
Schwarz reflection of the SciPy upper-bank value; SciPy ignores the zero sign.
Positive real K and exp(x)K are computed from the same 85-digit mpmath
reference and checked against independent SciPy kv/kve values. Any expected
output that underflows to zero or overflows in binary64 is omitted.
The Rust regression tests run without Python.
"""

import csv
import math
from pathlib import Path

import mpmath as mp
from scipy import special as sc

mp.mp.dps = 85
OUT = Path(__file__).with_name("native_reference.csv")
FUNCTIONS = (
    "jv", "yv", "iv", "kv", "jve", "yve", "ive", "kve",
    "hankel1", "hankel2", "hankel1e", "hankel2e",
)
SCALED = {"jve", "yve", "ive", "kve", "hankel1e", "hankel2e"}

# Coordinates chosen to expose distinct algorithmic regimes. Each complex
# point is used for all 12 public APIs, except enormous arguments (scaled only).
POINTS = [
    ("tiny_positive", 0.3, 1e-8, 2e-8),
    ("tiny_negative_order", -0.7, 3e-7, -5e-7),
    ("small_half", 0.5, 0.125, 0.25),
    ("small_negative_half", -0.5, 0.2, -0.15),
    ("quadrant_one", 1.3, 2.0, 3.0),
    ("quadrant_two", -2.25, -2.0, 3.0),
    ("quadrant_three", 3.25, -2.0, -3.0),
    ("quadrant_four", 4.5, 3.0, -2.0),
    ("imag_up", 2.5, 0.0, 11.0),
    ("imag_down", -1.5, 0.0, -9.0),
    ("pure_real", 6.5, 17.0, 0.0),
    ("cut_upper", 0.3, -2.0, 0.0),
    ("cut_lower", 0.3, -2.0, -0.0),
    ("cut_integer_upper", 2.0, -5.0, 0.0),
    ("cut_integer_lower", 2.0, -5.0, -0.0),
    ("cut_half_upper", -0.5, -2.0, 0.0),
    ("cut_half_lower", -0.5, -2.0, -0.0),
    ("cut_near_up", 0.5, -2.0, 1e-14),
    ("cut_near_down", 0.5, -2.0, -1e-14),
    ("cut_near_two_up", 3.00000001, -3.0, 1e-12),
    ("cut_near_two_down", 3.00000001, -3.0, -1e-12),
    ("near_zero_order_low", -2**-24, 1.25, 0.75),
    ("near_zero_order_high", 2**-24, 1.25, 0.75),
    ("near_one_order_low", 1.0 - 2**-24, -1.75, 0.4),
    ("near_one_order_high", 1.0 + 2**-24, -1.75, 0.4),
    ("near_negative_integer_low", -4.0 - 2**-24, 4.0, -1.0),
    ("near_negative_integer_high", -4.0 + 2**-24, 4.0, -1.0),
    ("near_half_low", 5.5 - 2**-24, 7.0, -3.0),
    ("near_half_high", 5.5 + 2**-24, 7.0, -3.0),
    ("near_negative_half_low", -2.5 - 2**-24, -3.0, 1.5),
    ("near_negative_half_high", -2.5 + 2**-24, -3.0, 1.5),
    ("moderate_large_order", 18.0, 8.0, 11.0),
    ("moderate_negative_order", -17.0 + 2**-24, 8.0, -11.0),
    ("turning_point_low", 30.0, 29.0, 0.15),
    ("turning_point_high", 30.0, 31.0, -0.15),
    ("turning_point_complex", 45.5, 46.0, 3.0),
    ("large_q_one", 4.5, 30.0, 15.0),
    ("large_q_two", 8.5, -30.0, 15.0),
    ("large_q_three", -0.5, -30.0, -15.0),
    ("large_q_four", 1.5, 30.0, -15.0),
    ("large_imag_up", 0.5, 2.0, 80.0),
    ("large_imag_down", -0.5, -2.0, -80.0),
    ("large_real_right", 6.0, 100.0, 0.25),
    ("large_real_left", 3.0, -100.0, -0.25),
    ("large_high_order", 75.5, 75.0, 5.0),
    ("large_negative_high_order", -35.5, 32.0, 5.0),
    ("extreme_real_right", 0.5, 750.0, 0.25),
    ("extreme_real_left", 3.5, -750.0, -0.25),
    ("extreme_imag_up", 1.5, 2.0, 750.0),
    ("extreme_imag_down", 2.5, -2.0, -750.0),
    ("extreme_q_one", 18.0, 600.0, 600.0),
    ("extreme_q_two", 18.0, -600.0, 600.0),
    ("extreme_q_three", 18.0, -600.0, -600.0),
    ("extreme_q_four", 18.0, 600.0, -600.0),
]

AIRY_X = (
    -100.0, -30.0, -20.0, -10.001, -10.0, -9.999,
    -8.3203353, -2.0, -1e-12, 0.0, 1e-12,
    2.0, 8.3203353, 8.33, 9.999, 10.0, 10.001,
    19.999, 20.0, 20.001, 30.0, 100.0,
)
AIRYE_X = (
    0.0, 1e-12, 1.0, 8.3203353, 9.999, 10.0, 10.001,
    19.999, 20.0, 20.001, 30.0, 100.0, 500.0, 1000.0,
)

# Real noninteger orders stress direct integration, near-integer reflection,
# high-order asymptotics and K_-v = K_v. Exact half-integers are deliberate
# control points, adjacent noninteger orders exercise the general path.
REAL_K_ORDERS = (
    0.01, 0.1, 0.3, 0.4999999403953552, 0.5, 0.5000000596046448,
    0.7, 0.9999999403953552, 1.0000000596046448, 1.3,
    2.25, 5.5, 10.25, 18.75, 30.25, 49.5, 75.25,
    100.25, 150.25, 250.25,
)
REAL_K_NEG_ORDERS = (-0.1, -0.5, -1.0000000596046448, -10.25, -30.25, -75.25)
REAL_K_ULTRATINY_ORDERS = (1e-8, 1e-6, 0.01, 0.1, 0.5, 1.01)
REAL_K_ULTRATINY_X = (1e-308, 1e-300)
# Use the actual binary64 values (rather than decimal approximations):
# 5e-324 rounds to 2^-1074, and 1e-320 rounds to a nearby subnormal.
REAL_K_SUBNORMAL_ORDERS = (1e-8, 0.01, 0.1, 0.5)
REAL_K_SUBNORMAL_X = (float.fromhex("0x0.0000000000001p-1022"), 1e-320, 1e-310)
REAL_K_X = (
    1e-12, 1e-9, 1e-6, 0.0001, 0.01, 0.1, 0.25,
    0.5, 0.99, 1.0, 1.01, 2.0, 3.0, 5.0, 8.0,
    10.0, 15.0, 20.0, 30.0, 50.0, 75.0, 100.0,
    200.0, 500.0, 700.0, 705.0, 710.0, 725.0,
    740.0, 750.0, 1000.0, 1e4, 1e6,
)


def real_k_inputs():
    # Base logarithmic/linear grid and order-dependent points near the
    # Debye transition x~v and the inverse-power asymptotic threshold.
    for v in REAL_K_ORDERS:
        positions = set(REAL_K_X)
        positions.update((v / 2, v, 2 * v, 0.5 * v * v, v * v, 2 * v * v))
        threshold = 100 * (abs(v) + 1)**2
        positions.update((threshold * (1 - 2**-20), threshold * (1 + 2**-20)))
        for x in sorted(x for x in positions if 0 < x <= 1e6):
            yield v, x
    for v in REAL_K_NEG_ORDERS:
        for x in (1e-6, 0.25, 1.0, 10.0, 100.0, 700.0, 1e4, 1e6):
            yield v, x
    # These inputs exercise the extreme-small-x behavior of positive-real
    # K without discarding representable results near order zero. In
    # particular K_0.01(1e-308) is ~6e4, despite the quadrature mode
    # exceeding 700; direct exponentiation of intermediate terms overflows.
    for v in REAL_K_ULTRATINY_ORDERS:
        for x in REAL_K_ULTRATINY_X:
            yield v, x
    for v in REAL_K_SUBNORMAL_ORDERS:
        for x in REAL_K_SUBNORMAL_X:
            yield v, x


def generate_real_k():
    maximum = mp.mpf(float.fromhex("0x1.fffffffffffffp+1023"))
    minimum = mp.mpf(2)**-1074  # minimum positive representable binary64
    worst = (0.0, "")
    checked = 0
    skipped = 0
    scipy_limited = 0
    rows = []
    for v, x in real_k_inputs():
        vv = mp.mpf(repr(v))
        # Float conversion preserves the *exact* binary64 argument on the
        # subnormal and ultratiny grids. Parsing decimal repr(5e-324)
        # instead would use 5.0e-324 rather than 2^-1074 and bias
        # K_0.5 by around 0.6%, well beyond the fixture tolerance.
        xx = mp.mpf(x) if x in REAL_K_SUBNORMAL_X + REAL_K_ULTRATINY_X else mp.mpf(repr(x))
        unscaled = mp.besselk(vv, xx)
        scaled = unscaled * mp.exp(xx)
        # Keep every *nonzero, finite, representable* expected output.
        # Filtering mpmath results before conversion avoids accidental
        # reference zeros/inf in the CSV, particularly beyond x~700.
        for api, value in (("kv", unscaled), ("kve", scaled)):
            if not (minimum <= value <= maximum):
                skipped += 1
                continue
            reference = float(value)
            scipy = float(getattr(sc, api)(v, x))
            if math.isfinite(scipy) and scipy > 0.0:
                disagreement = abs(scipy / reference - 1)
                if disagreement > worst[0]:
                    worst = (disagreement, f"{api}(v={v}, x={x})")
                checked += 1
            else:
                # AMOS may overflow for finite results near DBL_MAX, and
                # underflow for nonzero subnormal K at large x. mpmath
                # independently certifies these representable cases.
                scipy_limited += 1
            # Require 2e-10 relative accuracy; values near the subnormal
            # limit allow one absolute ulp for unavoidable rounding.
            if x in REAL_K_SUBNORMAL_X:
                tag = f"subnormal_v{v}_x{x}"
            elif x in REAL_K_ULTRATINY_X:
                tag = f"ultratiny_v{v}_x{x}"
            else:
                tag = f"v{v}_x{x}"
            entry = row("real_k", tag, api, v, complex(x, 0), complex(reference), 2e-10)
            entry[-1] = "5e-324"  # one binary64 subnormal ulp
            rows.append(entry)
    print(
        f"Real K: {len(rows)} representable references, {checked} SciPy "
        f"cross-checks, {scipy_limited} SciPy over/underflow cases, "
        f"{skipped} unrepresentable outputs omitted; "
        f"worst SciPy/mpmath discrepancy={worst}"
    )
    if worst[0] > 2e-8:
        raise ValueError("SciPy/mpmath real K discrepancy requires adjudication")
    return rows

def mp_bessel(name, v, z):
    base = name.removesuffix("e")
    if base == "jv":
        value = mp.besselj(v, z)
    elif base == "yv":
        value = mp.bessely(v, z)
    elif base == "iv":
        value = mp.besseli(v, z)
    elif base == "kv":
        value = mp.besselk(v, z)
    elif base == "hankel1":
        value = mp.hankel1(v, z)
    elif base == "hankel2":
        value = mp.hankel2(v, z)
    else:
        raise ValueError(name)
    if name in ("jve", "yve"):
        value *= mp.exp(-abs(mp.im(z)))
    elif name == "ive":
        value *= mp.exp(-abs(mp.re(z)))
    elif name == "kve":
        value *= mp.exp(z)
    elif name == "hankel1e":
        value *= mp.exp(-1j * z)
    elif name == "hankel2e":
        value *= mp.exp(1j * z)
    return complex(value)


def scipy_bessel(name, v, z):
    if z.real < 0.0 and z.imag == 0.0 and math.copysign(1, z.imag) < 0:
        # Schwarz reflection with exchange of Hankel kinds. AMOS itself
        # interprets either signed zero on the cut as +i0.
        partner = name.replace("hankel1", "hankelX").replace(
            "hankel2", "hankel1"
        ).replace("hankelX", "hankel2")
        return complex(getattr(sc, partner)(v, complex(z.real, 0.0))).conjugate()
    return complex(getattr(sc, name)(v, z))


def row(family, tag, api, v, z, value, tol):
    assert math.isfinite(value.real) and math.isfinite(value.imag), (
        family, tag, api, v, z, value
    )
    return [
        family, tag, api, repr(v), repr(z.real), repr(z.imag),
        format(value.real, ".17g"), format(value.imag, ".17g"),
        format(tol, ".3g"), "1e-300",
    ]


def main():
    data = []
    agreement_worst = (0.0, "")
    for tag, v, re, im in POINTS:
        z = complex(re, im)
        for name in FUNCTIONS:
            if tag.startswith("extreme_") and name not in SCALED:
                continue
            expect_scipy = scipy_bessel(name, v, z)
            if tag.startswith("extreme_"):
                value = expect_scipy
                tol = 3e-9
            elif im == 0.0 and math.copysign(1, im) < 0 and re < 0.0:
                # SciPy upper bank above is the independent reference.
                value = expect_scipy
                tol = 7e-10
            else:
                value = mp_bessel(name, mp.mpf(repr(v)), mp.mpc(re, im))
                tol = 6e-10 if abs(z) < 40 else 3e-9
                if "near_" in tag or "turning" in tag or "cut_" in tag:
                    tol = max(tol, 2e-9)
                if abs(value) > 1e-285:
                    error = abs(expect_scipy - value) / abs(value)
                    if math.isfinite(error) and error > agreement_worst[0]:
                        agreement_worst = (error, f"{name} {tag} {v} {z}")
            data.append(row("complex", tag, name, v, z, value, tol))

    for scaled, coords in ((False, AIRY_X), (True, AIRYE_X)):
        family = "airye" if scaled else "airy"
        for x in coords:
            with mp.workdps(85):
                ai = mp.airyai(x)
                aip = mp.airyai(x, derivative=1)
                bi = mp.airybi(x)
                bip = mp.airybi(x, derivative=1)
                if scaled:
                    zeta = mp.mpf(2) / 3 * mp.power(mp.mpf(x), mp.mpf("1.5"))
                    ai *= mp.exp(zeta)
                    aip *= mp.exp(zeta)
                    bi *= mp.exp(-zeta)
                    bip *= mp.exp(-zeta)
            for name, value in zip(("ai", "aip", "bi", "bip"), (ai, aip, bi, bip)):
                data.append(row(family, "x" + str(x), name, 0.0, complex(x, 0.0), complex(value), 8e-10))

    data.extend(generate_real_k())
    with OUT.open("w", newline="") as file:
        # Keep committed numerical fixtures LF-terminated on every platform.
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow((
            "family", "tag", "api", "order", "re", "im",
            "expected_re", "expected_im", "rel_tolerance", "abs_tolerance"
        ))
        writer.writerows(data)
    print(f"Wrote {len(data)} cases to {OUT.name}; max complex SciPy/mpmath disagreement: {agreement_worst}")
    if agreement_worst[0] > 2e-8:
        raise ValueError("SciPy/mpmath disagree; manually adjudicate the reported case")


if __name__ == "__main__":
    main()
