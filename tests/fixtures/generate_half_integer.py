#!/usr/bin/env -S uv run --script
# /// script
# dependencies = ["mpmath==1.4.1"]
# ///
"""Generate independent 100-digit half-integer J/Y and Sci derivative references."""

import csv
import math
from pathlib import Path

import mpmath as mp

mp.mp.dps = 100
out = Path(__file__).with_name("half_integer_reference.csv")

with out.open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(("api", "order", "x", "expected"))
    for a in (0.5, 1.5, 3.5, 6.5, 12.5, 32.5, 64.5):
        for v in (a, -a):
            for x in (a, a + 0.125, 2 * a + 0.35):
                for name, fn in (("jv", mp.besselj), ("yv", mp.bessely)):
                    writer.writerow((name, repr(v), repr(x), repr(float(fn(mp.mpf(v), mp.mpf(x))))))

    for i in range(32):
        x = 0.5 + 19.5 * i / 31
        for name, fn in (("jvp", mp.besselj), ("yvp", mp.bessely)):
            v = mp.mpf("1.5")
            z = mp.mpf(x)
            value = (fn(v - 2, z) - 2 * fn(v, z) + fn(v + 2, z)) / 4
            writer.writerow((name, "1.5", repr(x), repr(float(value))))

    for v in (1.5, -6.5, 6.5, 16.5, 32.5):
        for x in (abs(v) + 2.75, 2.5 * abs(v) + 0.125):
            for n in (1, 3, 4):
                for name, fn in (("jvp", mp.besselj), ("yvp", mp.bessely)):
                    order = mp.mpf(v)
                    z = mp.mpf(x)
                    value = sum(
                        (-1) ** k * mp.binomial(n, k) * fn(order - n + 2 * k, z)
                        for k in range(n + 1)
                    ) / 2**n
                    writer.writerow((f"{name}{n}", repr(v), repr(x), repr(float(value))))

    # Test cancellation at the nearest binary64 arguments to zeros of J and Y.
    for a in (0.5, 3.5, 16.5, 32.5, 64.5):
        for root in (mp.besseljzero(a, 1), mp.besselyzero(a, 1)):
            for x in (float(root), math.nextafter(float(root), math.inf)):
                for v in (a, -a):
                    for name, fn in (("jv", mp.besselj), ("yv", mp.bessely)):
                        writer.writerow((name, repr(v), repr(x), repr(float(fn(mp.mpf(v), mp.mpf(x))))))

    # Verify dispatch transitions and high-frequency oscillatory evaluation.
    for a in (1.5, 2.5, 8.5, 32.5, 64.5):
        for x in (math.nextafter(a, 0.0), a, math.nextafter(a, math.inf), 1e6, 1e12):
            for v in (a, -a):
                for name, fn in (("jv", mp.besselj), ("yv", mp.bessely)):
                    writer.writerow((name, repr(v), repr(x), repr(float(fn(mp.mpf(v), mp.mpf(x))))))

    for v in (1.5, -1.5, 8.5, -8.5, 32.5, -32.5):
        a = abs(v)
        for x in (math.nextafter(a, 0.0), a, math.nextafter(a, math.inf)):
            for n in (1, 3, 4):
                for name, fn in (("jvp", mp.besselj), ("yvp", mp.bessely)):
                    value = sum(
                        (-1) ** k * mp.binomial(n, k) * fn(mp.mpf(v - n + 2 * k), mp.mpf(x))
                        for k in range(n + 1)
                    ) / 2**n
                    writer.writerow((f"{name}{n}", repr(v), repr(x), repr(float(value))))
