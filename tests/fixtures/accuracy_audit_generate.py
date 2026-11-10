# /// script
# requires-python = ">=3.11"
# dependencies = ["mpmath==1.4.1"]
# ///
"""Regenerate audit regressions using 90-digit references, verified at 130 digits.

Run from the crate root: uv run --script tests/fixtures/accuracy_audit_generate.py
Inputs are converted from their exact binary64 values before evaluation.
"""
import csv
import math
from pathlib import Path

import mpmath as mp


def reference(api, order, argument):
    v, x = mp.mpf(order), mp.mpf(argument)
    if api in ("jv", "yv"):
        return {"jv": mp.besselj, "yv": mp.bessely}[api](v, x)
    if api == "struve":
        return mp.struveh(v, x)
    modified = api.startswith("si")
    function = mp.besseli if modified else mp.besselj
    prefactor = mp.sqrt(mp.pi / (2 * x))
    value = prefactor * function(v + mp.mpf("0.5"), x)
    if api.endswith("d"):
        value = prefactor * function(v - mp.mpf("0.5"), x) - (v + 1) / x * value
    return value


cases = [("jv", -265.9768637481902, 13.301565450454433),
         ("yv", -265.49999999999, 12.5)]
for v in (0.5, -0.5, 0.01, -0.01):
    for x in (5e-324, 1e-310):
        cases.append(("jv", v, x))
for n in (3.0, -3.0, 5.0):
    cases.append(("jv", n, -1e-10))
for api in ("jv", "yv"):
    for v in (50.01, 100.25, 200.1, 279.93588689722924, 283.8627249559925,
              295.60954413659204, 302.37886735015775, 499.9, 500.1):
        for x in (v * 0.22, v * 0.48, v):
            for order in (v, -v):
                cases.append((api, order, x))
# Exact reproductions from the seeded sweep.
cases += [("jv", 283.8627249559925, 135.55386283898832),
          ("yv", 279.93588689722924, 138.9798178044934),
          ("yv", 295.60954413659204, 65.40689047242715),
          ("jv", -302.37886735015775, 75.24771249176645)]
for api in ("si", "sid"):
    for n, x in ((1, 715.0), (2, 717.0), (256, 720.0),
                 (257, 30.0), (500, 100.0), (1000, 500.0), (1000, 1000.0),
                 (68, 0.0011338333600859067), (112, 0.1260233263455215)):
        cases.append((api, float(n), x))
for api in ("sjd", "sid"):
    for n, x in ((1, 1e-200), (2, 1e-200), (103, 0.08776684995814026),
                 (219, 6.65254869143548), (68, 0.0011338333600859067)):
        cases.append((api, float(n), x))
for x in (1e-200, 1e-10, 0.5, 10.0, 710.0, 1000.0, 1e20, 1e100, 1e308):
    cases.append(("struve", 0.5, x))
for k in (1, 2, 16, 113, 1000):
    root = float(2 * k * math.pi)
    for x in (math.nextafter(root, 0.0), root, math.nextafter(root, math.inf)):
        cases.append(("struve", 0.5, x))
# Independently confirm cases for which SciPy's extreme-input values are wrong.
cases.append(("yv", -0.5, 1e-310))

path = Path(__file__).with_name("accuracy_audit.csv")
with path.open("w", newline="") as output:
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(("api", "order", "x", "expected"))
    for api, order, x in cases:
        with mp.workdps(90):
            expected = float(reference(api, order, x))
        with mp.workdps(130):
            verified = float(reference(api, order, x))
        assert expected == verified, (api, order, x, expected, verified)
        assert math.isfinite(expected), (api, order, x, expected)
        writer.writerow((api, repr(order), repr(x), repr(expected)))
print(f"Wrote {len(cases)} references to {path}")
