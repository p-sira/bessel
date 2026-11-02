# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1"]
# ///
"""Independent numerical references for small, complex and derivative cases.

Run `uv run --offline --script tests/scipy/generate_tails.py` from the bessel
crate. The Python script imports only NumPy and SciPy. Inputs and expected
values are separate, with case IDs to prevent accidental row reordering.
All references are evaluated with SciPy 1.18.1; binary32 subnormal expected
values are explicitly rounded to np.float32 after quantizing each input.
"""
import csv
from pathlib import Path
import numpy as np
import scipy
from scipy import special

assert np.__version__ == "2.5.4", np.__version__
assert scipy.__version__ == "1.18.1", scipy.__version__
DEST = Path(__file__).resolve().parents[1] / "fixtures"


def write(name, headings, expected_columns, rows):
    with (DEST / f"{name}_inputs.csv").open("w", newline="") as input_file, (
        DEST / f"{name}_expected.csv"
    ).open("w", newline="") as expected_file:
        inputs = csv.writer(input_file, lineterminator="\n")
        outputs = csv.writer(expected_file, lineterminator="\n")
        inputs.writerow(["id", *headings])
        outputs.writerow(["id", *expected_columns])
        for i, (args, expected) in enumerate(rows):
            inputs.writerow([i, *args])
            outputs.writerow([i, *expected])
    print(f"{name}: {i + 1} cases")


# These are near the underflow/cancellation regions that a 2e-13 absolute
# tolerance permits to return 0. All expected values are finite and nonzero.
tiny_cases = [
    ("jv", 8.0, 0.01),
    ("jv", 10.0, 0.1),
    ("jv", 12.0, 0.2),
    ("spherical_jn", 10.0, 0.1),
    ("spherical_jn", 12.0, 0.3),
    ("spherical_jn", 15.0, 0.5),
    ("kv", 0.0, 30.0),
    ("kv", 0.0, 50.0),
    ("kv", 1.5, 45.0),
    ("k0", 0.0, 30.0),
    ("ker", 0.0, 50.0),
    ("kei", 0.0, 50.0),
    ("iv", 8.0, 0.01),
    ("ive", 8.0, 0.01),
]

def tiny_eval(name, order, x):
    if name == "spherical_jn":
        return special.spherical_jn(int(order), x)
    if name in {"k0", "ker", "kei"}:
        return getattr(special, name)(x)
    return getattr(special, name)(order, x)

tiny = []
for name, v, x in tiny_cases:
    want = float(tiny_eval(name, v, x))
    assert want != 0 and np.isfinite(want), (name, v, x, want)
    tiny.append(((name, v, x), (want,)))
write("tails", ["api", "order", "x"], ["value"], tiny)

# Only nonzero binary32 subnormals are selected. Exact input conversion to f32
# and independent SciPy evaluation at the quantized input agree with the
# expected public (binary64-core, then f32 rounding) interface.
f32_pairs = [
    (26.0, 0.5), (29.0, 1.0), (30.0, 1.0),
    (34.0, 2.0), (35.0, 2.0), (37.0, 2.0),
    (39.0, 3.0), (40.0, 3.0), (44.0, 4.0),
    (45.0, 4.0), (46.0, 4.0), (47.0, 6.0),
    (48.0, 6.0), (49.0, 5.0), (49.75, 5.0),
    (50.0, 5.0), (-50.0, 5.0), (50.0, 6.0),
    (51.0, 7.0), (52.0, 6.0), (53.0, 8.0),
]

subnormal = []
for v, x in f32_pairs:
    v, x = float(np.float32(v)), float(np.float32(x))
    for name in ("iv", "ive"):
        want = np.float32(getattr(special, name)(v, x))
        if not (0 < abs(float(want)) < np.finfo(np.float32).tiny):
            continue
        subnormal.append(((name, v, x), (repr(float(want)),)))
assert len(subnormal) >= 8, subnormal
write("subnormal", ["api", "order_f32", "x_f32"], ["value_f32"], subnormal)

# Complex principal branches sampled on both sides of the negative-real cut,
# on both sides of the imaginary axis, and at small / moderate / large radii.
# For each point, evaluate all twelve complex cylindrical Bessel/Hankel APIs.
complex_points = [
    (0.5, 0.4, 0.8),
    (1.3, 1.5, -0.7),
    (2.5, -1.4, 0.6),
    (2.5, -1.4, -0.6),
    (0.5, -2.0, 1e-8),
    (0.5, -2.0, -1e-8),
    (1.3, -0.5, 1.1),
    (1.3, -0.5, -1.1),
    (0.0, 0.1, 0.2),
    (3.25, 3.2, 4.1),
    (3.25, 3.2, -4.1),
    (0.5, 9.0, 7.0),
]
complex_apis = (
    "jv", "yv", "iv", "kv", "jve", "yve", "ive", "kve",
    "hankel1", "hankel2", "hankel1e", "hankel2e",
)
complex_rows = []
for v, re, im in complex_points:
    for name in complex_apis:
        w = complex(getattr(special, name)(v, complex(re, im)))
        assert np.isfinite(w.real) and np.isfinite(w.imag), (name, v, re, im, w)
        complex_rows.append(((name, v, re, im), (repr(w.real), repr(w.imag))))
write("complex", ["api", "order", "re", "im"], ["expected_re", "expected_im"], complex_rows)

# Independent SciPy derivatives with respect to *real* x, not differences of
# the Rust J/Y functions. n=0 is a normalization check, 1-3 exercise the
# public h1vp/h2vp implementation.
derivative_rows = []
for v in (0.0, 0.5, 2.25):
    for x in (0.25, 1.5, 5.0, 12.0):
        for n in (0, 1, 2, 3):
            for name in ("h1vp", "h2vp"):
                w = complex(getattr(special, name)(v, x, n))
                assert np.isfinite(w.real) and np.isfinite(w.imag)
                derivative_rows.append(((name, v, x, n), (repr(w.real), repr(w.imag))))
write("hankel_derivatives", ["api", "order", "x", "derivative"], ["expected_re", "expected_im"], derivative_rows)
