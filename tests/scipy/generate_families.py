# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1"]
# ///
"""Regenerate independent SciPy fixtures from committed families_inputs.csv.

Run with:
    uv run --script tests/scipy/generate_families.py
SciPy is the independent reference; Rust tests read the checked-in CSV only.
"""
import csv
from pathlib import Path

import numpy
import scipy
from scipy import special, optimize, integrate


def main():
    assert numpy.__version__ == "2.5.4", numpy.__version__
    assert scipy.__version__ == "1.18.1", scipy.__version__
    root = Path(__file__).resolve().parent
    with (root / "families_inputs.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    output = root.parent / "fixtures" / "families_expected.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["real", "imaginary"])
        for row in rows:
            name = row["function"]
            order = float(row["order"])
            x = float(row["argument"])
            if name in {"itj0y0", "it2j0y0", "iti0k0", "it2i0k0"}:
                result = tuple(getattr(special, name)(x))
                if name == "it2i0k0" and x >= 4.0:
                    # SciPy 1.18.1's specfun second component is inaccurate
                    # here. Validate the mathematically defined K0(t)/t tail
                    # independently using a positive shifted, scaled integral.
                    scaled, _ = integrate.quad(
                        lambda u: numpy.exp(-u) * special.k0e(x + u) / (x + u),
                        0.0, numpy.inf, epsabs=1e-14, epsrel=1e-13
                    )
                    # The original SciPy first component also loses accuracy
                    # at large x; verify the defining integral separately.
                    first, _ = integrate.quad(
                        lambda t: (special.i0(t) - 1.0) / t if t else 0.0,
                        0.0, x, epsrel=1e-12
                    )
                    result = (first, numpy.exp(-x) * scaled)
            elif name in {"itstruve0", "it2struve0", "itmodstruve0"}:
                result = getattr(special, name)(x)
                if name == "itstruve0":
                    # SciPy's old specfun integral loses many digits already at
                    # x=20 and is grossly wrong for x>=40; use the defining
                    # integral and independent SciPy Struve kernel.
                    result, _ = integrate.quad(
                        lambda t: special.struve(0, t), 0.0, x,
                        epsabs=1e-12, epsrel=1e-12, limit=1000
                    )
            elif name == "iv_ratio" or name == "iv_ratio_c":
                ratio = special.ive(order, x) / special.ive(order - 1.0, x)
                result = ratio if name == "iv_ratio" else 1.0 - ratio
            elif name == "iv_ratioinv":
                r = x
                func = lambda t: special.ive(order, t) / special.ive(order - 1.0, t) - r
                result = optimize.brentq(func, 1e-12, 1000.0, xtol=1e-13)
            elif name.startswith("besselpoly:"):
                lambd = float(name.split(":")[1])
                result = special.besselpoly(x, lambd, order)
            elif name.startswith("spherical_"):
                fn, _, derivative = name.partition(":")
                result = getattr(special, fn)(int(order), x, derivative=bool(derivative))
            elif name in {"hankel1", "hankel2", "hankel1e", "hankel2e"}:
                result = getattr(special, name)(order, x)
            elif ":" in name:
                fn, _, derivative = name.partition(":")
                result = getattr(special, fn)(order, x, int(derivative))
            elif name in {"ber", "bei", "ker", "kei", "berp", "beip", "kerp", "keip"}:
                result = getattr(special, name)(x)
            else:
                result = getattr(special, name)(order, x)
            if isinstance(result, tuple):
                writer.writerow([float(result[0]), float(result[1])])
            else:
                result = complex(result)
                writer.writerow([result.real, result.imag])
    print(f"wrote {len(rows)} SciPy cases to {output}")


if __name__ == "__main__":
    main()
