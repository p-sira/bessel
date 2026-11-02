# /// script
# requires-python = ">=3.13"
# dependencies = ["mpmath==1.4.1", "numpy==2.5.4", "scipy==1.18.1"]
# ///
"""Generate independent high-precision Struve H0-integral fixtures.

Run with:
    uv run --script tests/scipy/generate_itstruve0_large.py

The integrated Struve series is exactly
    x^2/pi * 2F3(1,1; 2,3/2,3/2; -x^2/4).
Evaluate this identity in mpmath at 50-digit precision; for arguments up to
5000 also verify against independent SciPy quadrature of special.struve.
"""

import csv
from pathlib import Path

import mpmath as mp
import numpy
import scipy
from scipy import integrate, special


def main():
    assert mp.__version__ == "1.4.1", mp.__version__
    assert numpy.__version__ == "2.5.4", numpy.__version__
    assert scipy.__version__ == "1.18.1", scipy.__version__
    mp.mp.dps = 50
    root = Path(__file__).resolve().parent
    with (root / "itstruve0_large_inputs.csv").open(newline="") as handle:
        inputs = list(csv.DictReader(handle))

    output = root.parent / "fixtures" / "itstruve0_large_expected.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["integral"])
        for row in inputs:
            # Convert the exact binary64 argument that Rust will receive.
            x = float(row["argument"])
            xx = mp.mpf(x)
            result = xx**2 / mp.pi * mp.hyper(
                [1, 1], [2, mp.mpf("1.5"), mp.mpf("1.5")], -xx**2 / 4
            )
            if 0 <= x <= 5000:
                reference, _ = integrate.quad(
                    lambda t: special.struve(0, t),
                    0.0,
                    x,
                    epsabs=1e-12,
                    epsrel=1e-12,
                    limit=2000,
                )
                assert abs(reference - float(result)) <= 4e-10, (x, reference, result)
            writer.writerow([mp.nstr(result, 45)])
    print(f"wrote {len(inputs)} mpmath cases to {output}")


if __name__ == "__main__":
    main()
