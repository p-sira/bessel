"""Generate real modified-Bessel K references with independent mpmath 1.3.0.

Run with Python containing mpmath, e.g.:
PYTHONPATH=... python3 tests/scipy/generate_real_k_small_order.py
"""
from pathlib import Path
import mpmath as mp
mp.mp.dps = 90
orders = ['0.0002', '0.02', '0.25', '0.45', '0.55', '0.8', '0.999', '1.0002', '1.75', '4.25']
arguments = ['0.01', '0.5', '1.99', '2.01', '5.0', '20.0']
dest = Path(__file__).resolve().parents[1] / 'fixtures' / 'real_k_small_order.csv'
with dest.open('w') as f:
    f.write('order,x,k_scaled\n')
    for v in orders:
        for x in arguments:
            val = mp.exp(mp.mpf(x)) * mp.besselk(mp.mpf(v), mp.mpf(x))
            f.write(f'{v},{x},{mp.nstr(val, 30)}\n')
print(dest, 'rows:', len(orders) * len(arguments))
