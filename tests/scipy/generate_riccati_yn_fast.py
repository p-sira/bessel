"""Independent 90-digit mpmath Riccati y and derivative references."""
from pathlib import Path
import mpmath as mp
mp.mp.dps = 90
out = Path(__file__).resolve().parents[1] / 'fixtures' / 'riccati_yn_fast.csv'
samples = {(n, 7.5) for n in range(33)}
samples.update((n, x) for n in [0,1,2,3,7,12,20,32] for x in [.05,.5,1,2,20])
with out.open('w') as f:
    f.write('n,x,value,derivative\n')
    for n, x in sorted(samples):
        order = n + mp.mpf('0.5')
        xv = mp.mpf(x)
        def riccati(t):
            return mp.sqrt(mp.pi*t/2) * mp.bessely(order,t)
        v = riccati(xv)
        derivative = mp.diff(riccati,xv)
        f.write(f'{n},{x!r},{mp.nstr(v,32)},{mp.nstr(derivative,32)}\n')
print(out, 'rows', len(samples))
