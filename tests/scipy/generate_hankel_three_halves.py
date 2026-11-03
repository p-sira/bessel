"""80-digit independent mpmath references for real positive Hankel 3/2."""
from pathlib import Path
import mpmath as mp
mp.mp.dps = 80
out = Path(__file__).resolve().parents[1] / 'fixtures' / 'hankel_three_halves.csv'
with out.open('w') as f:
    f.write('kind,x,re,im\n')
    for x in ['0.5','1.25','3','8','20','200']:
        z = mp.mpf(x)
        for kind in ['hankel1','hankel2','hankel1e','hankel2e']:
            h = mp.hankel1(mp.mpf('1.5'),z) if 'hankel1' in kind else mp.hankel2(mp.mpf('1.5'),z)
            if kind == 'hankel1e': h *= mp.exp(-1j*z)
            if kind == 'hankel2e': h *= mp.exp(1j*z)
            f.write(f'{kind},{x},{mp.nstr(mp.re(h),30)},{mp.nstr(mp.im(h),30)}\n')
print(out)
