# Reference fixtures

The Rust tests read committed input and expected-value files; they do not require
Python. The independent references were generated with SciPy 1.18.1 and NumPy
2.5.4 (and mpmath 1.4.1 for large Struve integrals), pinned in each script.

From the crate root, regenerate the fixtures with:

```bash
uv run --script tests/scipy/generate_families.py
uv run --script tests/scipy/generate_itstruve0_large.py
uv run --script tests/scipy/generate_bessel_iv.py
uv run --script tests/scipy/generate_tails.py
```

The `families_inputs.csv` suite contains 1,192 input cases; real I/IE regressions
cover another 1,259 cases. The SciPy 1.18.1 `it2i0k0` routine is inaccurate at
larger arguments, as is `itstruve0` from about `x = 20`, so their expected
values use independent quadrature of the defining integrals. The 21 large
`itstruve0` cases around and beyond `|x| = 512` use 50-digit mpmath ₂F₃
references, cross-checked against SciPy quadrature through 5,000. Additional
fixtures cover complex branch cuts, Hankel derivatives, small nonzero tails,
and binary32 subnormal results.

After regeneration, inspect input and expected-value changes and run:

```bash
cargo fmt --all -- --check
cargo test
cargo clippy --all-targets -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps
```
