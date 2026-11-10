# Reference fixtures

The Rust tests read committed input and expected-value files; they do not require
Python. The independent references were generated with SciPy 1.18.1 and NumPy
2.5.4, plus mpmath 1.4.1 where high precision is needed, pinned in each script.

From the crate root, regenerate the fixtures with:

```bash
uv run --script tests/scipy/generate_families.py
uv run --script tests/scipy/generate_itstruve0_large.py
uv run --script tests/scipy/generate_bessel_iv.py
uv run --script tests/scipy/generate_tails.py
uv run --script tests/fixtures/native_generate.py
uv run --script tests/fixtures/accuracy_audit_generate.py
```

The `families_inputs.csv` suite contains 1,192 input cases; real I/IE regressions
cover another 1,259 cases. The SciPy 1.18.1 `it2i0k0` routine is inaccurate at
larger arguments, as is `itstruve0` from about `x = 20`, so their expected
values use independent quadrature of the defining integrals. The 21 large
`itstruve0` cases around and beyond `|x| = 512` use 50-digit mpmath ₂F₃
references, cross-checked against SciPy quadrature through 5,000. Additional
fixtures cover complex branch cuts, Hankel derivatives, small nonzero tails,
and binary32 subnormal results.

`tests/fixtures/native_generate.py` generates `tests/fixtures/native_reference.csv`
for `cargo test --test native_accuracy` (2,299 cases). It checks the native
`backend::complex` cylindrical J/Y/I/K and Hankel evaluators (all twelve
scaled and unscaled public `bessel::complex` APIs), the real `airy`/`airye`
quartet implemented by `backend::airy`, and positive-real `kv`/`kve` from
`backend::real_k`. Most complex values use 85-digit mpmath and are cross-checked
against SciPy AMOS; extreme scaled values use SciPy directly. Signed negative-real
cut references account for SciPy's handling of negative imaginary zero. Real K
references include tiny and subnormal arguments, with finite representable values
retained even where SciPy overflows or underflows. The complex Airy kernel is
internal; `bessel::airy(x)` and `bessel::airye(x)` are the public real-input APIs.

After regeneration, inspect input and expected-value changes and run:

```bash
cargo fmt --all -- --check
cargo test
cargo clippy --all-targets -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps
```

The [numerical accuracy audit](../../docs/numerical-accuracy-audit.md) describes
178 additional high-precision regressions and the broader sweep methodology.
