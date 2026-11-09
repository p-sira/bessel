# Changelog

## 0.1

### 0.1.2

**Performance Improvements**
- Bound merged `jnjnp_zeros` searches to avoid unnecessary roots while retaining an unbounded fallback; fast-path the exact origin-only case ([#6](https://github.com/p-sira/bessel/pull/6)). 
- Reduce redundant work in Kelvin series, modified Bessel ratios/large-order K, spherical K derivatives, complex K integrals and integer-power integral recurrences ([#6](https://github.com/p-sira/bessel/pull/6)).

### 0.1.1

**Performance Improvements**
- `itstruve0`: Restore the moderate-large-argument Chebyshev/Clenshaw fast path to avoid expensive quadrature ([#4](https://github.com/p-sira/bessel/pull/4)).
- `kv`, `kve`: Extend the reduced-order K recurrence to well-conditioned orders below one, eliminating quadrature in common fractional-order cases ([#4](https://github.com/p-sira/bessel/pull/4)).
- `hankel1`, `hankel2`, `hankel1e`, `hankel2e`: Evaluate the exact order-3/2 closed form for positive real arguments within a bounded range ([#4](https://github.com/p-sira/bessel/pull/4)).
- `riccati_yn`: Evaluate the value and derivative sequences in a single O(n) upward recurrence instead of recomputing each order ([#4](https://github.com/p-sira/bessel/pull/4)).

**Tests**
- Add independent high-precision reference grids for Struve transitions, small noninteger K orders, half-integer Hankel functions, and Riccati–Neumann sequences.
- Verify behavior at algorithm boundaries and invalid inputs.

### 0.1.0

**New Features**
- Initial standalone extraction of real and complex Bessel, Hankel, Kelvin, Struve, spherical and Riccati–Bessel, and Airy functions from Sci.
- Add checked APIs and independent numerical regression fixtures.
