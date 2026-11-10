# Changelog

## 0.1.2

### Performance

- Speed up `jnjnp_zeros` by avoiding unnecessary root searches ([#6](https://github.com/p-sira/bessel/pull/6)).
- Reduce redundant computation in Kelvin functions, modified Bessel ratios and large-order K, spherical K derivatives, and Bessel integrals ([#6](https://github.com/p-sira/bessel/pull/6)).
- Speed up real K evaluation, spherical J/I derivatives, and complex Bessel series and recurrences  ([#6](https://github.com/p-sira/bessel/pull/6)).

## 0.1.1

### Performance

- Speed up `itstruve0` for moderate to large arguments ([#4](https://github.com/p-sira/bessel/pull/4)).
- Speed up `kv` and `kve` for well-conditioned fractional orders below one ([#4](https://github.com/p-sira/bessel/pull/4)).
- Speed up `hankel1`, `hankel2`, `hankel1e`, and `hankel2e` at order 3/2 over a bounded range of positive real arguments ([#4](https://github.com/p-sira/bessel/pull/4)).
- Compute `riccati_yn` values and derivatives together in O(n) time ([#4](https://github.com/p-sira/bessel/pull/4)).

## 0.1.0

### Added

- Initial standalone release of Bessel, Hankel, Kelvin, Struve, spherical Bessel, Riccati–Bessel, and Airy functions extracted from Sci, including checked APIs.
