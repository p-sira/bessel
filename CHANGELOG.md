# Changelog

## 0.1.3

### Performance

- Evaluate half-integer cylindrical J and Y in the oscillatory region using the stable spherical recurrence, reducing the cost of Bessel and Hankel derivatives while preserving the general kernels for other inputs.
- Add broad high-precision references for positive and negative half-integer orders, recurrence boundaries, and derivative orders.

## 0.1.2

### Bug Fixes

- Fix tiny-argument J underflow and parity, power-series cancellation, and fractional-order J/Y reflection overflow and precision loss ([#7](https://github.com/p-sira/bessel/pull/7)).
- Preserve finite spherical I values and derivatives when their scaled intermediates underflow ([#7](https://github.com/p-sira/bessel/pull/7)).
- Reduce premature underflow in spherical J/I derivatives at tiny arguments ([#7](https://github.com/p-sira/bessel/pull/7)).
- Fix half-order Struve H cancellation near its zeros and phase loss at very large arguments ([#7](https://github.com/p-sira/bessel/pull/7)).

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
