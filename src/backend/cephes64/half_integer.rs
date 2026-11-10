//! Elementary half-integer J/Y values with bounded recurrence depth.

/// Return J_v(x) and Y_v(x) for half-integer orders in the oscillatory region.
///
/// Recurrence from +/-1/2 avoids cancellation in the general J power series.
/// Forward recurrence of J is well-conditioned when x >= |v|; elsewhere the
/// general kernel computes the minimal solution using backward recurrence.
/// Limit the steps to bound work and accumulated roundoff at large orders.
pub(super) fn jy(v: f64, x: f64) -> Option<(f64, f64)> {
    if !x.is_finite() {
        return None;
    }
    let a = v.abs();
    if !(0.5..=64.5).contains(&a) || (a - 0.5).fract() != 0.0 || x < a {
        return None;
    }

    let n = (a - 0.5) as usize;
    let (sin, cos) = x.sin_cos();
    let (mut j_prev, mut j) = (cos, sin);
    let (mut y_prev, mut y) = (sin, -cos);
    for k in 0..n {
        let factor = (2 * k + 1) as f64 / x;
        let next_j = factor * j - j_prev;
        let next_y = factor * y - y_prev;
        (j_prev, j) = (j, next_j);
        (y_prev, y) = (y, next_y);
    }

    // J_-a = -sin(pi*a) Y_a, Y_-a = sin(pi*a) J_a.
    if v < 0.0 {
        let sign = if n.is_multiple_of(2) { 1.0 } else { -1.0 };
        (j, y) = (-sign * y, sign * j);
    }
    let scale = std::f64::consts::FRAC_2_PI.sqrt() / x.sqrt();
    Some((scale * j, scale * y))
}
