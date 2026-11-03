use bessel::{jn_zeros, jnjnp_zeros, jnp_zeros};

#[test]
fn merged_roots_agree_with_independently_computed_complete_sequences() {
    for nt in [1, 2, 10, 32, 64, 128] {
        let cutoff = ((nt as f64) * 8.0).sqrt() + 4.0;
        let mut expected: Vec<(f64, f64, f64, f64)> = vec![(0.0, 0.0, 0.0, 1.0)];
        for n in 0..=cutoff.ceil() as u32 {
            for (flag, roots) in [
                (0.0, jn_zeros(n, nt).unwrap()),
                (1.0, jnp_zeros(n, nt).unwrap()),
            ] {
                for (index, x) in roots.into_iter().enumerate() {
                    expected.push((x, f64::from(n), (index + 1) as f64, flag));
                }
            }
        }
        expected.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| b.3.total_cmp(&a.3)));
        let mut start = 0;
        while start < expected.len() {
            let mut end = start + 1;
            while end < expected.len() && (expected[end].0 - expected[start].0).abs() <= 1e-12 {
                end += 1;
            }
            expected[start..end].sort_by(|a, b| a.3.total_cmp(&b.3));
            start = end;
        }
        let actual = jnjnp_zeros(nt).unwrap();
        for (index, &(root, order, ordinal, flag)) in expected.iter().take(nt).enumerate() {
            assert!(
                (actual[0][index] - root).abs() < 2e-12,
                "nt={nt}, index={index}"
            );
            assert_eq!(
                (actual[1][index], actual[2][index], actual[3][index]),
                (order, ordinal, flag),
                "nt={nt}, index={index}"
            );
        }
    }
}
