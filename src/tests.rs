#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)] // The floating-point reference model is lossy by nature.

extern crate std;

use super::*;

mod reference;

use reference::REFERENCE;

/// Fraction bits of the reference table.
const REFERENCE_BITS: u32 = 100;

/// f64 carries 53 bits, so allow its own rounding on top of a bound.
const FLOAT_SLACK: f64 = 1e-13;

/// Calls `$check::<F>` once per listed fraction-bit count `F`.
macro_rules! for_fraction_bits {
    ([$($bits:literal),+], $check:ident $arguments:tt) => {
        $($check::<$bits> $arguments;)+
    };
}

/// Deterministic xorshift64* generator.
struct Rng(u64);

impl Rng {
    const fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    const fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound
    }
}

fn real(value: u64, bits: u32) -> f64 {
    value as f64 / (1_u64 << bits) as f64
}

fn samples() -> impl Iterator<Item = u64> {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let random = core::iter::from_fn(move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        Some(state >> (state % 64))
    })
    .take(20_000);
    let edges = (0..64).flat_map(|shift| {
        let power = 1_u64 << shift;
        [power.saturating_sub(1), power, power + 1]
    });
    (1..=4_096).chain(edges).chain(random).chain([u64::MAX])
}

/// `ceil(ln 2 * amount)` for an `amount` that is a multiple of `2^40`.
const fn times_ln_2(amount: u128) -> u128 {
    (((amount >> 40) * LN_2_Q64 as u128) >> 24) + 1
}

fn check_exact_and_limits<const F: u32>() {
    assert_eq!(log2_floor::<F>(0), 0);
    assert_eq!(log2_floor::<F>(1), 0);
    assert_eq!(log2_floor::<F>(2), 1 << F);
    for shift in 0..64_u64 {
        assert_eq!(log2_floor::<F>(1 << shift), shift << F);
    }
    assert!(u128::from(log2_floor::<F>(u64::MAX)) < 64_u128 << F);
    assert!(log2_floor::<F>(u64::MAX) >= 63 << F);
    assert!(F == 0 || log2_floor::<F>(u64::MAX) > 63 << F);
    assert_eq!(ln_floor::<F>(0), 0);
    assert_eq!(ln_floor::<F>(1), 0);
    assert_eq!(
        u128::from(ln_floor::<F>(2)),
        u128::from(LN_2_Q64) >> (64 - F)
    );
}

#[test]
fn log2_is_exact_for_powers_of_two_and_total_at_the_limits() {
    for_fraction_bits!([0, 1, 24, 48, 56, 58], check_exact_and_limits());
    assert_eq!(log2_floor::<58>(u64::MAX), u64::MAX);
    assert_eq!(log2_floor::<0>(u64::MAX), 63);
    assert_eq!(ln_floor::<0>(u64::MAX), 43);
}

#[test]
fn log2_and_ln_stay_inside_their_one_sided_bounds_against_f64() {
    for value in samples().filter(|value| *value > 0) {
        let float = value as f64;
        let log2_gap = float.log2() - real(log2_floor::<56>(value), 56);
        assert!(log2_gap > -FLOAT_SLACK, "log2({value}) overshoots");
        assert!(log2_gap < 2f64.powi(-55) + FLOAT_SLACK, "log2({value})");
        let ln_gap = float.ln() - real(ln_floor::<56>(value), 56);
        assert!(ln_gap > -FLOAT_SLACK, "ln({value}) overshoots");
        assert!(ln_gap < 2f64.powi(-55) + FLOAT_SLACK, "ln({value})");
    }
}

fn check_reference<const F: u32>(value: u64, log2: u128, ln: u128) {
    let unit = 1_u128 << (REFERENCE_BITS - F);
    let tail = 1_u128 << (REFERENCE_BITS - 60);
    let computed = u128::from(log2_floor::<F>(value)) << (REFERENCE_BITS - F);
    assert!(computed <= log2, "log2({value}) overshoots at {F} bits");
    assert!(log2 - computed < unit + tail, "log2({value}) at {F} bits");
    let computed = u128::from(ln_floor::<F>(value)) << (REFERENCE_BITS - F);
    let bound = times_ln_2(unit + tail) + 4 * tail + unit;
    assert!(computed <= ln, "ln({value}) overshoots at {F} bits");
    assert!(ln - computed < bound, "ln({value}) at {F} bits");
    if F <= 56 {
        assert!(ln - computed < 2 * unit, "ln({value}) at {F} bits");
    }
    if F == WORK_BITS {
        let slack = WORK_SLACK.cast_unsigned() * unit;
        assert!(
            ln - computed < slack,
            "ln({value}) exceeds the interval slack"
        );
    }
}

#[test]
fn logarithms_stay_inside_their_documented_bounds_against_the_reference() {
    for &(value, log2, ln) in &REFERENCE {
        for_fraction_bits!(
            [0, 1, 7, 24, 48, 55, 56, 57, 58],
            check_reference(value, log2, ln)
        );
    }
}

#[test]
fn fewer_fraction_bits_are_a_prefix_of_more() {
    for value in samples() {
        let full = log2_floor::<58>(value);
        assert_eq!(log2_floor::<56>(value), full >> 2);
        assert_eq!(log2_floor::<48>(value), full >> 10);
        assert_eq!(log2_floor::<1>(value), full >> 57);
        assert_eq!(log2_floor::<0>(value), full >> 58);
        assert_eq!(full >> 58, u64::from(value.max(1).ilog2()));
    }
}

#[test]
fn logarithms_are_monotonic() {
    let mut previous = (0, 0);
    for value in 1..=100_000_u64 {
        let current = (ln_floor::<56>(value), ln_floor::<58>(value));
        assert!(current.0 >= previous.0, "ln({value}) went backwards");
        assert!(current.1 >= previous.1, "ln({value}) went backwards");
        previous = current;
    }
    let mut previous = 0;
    for value in (u64::MAX - 100_000)..=u64::MAX {
        let current = log2_floor::<58>(value);
        assert!(current >= previous, "log2({value}) went backwards");
        previous = current;
    }
}

#[test]
fn logarithms_evaluate_in_const_context_to_pinned_values() {
    const LOG2_TEN: u64 = log2_floor::<56>(10);
    const LN_TEN: u64 = ln_floor::<56>(10);
    const LN_TWO: Interval = ln_ratio::<48>(&[2], &[]);
    assert_eq!(LOG2_TEN, 0x0352_69E1_2F34_6E2B);
    assert_eq!(LN_TEN, 0x024D_7637_76AA_A2AF);
    assert_eq!(ln_floor::<56>(2), LN_2_Q64 >> 8);
    assert_eq!(ln_floor::<56>(u64::MAX), 0x2C5C_85FD_F473_DE6A);
    assert_eq!(
        LN_TWO,
        Interval {
            lo: 0xB172_17F7_D1CF,
            hi: 0xB172_17F7_D1D0
        }
    );
}

#[test]
fn ln_two_constant_matches_an_independent_series() {
    // ln 2 = sum 1 / (k 2^k), in 120 fraction bits with every term rounded down.
    let sum: u128 = (1..=120u32)
        .map(|k| (1u128 << (120 - k)) / u128::from(k))
        .sum();
    assert_eq!(u128::from(LN_2_Q64), sum >> (120 - 64));
    let interval = ln_ratio::<48>(&[2, 1], &[1, 1]);
    let exact = i128::try_from(sum >> (120 - 48)).unwrap_or(0);
    assert!(interval.lo <= exact && exact < interval.hi);
    assert!(interval.hi - interval.lo <= 2);
}

#[test]
fn ratio_interval_brackets_the_floating_point_value() {
    let mut rng = Rng(0x5EED);
    let draw = |rng: &mut Rng| match rng.below(4) {
        0 => rng.below(1000) + 1,
        1 => rng.below(u64::from(u32::MAX)) + 1,
        2 => u64::MAX - rng.below(1000),
        _ => rng.next() | 1,
    };
    for _ in 0..20_000 {
        let (a, b, c, d) = (
            draw(&mut rng),
            draw(&mut rng),
            draw(&mut rng),
            draw(&mut rng),
        );
        let interval = ln_ratio::<48>(&[a, b], &[c, d]);
        let float = (a as f64).ln() + (b as f64).ln() - (c as f64).ln() - (d as f64).ln();
        assert!(interval.lo < interval.hi && interval.hi - interval.lo <= 2);
        let lo = interval.lo as f64 / (1_u64 << 48) as f64;
        assert!((lo - float).abs() < 1e-12, "{a} {b} {c} {d}");
    }
    assert_eq!(
        ln_ratio::<48>(&[1, 1], &[1, 1]),
        ln_ratio::<48>(&[7, 9], &[9, 7])
    );
    assert!(ln_ratio::<48>(&[u64::MAX, u64::MAX], &[1, 1]).hi < 89 << 48);
}

#[test]
fn ratio_interval_encloses_the_reference_value() {
    // Each table entry is rounded down, so a sum of four is within two units.
    let shift = REFERENCE_BITS - 48;
    let wide = REFERENCE_BITS - WORK_BITS;
    for &(a, _, ln_a) in &REFERENCE {
        for &(b, _, ln_b) in &REFERENCE {
            for &(c, _, ln_c) in &REFERENCE {
                for &(d, _, ln_d) in REFERENCE.iter().step_by(3) {
                    let exact = (ln_a + ln_b).cast_signed() - (ln_c + ln_d).cast_signed();
                    let interval = ln_ratio::<48>(&[a, b], &[c, d]);
                    assert!((interval.lo << shift) < exact + 2, "{a} {b} {c} {d}");
                    assert!((interval.hi << shift) > exact - 2, "{a} {b} {c} {d}");
                    assert!(interval.lo < interval.hi && interval.hi - interval.lo <= 2);
                    let full = ln_ratio::<WORK_BITS>(&[a, b], &[c, d]);
                    assert!((full.lo << wide) < exact + 2, "{a} {b} {c} {d}");
                    assert!((full.hi << wide) > exact - 2, "{a} {b} {c} {d}");
                    assert_eq!(full.hi - full.lo, 4 * WORK_SLACK);
                    assert_eq!(interval.lo, full.lo >> (WORK_BITS - 48));
                }
            }
        }
    }
}

#[test]
fn ratio_interval_holds_for_ratios_near_one() {
    let scale = (1_u64 << 48) as f64;
    let mut rng = Rng(0xD1CE);
    let near = (1..=5_000_u64)
        .chain((0..53).map(|shift| 1 << shift))
        .chain((0..5_000).map(|_| rng.below(1 << 53) + 1));
    for n in near {
        // ln((n + 1) / n) = ln_1p(1 / n), which f64 resolves far below 2^-48.
        let float = (1.0 / n as f64).ln_1p() * scale;
        let up = ln_ratio::<48>(&[n + 1], &[n]);
        assert!(up.lo as f64 <= float + 0.01 && float - 0.01 <= up.hi as f64);
        assert!(up.lo >= 0 && up.hi - up.lo <= 2, "{n}");
        let down = ln_ratio::<48>(&[n], &[n + 1]);
        assert_eq!((down.lo, down.hi), (-up.hi, -up.lo), "{n}");
    }
    for n in [u64::MAX - 1, (1 << 63) - 1, 1 << 62, (1 << 53) + 1] {
        // The exact value is positive and far below one unit of 2^-48.
        let up = ln_ratio::<48>(&[n + 1], &[n]);
        assert!(up.lo <= 0 && 0 < up.hi && up.hi - up.lo <= 2, "{n}");
    }
}

#[test]
fn ratio_interval_handles_empty_products_zeros_and_every_width() {
    assert_eq!(ln_ratio::<48>(&[], &[]), Interval { lo: 0, hi: 0 });
    assert_eq!(ln_ratio::<48>(&[0], &[]), ln_ratio::<48>(&[1], &[]));
    assert_eq!(ln_ratio::<48>(&[], &[0]), ln_ratio::<48>(&[], &[1]));
    let one_sided = ln_ratio::<58>(&[3], &[]);
    assert_eq!(one_sided.lo, i128::from(ln_floor::<58>(3)));
    assert_eq!(one_sided.hi, one_sided.lo + WORK_SLACK);
    let coarse = ln_ratio::<0>(&[1000, 1000, 1000], &[7]);
    assert_eq!((coarse.lo, coarse.hi), (18, 19));
    let negative = ln_ratio::<0>(&[7], &[1000, 1000, 1000]);
    assert_eq!((negative.lo, negative.hi), (-19, -18));
    let many = [u64::MAX; 64];
    let top = ln_ratio::<48>(&many, &[]);
    assert!(top.lo < top.hi && top.hi - top.lo <= 2);
    assert_eq!(top.lo >> 48, 2839);
}
