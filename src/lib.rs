#![no_std]
#![doc = include_str!("../README.md")]

/// Largest fraction-bit count a logarithm can be asked for.
///
/// `log2` of a `u64` is below 64, so with 58 fraction bits it still fits a
/// `u64`, and its product with [`LN_2_Q64`] still fits a `u128`.
pub const MAX_FRACTION_BITS: u32 = 58;

/// `floor(ln(2) * 2^64)`: the one `ln 2` every natural logarithm here uses.
pub const LN_2_Q64: u64 = 0xB172_17F7_D1CF_79AB;

/// Fraction bits [`ln_ratio`] works in before it rounds outwards.
const WORK_BITS: u32 = MAX_FRACTION_BITS;

/// Strict bound on `ln(n) - ln_floor::<WORK_BITS>(n)` in units of
/// `2^-WORK_BITS`: `ln 2 * (1 + 1/4) + 1 + 1 < 2.87`, see [`ln_floor`].
const WORK_SLACK: i128 = 3;

/// Binary logarithm of `value` with `F` fraction bits, rounded towards zero.
///
/// The result `r` stands for `r / 2^F`. Zero is treated as one, so the result
/// is `0`. `F` above [`MAX_FRACTION_BITS`] does not compile.
///
/// Error bound: `0 <= log2(value) - r / 2^F < 2^-F + 2^-60`. Exact for powers
/// of two. The bits come out most significant first, so a result with fewer
/// fraction bits is always a prefix of one with more.
///
/// The mantissa is normalised into `[1, 2)` as Q2.62 and squared once per
/// result bit. Each squaring truncates below `2^-62` relative, which the
/// following steps halve again and again, so the truncations sum to less than
/// `2^-60`; the remaining term is the cut after `F` bits. Every truncation is
/// downwards, so the error is one-sided.
#[must_use]
#[allow(clippy::cast_possible_truncation)] // Every narrowed value is below 2^64.
pub const fn log2_floor<const F: u32>(value: u64) -> u64 {
    const { assert!(F <= MAX_FRACTION_BITS) };
    if value <= 1 {
        return 0;
    }
    let exponent = value.ilog2();
    let mut mantissa = if exponent == 63 {
        value >> 1
    } else {
        value << (62 - exponent)
    };
    let mut result = (exponent as u64) << F;
    let mut bit = (1_u64 << F) >> 1;
    while bit != 0 {
        let wide = (mantissa as u128) * (mantissa as u128);
        let squared = (wide >> 62) as u64;
        if squared >= 1 << 63 {
            result |= bit;
            mantissa = squared >> 1;
        } else {
            mantissa = squared;
        }
        bit >>= 1;
    }
    result
}

/// Natural logarithm of `value` with `F` fraction bits, rounded towards zero.
///
/// The result `r` stands for `r / 2^F`. Zero is treated as one, so the result
/// is `0`. `F` above [`MAX_FRACTION_BITS`] does not compile.
///
/// Error bound, one-sided:
/// `0 <= ln(value) - r / 2^F < ln 2 * (2^-F + 2^-60) + 2^-58 + 2^-F`.
/// That is below `2^-55` for `F = 56`, below `3 * 2^-58` for `F = 58`, and
/// below `2^(1 - F)` for every `F <= 56`.
///
/// Computed as `floor(log2_floor::<F>(value) * LN_2_Q64 / 2^64)`. The three
/// terms of the bound are the error of [`log2_floor`] scaled by `ln 2`, the
/// truncated [`LN_2_Q64`] times a `log2` below 64, and the final cut.
#[must_use]
#[allow(clippy::cast_possible_truncation)] // Every narrowed value is below 2^64.
pub const fn ln_floor<const F: u32>(value: u64) -> u64 {
    let wide = (log2_floor::<F>(value) as u128) * (LN_2_Q64 as u128);
    (wide >> 64) as u64
}

/// An interval guaranteed to contain an exact real value.
///
/// Both ends are in units of `2^-F`, for the `F` the interval was asked for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    /// Not above the exact value.
    pub lo: i128,
    /// Not below the exact value.
    pub hi: i128,
}

/// Interval around the natural logarithm of a ratio of products,
/// `ln(n_1 * n_2 * ... / (d_1 * d_2 * ...))`, in units of `2^-F`.
///
/// The interval always contains the exact logarithm. An empty slice is the
/// empty product, one; a zero factor is treated as one, like [`ln_floor`].
/// `F` above [`MAX_FRACTION_BITS`] does not compile.
///
/// Each factor's logarithm is taken with [`MAX_FRACTION_BITS`] fraction bits,
/// where it falls short by less than `3 * 2^-58`. Numerators can therefore
/// only pull the sum below the exact value and denominators above it; the sum
/// is widened by that much on each side and then rounded outwards to `F`
/// bits. With `k` factors in total the working value is within
/// `3 * k * 2^-58` of the exact logarithm, and for `3 * k < 2^(58 - F)` the
/// interval is one or two units of `2^-F` wide.
#[must_use]
pub const fn ln_ratio<const F: u32>(numerators: &[u64], denominators: &[u64]) -> Interval {
    const { assert!(F <= MAX_FRACTION_BITS) };
    let (above, raised) = sum_ln(numerators);
    let (below, lowered) = sum_ln(denominators);
    let value = above - below;
    let drop = WORK_BITS - F;
    Interval {
        lo: (value - lowered) >> drop,
        hi: -((-(value + raised)) >> drop),
    }
}

/// Sum of the working-precision logarithms of `factors`, and a strict bound
/// on how far that sum falls short of the exact one.
const fn sum_ln(factors: &[u64]) -> (i128, i128) {
    let mut sum = 0_i128;
    let mut slack = 0_i128;
    let mut index = 0;
    while index < factors.len() {
        sum += ln_floor::<WORK_BITS>(factors[index]) as i128;
        slack += WORK_SLACK;
        index += 1;
    }
    (sum, slack)
}

#[cfg(test)]
mod tests;
