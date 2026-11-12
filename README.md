# sliderule-machines

> **AMS Napier**, NCC-4235, *Napier*-class computation tender.

Engine-room flight manual: take the logarithm of a 64-bit integer using
nothing but integers, and say exactly how far off the answer can be.

John Napier published the first table of logarithms in 1614, after twenty
years of computing it by hand, so that nobody else would have to. The tender
carries the same cargo for the fleet: one table, computed once, that every
ship reads the same way. Two ships had each started carving their own. They
agreed to about fifty-five bits and disagreed after that, which is the kind of
thing nobody notices until two machines rank the same arms differently.

Allocation-free `no_std`, no dependencies, no floating point, no `unsafe`.
Every function is a `const fn` and gives the same bits on every target.

```rust
use sliderule_machines::{Interval, ln_floor, ln_ratio, log2_floor};

// Q8.56: the integer r stands for r / 2^56. Never above the true value.
assert_eq!(log2_floor::<56>(8), 3 << 56);
assert_eq!(ln_floor::<56>(2), 0x00B1_7217_F7D1_CF79);

// ln(3 * 5 / 7) with 48 fraction bits, as an interval that contains it.
let Interval { lo, hi } = ln_ratio::<48>(&[3, 5], &[7]);
assert!(lo < hi && hi - lo <= 2);

// Usable in constants.
const LN_TEN: u64 = ln_floor::<56>(10);
assert_eq!(LN_TEN >> 56, 2);
```

## Contracts

Fraction bits are a const parameter `F`, at most `MAX_FRACTION_BITS = 58`. A
result `r` stands for `r / 2^F`. Asking for more does not compile. Zero has no
logarithm; it is treated as one and yields zero.

**`log2_floor::<F>(value)`** rounds towards zero.
`0 <= log2(value) - r / 2^F < 2^-F + 2^-60`. Exact for powers of two. A result
with fewer fraction bits is a prefix of one with more.

**`ln_floor::<F>(value)`** rounds towards zero.
`0 <= ln(value) - r / 2^F < ln 2 * (2^-F + 2^-60) + 2^-58 + 2^-F`: below
`2^-55` at `F = 56`, below `3 * 2^-58` at `F = 58`. The error is one-sided:
the result is never above the true logarithm.

**`ln_ratio::<F>(numerators, denominators)`** returns an `Interval { lo, hi }`
around `ln` of the product of the numerators over the product of the
denominators, in units of `2^-F`. The interval always contains the exact
value. Each factor is taken at 58 fraction bits, the sum is widened by
`3 * 2^-58` per factor on the side that factor can be wrong on, and the ends
are rounded outwards. Four factors at `F = 48` are within `2^-54` of the exact
logarithm before rounding, and the interval is one or two units wide.

## Arithmetic

One kernel. The mantissa is normalised into `[1, 2)` as Q2.62 and squared once
per result bit; every truncation is downwards, which is where the one-sided
bound comes from. The natural logarithm is that result times `LN_2_Q64`,
`floor(ln 2 * 2^64)`, the only `ln 2` in the crate. No series, no division,
no table.

The bounds are checked against `f64`, against 90-digit reference values in
exact integer arithmetic, and at the edges: one, two, every power of two and
its neighbours, `u64::MAX`, and ratios a hair from one.

## Scope

Integer logarithms of `u64` values. That is all of it. No `exp`, no `sqrt`,
no signed or fractional inputs, no arbitrary precision. A ratio goes in as
its integer factors.

`suspect_machines` keeps its own `ln`, `exp` and `sqrt`, and that is on
purpose. Its suspicion level is an `f64`, so its routines are deterministic
IEEE-754 sequences with relative error bounds. The ones here work on
fixed-point integers with absolute, one-sided bounds. They share a name and
nothing else; one crate holding both would be two crates in a trench coat.

## License

MIT License. See [LICENSE.txt](LICENSE.txt) for details.
