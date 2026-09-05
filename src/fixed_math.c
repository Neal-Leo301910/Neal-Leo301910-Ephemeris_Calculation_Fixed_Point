/**
 * @file fixed_math.c
 * @brief Implements fixed-point arithmetic functions for ephemeris calculations.
 * Provides functions for fixed-point multiplication, division, and square
 * root, all using pure integer arithmetic.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 */
#include "fixed_math.h"
#include <stdio.h>

/*
 * BUGFIX (naming): the original function was declared/defined as `quml`,
 * a typo for `qmul`. It compiled fine (header and .c agreed on the typo),
 * but the name is confusing and doesn't match the naming used everywhere
 * else in the project (design doc, ephemeris code, etc). Renamed to `qmul`.
 */
fixed32_t qmul(fixed32_t a, fixed32_t b, int q) {
    /* Perform fixed-point multiplication with specified Q format.
     * a and b are fixed-point numbers in Q format with q fractional bits.
     * The int64 intermediate is required: two int32 values near their max
     * magnitude can multiply to ~1e18, which overflows int32/int64-shift-back
     * silently if done in 32-bit arithmetic. */
    return (fixed32_t)(((fixed64_t)a * b) >> q);
}

fixed32_t qmul_ab(fixed32_t a, fixed32_t b, int qa, int qb, int qc) {
    /* Perform fixed-point multiplication with specified Q formats for a, b, and result.
     * a is in Q format with qa fractional bits.
     * b is in Q format with qb fractional bits.
     * qc is the result in Q format with qc fractional bits. */
    return (fixed32_t)(((fixed64_t)a * b) >> (qa + qb - qc));
}

fixed32_t qdiv(fixed32_t a, fixed32_t b, int q) {
    /* Perform fixed-point division with specified Q format.
     * a and b are fixed-point numbers in Q format with q fractional bits. */
    if (b == 0) {
        fprintf(stderr, "Warning: qdiv() division by zero.\n");
        return 0;
    }
    return (fixed32_t)(((fixed64_t)a << q) / b);
}

fixed32_t fixed_sqrt(fixed32_t a, int q) {
    /* Perform fixed-point square root with specified Q format.
     * a is a fixed-point number in Q format with q fractional bits.
     * The result is in the same Q format.
     */
    if (a < 0) {
        fprintf(stderr, "Warning: fixed_sqrt() called with a negative number.\n");
        return 0;
    }
    if (a == 0) {
        return 0;
    }

    /* Q-format sqrt identity: sqrt_fixed(a_Q) = sqrt(a_Q * 2^Q), so the
     * result comes back out in the same Q format instead of losing half
     * the scale. */
    fixed64_t s = ((fixed64_t)a) << q;

    /* Find the largest power of 4 that is <= s, as the starting bit for the
     * digit-by-digit algorithm (must start on an even bit position). */
    fixed64_t bit = (fixed64_t)1 << 62;
    while (bit > s) {
        bit >>= 2;
    }

    fixed64_t result = 0;
    fixed64_t rem = s;
    while (bit != 0) {
        if (rem >= result + bit) {
            rem -= result + bit;
            result = (result >> 1) + bit;
        } else {
            result >>= 1;
        }
        bit >>= 2;
    }

    return (fixed32_t)result;
}

fixed32_t fixed_pow(fixed32_t a, int n, int q) {
    /* n=0 -> 1.0 in Q format; otherwise repeated qmul(). See the overflow
     * warning in fixed_math.h -- this is only safe for small-magnitude a
     * and/or small n. */
    if (n <= 0) {
        return (fixed64_t)1 << q;
    }
    fixed64_t result = a;
    for (int i = 1; i < n; i++) {
        result = qmul(result, a, q);
    }
    return (fixed32_t)result;
}



double fixed32_to_double(fixed32_t value, int q) {
    return (double)value / (1 << q);
}

fixed32_t double_to_fixed32(double value, int q) {
    return (fixed32_t)llround(value * (double)(1LL << q));  // Convert double to fixed-point with Q format
}