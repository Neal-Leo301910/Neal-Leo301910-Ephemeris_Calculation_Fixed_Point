/**
 * @file fixed_math.h
 * @brief Provides fixed-point arithmetic functions for ephemeris calculations.
 * This file defines functions for performing fixed-point multiplication,
 * division, and square root, all using pure integer arithmetic (int64
 * intermediates, no floating point), matching the target int32 fixed-point
 * pipeline used by the ephemeris_fixed implementation.
 * @note All fixed-point values use the Q format specified in common_types.h.
 * @author Qiwei Yang
 * @date 2024-06-05
 */

#pragma once
#include "common_types.h"
#if defined(__cplusplus)
extern "C" {
#endif

/**
 * Multiply two fixed-point numbers that share the same Q format.
 * @param a,b  Operands in Q format with q fractional bits.
 * @param q    Fractional bit count shared by a, b and the result.
 * @return     a*b in the same Q format.
 */
fixed32_t qmul(fixed32_t a, fixed32_t b, int q);

/**
 * Multiply two fixed-point numbers in (possibly) different Q formats and
 * produce a result in a third Q format.
 * @param a   Operand in Q format with qa fractional bits.
 * @param b   Operand in Q format with qb fractional bits.
 * @param qa  Fractional bit count of a.
 * @param qb  Fractional bit count of b.
 * @param qc  Desired fractional bit count of the result.
 * @return    a*b expressed in Q(qc) format.
 */
fixed32_t qmul_ab(fixed32_t a, fixed32_t b, int qa, int qb, int qc);

/**
 * Divide two fixed-point numbers that share the same Q format.
 * @param a,b  Operands in Q format with q fractional bits.
 * @param q    Fractional bit count shared by a, b and the result.
 * @return     a/b in the same Q format. Returns 0 and logs a warning if b == 0.
 */
fixed32_t qdiv(fixed32_t a, fixed32_t b, int q);

/**
 * Integer square root of a fixed-point number, staying in the same Q format.
 * Uses a pure-integer digit-by-digit algorithm (no library sqrt()/floating
 * point), so it is safe to use on targets without hardware floating point.
 * @param a  Non-negative operand in Q format with q fractional bits.
 * @param q  Fractional bit count of a and of the result.
 * @return   sqrt(a) in the same Q format. Returns 0 and logs a warning if a < 0.
 */
fixed32_t fixed_sqrt(fixed32_t a, int q);

/**
 * Raise a fixed-point number to a small non-negative integer power, staying
 * in the same Q format (implemented as repeated qmul()).
 *
 * WARNING: intermediate products grow fast. For a value with magnitude M in
 * Q format q, computing a^n involves an intermediate around M^n * 2^q, which
 * can overflow the int64 multiply step even when the final Q-format result
 * would fit -- e.g. cubing a distance-scale value like a GPS semi-major axis
 * (~2.6e7 m) overflows int64 well before the cube completes. Prefer
 * reformulating the math to avoid large powers of large-magnitude values
 * (e.g. sqrt(mu/A^3) as sqrt(mu/A)/A) instead of calling fixed_pow(a, 3, q)
 * directly on such values.
 *
 * @param a  Operand in Q format with q fractional bits.
 * @param n  Non-negative integer exponent.
 * @param q  Fractional bit count of a and of the result.
 * @return   a^n in the same Q format.
 */
fixed32_t fixed_pow(fixed32_t a, int n, int q);

/**
 * Convert a fixed-point number to a double-precision floating-point number.
 * @param value  Operand in Q format with q fractional bits.
 * @param q      Fractional bit count of the operand.
 * @return       The double-precision representation of the fixed-point value.
 */
double fixed32_to_double(fixed32_t value, int q);


/**
 * Convert a double-precision floating-point number to a fixed-point number.
 * @param value  Operand as a double-precision floating-point number.
 * @param q      Desired fractional bit count of the result.
 * @return       The fixed-point representation of the double-precision value in Q format with q fractional bits.
 */
fixed32_t double_to_fixed32(double value, int q);

#if defined(__cplusplus)
}
#endif