/**
 * @file cordic.h
 * @brief CORDIC-based fixed-point trigonometric functions.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 */
#pragma once
#include <stdint.h>
#include "common_types.h"

#if defined(__cplusplus)
extern "C" {
#endif

#define CORDIC_ITERATIONS 28
#define K_CORDIC_Q30 0x26DD3B6A // CORDIC gain in Q30 format

/**
 * Compute sin/cos of an angle using the CORDIC algorithm (rotation mode).
 * @param z          Angle in Q28 format (radians). Must be pre-normalized by
 *                    the caller into roughly [-pi/2, pi/2] -- CORDIC only
 *                    converges for |angle| <~ 99.9 degrees.
 * @param sinAngle   Output: sin(z), in Q30 format.
 * @param cosAngle   Output: cos(z), in Q30 format.
 * @return           Residual angle after CORDIC_ITERATIONS steps (should be
 *                    close to 0 for a valid input angle).
 */
fixed32_t cordic_sincos_q28(fixed32_t z, fixed32_t *sinAngle, fixed32_t *cosAngle);

/**
 * Compute atan2(y0, x0) using the CORDIC algorithm (vector mode).
 * @param x0  x component, Q30 format.
 * @param y0  y component, Q30 format.
 * @return    atan2(y0, x0) in Q28 format (radians).
 */
fixed32_t cordic_atan_q28(fixed32_t x0, fixed32_t y0);

#if defined(__cplusplus)
}
#endif