/**
 * @file cordic.c
 * @brief Implements CORDIC algorithm for fixed-point arithmetic.
 * Provides functions for trigonometric calculations using fixed-point numbers.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 */

#include "cordic_table.h"
#include "cordic.h"
#include "common_types.h"
#include "fixed_math.h"
#include <stdio.h>


/* Sum of atan(2^-i) for i=0..CORDIC_ITERATIONS-1 is ~1.7433 rad (~99.88 deg,
 * the well-known CORDIC convergence limit for a single elementary rotation
 * sequence). PI_Q28/HALF_PI_Q28 below fold any input in [-pi, pi] down into
 * [-pi/2, pi/2] first (comfortably inside that 99.88 deg limit), so the
 * caller only needs to have done the *coarser* [-pi, pi] normalization
 * (e.g. via norm_angle_q28()) before calling this function.
 */
#define PI_Q28 843314857
#define HALF_PI_Q28 (PI_Q28 / 2)
#define CORDIC_MAX_ANGLE_Q28 467959937 /* ~1.7433 rad in Q28 */

fixed32_t cordic_sincos_q28(fixed32_t z, fixed32_t *sinAngle, fixed32_t *cosAngle){
    if (z > PI_Q28 || z < -PI_Q28) {
        fprintf(stderr,
                "Warning: cordic_sincos_q28() called with |angle| > 180 deg "
                "(z=%d); the caller must normalize the angle into [-pi, pi] "
                "first (e.g. via norm_angle_q28()), result will not converge "
                "correctly.\n", z);
    }

    int negate = 0;
    if (z > HALF_PI_Q28) {
        z -= PI_Q28;
        negate = 1;
    } else if (z < -HALF_PI_Q28) {
        z += PI_Q28;
        negate = 1;
    }

    fixed32_t x0 = ONE_Q30; // 1.0 in Q30 format
    fixed32_t y0 = 0;
    fixed32_t z0 = z;

    for (int i = 0; i < CORDIC_ITERATIONS; i++) {
        int sign = (z0 >= 0) ? 1 : -1;
        fixed32_t x1 = x0 - sign * (y0 >> i); // CORDIC rotation step for x
        fixed32_t y1 = y0 + sign * (x0 >> i); // CORDIC rotation step for y
        fixed32_t z1 = z0 - sign * CORDIC_TABLE_Q28[i]; // Update the remaining angle

        x0 = x1;
        y0 = y1;
        z0 = z1;
    }

    fixed32_t cosResult = qmul(x0, K_CORDIC_Q30, 30);
    fixed32_t sinResult = qmul(y0, K_CORDIC_Q30, 30);
    if (negate) {
        cosResult = -cosResult;
        sinResult = -sinResult;
    }
    *cosAngle = cosResult;
    *sinAngle = sinResult;

    return z0;
}


fixed32_t cordic_atan_q28(fixed32_t x0, fixed32_t y0){
/* Compute atan2(y0, x0) using the CORDIC algorithm (vector mode). */

    int reflected = 0;
    if (x0 < 0) {
        x0 = -x0;
        y0 = -y0;
        reflected = 1;
    }

    fixed32_t z0 = 0;
    for (int i = 0; i < CORDIC_ITERATIONS; i++) {
        fixed32_t sign = (y0 >= 0) ? 1 : -1;
        fixed32_t x1 = x0 + sign * (y0 >> i);
        fixed32_t y1 = y0 - sign * (x0 >> i);
        fixed32_t z1 = z0 + sign * CORDIC_TABLE_Q28[i];

        x0 = x1;
        y0 = y1;
        z0 = z1;
    }

    if (reflected) {
        /* y0's ORIGINAL sign (before negation above) decided which way to
         * fold; we already overwrote y0, so use z0's sign as a stand-in --
         * z0 here is atan2(-y0_orig, -x0_orig), which has the same sign as
         * -y0_orig, i.e. the opposite sign of y0_orig. */
        z0 = (z0 <= 0) ? (z0 + PI_Q28) : (z0 - PI_Q28);
    }

    return z0;
}