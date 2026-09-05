/**
 * @file ephemeris_fixed.h
 * @brief Fixed-point satellite ephemeris position calculation.
 * @details This header file declares functions for calculating satellite positions using fixed-point arithmetic.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 */
#pragma once
#include "common_types.h"

#if defined(__cplusplus)
extern "C" {
#endif

/** Normalize an angle (Q28/ANG_Q format) into [-pi, pi]. */
fixed32_t norm_angle_q28(fixed32_t angle_q28);

/**
 * Compute satellite position at time t from Keplerian ephemeris parameters,
 * using pure int32/int64 fixed-point arithmetic (see fixed_math.h, cordic.h).
 *
 * @param t_q11      Time of interest (same time base as eph_fixed.toe_32_q11),
 *                    in fixed-point Q11 (TIME_Q) format.
 * @param eph_fixed  Ephemeris parameter set (see Eph_fixed in common_types.h).
 * @return           Satellite position in the Earth-fixed frame, in Q5
 *                    (DIST_Q) fixed-point meters. If eph_fixed.A_32_q5 <= 0
 *                    (degenerate/invalid ephemeris), returns {0,0,0}.
 */
Vec3 calc_sat_pos_fixed(fixed32_t t_q11, Eph_fixed eph_fixed);

#if defined(__cplusplus)
}
#endif