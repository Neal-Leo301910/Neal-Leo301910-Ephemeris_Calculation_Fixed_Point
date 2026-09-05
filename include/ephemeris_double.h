/**
 * @file ephemeris_double.h
 * @brief Double-precision satellite ephemeris position calculation (reference implementation).
 */
#pragma once
#include "common_types.h"

#if defined(__cplusplus)
extern "C" {
#endif

/** Normalize an angle (radians) into [-PI, PI]. Safe for NaN/Inf input. */
double norm_angle(double angle);

/**
 * Compute satellite position at time t from Keplerian ephemeris parameters.
 *
 * @param t    Time of interest (same time base as eph.toe), seconds.
 * @param eph  Ephemeris parameter set.
 * @return     Satellite position in the Earth-fixed frame, in meters.
 *             If eph.A <= 0 (degenerate/invalid ephemeris), returns {0,0,0}.
 */
Vec3d calc_sat_pos_double(double t, Eph eph);

#if defined(__cplusplus)
}
#endif