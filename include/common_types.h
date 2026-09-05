/**
 * @file common_types.h
 * @brief Defines common types and fixed-point Q format constants for ephemeris calculations.
 * Defines common fixed-point types and Q format constants.
 * Provides typedefs for fixed-point and double-precision ephemeris structures, as well as a 3D vector structure.
 * Provides a consistent way to represent and manipulate ephemeris data in both fixed-point and double-precision formats.
 * @note This file is intended to be included in both C and C++ projects.
 * @note All fixed-point values use the Q format specified by the corresponding constants.
 * @author Qiwei Yang
 * @date 2024-06-05
 */

#pragma once
#include <stdint.h>
#include <math.h>

#if defined(__cplusplus)
extern "C" {    
#endif

/**
 * @brief Q format constants for fixed-point representations.
 */
#define ANG_Q 28        // Q format for angles: M, E, v, u, omega, i0, Omega_e, Omega_0. 
#define TRIG_Q 30        // Q format for trigonometric values: sin, cos, tan
#define DIST_Q 5        // Q format for distances: A, r, x, y, z
#define TIME_Q 11        // Q format for time: tk
#define RATE_Q 31        // Q format for rates: vx, vy, vz

#define ONE_Q30 (1 << 30)

/**
 * @brief Fixed-point type definitions.
 */
typedef int32_t fixed32_t;
typedef int64_t fixed64_t;

/**
 * @brief Ephemeris structure in fixed-point representation.
 */
typedef struct Eph_fixed {
    fixed32_t toe_32_q11;    // Time of ephemeris, Q11 (TIME_Q)
    fixed32_t A_32_q5;       // Semi-major axis, Q5 (DIST_Q)
    fixed32_t e_32_q30;      // Eccentricity, Q30 (TRIG_Q)
    fixed32_t M0_32_q28;     // Mean anomaly at reference time, Q28 (ANG_Q)
    fixed32_t delta_n_32_q28; // Mean motion difference, Q28 (ANG_Q)
    fixed32_t omega_32_q28;  // Argument of perigee, Q28 (ANG_Q)
    fixed32_t Omega_0_32_q28; // Longitude of ascending node at reference time, Q28 (ANG_Q)
    fixed32_t i0_32_q28;     // Inclination angle at reference time, Q28 (ANG_Q)
    fixed32_t mu_32_q28;       // Gravitational parameter, Q5 (DIST_Q)
    fixed32_t Omega_e_32_q31; // Earth rotation rate, Q31 (RATE_Q). 0 => use OMEGA_E_DEFAULT.
} Eph_fixed;


/**
 * @brief Ephemeris structure in double-precision representation.
 */
typedef struct Eph {
    double toe; // Time of ephemeris
    double A;    // Semi-major axis
    double e;    // Eccentricity
    double M0;   // Mean anomaly at reference time
    double delta_n; // Mean motion difference
    double omega; // Argument of perigee
    double Omega_0; // Longitude of ascending node at reference time
    double i0;   // Inclination angle at reference time
    double mu;   // Gravitational parameter
    double Omega_e; // Earth rotation rate
} Eph;

/**
 * @brief 3D vector structure in fixed-point representation.
 */
typedef struct Vec3 {
    fixed32_t x;
    fixed32_t y;
    fixed32_t z;
} Vec3;

/**
 * @brief 3D vector structure in double-precision representation.
 */
typedef struct Vec3d {
    double x;
    double y;
    double z;
} Vec3d;

#if defined(__cplusplus)
}
#endif