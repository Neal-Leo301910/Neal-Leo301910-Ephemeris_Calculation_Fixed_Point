/**
 * @brief Implements ephemeris calculations using double-precision floating-point arithmetic.
 * Provides functions to compute satellite positions based on ephemeris data.
 * @details This module provides functions to calculate satellite positions using ephemeris data in double-precision floating-point format.
 * @note Requires the ephemeris_double.h and kepler_solver_double.h headers.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 * @remark This module is part of a larger system for satellite navigation and positioning.
 */

#include <math.h>
#include <stdio.h>
#include "ephemeris_double.h"
#include "kepler_solver_double.h"




/* Fallback constants: used only if a per-satellite eph.mu / eph.Omega_e is
 * not provided (i.e. left as 0). Prefer eph.mu / eph.Omega_e, since the
 * struct exists precisely so these can vary per constellation (GPS/BDS/etc).
 */
#define MU_DEFAULT 3.986004418e14      /* Earth's gravitational constant, m^3/s^2 */
#define OMEGA_E_DEFAULT 7.2921151467e-5 /* Earth's rotation rate, rad/s */
#define PI 3.14159265358979323846

#ifndef EPH_DEBUG_PRINT
#define EPH_DEBUG_PRINT 0   /* set to 1 while debugging to re-enable printf traces */
#endif

/*
 * use fmod() so the loop always terminates in O(1) instead of O(angle/2PI) iterations, 
 * and explicitly guard against non-finite input so a bad upstream value produces a visible NaN/Inf result instead of a hang.
 */
double norm_angle(double angle) {
    if (!isfinite(angle)) {
        return angle; /* propagate NaN/Inf visibly instead of looping forever */
    }
    angle = fmod(angle, 2.0 * PI);
    if (angle > PI) {
        angle -= 2.0 * PI;
    } else if (angle < -PI) {
        angle += 2.0 * PI;
    }
    return angle;
}

Vec3d calc_sat_pos_double(double t, Eph eph) {
    Vec3d sat_pos = {0.0, 0.0, 0.0};

    /* Guard against degenerate ephemeris (e.g. an unfilled/zeroed struct in a
     * test). A semi-major axis of 0 makes n0 = sqrt(MU/A^3) blow up to infinity,
     * which used to hang norm_angle(). */
    if (eph.A <= 0.0) {
        fprintf(stderr, "calc_sat_pos_double: invalid ephemeris, A=%f <= 0\n", eph.A);
        return sat_pos;
    }

    double mu = (eph.mu > 0.0) ? eph.mu : MU_DEFAULT;
    double omega_e = (eph.Omega_e != 0.0) ? eph.Omega_e : OMEGA_E_DEFAULT;

    double tk = t - eph.toe; /* Time from ephemeris reference epoch */

    /* Normalize tk to within half a week (handle GPS week rollover) */
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;

    double n0 = sqrt(mu / (eph.A * eph.A * eph.A)); /* Computed mean motion */
    double n = n0 + eph.delta_n;                    /* Corrected mean motion */
#if EPH_DEBUG_PRINT
    printf("Corrected mean motion: %lf\n", n);
#endif

    double Mk = eph.M0 + n * tk; /* Mean anomaly */
    Mk = norm_angle(Mk);
#if EPH_DEBUG_PRINT
    printf("Mean anomaly: %lf\n", Mk);
#endif

    double E = solve_kepler_double(Mk, eph.e);
    E = norm_angle(E);
#if EPH_DEBUG_PRINT
    printf("Eccentric anomaly: %lf\n", E);
#endif

    double v = atan2(sqrt(1.0 - eph.e * eph.e) * sin(E), cos(E) - eph.e); /* True anomaly */
    v = norm_angle(v);
#if EPH_DEBUG_PRINT
    printf("True anomaly: %lf\n", v);
#endif

    double u = v + eph.omega; /* Argument of latitude */
    u = norm_angle(u);
#if EPH_DEBUG_PRINT
    printf("Argument of latitude: %lf\n", u);
#endif

    double r = eph.A * (1.0 - eph.e * cos(E)); /* Radius */
#if EPH_DEBUG_PRINT
    printf("Radius: %lf\n", r);
#endif


    double i = eph.i0;
#if EPH_DEBUG_PRINT
    printf("Inclination: %lf\n", i);
#endif

    double x_prime = r * cos(u); /* Orbital plane coordinates x */
    double y_prime = r * sin(u); /* Orbital plane coordinates y */
#if EPH_DEBUG_PRINT
    printf("cos(u): %lf, sin(u): %lf\n", cos(u), sin(u));
    printf("x_prime: %lf, y_prime: %lf\n", x_prime, y_prime);
#endif


    double Omega = eph.Omega_0 - omega_e * tk; /* Longitude of ascending node */
    Omega = norm_angle(Omega);
#if EPH_DEBUG_PRINT
    printf("Longitude of ascending node: %lf\n", Omega);
#endif

    double cosI_sineOmega = cos(i) * sin(Omega);
    double cosI_cosOmega = cos(i) * cos(Omega);
#if EPH_DEBUG_PRINT
    printf("cosI_sineOmega: %lf\n", cosI_sineOmega);
    printf("cosI_cosOmega: %lf\n", cosI_cosOmega);
#endif

    /* Transform from orbital plane coordinates to ECEF coordinates */
    sat_pos.x = x_prime * cos(Omega) - y_prime * cosI_sineOmega; 
    sat_pos.y = x_prime * sin(Omega) + y_prime * cosI_cosOmega;
    sat_pos.z = y_prime * sin(i);
    
    
    return sat_pos;
}