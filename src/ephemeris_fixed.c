/**
 * @file ephemeris_fixed.c
 * @brief Fixed-point ephemeris calculations.
 * @details This source file provides functions for calculating satellite positions using fixed-point arithmetic.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 */

#include <stdio.h>
#include <stdint.h>
#include "ephemeris_fixed.h"
#include "kepler_solver_fixed.h"
#include "common_types.h"
#include "cordic.h"
#include "cordic_table.h"
#include "fixed_math.h"

#define MU_DEFAULT 3.986004418e14 /* Earth's gravitational constant, m^3/s^2 */
#define MU_TIMES_2DISTQ 12755214137600000LL /* round(MU_DEFAULT * (1<<DIST_Q)) */
#define OMEGA_E_DEFAULT_Q31 156543 /* round(7.2921151467e-5 * (1<<31)) */

/* pi in Q28, precomputed (matches cordic.c's PI_Q28) -- avoids doing float
 * multiplication on every norm_angle_q28() call like the old code did. */
#define PI_Q28 843314857
#define TWO_PI_Q28 (2 * PI_Q28)

/* Half/full GPS week, in Q11 (TIME_Q) format, for tk rollover handling.
 * BUGFIX: the old code compared tk against `(1 << TIME_Q) / 2`, which is
 * 1024 -- that's the Q11 encoding of 0.5 SECONDS, not half a week. It would
 * incorrectly "roll over" tk for basically any normal call. */
#define HALF_WEEK_Q11 619315200 /* round(302400.0 * (1<<TIME_Q)) */
#define WEEK_Q11 1238630400     /* round(604800.0 * (1<<TIME_Q)) */

#ifndef EPH_DEBUG_PRINT
#define EPH_DEBUG_PRINT 0   /* set to 1 while debugging to re-enable printf traces */
#endif

/**
 * @brief Normalize an angle in Q28 format to the range [-pi, pi].
 * @param angle_q28 Angle in Q28 format.
 * @return Normalized angle in Q28 format.
 */
fixed32_t norm_angle_q28(fixed32_t angle_q28) {
    angle_q28 = angle_q28 % TWO_PI_Q28;
    if (angle_q28 > PI_Q28) {
        angle_q28 -= TWO_PI_Q28;
    } else if (angle_q28 < -PI_Q28) {
        angle_q28 += TWO_PI_Q28;
    }
    return angle_q28;
}


Vec3 calc_sat_pos_fixed(fixed32_t t_q11, Eph_fixed eph_fixed) {
    Vec3 sat_pos = {0, 0, 0};

    /* Use default Earth rotation rate if not provided in ephemeris */
    fixed32_t omega_e_q31 = (eph_fixed.Omega_e_32_q31 != 0)
                                 ? eph_fixed.Omega_e_32_q31
                                 : OMEGA_E_DEFAULT_Q31;
    
    /* Time from ephemeris reference epoch, with week rollover handling */
    fixed32_t tk_q11 = t_q11 - eph_fixed.toe_32_q11;
    if (tk_q11 > HALF_WEEK_Q11) tk_q11 -= WEEK_Q11;
    if (tk_q11 < -HALF_WEEK_Q11) tk_q11 += WEEK_Q11;

    if (eph_fixed.A_32_q5 <= 0) {
        fprintf(stderr, "calc_sat_pos_fixed: invalid ephemeris, A_32_q5=%d <= 0\n",
                eph_fixed.A_32_q5);
        return sat_pos;
    }

    /*
     * Precision upgrade: MOA_Q used to be a hardcoded 4, leaving a lot of
     * headroom unused for typical MEO/GEO orbits (A ~ 2.6e7 m), which was
     * the single biggest contributor to the overall fixed-vs-double error
     * budget (~475m RMSE end to end). mu_over_A_raw (~1.5e7 for a GPS-like
     * orbit) has plenty of room to be shifted further left before it would
     * overflow int32 -- but exactly how much room depends on A (a smaller,
     * low-orbit A gives a larger mu/A and less headroom), so MOA_Q is now
     * picked at runtime as the largest safe value instead of a fixed
     * constant, capped at 8 (beyond which the gain in end-to-end accuracy
     * flattens out, per empirical testing).
     */
    fixed64_t mu_over_A_raw = MU_TIMES_2DISTQ / (fixed64_t)eph_fixed.A_32_q5; /* mu/A, Q0 */
    int moa_q = 0;
    while (moa_q < 8 && (mu_over_A_raw << (moa_q + 1)) <= INT32_MAX) {
        moa_q++;
    }
    fixed32_t mu_over_A_qN = (fixed32_t)(mu_over_A_raw << moa_q);
    fixed32_t sqrt_moa_qN = fixed_sqrt(mu_over_A_qN, moa_q);
    fixed64_t n0_raw = ((fixed64_t)sqrt_moa_qN << (DIST_Q - moa_q + ANG_Q)) / eph_fixed.A_32_q5;
    fixed32_t n0_q28 = (fixed32_t)n0_raw;

    /* Corrected mean motion calculation */
    fixed32_t n_q28 = n0_q28 + eph_fixed.delta_n_32_q28;
#if EPH_DEBUG_PRINT
    printf("Corrected mean motion (Q28 raw): %d\n", n_q28);
#endif

    /* Mean anomaly calculation */
    fixed32_t Mk_q28 = eph_fixed.M0_32_q28 + qmul_ab(n_q28, tk_q11, ANG_Q, TIME_Q, ANG_Q);
    Mk_q28 = norm_angle_q28(Mk_q28);
#if EPH_DEBUG_PRINT
    printf("Mean anomaly (Q28 raw): %d\n", Mk_q28);
#endif

    /* Eccentric anomaly calculation */
    fixed32_t Ek_q30 = kepler_solve_fixed(Mk_q28, eph_fixed.e_32_q30);
    fixed32_t Ek_q28 = norm_angle_q28(Ek_q30 >> (TRIG_Q - ANG_Q));
#if EPH_DEBUG_PRINT
    printf("Eccentric anomaly (Q28 raw): %d\n", Ek_q28);
#endif

    /* True anomaly calculation */
    fixed32_t sin_Ek_q30, cos_Ek_q30;
    cordic_sincos_q28(Ek_q28, &sin_Ek_q30, &cos_Ek_q30);
    fixed32_t sin_Ek_q28 = sin_Ek_q30 >> (TRIG_Q - ANG_Q);
    fixed32_t cos_Ek_q28 = cos_Ek_q30 >> (TRIG_Q - ANG_Q);

    fixed32_t e_32_shifted_q28 = eph_fixed.e_32_q30 >> (TRIG_Q - ANG_Q);

    fixed32_t sin_v_term_q28 = qmul(
        fixed_sqrt((1 << ANG_Q) - qmul(e_32_shifted_q28, e_32_shifted_q28, ANG_Q), ANG_Q),
        sin_Ek_q28, ANG_Q); /* sqrt(1-e^2) * sinE */
    fixed32_t cos_v_term_q28 = cos_Ek_q28 - e_32_shifted_q28; /* cosE - e */

    fixed32_t v_q28 = cordic_atan_q28(cos_v_term_q28, sin_v_term_q28);
    v_q28 = norm_angle_q28(v_q28);
#if EPH_DEBUG_PRINT
    printf("True anomaly (Q28 raw): %d\n", v_q28);
#endif

    /* Argument of latitude calculation */
    fixed32_t u_q28 = norm_angle_q28(v_q28 + eph_fixed.omega_32_q28);
#if EPH_DEBUG_PRINT
    printf("Argument of latitude (Q28 raw): %d\n", u_q28);
#endif

    /* Radius calculation */
    fixed32_t radius_factor_q28 = (1 << ANG_Q) - qmul(e_32_shifted_q28, cos_Ek_q28, ANG_Q);
    fixed32_t r_q5 = qmul_ab(eph_fixed.A_32_q5, radius_factor_q28, DIST_Q, ANG_Q, DIST_Q);
#if EPH_DEBUG_PRINT
    printf("Radius (Q5 raw): %d\n", r_q5);
#endif

    /* Argument of latitude, radius, and inclination calculations follow */
    fixed32_t sin_u_q30, cos_u_q30;
    cordic_sincos_q28(u_q28, &sin_u_q30, &cos_u_q30);
    fixed32_t sin_u_q28 = sin_u_q30 >> (TRIG_Q - ANG_Q);
    fixed32_t cos_u_q28 = cos_u_q30 >> (TRIG_Q - ANG_Q);

    /* x'/y' also kept in Q5 (same reasoning as r above). */
    fixed32_t x_prime_q5 = qmul_ab(r_q5, cos_u_q28, DIST_Q, ANG_Q, DIST_Q);
    fixed32_t y_prime_q5 = qmul_ab(r_q5, sin_u_q28, DIST_Q, ANG_Q, DIST_Q);
#if EPH_DEBUG_PRINT
    printf("x' (Q5 raw): %d, y' (Q5 raw): %d\n", x_prime_q5, y_prime_q5);
#endif

    /* Longitude of ascending node */
    fixed32_t Omega_q28 = eph_fixed.Omega_0_32_q28
                           - qmul_ab(omega_e_q31, tk_q11, RATE_Q, TIME_Q, ANG_Q);
    Omega_q28 = norm_angle_q28(Omega_q28);
#if EPH_DEBUG_PRINT
    printf("Longitude of ascending node (Q28 raw): %d\n", Omega_q28);
#endif

    /* Precompute sine and cosine of the longitude of ascending node for rotation */
    fixed32_t sin_Omega_q30, cos_Omega_q30;
    cordic_sincos_q28(Omega_q28, &sin_Omega_q30, &cos_Omega_q30);
    fixed32_t sin_Omega_q28 = sin_Omega_q30 >> (TRIG_Q - ANG_Q);
    fixed32_t cos_Omega_q28 = cos_Omega_q30 >> (TRIG_Q - ANG_Q);
    
    /* Inclination */
    fixed32_t i_q28 = eph_fixed.i0_32_q28;
    fixed32_t sin_i_q30, cos_i_q30;
    cordic_sincos_q28(i_q28, &sin_i_q30, &cos_i_q30);
    fixed32_t sin_i_q28 = sin_i_q30 >> (TRIG_Q - ANG_Q);
    fixed32_t cos_i_q28 = cos_i_q30 >> (TRIG_Q - ANG_Q);

    fixed32_t cosI_sinOmega_q28 = qmul(cos_i_q28, sin_Omega_q28, ANG_Q);
    fixed32_t cosI_cosOmega_q28 = qmul(cos_i_q28, cos_Omega_q28, ANG_Q);

    /* Compute satellite position in ECEF coordinates */
    sat_pos.x = qmul_ab(x_prime_q5, cos_Omega_q28, DIST_Q, ANG_Q, DIST_Q)
                - qmul_ab(y_prime_q5, cosI_sinOmega_q28, DIST_Q, ANG_Q, DIST_Q);
    sat_pos.y = qmul_ab(x_prime_q5, sin_Omega_q28, DIST_Q, ANG_Q, DIST_Q)
                + qmul_ab(y_prime_q5, cosI_cosOmega_q28, DIST_Q, ANG_Q, DIST_Q);
    sat_pos.z = qmul_ab(y_prime_q5, sin_i_q28, DIST_Q, ANG_Q, DIST_Q);

    return sat_pos;
}