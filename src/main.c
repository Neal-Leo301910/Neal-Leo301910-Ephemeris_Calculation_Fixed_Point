/**
 * @brief Main program for testing ephemeris calculations in both double and fixed-point.
 * 
 * Compares the satellite position calculations using double precision and fixed-point arithmetic.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @details This program tests and compares satellite position calculations using both double precision and fixed-point arithmetic.
 * @note This program requires the ephemeris_double.c and ephemeris_fixed.c implementations.
 * @debug Controlled by the EPH_DEBUG_PRINT macro.
 * @version
 * @output Prints the sample ephemeris values for both double and fixed-point representations.
 */

#include "cordic.h"
#include "common_types.h"
#include "cordic_table.h"
#include "fixed_math.h"
#include "ephemeris_double.h"
#include "ephemeris_fixed.h"
#include "error_eval.h"
#include "kepler_solver_double.h"
#include "kepler_solver_fixed.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Enable debug printing for ephemeris calculations (0 = off, 1 = on)
#define EPH_DEBUG_PRINT 0

/**
 * Creates a sample ephemeris structure with predefined values.
 * Used for testing and comparison between double and fixed-point calculations.
 */
static Eph make_sample_eph(void) {
    Eph eph;
    eph.toe     = 352780.0;                 /* reference epoch, s */
    eph.A       = 26560623.69;          /* semi-major axis, m (~MEO) */
    eph.e       = 0.009451;                /* eccentricity */
    eph.M0      = 2.0809;                 /* mean anomaly at toe, rad */
    eph.delta_n = 4.3363e-09;                 /* mean motion correction, rad/s */
    eph.omega   = 0.63488;                 /* argument of perigee, rad */
    eph.Omega_0 = 1.2866;                 /* RAAN at toe, rad */
    eph.i0      = 0.97551;  /* inclination, rad (~55 deg) */
    eph.mu      = 3.986004418e14;      /* Earth gravitational constant */
    eph.Omega_e = 7.2921151467e-5;     /* Earth rotation rate, rad/s */
    return eph;
}

/**
 * Creates a sample fixed-point ephemeris structure.
 * Converts the same values used by make_sample_eph() into fixed-point representation.
 */
static Eph_fixed make_sample_eph_fixed(void) {
    Eph_fixed eph_fixed_32;
    eph_fixed_32.toe_32_q11     = double_to_fixed32(352780.0, TIME_Q);
    eph_fixed_32.A_32_q5        = double_to_fixed32(26560623.69, DIST_Q);
    eph_fixed_32.e_32_q30       = double_to_fixed32(0.009451, TRIG_Q);
    eph_fixed_32.M0_32_q28      = double_to_fixed32(2.0809, ANG_Q);
    eph_fixed_32.delta_n_32_q28 = double_to_fixed32(4.3363e-09, ANG_Q);
    eph_fixed_32.omega_32_q28   = double_to_fixed32(0.63488, ANG_Q);
    eph_fixed_32.Omega_0_32_q28 = double_to_fixed32(1.2866, ANG_Q);
    eph_fixed_32.i0_32_q28      = double_to_fixed32(0.97551, ANG_Q);
    eph_fixed_32.Omega_e_32_q31 = double_to_fixed32(7.2921151467e-5, RATE_Q);
    return eph_fixed_32;
}

int main(void) {
    /* Create sample ephemeris structures */
    Eph sample_eph = make_sample_eph();
    Eph_fixed sample_eph_fixed_32 = make_sample_eph_fixed();

    /* Print sample ephemeris in double-precision format */
    printf("Sample Eph:\n");
    printf("toe: %f\n", sample_eph.toe);
    printf("A: %f\n", sample_eph.A);
    printf("e: %f\n", sample_eph.e);
    printf("M0: %f\n", sample_eph.M0);
    printf("delta_n: %f\n", sample_eph.delta_n);
    printf("omega: %f\n", sample_eph.omega);
    printf("Omega_0: %f\n", sample_eph.Omega_0);
    printf("i0: %f\n", sample_eph.i0);
    printf("mu: %f\n", sample_eph.mu);
    printf("Omega_e: %f\n", sample_eph.Omega_e);

    printf("==================================\n");
    /* Print sample ephemeris in fixed-point format */
    printf("\nSample Eph Fixed 32:\n");
    printf("toe: %d\n", sample_eph_fixed_32.toe_32_q11);
    printf("A: %d\n", sample_eph_fixed_32.A_32_q5);
    printf("e: %d\n", sample_eph_fixed_32.e_32_q30);
    printf("M0: %d\n", sample_eph_fixed_32.M0_32_q28);
    printf("delta_n: %d\n", sample_eph_fixed_32.delta_n_32_q28);
    printf("omega: %d\n", sample_eph_fixed_32.omega_32_q28);
    printf("Omega_0: %d\n", sample_eph_fixed_32.Omega_0_32_q28);
    printf("i0: %d\n", sample_eph_fixed_32.i0_32_q28);
    /* NOTE: 
     * mu is intentionally not printed here -- Eph_fixed has no mu field (see common_types.h);
     * calc_sat_pos_fixed() uses a compile-time MU_DEFAULT constant instead. 
     */
    printf("Omega_e: %d\n", sample_eph_fixed_32.Omega_e_32_q31);

    printf("\n");
    printf("End of sample output.\n");

    printf("==================================\n");
    /* Prompt the user for the time of interest */
    printf("Please enter the time of interest: \n");
    double t;
    if (scanf("%lf", &t) != 1) {
        fprintf(stderr, "Invalid input.\n");
        return 1;
    }
    fixed32_t t_fixed_q11 = double_to_fixed32(t, TIME_Q);
    printf("Time of interest (fixed-point Q11): %d\n", t_fixed_q11);
    printf("==================================\n");
    printf("double precision calculation:\n");
    Vec3d position = calc_sat_pos_double(t, sample_eph);
    printf("Satellite position at time of interest: x=%f, y=%f, z=%f\n",
           position.x, position.y, position.z);

    printf("==================================\n");
    printf("fixed-point calculation:\n");
    Vec3 position_fixed = calc_sat_pos_fixed(t_fixed_q11, sample_eph_fixed_32);
    printf("Satellite position at time of interest (fixed-point Q%d): x=%d, y=%d, z=%d\n",
           DIST_Q, position_fixed.x, position_fixed.y, position_fixed.z);
    printf("  (converted back to meters: x=%f, y=%f, z=%f)\n",
           fixed32_to_double(position_fixed.x, DIST_Q),
           fixed32_to_double(position_fixed.y, DIST_Q),
           fixed32_to_double(position_fixed.z, DIST_Q));

    printf("\n");

    printf("==================================\n");
    printf("Quantify the error between double and fixed-point positions:\n");

    show_error(sample_eph, sample_eph_fixed_32, DIST_Q);

    printf("Program finished successfully.\n");

    return 0;
}