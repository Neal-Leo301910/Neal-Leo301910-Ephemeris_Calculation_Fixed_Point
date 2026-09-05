/**
 * @file error_eval.c
 * @brief Evaluate errors between fixed-point and double-precision satellite positions.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 */
#include <math.h>
#include "common_types.h"
#include "ephemeris_double.h"
#include "ephemeris_fixed.h"
#include "error_eval.h"
#include "fixed_math.h"


#include <math.h>
#include <stdio.h>



double evaluate_error(Vec3 sat_pos_fixed_32, Vec3d sat_pos_double, int q) {
    /**
     * Convert the fixed-point satellite position to double-precision using the specified Q format,
     * then compute the Euclidean distance to the double-precision reference position.
     */
    
    double error_x = fabs(fixed32_to_double(sat_pos_fixed_32.x, q) - sat_pos_double.x);
    double error_y = fabs(fixed32_to_double(sat_pos_fixed_32.y, q) - sat_pos_double.y);
    double error_z = fabs(fixed32_to_double(sat_pos_fixed_32.z, q) - sat_pos_double.z);
    return sqrt(error_x * error_x + error_y * error_y + error_z * error_z);
}

int show_error(Eph eph, Eph_fixed eph_fixed, int q) {
    /**
     * Sweep tk from MIN_TK to MAX_TK, compute the satellite position errors,
     * and accumulate statistics (max, mean, RMSE).
     */
    
    double max_error = 0.0;
    double mean_error = 0.0;
    double rmse = 0.0;

    long count = (long)MAX_TK - (long)MIN_TK + 1;

    for (long tk_s = MIN_TK; tk_s <= MAX_TK; tk_s++) {
        fixed32_t tk_q11 = double_to_fixed32((double)tk_s, TIME_Q);
        Vec3 sat_pos_fixed_32 = calc_sat_pos_fixed(tk_q11, eph_fixed);
        Vec3d sat_pos_double = calc_sat_pos_double((double)tk_s, eph);
        double error = evaluate_error(sat_pos_fixed_32, sat_pos_double, q);
        if (error > max_error) {
            max_error = error;
        }
        mean_error += error;
        rmse += error * error;
    }
    mean_error /= (double)count;
    rmse = sqrt(rmse / (double)count);

    printf("Max Error: %f m\n", max_error);
    printf("Mean Error: %f m\n", mean_error);
    printf("RMSE: %f m\n", rmse);

    return 0;
}