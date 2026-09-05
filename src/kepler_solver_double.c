/**
 * @file kepler_solver_double.c
 * @brief Double-precision Kepler equation solver (reference implementation).
 * @details This source file provides the implementation for solving Kepler's equation in double precision.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 */

#include "kepler_solver_double.h"
#include <math.h>

double solve_kepler_double(double Mk, double e) {
    /* use Newton-Raphson iteration to solve Kepler's equation */
    double E = Mk;
    for (int i = 0; i < 10; i++) {
        double f = E - e * sin(E) - Mk;
        double f_prime = 1.0 - e * cos(E);
        double delta = f / f_prime;
        E = E - delta;
        /* Early exit once converged; also avoids extra sin/cos calls once
         * the update is already below double precision's useful resolution. */
        if (fabs(delta) < 1e-14) {
            break;
        }
    }
    return E;
}