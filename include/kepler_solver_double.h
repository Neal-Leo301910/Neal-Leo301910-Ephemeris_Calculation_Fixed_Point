/**
 * @file kepler_solver_double.h
 * @brief Double-precision Kepler equation solver (reference implementation).
 * @details This header provides the interface for solving Kepler's equation in double precision.
 * @author Qiwei Yang
 * @date 2024-06-05
 * @version 1.0
 */
#pragma once
#include "common_types.h"

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * Solve Kepler's equation E = Mk + e*sin(E) for the eccentric anomaly E,
 * using Newton-Raphson iteration.
 *
 * @param Mk  Mean anomaly (radians), should already be normalized.
 * @param e   Eccentricity, expected in [0, 1).
 * @return    Eccentric anomaly E (radians), NOT normalized by this function.
 */
double solve_kepler_double(double Mk, double e);

#if defined(__cplusplus)
}
#endif