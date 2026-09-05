/**
 * @file kepler_solver_fixed.h
 * @brief Fixed-point Kepler equation solver (reference implementation).
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
fixed32_t kepler_solve_fixed(fixed32_t Mk_q28, fixed32_t e_32_q30);

#if defined(__cplusplus)
}
#endif