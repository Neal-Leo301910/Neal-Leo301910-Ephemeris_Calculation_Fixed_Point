/**
 * @file error_eval.h
 * @brief Compares calc_sat_pos_fixed() against calc_sat_pos_double() over a
 * range of tk and reports position-error statistics (max / mean / RMS).
 */
#pragma once

#include "common_types.h"
#include "ephemeris_double.h"
#include "ephemeris_fixed.h"

#include <math.h>
#include <stdio.h>

/* tk sweep range, seconds (matches the design doc's suggested +-7200s scan). */
#define MIN_TK -7200
#define MAX_TK 7200

/**
 * 3D Euclidean distance between a fixed-point position (Q(q) format, e.g.
 * DIST_Q) and a double-precision reference position, in meters.
 * @param sat_pos_fixed_32  Fixed-point satellite position in Q(q) format.
 * @param sat_pos_double    Double-precision reference satellite position.
 * @param q                 Q format of the fixed-point position.
 * @return                  Euclidean distance between the fixed-point and double-precision positions, in meters.
 * 
 */
double evaluate_error(Vec3 sat_pos_fixed_32, Vec3d sat_pos_double, int q);

/**
 * Sweep tk from MIN_TK to MAX_TK (one point per second), compute
 * calc_sat_pos_fixed() vs calc_sat_pos_double() at each point, and print the
 * max/mean/RMS position error.
 * @param eph          Double-precision ephemeris data.
 * @param eph_fixed    Fixed-point ephemeris data.
 * @param q  Q format the fixed-point position (Vec3) is stored in -- this is
 *           DIST_Q, NOT the time format tk is expressed in.
 * @return                  0 on success, non-zero on failure.
 */
int show_error(Eph eph, Eph_fixed eph_fixed, int q);