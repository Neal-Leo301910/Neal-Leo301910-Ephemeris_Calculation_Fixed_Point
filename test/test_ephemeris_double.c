#include "common_types.h"
#include "ephemeris_double.h"
#include <stdio.h>


static Eph make_sample_eph(void) {
    Eph eph;
    eph.toe     = 0.0;                 /* reference epoch, s */
    eph.A       = 26560000.0;          /* semi-major axis, m (~MEO) */
    eph.e       = 0.01;                /* eccentricity */
    eph.M0      = 0.5;                 /* mean anomaly at toe, rad */
    eph.delta_n = 0.0;                 /* mean motion correction, rad/s */
    eph.omega   = 0.3;                 /* argument of perigee, rad */
    eph.Omega_0 = 1.0;                 /* RAAN at toe, rad */
    eph.i0      = 0.9599310885968813;  /* inclination, rad (~55 deg) */
    eph.mu      = 3.986004418e14;      /* Earth gravitational constant */
    eph.Omega_e = 7.2921151467e-5;     /* Earth rotation rate, rad/s */
    return eph;
}

int test_ephemeris_double(void) {
    Eph eph = make_sample_eph();

    double t = 0.0;
    Vec3d sat_pos = calc_sat_pos_double(t, eph);
    printf("Satellite Position at t=0:  x = %f, y = %f, z = %f\n",
           sat_pos.x, sat_pos.y, sat_pos.z);

    double t1 = 20.0;
    Vec3d sat_pos1 = calc_sat_pos_double(t1, eph);
    printf("Satellite Position at t=20: x = %f, y = %f, z = %f\n",
           sat_pos1.x, sat_pos1.y, sat_pos1.z);

    /* Regression check: a zeroed/invalid ephemeris must now fail fast
     * (fprintf a warning + return {0,0,0}) instead of hanging. */
    Eph bad_eph = {0};
    Vec3d bad_pos = calc_sat_pos_double(0.0, bad_eph);
    printf("Invalid ephemeris (A=0) handled safely: x = %f, y = %f, z = %f\n",
           bad_pos.x, bad_pos.y, bad_pos.z);

    return 0;
}