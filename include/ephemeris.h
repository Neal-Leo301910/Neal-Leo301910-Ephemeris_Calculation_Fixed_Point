#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <stdint.h>

/*
 * All values use signed Q16.16 fixed-point notation.  Angles are degrees,
 * and distances are astronomical units.
 */
typedef int32_t ephemeris_fixed_t;

#define EPHEMERIS_FIXED_ONE ((ephemeris_fixed_t)65536)

typedef struct {
    ephemeris_fixed_t ecliptic_longitude_deg;
    ephemeris_fixed_t ecliptic_latitude_deg;
    ephemeris_fixed_t radius_au;
} ephemeris_sun_position_t;

enum {
    EPHEMERIS_OK = 0,
    EPHEMERIS_INVALID_ARGUMENT = -1,
    EPHEMERIS_OUT_OF_RANGE = -2
};

/*
 * Converts a Gregorian UTC date and time to days from J2000.0
 * (2000-01-01 12:00 UTC).  The representable Q16.16 range is enforced.
 */
int ephemeris_days_since_j2000(
    int year, unsigned month, unsigned day,
    unsigned hour, unsigned minute, unsigned second,
    ephemeris_fixed_t *days);

/* Calculates the geocentric ecliptic position of the Sun. */
int ephemeris_sun_position(
    ephemeris_fixed_t days_since_j2000,
    ephemeris_sun_position_t *position);

#endif
