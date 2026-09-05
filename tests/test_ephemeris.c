#include "ephemeris.h"

#include <assert.h>
#include <stdio.h>

static int absolute_difference(ephemeris_fixed_t left, ephemeris_fixed_t right)
{
    int difference = left - right;
    return difference < 0 ? -difference : difference;
}

int main(void)
{
    ephemeris_fixed_t days;
    ephemeris_sun_position_t position;

    assert(ephemeris_days_since_j2000(2000, 1, 1, 12, 0, 0, &days) == EPHEMERIS_OK);
    assert(days == 0);
    assert(ephemeris_days_since_j2000(2000, 2, 29, 12, 0, 0, &days) == EPHEMERIS_OK);
    assert(days == 59 * EPHEMERIS_FIXED_ONE);
    assert(ephemeris_days_since_j2000(2001, 2, 29, 12, 0, 0, &days) ==
           EPHEMERIS_INVALID_ARGUMENT);
    assert(ephemeris_days_since_j2000(1800, 1, 1, 0, 0, 0, &days) ==
           EPHEMERIS_OUT_OF_RANGE);

    assert(ephemeris_sun_position(0, &position) == EPHEMERIS_OK);
    assert(absolute_difference(position.ecliptic_longitude_deg,
           18374984) < 3277);
    assert(absolute_difference(position.radius_au,
           64442) < 66);
    assert(position.ecliptic_latitude_deg == 0);
    assert(ephemeris_sun_position(0, NULL) == EPHEMERIS_INVALID_ARGUMENT);

    puts("ephemeris tests passed");
    return 0;
}
