#include "ephemeris.h"

#include <limits.h>
#include <stddef.h>

static ephemeris_fixed_t fixed_mul(ephemeris_fixed_t left, ephemeris_fixed_t right)
{
    return (ephemeris_fixed_t)(((int64_t)left * right) >> 16);
}

static ephemeris_fixed_t normalize_degrees(ephemeris_fixed_t degrees)
{
    const ephemeris_fixed_t circle = 360 * EPHEMERIS_FIXED_ONE;

    degrees %= circle;
    return degrees < 0 ? degrees + circle : degrees;
}

static ephemeris_fixed_t sine_degrees(ephemeris_fixed_t degrees)
{
    static const int32_t atan_q30[] = {
        843314857, 497837829, 263043837, 133525159, 67021687,
        33543516, 16775851, 8388437, 4194283, 2097149,
        1048575, 524288, 262144, 131072, 65536,
        32768, 16384, 8192, 4096, 2048,
        1024, 512, 256, 128, 64,
        32, 16, 8, 4, 2
    };
    const int32_t half_circle = 180 * EPHEMERIS_FIXED_ONE;
    const int32_t quadrant = 90 * EPHEMERIS_FIXED_ONE;
    int32_t angle;
    int32_t x = 652032874;
    int32_t y = 0;
    int32_t sign = 1;
    unsigned int i;

    degrees = normalize_degrees(degrees);
    if (degrees > half_circle) {
        degrees -= half_circle;
        sign = -1;
    }
    if (degrees > quadrant) {
        degrees = half_circle - degrees;
    }

    angle = (int32_t)(((int64_t)degrees * 18740330) >> 16);
    for (i = 0; i < sizeof(atan_q30) / sizeof(atan_q30[0]); ++i) {
        int32_t x_shift = x >> i;
        int32_t y_shift = y >> i;

        if (angle >= 0) {
            x -= y_shift;
            y += x_shift;
            angle -= atan_q30[i];
        } else {
            x += y_shift;
            y -= x_shift;
            angle += atan_q30[i];
        }
    }
    return (ephemeris_fixed_t)((sign * (int64_t)y) >> 14);
}

static int is_leap_year(int year)
{
    return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

int ephemeris_days_since_j2000(
    int year, unsigned month, unsigned day,
    unsigned hour, unsigned minute, unsigned second,
    ephemeris_fixed_t *days)
{
    static const unsigned char month_lengths[] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };
    int adjusted_year;
    unsigned adjusted_month;
    int64_t julian_day;
    int64_t seconds_since_noon;
    int64_t result;

    if (days == NULL || month < 1 || month > 12 || hour > 23 ||
        minute > 59 || second > 59 ||
        day < 1 || day > (unsigned int)(month_lengths[month - 1] +
            (month == 2 && is_leap_year(year)))) {
        return EPHEMERIS_INVALID_ARGUMENT;
    }

    adjusted_year = year - (month <= 2);
    adjusted_month = month + (month <= 2 ? 9 : -3);
    julian_day = (int64_t)(365 * adjusted_year) + adjusted_year / 4 -
        adjusted_year / 100 + adjusted_year / 400 +
        (153 * adjusted_month + 2) / 5 + day - 1 + 1721120;
    seconds_since_noon = ((int64_t)hour - 12) * 3600 +
        (int64_t)minute * 60 + second;
    result = (julian_day - 2451545) * EPHEMERIS_FIXED_ONE +
        (seconds_since_noon * EPHEMERIS_FIXED_ONE) / 86400;

    if (result < INT32_MIN || result > INT32_MAX) {
        return EPHEMERIS_OUT_OF_RANGE;
    }
    *days = (ephemeris_fixed_t)result;
    return EPHEMERIS_OK;
}

int ephemeris_sun_position(
    ephemeris_fixed_t days_since_j2000,
    ephemeris_sun_position_t *position)
{
    ephemeris_fixed_t mean_longitude;
    ephemeris_fixed_t mean_anomaly;
    ephemeris_fixed_t sin_anomaly;
    ephemeris_fixed_t longitude;

    if (position == NULL) {
        return EPHEMERIS_INVALID_ARGUMENT;
    }

    mean_longitude = normalize_degrees(
        18380650 + fixed_mul(64595, days_since_j2000));
    mean_anomaly = normalize_degrees(
        23431028 + fixed_mul(64592, days_since_j2000));
    sin_anomaly = sine_degrees(mean_anomaly);
    longitude = mean_longitude +
        fixed_mul(125475, sin_anomaly) +
        fixed_mul(1310, sine_degrees(2 * mean_anomaly)) +
        fixed_mul(19, sine_degrees(3 * mean_anomaly));

    position->ecliptic_longitude_deg = normalize_degrees(longitude);
    position->ecliptic_latitude_deg = 0;
    position->radius_au = 65545 -
        fixed_mul(1095, sine_degrees(mean_anomaly + 90 * EPHEMERIS_FIXED_ONE)) -
        fixed_mul(9, sine_degrees(2 * mean_anomaly + 90 * EPHEMERIS_FIXED_ONE));
    return EPHEMERIS_OK;
}
