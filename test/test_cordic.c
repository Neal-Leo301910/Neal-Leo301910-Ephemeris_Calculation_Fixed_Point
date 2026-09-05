#include "cordic.h"
#include "common_types.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static int g_failures = 0;
#define ANG_Q_TEST 28
#define TRIG_Q_TEST 30

static double q_to_double(fixed32_t v, int q) {
    return (double)v / (double)(1LL << q);
}
static fixed32_t double_to_q(double v, int q) {
    return (fixed32_t)llround(v * (double)(1LL << q));
}

static void check_close(const char *name, double got, double expected, double tol) {
    double err = fabs(got - expected);
    if (err > tol) {
        printf("[FAIL] %s: got %.8f, expected %.8f, err=%.2e (tol=%.2e)\n",
               name, got, expected, err, tol);
        g_failures++;
    } else {
        printf("[ OK ] %s: got %.8f, expected %.8f, err=%.2e\n",
               name, got, expected, err);
    }
}

int test_cordic(void) {
    printf("=== cordic_sincos_q28 vs libm sin/cos ===\n");
    double test_angles_deg[] = {0, 1, 15, 30, 45, 60, 75, 89, -30, -60, -89};
    int n = sizeof(test_angles_deg) / sizeof(test_angles_deg[0]);
    for (int i = 0; i < n; i++) {
        double rad = test_angles_deg[i] * M_PI / 180.0;
        fixed32_t z_q28 = double_to_q(rad, ANG_Q_TEST);
        fixed32_t sinQ, cosQ;
        cordic_sincos_q28(z_q28, &sinQ, &cosQ);
        double sin_got = q_to_double(sinQ, TRIG_Q_TEST);
        double cos_got = q_to_double(cosQ, TRIG_Q_TEST);
        char label[64];
        snprintf(label, sizeof(label), "sin(%g deg)", test_angles_deg[i]);
        check_close(label, sin_got, sin(rad), 1e-6);
        snprintf(label, sizeof(label), "cos(%g deg)", test_angles_deg[i]);
        check_close(label, cos_got, cos(rad), 1e-6);
    }

    printf("\n=== cordic_atan_q28 vs libm atan2 ===\n");
    /* (x, y) pairs spanning all four quadrants, magnitudes < 1 to match Q30 usage */
    double pairs[][2] = {
        {1.0, 0.0}, {0.7071, 0.7071}, {0.0, 1.0}, {-0.7071, 0.7071},
        {-1.0, 0.0}, {-0.7071, -0.7071}, {0.0, -1.0}, {0.7071, -0.7071},
        {0.5, 0.1}, {0.9, -0.2},
    };
    int np = sizeof(pairs) / sizeof(pairs[0]);
    for (int i = 0; i < np; i++) {
        fixed32_t x_q30 = double_to_q(pairs[i][0], TRIG_Q_TEST);
        fixed32_t y_q30 = double_to_q(pairs[i][1], TRIG_Q_TEST);
        fixed32_t z_q28 = cordic_atan_q28(x_q30, y_q30);
        double got = q_to_double(z_q28, ANG_Q_TEST);
        double expected = atan2(pairs[i][1], pairs[i][0]);
        char label[64];
        snprintf(label, sizeof(label), "atan2(%.4f, %.4f)", pairs[i][1], pairs[i][0]);
        check_close(label, got, expected, 1e-5);
    }

    printf("\n%d check(s) failed.\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}