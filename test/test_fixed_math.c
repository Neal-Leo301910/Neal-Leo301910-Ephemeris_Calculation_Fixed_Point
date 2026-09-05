#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fixed_math.h"


static int g_failures = 0;

static void check_eq_i32(const char *name, fixed32_t got, fixed32_t expected) {
    if (got != expected) {
        printf("[FAIL] %s: got %d, expected %d\n", name, got, expected);
        g_failures++;
    } else {
        printf("[ OK ] %s: %d\n", name, got);
    }
}

int test_fixed_math(void) {
    /* --- Basic sanity checks (same as original test) --- */
    fixed32_t a = 1 << 16; /* 1.0 in Q16 */
    fixed32_t b = 2 << 16; /* 2.0 in Q16 */
    int q = 16;
    int qc = 20;

    check_eq_i32("qmul(1.0, 2.0) Q16 -> 2.0",
                 qmul(a, b, q), 2 << 16);
    check_eq_i32("qmul_ab(1.0 Q16, 2.0 Q16 -> Q20) -> 2.0",
                 qmul_ab(a, b, q, q, qc), 2 << 20);
    check_eq_i32("qdiv(1.0/2.0) Q16 -> 0.5",
                 qdiv(a, b, q), 1 << 15);
    check_eq_i32("fixed_sqrt(1.0) Q16 -> 1.0",
                 fixed_sqrt(a, q), 1 << 16);

    /* --- From the design doc worked example: sqrt(2.25) in Q4 -> 1.5 --- */
    check_eq_i32("fixed_sqrt(2.25) Q4 -> 1.5",
                 fixed_sqrt(36, 4), 24);

    /* --- Negative operands (qmul/qmul_ab/qdiv must handle sign correctly) --- */
    check_eq_i32("qmul(-1.5, 2.0) Q16 -> -3.0",
                 qmul((fixed32_t)(-1.5 * (1 << 16)), b, q),
                 (fixed32_t)(-3.0 * (1 << 16)));
    check_eq_i32("qdiv(-1.0/2.0) Q16 -> -0.5",
                 qdiv((fixed32_t)(-1.0 * (1 << 16)), b, q),
                 (fixed32_t)(-0.5 * (1 << 16)));

    /* --- Overflow guard: two large int32 Q30 values whose raw product would
     * overflow a 32-bit multiply (this is the a=b=0.8 example from the design
     * doc, ~7.38e17 intermediate). Must go through the int64 path correctly. */
    fixed32_t big_q30 = (fixed32_t)(0.8 * (1 << 30)); /* ~0.8 in Q30 */
    fixed32_t big_result = qmul(big_q30, big_q30, 30); /* expect ~0.64 in Q30 */
    double big_result_d = (double)big_result / (1 << 30);
    printf("[INFO] qmul(0.8,0.8) Q30 -> %d (%.6f), expected ~0.64\n",
           big_result, big_result_d);
    if (big_result_d < 0.635 || big_result_d > 0.645) {
        printf("[FAIL] qmul Q30 overflow case out of tolerance\n");
        g_failures++;
    } else {
        printf("[ OK ] qmul Q30 overflow case within tolerance\n");
    }

    /* --- Distance-scale sqrt, matching DIST_Q usage (Q5, ~2.6e7 m range) --- */
    fixed32_t r2_q5 = (fixed32_t)(26000000.0 * (1 << 5)); /* (26,000,000 m)^0-ish placeholder magnitude */
    fixed32_t sqrt_r2 = fixed_sqrt(r2_q5, 5);
    printf("[INFO] fixed_sqrt(26000000) Q5 -> %d (%.3f)\n",
           sqrt_r2, (double)sqrt_r2 / (1 << 5));

    /* --- Error-path checks: must not crash --- */
    fixed32_t sqrt_neg = fixed_sqrt(-100, q);
    check_eq_i32("fixed_sqrt(negative) returns 0", sqrt_neg, 0);

    fixed32_t div_by_zero = qdiv(a, 0, q);
    check_eq_i32("qdiv(x, 0) returns 0", div_by_zero, 0);

    printf("\n%d check(s) failed.\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}