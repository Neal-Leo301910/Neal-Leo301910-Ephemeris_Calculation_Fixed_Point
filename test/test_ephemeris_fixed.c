#include "common_types.h"
#include "kepler_solver_fixed.h"
#include "ephemeris_fixed.h"
#include "fixed_math.h"
#include "cordic.h"
#include "cordic_table.h"
#include <stdio.h>
#include <math.h>



static int g_failures = 0;

static void check(const char *name, int passed, const char *detail) {
    if (passed) {
        printf("[ OK ] %s%s%s\n", name, detail ? " -- " : "", detail ? detail : "");
    } else {
        printf("[FAIL] %s%s%s\n", name, detail ? " -- " : "", detail ? detail : "");
        g_failures++;
    }
}

static Eph_fixed make_sample_eph_fixed(void) {
    Eph_fixed f;
    f.toe_32_q11     = double_to_fixed32(0.0, TIME_Q);
    f.A_32_q5        = double_to_fixed32(26560000.0, DIST_Q);
    f.e_32_q30       = double_to_fixed32(0.01, TRIG_Q);
    f.M0_32_q28      = double_to_fixed32(0.5, ANG_Q);
    f.delta_n_32_q28 = double_to_fixed32(0.0, ANG_Q);
    f.omega_32_q28   = double_to_fixed32(0.3, ANG_Q);
    f.Omega_0_32_q28 = double_to_fixed32(1.0, ANG_Q);
    f.i0_32_q28      = double_to_fixed32(0.9599310885968813, ANG_Q);
    f.Omega_e_32_q31 = double_to_fixed32(7.2921151467e-5, RATE_Q);
    return f;
}

static double vec3_norm_m(Vec3 v) {
    double x = fixed32_to_double(v.x, DIST_Q);
    double y = fixed32_to_double(v.y, DIST_Q);
    double z = fixed32_to_double(v.z, DIST_Q);
    return sqrt(x * x + y * y + z * z);
}

int test_ephemeris_fixed(void) {
    Eph_fixed eph = make_sample_eph_fixed();
    double A_real = 26560000.0;

    /* --- basic sanity: a valid ephemeris must produce a non-zero, ---
     * --- plausible-magnitude position                              --- */
    Vec3 pos0 = calc_sat_pos_fixed(double_to_fixed32(0.0, TIME_Q), eph);
    double r0 = vec3_norm_m(pos0);
    check("t=0 position is non-zero", pos0.x != 0 || pos0.y != 0 || pos0.z != 0, NULL);

    char detail[128];
    snprintf(detail, sizeof(detail), "r=%.1f m, A=%.1f m", r0, A_real);
    /* orbit is nearly circular (e=0.01), so |r - A| / A should be small
     * (a few percent at most). This mainly catches gross Q-format /
     * overflow regressions, not tight precision. */
    check("t=0 orbital radius close to A", fabs(r0 - A_real) / A_real < 0.05, detail);

    /* --- invalid ephemeris (A<=0) must fail safely, not crash --- */
    Eph_fixed bad_eph = {0};
    Vec3 bad_pos = calc_sat_pos_fixed(0, bad_eph);
    check("invalid ephemeris (A<=0) returns {0,0,0}",
          bad_pos.x == 0 && bad_pos.y == 0 && bad_pos.z == 0, NULL);

    /* --- position must actually move over time (not stuck/degenerate) --- */
    Vec3 pos_600 = calc_sat_pos_fixed(double_to_fixed32(600.0, TIME_Q), eph);
    check("position changes between t=0 and t=600",
          pos0.x != pos_600.x || pos0.y != pos_600.y || pos0.z != pos_600.z, NULL);

    /* --- sweep a wide tk range: regression test for the CORDIC 90-degree
     * folding bug, the half-week rollover bug, and general hangs/crashes.
     * Every sample's radius should stay close to A the whole way through. */
    int sweep_ok = 1;
    double max_radius_err_pct = 0.0;
    for (int t = -7200; t <= 7200; t += 30) {
        Vec3 p = calc_sat_pos_fixed(double_to_fixed32((double)t, TIME_Q), eph);
        double r = vec3_norm_m(p);
        double err_pct = fabs(r - A_real) / A_real * 100.0;
        if (err_pct > max_radius_err_pct) max_radius_err_pct = err_pct;
        if (err_pct > 5.0) { /* 5% is a generous bound; a folding/rollover
                                 bug tends to blow this up by orders of
                                 magnitude, not a few percent */
            sweep_ok = 0;
        }
    }
    snprintf(detail, sizeof(detail), "max radius error over sweep: %.3f%%", max_radius_err_pct);
    check("tk sweep [-7200, 7200]s stays within 5% of A", sweep_ok, detail);

    printf("\n%d check(s) failed.\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}