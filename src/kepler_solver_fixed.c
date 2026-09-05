#include "kepler_solver_fixed.h"
#include "fixed_math.h"
#include "cordic.h"
#include "common_types.h"

fixed32_t kepler_solve_fixed(fixed32_t Mk_q28, fixed32_t e_32_q30) {
    fixed32_t Mk_q30 = Mk_q28 << (TRIG_Q - ANG_Q); /* Mk converted once, Q28 -> Q30 */
    fixed32_t E_q30 = Mk_q30;

    for (int i = 0; i < 10; i++) {
        fixed32_t sin_E_q30 = 0;
        fixed32_t cos_E_q30 = 0;
        
        cordic_sincos_q28(E_q30 >> (TRIG_Q - ANG_Q), &sin_E_q30, &cos_E_q30);

        fixed32_t f_q30 = E_q30 - qmul_ab(e_32_q30, sin_E_q30, TRIG_Q, TRIG_Q, TRIG_Q) - Mk_q30;
        fixed32_t f_prime_q30 = (1 << TRIG_Q) - qmul_ab(e_32_q30, cos_E_q30, TRIG_Q, TRIG_Q, TRIG_Q);


        fixed32_t delta_q30 = qdiv(f_q30, f_prime_q30, TRIG_Q);
        E_q30 = E_q30 - delta_q30;

        if (delta_q30 == 0) {
            break;
        }
    }
    return E_q30;
}