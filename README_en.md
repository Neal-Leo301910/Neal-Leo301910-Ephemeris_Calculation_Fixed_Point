# Fixed-Point Ephemeris Computation

A pure fixed-point (32-bit integer) implementation of GNSS satellite ephemeris computation, including Kepler equation solving, CORDIC-based trigonometric functions, and full coordinate transformations. This project provides both a floating-point reference implementation and a fixed-point production implementation for resource-constrained embedded systems.

## Table of Contents

- [Overview](#overview)
- [Data Types](#data-types)
- [Fixed-Point Format](#fixed-point-format)
- [Fixed-Point Arithmetic Library](#fixed-point-arithmetic-library)
- [Floating-Point Ephemeris Computation](#floating-point-ephemeris-computation)
- [CORDIC Trigonometric Functions](#cordic-trigonometric-functions)
- [Fixed-Point Ephemeris Computation](#fixed-point-ephemeris-computation)
- [Error Analysis](#error-analysis)

---

## Overview

This project implements a complete satellite ephemeris computation pipeline using only `int32_t` fixed-point arithmetic. The algorithm computes Earth-centered Earth-fixed (ECEF) coordinates $(x, y, z)$ from broadcast ephemeris parameters at a given time $t$.

Key components:

- **Fixed-point math library**: multiplication, division, square root, and power operations with configurable Q-format.
- **CORDIC engine**: sine, cosine, and `atan2` computed via iterative integer shifts and additions.
- **Kepler solver**: Newton-Raphson iteration for eccentric anomaly.
- **Coordinate transformation**: orbit-plane to ECEF via three-axis rotation.
- **Error evaluator**: RMS/MAX error metrics against double-precision reference.

---

## Data Types

Input and output parameters are defined in `include/common_types.h` for both fixed-point (`int32_t`) and floating-point (`double`) domains.

### Input Parameters

**Ephemeris parameters** (`struct Eph_fixed` / `struct Eph`):

| Symbol | Description | Unit / Note | Fixed-Point Member | Float Member |
|--------|-------------|-------------|-------------------|--------------|
| $t_{oe}$ | Ephemeris reference time | seconds | `toe_32_q11` | `toe` |
| $A$ | Semi-major axis | meters | `A_32_q5` | `A` |
| $e$ | Eccentricity | dimensionless | `e_32_q30` | `e` |
| $M_0$ | Mean anomaly at reference time | radians | `M0_32_q28` | `M0` |
| $\Delta n$ | Mean motion correction | rad/s | `delta_n_32_q28` | `delta_n` |
| $\omega$ | Argument of perigee | radians | `omega_32_q28` | `omega` |
| $\Omega_0$ | Right ascension of ascending node at reference time | radians | `Omega_0_32_q28` | `Omega_0` |
| $i_0$ | Inclination angle at reference time | radians | `i0_32_q28` | `i0` |
| $\mu$ | Earth's gravitational constant | $3.986004418 \times 10^{14}\ \text{m}^3/\text{s}^2$ | `mu_32_q28` | `mu` |
| $\Omega_e$ | Earth's rotation rate | $7.292115 \times 10^{-5}\ \text{rad/s}$ | `Omega_e_32_q31` | `Omega_e` |

Type aliases:

```c
typedef int32_t fixed32_t;
typedef int64_t fixed64_t;
```

**Observation time**:

- Fixed-point: `t_fixed_q11` (`fixed32_t`)
- Floating-point: `t` (`double`)

### Output Parameters

**ECEF position vector** (`struct Vec3` / `struct Vec3d`):

| Domain | Type | Fields |
|--------|------|--------|
| Fixed-point | `struct Vec3` | `x`, `y`, `z` (`fixed32_t`) |
| Floating-point | `struct Vec3d` | `x`, `y`, `z` (`double`) |

---

## Fixed-Point Format

All fixed-point values use Q-format notation: **Qx.y** denotes $x$ integer bits (including sign), $y$ fractional bits, and $x + y = 31$ for signed 32-bit values.

| Category | Variables | Q-Format | Range / Precision | Macro |
|----------|-----------|----------|-------------------|-------|
| Angles | $M, E, v, u, \Omega_0, i$ | Q3.28 | $\pm 4\pi$, $\approx 3.7 \times 10^{-9}\ \text{rad}$ | `ANG_Q` (28) |
| Trigonometric | $\sin, \cos$ | Q1.30 | $[-1, 1]$, $\approx 9.3 \times 10^{-10}$ | `TRIG_Q` (30) |
| Eccentricity | $e$ | Q1.30 | $[-1, 1]$ | `ONE_Q30` (`1 << 30`) |
| Distance | $A, r, x, y, z$ | Q26.5 | $\pm 2.68 \times 10^8\ \text{m}$, $\approx 0.03\ \text{m}$ | `DIST_Q` (5) |
| Time | $t_k$ | Q20.11 | — | `TIME_Q` (11) |
| Angular rate | $n, \Omega_e$ | Q0.31 | — | `RATE_Q` (31) |
| Velocity | $v_x, v_y, v_z$ | Q15.16 | — | — |

---

## Fixed-Point Arithmetic Library

The fixed-point arithmetic primitives are declared in `include/fixed_math.h` and implemented in `src/fixed_math.c`.

### API

```c
fixed32_t qmul(fixed32_t a, fixed32_t b, int q);                          // same-Q multiply
fixed32_t qmul_ab(fixed32_t a, fixed32_t b, int qa, int qb, int qc);      // cross-Q multiply
fixed32_t qdiv(fixed32_t a, fixed32_t b, int q);                          // same-Q divide
fixed32_t fixed_sqrt(fixed32_t a, int q);                                 // integer square root
fixed32_t fixed_pow(fixed32_t a, int n, int q);                           // integer power
fixed32_t double_to_fixed32(double value, int q);                         // float → fixed
double    fixed32_to_double(fixed32_t value, int q);                      // fixed → float
```

### Arithmetic Definitions

Same-Q multiplication:
$$c_{\text{fixed-}Q} = \frac{(a \cdot 2^Q) \times (b \cdot 2^Q)}{2^Q}$$

Same-Q division:
$$c_{\text{fixed-}Q} = \frac{(a \cdot 2^Q)}{(b \cdot 2^Q)} \cdot 2^Q$$

Square root:
$$c_{\text{fixed-}Q} = \sqrt{(a \cdot 2^Q) \cdot 2^Q} = \sqrt{a} \cdot 2^Q$$

Integer power:
$$c_{\text{fixed-}Q} = (a \cdot 2^Q)^n = a^n \cdot 2^{Qn}$$

Cross-Q multiplication:
$$c_{\text{fixed-}Q_c} = \frac{(a \cdot 2^{Q_a}) \times (b \cdot 2^{Q_b})}{2^{Q_a + Q_b - Q_c}}$$

### Implementation Notes

1. **Overflow protection**: All intermediate products are promoted to `int64_t` before scaling to prevent `int32_t` overflow.
2. **Domain checks**: Division includes divide-by-zero protection; square root includes non-negative input validation.

### Integer Square Root Algorithm

The square root operates entirely on integers without floating-point conversion.

**Algorithm (bit-by-bit restoration)**:

Given input $S = a \cdot 2^{2Q}$ (pre-shifted), the algorithm finds the largest integer $R$ such that $R^2 \leq S$.

```c
/* Pseudocode */
fixed64_t bit = (fixed64_t)1 << 62;   /* highest 4-power <= 2^63 */
fixed64_t result = 0;
fixed64_t remainder = S;

while (bit != 0) {
    fixed64_t trial = result + bit;
    if (remainder >= trial) {
        remainder -= trial;
        result = (result >> 1) + bit;
    } else {
        result >>= 1;
    }
    bit >>= 2;
}
/* result = floor(sqrt(S)) */
```

**Key insight**: Because the operation is a square root, the trial bit is shifted by 2 positions per iteration (one bit position in the result corresponds to two bit positions in the square). The recurrence uses the identity $(R + b)^2 = R^2 + 2Rb + b^2$ with dynamic rescaling to avoid overflow.

### Floating-Point to Fixed-Point Conversion

For values requiring high fractional precision (e.g., $\Omega_e$ in Q31), direct computation of `1LL << q` must be performed in 64-bit (`long long`) to avoid sign-bit overflow in 32-bit `int`.

---

## Floating-Point Ephemeris Computation

The reference implementation is declared in `include/solve_kepler_double.h` and `include/ephemeris_double.h`, with implementations in `src/solve_kepler_double.c` and `src/ephemeris_double.c`.

### Computation Pipeline

```
Inputs: (toe, A, e, M0, Δn, ω, Ω0, i0, μ, Ωe) + time t
        ↓
   tk = t - toe  (normalized to ±302400 s)
        ↓
   n = sqrt(μ / A³) + Δn
        ↓
   Mk = M0 + n · tk
        ↓
   Ek  (Kepler equation via Newton-Raphson)
        ↓
   v = atan2( sqrt(1-e²)·sin(Ek), cos(Ek)-e )
   u = v + ω
        ↓
   r = A · (1 - e·cos(Ek))
   x' = r·cos(u),  y' = r·sin(u)
        ↓
   i  = i0
   Ω  = Ω0 + (Ω_dot - Ωe)·tk - Ωe·toe
        ↓
   x = x'·cos(Ω) - y'·cos(i)·sin(Ω)
   y = x'·sin(Ω) + y'·cos(i)·cos(Ω)
   z = y'·sin(i)
```

### Step-by-Step Algorithm

1. **Time difference**: $t_k = t - t_{oe}$. Normalize across week boundary:
   $$t_k = \begin{cases} t_k - 604800 & t_k > 302400 \\ t_k + 604800 & t_k < -302400 \end{cases}$$

2. **Mean motion**: $n_0 = \sqrt{\mu / A^3}$, then $n = n_0 + \Delta n$.

3. **Mean anomaly**: $M_k = M_0 + n \cdot t_k$.

4. **Kepler equation** (Newton-Raphson iteration for eccentric anomaly $E_k$):
   $$\begin{aligned}
   E_k^{(0)} &= M_k \\
   E_k^{(i)} &= E_k^{(i-1)} - \frac{E_k^{(i-1)} - e \cdot \sin(E_k^{(i-1)}) - M_k}{1 - e \cdot \cos(E_k^{(i-1)})}
   \end{aligned}$$

5. **True anomaly**: $v = \text{atan2}\big(\sqrt{1-e^2} \cdot \sin(E_k),\; \cos(E_k) - e\big)$.

6. **Argument of latitude**: $u = v + \omega$.

7. **Orbital radius and plane coordinates**:
   $$r = A \cdot \big(1 - e \cdot \cos(E_k)\big),\quad x' = r \cos u,\quad y' = r \sin u$$

8. **Orbital orientation**:
   $$i = i_0,\quad \Omega = \Omega_0 + (\dot\Omega - \Omega_e) \cdot t_k - \Omega_e \cdot t_{oe}$$

9. **ECEF transformation**:
   $$\begin{aligned}
   x &= x' \cos\Omega - y' \cos i \,\sin\Omega \\
   y &= x' \sin\Omega + y' \cos i \,\cos\Omega \\
   z &= y' \sin i
   \end{aligned}$$

**Normalization requirements**: All angular quantities ($M_k$, $E_k$, $v$, $u$, $\Omega$) must be normalized to $(-\pi, \pi]$. Time difference $t_k$ is normalized to one week.

---

## CORDIC Trigonometric Functions

CORDIC (COordinate Rotation DIgital Computer) computes $\sin$, $\cos$, and $\text{atan2}$ using only integer addition, subtraction, shifts, and a small lookup table.

### Rotation Principle

A vector $(x_i, y_i)$ rotated by angle $\theta_i$ becomes:
$$\begin{aligned}
x_{i+1} &= x_i \cos\theta_i - y_i \sin\theta_i = K_i \,(x_i - \sigma_i \, y_i \, 2^{-i}) \\
y_{i+1} &= y_i \cos\theta_i + x_i \sin\theta_i = K_i \,(y_i + \sigma_i \, x_i \, 2^{-i}) \\
K_i &= \cos\big(\arctan(2^{-i})\big),\quad \sigma_i = \pm 1
\end{aligned}$$

The aggregate scaling factor $K = \prod_i K_i$ is precomputed and applied once at the end.

### sin / cos Mode

Input angle $z$ (normalized to $[-\pi/2, \pi/2]$); the algorithm rotates from $0^\circ$ to target $z$.

```c
/* Pseudocode */
x0 = 1 / K;   y0 = 0;   z0 = z;
for (i = 0; i < N; i++) {
    sigma = (z_i >= 0) ? +1 : -1;
    x_{i+1} = x_i - sigma * (y_i >> i);
    y_{i+1} = y_i + sigma * (x_i >> i);
    z_{i+1} = z_i - sigma * atan_table[i];
}
cos(z) ≈ x_N * K
sin(z) ≈ y_N * K
```

### atan2 Mode

Input vector $(x_0, y_0)$; the algorithm rotates $y$ toward $0$ and accumulates the rotation angle.

```c
/* Pseudocode */
z0 = 0;
for (i = 0; i < N; i++) {
    sigma = (y_i < 0) ? +1 : -1;   /* rotate y toward 0 */
    x_{i+1} = x_i - sigma * (y_i >> i);
    y_{i+1} = y_i + sigma * (x_i >> i);
    z_{i+1} = z_i - sigma * atan_table[i];
}
atan2(y0, x0) ≈ z_N
```

### Implementation Notes

1. **Input range**: Ephemeris angles use Q28. Input range is validated to $[-\pi, \pi]$ (scaled to Q28).
2. **Scaling factor**: $K$ exceeds Q28 range; internal CORDIC computation uses Q30 for $x_0$ and $y_0$.
3. **Convergence range**: Standard CORDIC converges in quadrants I and IV only. For quadrants II and III, the input vector is pre-rotated by $\pi$ (180°) into the convergence domain; the recorded rotation is subtracted from the final result.
4. **Lookup table generation**: `include/cordic_table.h` is auto-generated by `tools/atan_table_gen.py` to eliminate manual calculation error and ensure consistency. The table is declared as `static const fixed32_t` to limit scope to the compilation unit and avoid global namespace pollution.

---

## Fixed-Point Ephemeris Computation

The fixed-point implementation follows the same algorithmic pipeline as the floating-point reference, with explicit Q-format management at each stage.

Declarations: `include/solve_kepler_fixed.h`, `include/ephemeris_fixed.h`  
Implementations: `src/solve_kepler_fixed.c`, `src/ephemeris_fixed.c`

### Q-Format Flow

```
Inputs: (toe_q11, A_q5, e_q30, M0_q28, Δn_q28, ω_q28, Ω0_q28, i0_q28, mu_qN, Ωe_q31) + t_q11
        ↓
   tk_q11
        ↓
   n_q28      (dynamic precision scaling, intermediate in Q5)
        ↓
   Mk_q28
        ↓
   Ek_q30     (Kepler iteration in Q30)
   Ek_q28     (down-converted for subsequent angle ops)
        ↓
   v_q28 → u_q28
        ↓
   r_q5
        ↓
   x'_q5, y'_q5
        ↓
   x_q5, y_q5, z_q5   (ECEF output)
```

### Dynamic-Precision Mean Motion

The computation $n_0 = \sqrt{\mu / A^3}$ involves a wide dynamic range. To maximize precision without overflow:

1. Compute initial ratio $\mu / A$ in Q0.
2. Dynamically left-shift the ratio by $n$ bits (where $n < 8$ and $(\mu/A) \ll n < \text{INT32\_MAX}$) to adaptive Q-format $Q_n$.
3. Compute fixed-point square root $\sqrt{\mu / A}$ in $Q_n$.
4. Compute fixed-point division $\sqrt{\mu / A} \,/\, A$ to obtain final $n_0$.

### Implementation Notes

1. Cross-Q operations must use `qmul_ab()` from `fixed_math.c`.
2. CORDIC `sin`/`cos` operates in Q30 internally; results are down-converted to Q28 for angle-domain compatibility.
3. Eccentricity $e$ remains in Q30 for Kepler iteration; down-converted to Q28 only when combined with angle quantities.

---

## Error Analysis

The error evaluation module is declared in `include/error_eval.h` and implemented in `src/error_eval.c`.

### Error Metric

$$\text{error} = \sqrt{(x_{\text{fixed}} - x_{\text{double}})^2 + (y_{\text{fixed}} - y_{\text{double}})^2 + (z_{\text{fixed}} - z_{\text{double}})^2}$$

### Test Conditions

- Time span: $t_k \in [-7200, +7200]$ seconds
- Sampling interval: 10 / 30 / 60 seconds

### Results

| Metric | Absolute Error | Relative Error ($/ A$) |
|--------|---------------|------------------------|
| Maximum | 105.998 m | $3.99 \times 10^{-6}$ |
| Mean | 52.669 m | $1.98 \times 10^{-6}$ |
| RMS | 60.816 m | $2.29 \times 10^{-6}$ |

Test semi-major axis reference: $A = 26\,560\,000$ m.

---
