# 星历定点化计算

面向嵌入式平台的卫星星历计算项目，提供双精度浮点实现与纯整数定点实现，并通过误差扫描和模块化测试验证定点计算的精度与稳定性。

----
## 目录

- [项目概述](#1. 项目概述)
- [输入与输出](#2. 输入与输出)
- [定点数格式](#3. 定点数格式)
- [计算流程](#4. 计算流程)
- [CORDIC 实现](#5. CORDIC 实现)
- [误差评估](#6. 误差评估)
- [项目结构](#7. 项目结构)
- [构建与运行](#8. 构建与运行)
- [开发约定](#9. 开发约定)

----
## 1. 项目概述

本项目实现基于开普勒轨道模型的卫星位置计算，支持：

- 双精度浮点星历计算，作为参考实现；
- 基于 `int32_t` / `int64_t` 的定点数计算；
- 基于 CORDIC 的 `sin`、`cos` 和 `atan2` 计算；
- 基于纯整数算法的平方根计算；
- 定点实现与浮点参考实现之间的位置误差评估；
- 定点数学库、CORDIC、开普勒方程求解器和星历计算模块的独立测试。

> 当前 README 记录算法原理、数据格式、接口约定、目录结构和构建方式。具体接口以源代码中的头文件声明为准。

----
## 2. 输入与输出

输入和输出数据类型定义于 `include/common_types.h`。

### 2.1 星历参数

项目同时提供浮点结构体 `Eph` 和定点结构体 `Eph_fixed`。

| 符号 | 含义 | 浮点字段 | 定点字段 |
|---|---|---|---|
| `toe` | 星历参考时刻 | `toe` | `toe_32_q11` |
| `A` | 轨道半长轴 | `A` | `A_32_q5` |
| `e` | 偏心率 | `e` | `e_32_q30` |
| `M0` | 参考时刻平近点角 | `M0` | `M0_32_q28` |
| `delta_n` | 平均角速度修正项 | `delta_n` | `delta_n_32_q28` |
| `omega` | 近地点幅角 | `omega` | `omega_32_q28` |
| `Omega_0` | 参考时刻升交点赤经 | `Omega_0` | `Omega_0_32_q28` |
| `i0` | 参考时刻轨道倾角 | `i0` | `i0_32_q28` |
| `mu` | 地球引力常数 | `mu` | `mu_32_q28` |
| `Omega_e` | 地球自转角速度 | `Omega_e` | `Omega_e_32_q31` |

常量：

- `mu = 3.986004418e14 m^3/s^2`
- `Omega_e = 7.292115e-5 rad/s`

关注时刻：

- 浮点版本：`double t`
- 定点版本：`fixed32_t t_fixed_q11`

输出为地球坐标系下的卫星位置：

- 浮点版本：`Vec3d`
- 定点版本：`Vec3`

两种结构均包含 `x`、`y`、`z` 三个坐标分量。

----
## 3. 定点数格式

项目使用 Q 格式表示定点数。对于 32 位有符号定点数，`Qx.y` 表示整数部分使用 `x` 位、小数部分使用 `y` 位，另有 1 位符号位，因此 `x + y = 31`。

| 变量类别  | 变量                      |   Q 格式 | 代码宏           |
| ----- | ----------------------- | -----: | ------------- |
| 角度    | `M`、`E`、`v`、`u`、`Omega` |  Q3.28 | `ANG_Q = 28`  |
| 三角函数值 | `sin`、`cos`             |  Q1.30 | `TRIG_Q = 30` |
| 偏心率   | `e`                     |  Q1.30 | `ONE_Q30`     |
| 距离    | `A`、`r`、`x`、`y`、`z`     |  Q26.5 | `DIST_Q = 5`  |
| 时间    | `tk`                    | Q20.11 | `TIME_Q = 11` |
| 角速度   | `n`、`Omega_e`           |  Q0.31 | `RATE_Q = 31` |
| 速度    | `vx`、`vy`、`vz`          | Q15.16 | —             |

定点类型别名：

```c
typedef int32_t fixed32_t;
typedef int64_t fixed64_t;
```

### 3.1 定点运算约定

相关接口定义于 `include/fixed_math.h`，实现位于 `src/fixed_math.c`：

```c
fixed32_t qmul(fixed32_t a, fixed32_t b, int q);
fixed32_t qmul_ab(fixed32_t a, fixed32_t b, int qa, int qb, int qc);
fixed32_t qdiv(fixed32_t a, fixed32_t b, int q);
fixed32_t fixed_sqrt(fixed32_t a, int q);
fixed32_t fixed_pow(fixed32_t a, int n, int q);
double fixed32_to_double(fixed32_t value, int q);
fixed32_t double_to_fixed32(double value, int q);
```

定点数在相同 Q 格式下：
- 乘法：
$$
c_{Q} = (a \cdot 2^Q) \times (b \cdot 2^Q) / 2^Q
$$
- 除法：
$$
c_{Q} = (a \cdot 2^Q) / (b \cdot 2^Q) \cdot 2^Q
$$
- 开方：
$$
c_{Q} = \sqrt{(a\cdot2^Q)\cdot2^Q} = \sqrt a\cdot 2^Q
$$
- 平方/幂：
$$
c_{Q} = (a\cdot 2^Q)^n=a^n \cdot 2^{Qn}
$$
```text
mul(a, b, Q) = (a * b) >> Q
div(a, b, Q) = (a << Q) / b
sqrt(a, Q)   = sqrt(a << Q)
```

不同 Q 格式下，目标格式为 `Qc`的定点数的乘法：
$$
c_{Qc} = (a \cdot 2^{Qa}) \times (b \cdot 2^{Qb}) / 2^{Qa+Qb-Qc}
$$

```text
c = (a * b) / 2^(Qa + Qb - Qc)
```

> 注意：
>
> 1. 中间乘法和左移操作可能溢出 32 位整数，应使用 `fixed64_t` 执行中间计算。
> 2. 除法必须进行除零检查。
> 3. 平方根必须检查输入是否为非负数。
> 4. 不同 Q 格式之间的运算应统一通过 `qmul_ab()` 等接口完成，避免隐式缩放错误。

----
## 4. 计算流程

### 4.1 浮点参考实现

浮点实现位于：

- `include/solve_kepler_double.h`
- `include/ephemeris_double.h`
- `src/solve_kepler_double.c`
- `src/ephemeris_double.c`

计算流程如下：

```text
输入星历参数和目标时刻
        ↓
计算时间差 tk，并进行跨周归一化
        ↓
计算平均角速度 n
        ↓
计算平近点角 Mk
        ↓
使用牛顿迭代求解偏近点角 Ek
        ↓
计算真近点角 v 和升交点角距 u
        ↓
计算轨道半径 r
        ↓
计算轨道面坐标 x'、y'
        ↓
进行轨道倾角和升交点赤经旋转
        ↓
输出地球坐标系位置 x、y、z
```

主要计算步骤：

1. `tk = t - toe`，并将其归一化到一周范围内；
2. `n0 = sqrt(mu / A^3)`，`n = n0 + delta_n`；
3. `Mk = M0 + n * tk`；
4. 使用牛顿法求解 `Ek - e sin(Ek) = Mk`；
5. `v = atan2(sqrt(1-e^2) sin(Ek), cos(Ek)-e)`；
6. `u = v + omega`；
7. `r = A(1-e cos(Ek))`；
8. 计算轨道面坐标并转换至地球坐标系。

角度统一归一化到 `(-pi, pi]`；时间差按一周 `604800 s` 进行跨周处理。

### 4.2 定点实现

定点实现位于：

- `include/solve_kepler_fixed.h`
- `include/ephemeris_fixed.h`
- `src/solve_kepler_fixed.c`
- `src/ephemeris_fixed.c`

定点实现与浮点实现采用相同的轨道计算流程，但在模块边界处必须明确输入和输出 Q 格式。

典型数据流：

```text
toe(Q11), A(Q5), e(Q30), M0(Q28), delta_n(Q28)
        ↓
tk(Q11)
        ↓
n(Q28)
        ↓
Mk(Q28)
        ↓
Ek(Q30) → Ek(Q28)
        ↓
CORDIC sin/cos
        ↓
v(Q28), u(Q28)
        ↓
r(Q5)
        ↓
x'(Q5), y'(Q5)
        ↓
坐标旋转
        ↓
x(Q5), y(Q5), z(Q5)
```

平均角速度计算采用动态定标：

1. 计算 `mu / A`；
2. 在不超过 32 位有符号整数范围的前提下动态放大中间结果；
3. 在动态 Q 格式下执行平方根；
4. 通过定点除法得到平均角速度；
5. 将结果转换为后续计算所需的 Q 格式。

----
## 5. CORDIC 实现

CORDIC 模块使用整数加法、减法、移位和查表操作计算三角函数，避免依赖浮点数学库。

相关实现位于：

- `include/cordic.h`
- `src/cordic.c`
- `include/cordic_table.h`

### 5.1 `sin` / `cos`

旋转模式下，CORDIC 通过逐次旋转使残余角度趋近于零：

```text
x(i+1) = x(i) - sigma * y(i) >> i
y(i+1) = y(i) + sigma * x(i) >> i
z(i+1) = z(i) - sigma * atan_table[i]
```

其中 `sigma` 由当前残余角度的符号决定。迭代完成后，对结果进行 CORDIC 增益补偿。

### 5.2 `atan2`

向量模式下，CORDIC 通过逐次旋转将输入向量的 `y` 分量逼近零，累计旋转角度即为 `atan2(y, x)`。

由于传统 CORDIC 的收敛范围有限，输入向量在进入迭代前需要进行象限处理。若向量位于第二或第三象限，则先进行 180° 对称变换，并在迭代结束后恢复象限。

### 5.3 查表生成

CORDIC 角度表由 Python 脚本自动生成：

```text
tools/atan_table_gen.py
```

脚本根据 `atan(2^-i)` 生成定点常量，输出至：

```text
include/cordic_table.h
```

使用 `static const` 保存查表数据，以限制符号作用域并避免重复分配。

----
## 6. 误差评估

误差评估模块位于：

- `include/error_eval.h`
- `src/error_eval.c`

位置误差定义为：

```text
error = sqrt(
    (x_fixed - x_double)^2 +
    (y_fixed - y_double)^2 +
    (z_fixed - z_double)^2
)
```

测试程序在 `tk = -7200 s ... +7200 s` 范围内进行扫描，并按照指定时间间隔采样。

当前记录的测试结果：

- 最大误差：`105.998050 m`
- 平均误差：`52.668958 m`
- 最小均方误差：`60.816272 m`
- 最大相对误差约为：`3.9909e-6`
- 平均相对误差约为：`1.9830e-6`

上述结果仅适用于对应测试参数、采样范围和当前实现版本。修改 Q 格式、迭代次数、CORDIC 表或输入星历后，应重新执行误差评估。

----
## 7. 项目结构

```text
ephemeris_fixedpoint/
├── Makefile
├── README.md
├── include/
│   ├── common_types.h
│   ├── fixed_math.h
│   ├── cordic.h
│   ├── cordic_table.h
│   ├── kepler_solver.h
│   ├── ephemeris_double.h
│   ├── ephemeris_fixed.h
│   └── error_eval.h
├── src/
│   ├── fixed_math.c
│   ├── cordic.c
│   ├── kepler_solver_double.c
│   ├── kepler_solver_fixed.c
│   ├── ephemeris_double.c
│   ├── ephemeris_fixed.c
│   └── error_eval.c
├── test/
│   ├── main.c
│   ├── test_fixed_math.c
│   ├── test_cordic.c
│   ├── test_kepler.c
│   └── test_data/
│       └── sample_ephemeris.csv
└── tools/
    ├── sample_ephemeris_generate.py
    └── atan_table_gen.py
    
```

----
## 8. 构建与运行

### 构建全部目标

```bash
make
# 或
make all
```

### 运行全部测试

```bash
make run
```

### 清理构建产物

```bash
make clean
```

该命令清理 `build/` 目录，并删除自动生成的 `include/cordic_table.h`。

### 构建指定目标

```bash
make build/main
make build/test_fixed_math
make build/test_cordic
make build/test_ephemeris
make build/test_ephemeris_fixed
```

### 运行主程序

```bash
echo "xxxx" | ./build/main
```

主程序读取目标时刻并输出计算结果。实际输入格式以 `test/main.c` 中的实现为准。

### 测试真实星历

本项目使用`sample_ephemeris_generate.py`脚本，加载真实星历的sample data到主程序中。

项目设定默认条目为：`meo_representative`

可自行添加条目后，将指定条目名替换主程序进行编译，以条目名`gps_prn27_20210624`为例：
```Shell
make DATA_TARGET=gps_prn27_20210624
```

`sample_ephemeris.csv`的存储规则：

,,,,,,,,,

| name                 | source                                                                                                                                                                                                             |    `toe` | `A`         | `e`      | `M0`   | `delta_n`  | `omega` | `Omega_0` | `i0`               | mu             | Omega_e         |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------: | ----------- | -------- | ------ | ---------- | ------- | --------- | ------------------ | -------------- | --------------- |
| `meo_representative` | "Synthetic representative MEO orbit (NOT from a real broadcast ephemeris) -- used throughout this project's earlier test files (test_ephemeris_double.c / test_ephemeris_fixed.c) as the baseline regression case" |      0.0 | 26560000.0  | 0.01     | 0.5    | 0.0        | 0.3     | 1.0       | 0.9599310885968813 | 3.986004418e14 | 7.2921151467e-5 |
| `gps_prn27_20210624` | "Real GPS broadcast ephemeris, PRN27, RINEX nav file 2021-06-24T01:59:44Z (GPS week 2163), via MathWorks rinexread() documentation: https://www.mathworks.com/help/nav/ref/rinexread.html",                        | 352780.0 | 26560623.69 | 0.009451 | 2.0809 | 4.3363e-09 | 0.63488 | 1.2866    | 0.97551            | 3.986004418e14 | 7.2921151467e-5 |


----
## 9. 开发约定

- 头文件只声明接口、类型和必要的常量；
- 源文件负责具体实现；
- 定点计算不得隐式混用不同 Q 格式；
- 所有可能溢出的中间计算优先使用 `int64_t`；
- 新增算法模块时，应同步增加对应单元测试；
- 修改定点格式、CORDIC 迭代次数或星历模型后，应重新生成查表并执行误差扫描；
- 浮点实现作为参考基准，不应与定点实现共享可能掩盖定点误差的计算路径。
