import csv
import os
import re

def update_c_file(target_name="meo_representative"):
    # 1. 动态获取当前脚本 (sample_ephemeris_generate.py) 的绝对路径
    script_path = os.path.abspath(__file__)
    # 2. 获取脚本所在的本地目录 (即 tools 文件夹的绝对路径)
    tools_dir = os.path.dirname(script_path)
    # 3. 向上推一级，强行锁定项目根目录的绝对路径
    project_root = os.path.dirname(tools_dir)

    # 4. 基于绝对的项目根目录，精准拼接出 CSV 和 main.c 的绝对路径
    csv_path = os.path.normpath(os.path.join(project_root, "test", "test_data", "sample_ephemeris.csv"))
    c_path = os.path.normpath(os.path.join(project_root, "src", "main.c"))

    # 5. 从 CSV 文件中提取目标行的数据
    ephemeris_data = {}
    if not os.path.exists(csv_path):
        print(f"【错误】: 找不到 CSV 文件。")
        print(f"  -> 尝试访问的绝对路径为: {csv_path}")
        return

    with open(csv_path, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        for row in reader:
            if row['name'] == target_name:
                ephemeris_data = {k: row[k].strip() for k in fieldnames if k not in ['name', 'source']}
                break
    
    if not ephemeris_data:
        print(f"【错误】: 在 CSV 中未找到 name 为 '{target_name}' 的数据")
        return

    # 6. 读取 C 文件
    if not os.path.exists(c_path):
        print(f"【错误】: 找不到 C 文件。")
        print(f"  -> 尝试访问的绝对路径为: {c_path}")
        return

    with open(c_path, mode='r', encoding='utf-8') as f:
        c_content = f.read()

    # 7. 构造要替换进去的全新标准 C 函数文本
    new_functions = f"""static Eph make_sample_eph(void) {{
    Eph eph;
    eph.toe     = {ephemeris_data['toe']};                 /* reference epoch, s */
    eph.A       = {ephemeris_data['A']};          /* semi-major axis, m (~MEO) */
    eph.e       = {ephemeris_data['e']};                /* eccentricity */
    eph.M0      = {ephemeris_data['M0']};                 /* mean anomaly at toe, rad */
    eph.delta_n = {ephemeris_data['delta_n']};                 /* mean motion correction, rad/s */
    eph.omega   = {ephemeris_data['omega']};                 /* argument of perigee, rad */
    eph.Omega_0 = {ephemeris_data['Omega_0']};                 /* RAAN at toe, rad */
    eph.i0      = {ephemeris_data['i0']};  /* inclination, rad (~55 deg) */
    eph.mu      = {ephemeris_data['mu']};      /* Earth gravitational constant */
    eph.Omega_e = {ephemeris_data['Omega_e']};     /* Earth rotation rate, rad/s */
    return eph;
}}

/**
 * Creates a sample fixed-point ephemeris structure.
 * Converts the same values used by make_sample_eph() into fixed-point representation.
 */
static Eph_fixed make_sample_eph_fixed(void) {{
    Eph_fixed eph_fixed_32;
    eph_fixed_32.toe_32_q11     = double_to_fixed32({ephemeris_data['toe']}, TIME_Q);
    eph_fixed_32.A_32_q5        = double_to_fixed32({ephemeris_data['A']}, DIST_Q);
    eph_fixed_32.e_32_q30       = double_to_fixed32({ephemeris_data['e']}, TRIG_Q);
    eph_fixed_32.M0_32_q28      = double_to_fixed32({ephemeris_data['M0']}, ANG_Q);
    eph_fixed_32.delta_n_32_q28 = double_to_fixed32({ephemeris_data['delta_n']}, ANG_Q);
    eph_fixed_32.omega_32_q28   = double_to_fixed32({ephemeris_data['omega']}, ANG_Q);
    eph_fixed_32.Omega_0_32_q28 = double_to_fixed32({ephemeris_data['Omega_0']}, ANG_Q);
    eph_fixed_32.i0_32_q28      = double_to_fixed32({ephemeris_data['i0']}, ANG_Q);
    eph_fixed_32.Omega_e_32_q31 = double_to_fixed32({ephemeris_data['Omega_e']}, RATE_Q);
    return eph_fixed_32;
}}"""

    # 8. 精准块替换算法 (跨行匹配)
    pattern = r"static\s+Eph\s+make_sample_eph\s*\(\s*void\s*\)\s*\{.*?return\s+eph_fixed_32\s*;\s*\}"
    modified_content, count = re.subn(pattern, new_functions, c_content, flags=re.DOTALL)

    if count == 0:
        print("【严重错误】: 脚本未能从 main.c 中识别到目标函数结构，未做任何修改！")
        return

    # 9. 写回文件
    with open(c_path, mode='w', encoding='utf-8') as f:
        f.write(modified_content)
        
    print(f"【成功】: 已精准定位并自动刷新了 {c_path} 内部的代码！")

import sys

if __name__ == "__main__":
    # 如果命令行传了参数，就用参数作为 target_name，否则默认使用 meo_representative
    target = sys.argv[1] if len(sys.argv) > 1 else "meo_representative"
    update_c_file(target_name=target)


