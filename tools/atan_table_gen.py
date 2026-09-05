import math
import os

def generate_cordic_h_file(filename="include/cordic_table.h", iterations=28):
    # 确保 include 文件夹存在
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("// 本文件由 Python 脚本自动生成，请勿手动修改\n")
        f.write("#pragma once\n")
        f.write("#ifndef CORDIC_TABLE_H\n")
        f.write("#define CORDIC_TABLE_H\n\n")
        f.write("#include \"common_types.h\"\n\n")
        
        # 使用 static const，这样多个 .c 文件 #include 它时才不会报重复定义错误
        f.write("static const fixed32_t CORDIC_TABLE_Q28[28] = {\n")
        
        for i in range(iterations):
            val = round(math.atan(2 ** -i) * (2 ** 28))
            hex_str = f"0x{val:08X}"
            
            if i == iterations - 1:
                f.write(f"    {hex_str}   // i = {i:2d} (Dec: {val})\n")
            else:
                f.write(f"    {hex_str},  // i = {i:2d} (Dec: {val})\n")
                
        f.write("};\n\n")
        f.write("#endif // CORDIC_TABLE_H\n")
    print(f"成功生成文件: {filename}")

if __name__ == "__main__":
    generate_cordic_h_file()
