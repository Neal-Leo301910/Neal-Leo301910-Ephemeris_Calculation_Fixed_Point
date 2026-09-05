CC = gcc
# 加上 -fcommon 彻底杜绝全局符号合并可能引发的重定义报错
CFLAGS = -Wall -Wextra -O2 -Iinclude -fcommon

# 自动获取当前 src 下的所有 .c 文件（此时包含 src/main.c）
SRC = $(wildcard src/*.c)

# 定义星历数据路径,路径完全匹配同级父目录结构
CSV_DATA      = test/test_data/sample_ephemeris.csv
UPDATE_SCRIPT = tools/sample_ephemeris_generate.py

# ==========================================
# 定义默认的数据条目变量，后续添加后，可修改，比如改为 gps_prn27_20210624
# ==========================================
DATA_TARGET ?= meo_representative

all: .check_csv build/test_fixed_math build/test_ephemeris build/test_cordic build/test_ephemeris_fixed build/main

# 1. 强制数据检查规则（使用一个隐藏文件记录上一次运行的状态）
.check_csv: $(CSV_DATA) $(UPDATE_SCRIPT)
	@echo "正在准备星历数据，当前选择的条目: $(DATA_TARGET) ..."
	python3 $(UPDATE_SCRIPT) $(DATA_TARGET)
	@touch src/main.c
	@echo "数据已刷新，正在清理旧的编译缓存以确保彻底重编..."
	@rm -rf build/
	@echo "$(DATA_TARGET)" > .check_csv



# 2. 指定头文件的生成规则
include/cordic_table.h: tools/atan_table_gen.py
	@mkdir -p include
	@echo "正在通过 Python 脚本生成 CORDIC 查找表..."
	python3 tools/atan_table_gen.py

# 3. 编译主程序
build/main: include/cordic_table.h $(SRC)
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRC) -lm -o $@

# 4. 编译测试程序
build/test_fixed_math: include/cordic_table.h $(SRC) test/test_fixed_math.c
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRC) test/test_fixed_math.c -lm -o $@

build/test_cordic: include/cordic_table.h $(SRC) test/test_cordic.c
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRC) test/test_cordic.c -lm -o $@

build/test_ephemeris: include/cordic_table.h $(SRC) test/test_ephemeris_double.c
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRC) test/test_ephemeris_double.c -lm -o $@

build/test_ephemeris_fixed: include/cordic_table.h $(SRC) test/test_ephemeris_fixed.c
	@mkdir -p build
	$(CC) $(CFLAGS) $(SRC) test/test_ephemeris_fixed.c -lm -o $@

run: all
	./build/test_fixed_math
	./build/test_ephemeris
	./build/test_cordic
	./build/test_ephemeris_fixed
	./build/main

clean:
	rm -rf build
	rm -f include/cordic_table.h .check_csv

.PHONY: all run clean
