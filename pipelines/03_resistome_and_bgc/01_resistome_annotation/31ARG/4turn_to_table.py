#!/usr/bin/python3
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 01_resistome_annotation
# ARG annotation: RGI/CARD and diamond against SARG
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:31ARG/4turn_to_table.py
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

import sys
import pandas as pd
from collections import defaultdict

# 检查参数是否正确
if len(sys.argv) != 2:
    print("Usage: python script.py <input_file>")
    sys.exit(1)

input_file = sys.argv[1]

# 从输入文件名中提取文件名（不包括路径和扩展名）
input_filename = input_file.split('/')[-1].split('.')[0]


# 使用defaultdict来存储数据，便于处理不存在的项
data = defaultdict(lambda: defaultdict(int))

# 读取文件并处理每一行
with open(input_file, 'r') as f:
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) >= 2:
            # 提取MAG名字所在列和基因注释所在列作为行名和列名
            row_name, col_name = cols[-2], cols[-1]
            data[row_name][col_name] += 1  # 填充数字，这里简单计数

# 转换为pandas DataFrame进行进一步处理
df = pd.DataFrame(data).fillna(0).astype(int)
df = df.T  # 转置，使第一列字符串去冗余后作为列名，第二列字符串作为行名

# 输出结果到CSV文件，文件名根据输入文件名动态生成
output_file = input_filename + '_type.csv'
df.to_csv(output_file, sep='\t')
print(f"Output saved to {output_file}")
