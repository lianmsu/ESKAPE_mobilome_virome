
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 01_resistome_annotation
# ARG annotation: RGI/CARD and diamond against SARG
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:31ARG/3diamond_res_pretreatment.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

#
cd ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/42ARG/1ARG
awk 'BEGIN{OFS="\t"} {
    # 从第一列提取 "bin.数字" 之前的所有内容
    extracted = "NA"  # 默认值
    if (match($1, /.*bin\.[0-9]+/)) {
        extracted = substr($1, 1, RLENGTH)
    }
    # 在原行末尾追加提取的内容
    print $0, extracted
}' diamond_ARG_out > diamond_ARG_out2
awk 'BEGIN{OFS="\t"} {print $0, "ARG"}' diamond_ARG_out2 > diamond_ARG_out3
python3 ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/42ARG/4turn_to_table.py diamond_ARG_out3

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG/2VFG
awk 'BEGIN{OFS="\t"} {
    # 提取第一个 "-" 之前的内容，如果没有 "-" 则保留原值
    extracted = $1  # 默认取整个第一列
    idx = index($1, "-")
    if (idx > 0) {
        extracted = substr($1, 1, idx-1)
    }
    # 在原行末尾追加提取的内容
    print $0, extracted
}' diamond_VFG_out > diamond_VFG_out2
awk 'BEGIN{OFS="\t"} {print $0, "VFG"}' diamond_VFG_out2 > diamond_VFG_out3
python3 ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG/4turn_to_table.py diamond_VFG_out3

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG/3MRG
awk 'BEGIN{OFS="\t"} {
    # 提取第一个 "-" 之前的内容，如果没有 "-" 则保留原值
    extracted = $1  # 默认取整个第一列
    idx = index($1, "-")
    if (idx > 0) {
        extracted = substr($1, 1, idx-1)
    }
    # 在原行末尾追加提取的内容
    print $0, extracted
}' diamond_MRG_out > diamond_MRG_out2
awk 'BEGIN{OFS="\t"} {print $0, "MRG"}' diamond_MRG_out2 > diamond_MRG_out3
python3 ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG/4turn_to_table.py diamond_MRG_out3

#################################################################
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG/0ARG_RGI
cut -f1,15,16,17,28 rgiargs.txt > rgiargs_2.txt

# 第一步：清理文件（去除^M）
sed "s/\r$//" rgiargs_2.txt > rgiargs_2_clean.txt
awk -F"\t" 'BEGIN{OFS="\t"} {
    # 提取第一个 "-" 之前的内容，如果没有 "-" 则保留原值
    extracted = $1  # 默认取整个第一列
    idx = index($1, "-")
    if (idx > 0) {
        extracted = substr($1, 1, idx-1)
    }
    # 在原行末尾追加提取的内容
    print $0, extracted
}' rgiargs_2_clean.txt > rgiargs_3.txt
awk 'BEGIN{OFS="\t"} {print $0, "ARG"}' rgiargs_3.txt > rgiargs_4.txt
python3 ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG/4turn_to_table.py rgiargs_4.txt

