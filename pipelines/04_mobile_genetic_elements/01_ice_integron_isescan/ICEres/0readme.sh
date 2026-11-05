
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 01_ice_integron_isescan
# ICE detection with ICEfinder2 (batches of 1000 and the full set)
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:5ICEfinder/ICEres/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

ICE_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder

#从下面文件夹中及其子文件夹中提取所有以.fa结尾的文件到${HOME_ROOT:-/path/to/home}/03GW_prokaryotes/11ICE_res/allICE_sequence_of19847MAG
${ICE_path}/ICEfinder2_1000/result/
${ICE_path}/ICEfinder2_2000/result/
${ICE_path}/ICEfinder2_3000/result/
${ICE_path}/ICEfinder2_4000/result/
${ICE_path}/ICEfinder2_5000/result/
${ICE_path}/ICEfinder2_5455/result/

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres
find "${ICE_path}/ICEfinder2_1000/result/" -type f -name "*.fa" -exec cp {} ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence \;

find "${ICE_path}/ICEfinder2_2000/result/" -type f -name "*.fa" -exec cp {} ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence \;

find "${ICE_path}/ICEfinder2_3000/result/" -type f -name "*.fa" -exec cp {} ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence \;

find "${ICE_path}/ICEfinder2_4538/result/" -type f -name "*.fa" -exec cp {} ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence \;



#对文件夹${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEres/all5455MAG_ICE_sequence下的.fa文件中的序列名进行更改。把文件名（不包括.fa后缀）加入到原本的序列名之前。

for file in ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence/*.fa; do
    # 获取不带路径和后缀的文件名
    filename=$(basename "$file" .fa)
    # 使用 sed 把序列名中的... 替换为_
    sed -i "s/\.\.\. /_ICE_/g" "$file"
done

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence/*.fa > all_combined_ICEs.fa

grep '>' ./all_combined_ICEs.fa | wc -l
#4799个ICE/IME

defense-finder run ./all_combined_ICEs.fa --preserve-raw -o ${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEres/defensefinder2.0.0

#606个ICE/IME中共发现84个DSs（210个DGs）

#所有ICE的序列名
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres
grep '>' ./all_combined_ICEs.fa > all_ICE_contigid.txt

# 提取信息并计算长度
awk -F"[>_ -]" '{
    contig = $2 "_" $3 "_" $4
    start = $(NF-1)
    end = $NF
    len = end - start + 1
    print contig "\t" start "\t" end "\t" len
}' ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all_ICE_contigid.txt > ICE_contig_info.tsv

# 使用sed在第一行前插入标题
sed -i "1icontig\tstart\tend\tice_length" ICE_contig_info.tsv

#################################################################################
#对文件夹${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEres/all5455MAG_ICE_sequence下的.fa文件中的序列名进行更改。把文件名（不包括.fa后缀）加入到原本的序列名之前。

for file in ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence/*.fa; do
    # 获取不带路径和后缀的文件名
    filename=$(basename "$file" .fa)
    # 使用 sed 在序列名前添加文件名（格式：>文件名_原序列名）
    sed -i "s/^>/>${filename}_/" "$file"
done

for file in ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence/*.fa; do
    # 获取不带路径和后缀的文件名
    filename=$(basename "$file" .fa)
    # 使用 sed 把序列名中的... 替换为_
    sed -i "s/\.\.\. /_/g" "$file"
done

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres
cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all4538genome_ICE_sequence/*.fa > all_combined_ICEs.fa

grep '>' ./all_combined_ICEs.fa | wc -l
#4799个ICE/IME


cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres
defense-finder run ./all_combined_ICEs.fa --preserve-raw -o ${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEres/defensefinder2.0.0

#606个ICE/IME中共发现84个DSs（210个DGs）

#所有ICE的序列名
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres
grep '>' ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/all_combined_ICEs.fa > all_ICE_contigid.txt

awk -F'_ICE[0-9]_' '{
    split($1, a, "_contig_"); 
    contig_id = a[1]"_contig_"a[2];
    range = $2; 
    sub(/^[^-]*-/, "", range);  # 删除开头的非"-"部分（如CJ_14-18JJ_SAB-_）
    split(range, coords, "-");
    print contig_id, coords[1], coords[2], coords[2]-coords[1]+1
}' all_ICE_contigid.txt

#后续用${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEres/process_iceid.py做进一步处理，得到${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEres/processed_ice_sequences_ID.txt

#把文件${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEres/processed_ice_sequences_ID.txt中contig_替换为contig-
sed -i 's/_renamed_contig_/-/g' ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEres/processed_ice_sequences_ID.txt
