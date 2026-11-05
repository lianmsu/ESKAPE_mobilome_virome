# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 01_mge_crispr_spacers
# Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets)
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:35MGE_spacer_competition/3spacer2db/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/3spacer2db

spacer_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/1MGE_spacer
db_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/2blastn_db
#spacer of plasmid to ICEdb
blastn -db "${db_path}/ICE.db" \
       -query "${spacer_path}/Plasmid_spacers.fa" \
       -out "./Plasmid_spacers2ICE.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Plasmid_spacers2ICE.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Plasmid_spacers2ICE.blastn.1mis95cov.out
cut -f1 Plasmid_spacers2ICE.blastn.1mis95cov.out | sort | uniq | wc -l
#315

cat ./Plasmid_spacers2ICE.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Plasmid_spacers2ICE.blastn.1mis100cov.out
cut -f1 Plasmid_spacers2ICE.blastn.1mis100cov.out | sort | uniq | wc -l
#301

#spacer of plasmid to Prophagedb
blastn -db "${db_path}/Prophage.db" \
       -query "${spacer_path}/Plasmid_spacers.fa" \
       -out "./Plasmid_spacers2Prophage.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Plasmid_spacers2Prophage.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Plasmid_spacers2Prophage.blastn.1mis95cov.out
cut -f1 Plasmid_spacers2Prophage.blastn.1mis95cov.out | sort | uniq | wc -l
#992

cat ./Plasmid_spacers2Prophage.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Plasmid_spacers2Prophage.blastn.1mis100cov.out
cut -f1 Plasmid_spacers2Prophage.blastn.1mis100cov.out | sort | uniq | wc -l
#939

#spacer of plasmid to Integrondb
blastn -db "${db_path}/Integron.db" \
       -query "${spacer_path}/Plasmid_spacers.fa" \
       -out "./Plasmid_spacers2Integron.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Plasmid_spacers2Integron.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Plasmid_spacers2Integron.blastn.1mis95cov.out
cut -f1 Plasmid_spacers2Integron.blastn.1mis95cov.out | sort | uniq | wc -l
#784

cat ./Plasmid_spacers2Integron.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Plasmid_spacers2Integron.blastn.1mis100cov.out
cut -f1 Plasmid_spacers2Integron.blastn.1mis100cov.out | sort | uniq | wc -l
#770

#####################################################
#spacer of Prophage to ICEdb
blastn -db "${db_path}/ICE.db" \
       -query "${spacer_path}/Prophage_spacers.fa" \
       -out "./Prophage_spacers2ICE.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Prophage_spacers2ICE.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Prophage_spacers2ICE.blastn.1mis95cov.out
cut -f1 Prophage_spacers2ICE.blastn.1mis95cov.out | sort | uniq | wc -l
#60

cat ./Prophage_spacers2ICE.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Prophage_spacers2ICE.blastn.1mis100cov.out
cut -f1 Prophage_spacers2ICE.blastn.1mis100cov.out | sort | uniq | wc -l
#58

#spacer of Prophage to Plasmiddb
blastn -db "${db_path}/Plasmid.db" \
       -query "${spacer_path}/Prophage_spacers.fa" \
       -out "./Prophage_spacers2Plasmid.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Prophage_spacers2Plasmid.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Prophage_spacers2Plasmid.blastn.1mis95cov.out
cut -f1 Prophage_spacers2Plasmid.blastn.1mis95cov.out | sort | uniq | wc -l
#263

cat ./Prophage_spacers2Plasmid.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Prophage_spacers2Plasmid.blastn.1mis100cov.out
cut -f1 Prophage_spacers2Plasmid.blastn.1mis100cov.out | sort | uniq | wc -l
#259

#spacer of Prophage to Integrondb
blastn -db "${db_path}/Integron.db" \
       -query "${spacer_path}/Prophage_spacers.fa" \
       -out "./Prophage_spacers2Integron.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Prophage_spacers2Integron.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Prophage_spacers2Integron.blastn.1mis95cov.out
cut -f1 Prophage_spacers2Integron.blastn.1mis95cov.out | sort | uniq | wc -l
#139

cat ./Prophage_spacers2Integron.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Prophage_spacers2Integron.blastn.1mis100cov.out
cut -f1 Prophage_spacers2Integron.blastn.1mis100cov.out | sort | uniq | wc -l
#138

####################################################
#spacer of ICE to Integrondb
blastn -db "${db_path}/Integron.db" \
       -query "${spacer_path}/ICE_spacers.fa" \
       -out "./ICE_spacers2Integron.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./ICE_spacers2Integron.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > ICE_spacers2Integron.blastn.1mis95cov.out
cut -f1 ICE_spacers2Integron.blastn.1mis95cov.out | sort | uniq | wc -l
#32

cat ./ICE_spacers2Integron.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > ICE_spacers2Integron.blastn.1mis100cov.out
cut -f1 ICE_spacers2Integron.blastn.1mis100cov.out | sort | uniq | wc -l
#31

#spacer of ICE to Plasmiddb
blastn -db "${db_path}/Plasmid.db" \
       -query "${spacer_path}/ICE_spacers.fa" \
       -out "./ICE_spacers2Plasmid.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./ICE_spacers2Plasmid.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > ICE_spacers2Plasmid.blastn.1mis95cov.out
cut -f1 ICE_spacers2Plasmid.blastn.1mis95cov.out | sort | uniq | wc -l
#17

cat ./ICE_spacers2Plasmid.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > ICE_spacers2Plasmid.blastn.1mis100cov.out
cut -f1 ICE_spacers2Plasmid.blastn.1mis100cov.out | sort | uniq | wc -l
#16

#spacer of ICE to Prophagedb
blastn -db "${db_path}/Prophage.db" \
       -query "${spacer_path}/ICE_spacers.fa" \
       -out "./ICE_spacers2Prophage.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./ICE_spacers2Prophage.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > ICE_spacers2Prophage.blastn.1mis95cov.out
cut -f1 ICE_spacers2Prophage.blastn.1mis95cov.out | sort | uniq | wc -l
#44

cat ./ICE_spacers2Prophage.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > ICE_spacers2Prophage.blastn.1mis100cov.out
cut -f1 ICE_spacers2Prophage.blastn.1mis100cov.out | sort | uniq | wc -l
#43

###########################################################################
#spacer of Integron to Prophagedb
blastn -db "${db_path}/Prophage.db" \
       -query "${spacer_path}/Integron_spacers.fa" \
       -out "./Integron_spacers2Prophage.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Integron_spacers2Prophage.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Integron_spacers2Prophage.blastn.1mis95cov.out
cut -f1 Integron_spacers2Prophage.blastn.1mis95cov.out | sort | uniq | wc -l
#310

cat ./Integron_spacers2Prophage.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Integron_spacers2Prophage.blastn.1mis100cov.out
cut -f1 Integron_spacers2Prophage.blastn.1mis100cov.out | sort | uniq | wc -l
#309

#spacer of Integron to Plasmiddb
blastn -db "${db_path}/Plasmid.db" \
       -query "${spacer_path}/Integron_spacers.fa" \
       -out "./Integron_spacers2Plasmid.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Integron_spacers2Plasmid.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Integron_spacers2Plasmid.blastn.1mis95cov.out
cut -f1 Integron_spacers2Plasmid.blastn.1mis95cov.out | sort | uniq | wc -l
#81

cat ./Integron_spacers2Plasmid.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Integron_spacers2Plasmid.blastn.1mis100cov.out
cut -f1 Integron_spacers2Plasmid.blastn.1mis100cov.out | sort | uniq | wc -l
#78

#spacer of Integron to ICEdb
blastn -db "${db_path}/ICE.db" \
       -query "${spacer_path}/Integron_spacers.fa" \
       -out "./Integron_spacers2ICE.blastn.out" \
       -outfmt "6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send" \
       -num_threads 12 \
       -mt_mode 1 \
       -evalue 1e-5 \
       -word_size 8 \
       -task blastn-short

cat ./Integron_spacers2ICE.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > Integron_spacers2ICE.blastn.1mis95cov.out
cut -f1 Integron_spacers2ICE.blastn.1mis95cov.out | sort | uniq | wc -l
#85

cat ./Integron_spacers2ICE.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > Integron_spacers2ICE.blastn.1mis100cov.out
cut -f1 Integron_spacers2ICE.blastn.1mis100cov.out | sort | uniq | wc -l
#79

