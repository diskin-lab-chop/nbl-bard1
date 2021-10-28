#!/bin/sh

# Adapted from https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/copy_number_consensus_call

# # hg38
# curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
# mv genomicSuperDups.txt.gz genomicSuperDups_hg38.txt.gz
# gunzip -c genomicSuperDups_hg38.txt.gz \
# | awk -v OFS='\t' '{if($27>0.95){print $2, $3, $4, $27}}' \
# | sort -k1,1 -k2,2n \
# | bedtools merge > segmental_dups_hg38.bed

# hg19
curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
mv genomicSuperDups.txt.gz genomicSuperDups_hg19.txt.gz
# gunzip -c genomicSuperDups_hg19.txt.gz | head > genomicSuperDups_hg19.head.txt
gunzip -c genomicSuperDups_hg19.txt.gz \
| awk -v OFS='\t' '{if($27>0.95){print $2, $3, $4, $27}}' \
| sort -k1,1 -k2,2n \
| bedtools merge > segmental_dups_hg19.bed
