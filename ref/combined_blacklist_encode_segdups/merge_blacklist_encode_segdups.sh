#!/bin/sh

module load bedtools
# bedtools merge requires that you presort your data by chromosome and then by start position (e.g., sort -k1,1 -k2,2n in.bed > in.sorted.bed for BED files).
cat  ../encode_blacklist/hg19-blacklist.v2.bed ../segmental_duplications/segmental_dups_hg19.bed | cut -f 1-3 | sort -k1,1 -k2,2n | bedtools merge > combined_blacklist_encode_segdups_hg19.bed







