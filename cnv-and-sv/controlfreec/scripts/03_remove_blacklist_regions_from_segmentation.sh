#!/bin/sh

# Add unique segment ID 
awk -F "\t" '{if (NR==1) print $0 "\tSegmentID"; else print $0 "\tSegment" NR-1}' intermediate_files/IMR5-merge_segmentation.txt > intermediate_files/IMR5-merge_segmentation_withSegmentID.txt

# Reformat to BED file; drop all info other than position and segment ID
tail -n+2 intermediate_files/IMR5-merge_segmentation_withSegmentID.txt | awk -F "\t" '{print "chr"$1 "\t" $2 "\t" $3 "\t" $8}' - > intermediate_files/IMR5-merge_segmentation_withSegmentID.bed

# Filter out segments that overlap 50% or more with exclude list (ENCODE blacklist merged w/ segmental duplications list)
module load BEDTools/2.29.2-GCC-9.3.0
bedtools subtract \
	-A \
	-a intermediate_files/IMR5-merge_segmentation_withSegmentID.bed \
	-b ../../ref/combined_blacklist_encode_segdups/combined_blacklist_encode_segdups_hg19.bed \
	-f 0.5 > intermediate_files/IMR5-merge_segmentation_withSegmentID_rmBlacklistAndSegDup.bed
