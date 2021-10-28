### Create segmentation file from ratio file

# Ignore the BAF and Subclone info because BAF calculation did not run properly

# Loop through the file bin by bin. For the current bin:
	# Check whether this bin is a continuation of the same segment or if it's the beginning of a new segment
		# (Do the Chromosome and MedianRatio for this bin match the segment?)
	# If it's part of the same segment:
		# Update segment_End and segment_BinCount
	# If it's the start of a new segment:
		# First write out the previous segment's info 
		# Then re-define the segment info

infile = open("results/IMR5-merge_ratio.txt", "r")
outfile = open("intermediate_files/IMR5-merge_segmentation.txt", "w")

next(infile)	# Skip header

segment_Chrom = None
segment_Start = None
segment_End = None
segment_MedianRatio = None
segment_CopyNumber = None
segment_Sample = None
segment_BinCount = 0

binSize = 1123	# For this particular dataset, the bins are either 1123 or 1124 bases long

first_bin = True

for line in infile:

	# Extract current bin info
	lineSplit = line.rstrip().split("\t")
	bin_Chrom = lineSplit[0]
	bin_Start = lineSplit[1]
	bin_MedianRatio = lineSplit[3]
	bin_CopyNumber = lineSplit[4]
	bin_Sample = lineSplit[11]

	# Is this bin is a continuation of the same segment?
	if bin_Chrom == segment_Chrom and bin_MedianRatio == segment_MedianRatio: 
		segment_End = int(bin_Start) + binSize
		segment_BinCount = segment_BinCount + 1

	# If not, this bin must be the beginning of a new segment
	else: 
		
		# Write out info for the old segment
		# (unless this is the very first bin, then just write out the file header and move on to defining the new segment)
		if first_bin == True:
			outfile.write("Chromosome\tStart\tEnd\tMedianRatio\tCopyNumber\tSample\tBinCount\n")
			first_bin = False
		else:
			outfile.write("\t".join([segment_Chrom, segment_Start, str(segment_End), segment_MedianRatio, 
				segment_CopyNumber, segment_Sample, str(segment_BinCount)]) + "\n")

		# Define info for the new segment
		segment_Chrom = bin_Chrom
		segment_Start = bin_Start
		segment_End = int(bin_Start) + binSize
		segment_MedianRatio = bin_MedianRatio
		segment_CopyNumber = bin_CopyNumber
		segment_Sample = bin_Sample
		segment_BinCount = 1

# Write out info for the final segment after the last bin is read
outfile.write("\t".join([segment_Chrom, segment_Start, str(segment_End), segment_MedianRatio, 
	segment_CopyNumber, segment_Sample, str(segment_BinCount)]) + "\n")

infile.close()
outfile.close()
