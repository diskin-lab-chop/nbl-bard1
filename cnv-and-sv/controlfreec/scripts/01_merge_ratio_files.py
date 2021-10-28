### Merge ratio files from each sample

import glob

outfile = open("results/IMR5-merge_ratio.txt", "w")
path = "output_from_cavatica/IMR5-*.downsampled.bam_ratio.txt"

first_file = True
for filename in glob.glob(path):
	file = open(filename, "r")
	sample = filename.replace("output_from_cavatica/", "").replace(".downsampled.bam_ratio.txt", "")
	
	if first_file==True:
		header = file.readline().rstrip()
		header = header + "\tSample\n"
		outfile.write(header)
	else:
		next(file)
	
	for line in file:
		line = line.rstrip() + "\t" + sample + "\n"
		outfile.write(line)

	file.close()
	first_file=False

outfile.close()