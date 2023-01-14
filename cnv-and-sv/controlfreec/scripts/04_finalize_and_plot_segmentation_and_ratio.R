library(data.table)
library(ggplot2)


##### Set up segmentation data

# Read in segmentation data and filtered bed file
seg <- read.table("intermediate_files/IMR5-merge_segmentation_withSegmentID.txt", stringsAsFactors = F, head=T, sep="\t")
seg_bed_filtered <- read.table("intermediate_files/IMR5-merge_segmentation_withSegmentID_rmBlacklistAndSegDup.bed", stringsAsFactors = F, head=F, sep="\t")
names(seg_bed_filtered) <- c("Chromosome", "Start", "End", "SegmentID")

# Remove segments that overlap 50% or more with exclude regions (keep only the segments listed in the filtered BED)
filter_name <- "rmBlacklistAndSegDup"
seg <- subset(seg, SegmentID %in% seg_bed_filtered$SegmentID)

# Remove segments with MedianRatio = -1 (this is ControlFREEC's annotation for NA)
seg <- subset(seg, MedianRatio != -1)

# Multiply MedianRatio in sex chromosomes by 2 
# (Not sure why it didn't normalize to the reference)
sex_chrom_logical <- seg$Chromosome %in% c("X", "Y")
seg[sex_chrom_logical, "MedianRatio"] <- seg[sex_chrom_logical, "MedianRatio"] * 2

# Annotate segments with gain or loss status based on new MedianRatio thresholds
upperThreshold <- 1.2
lowerThreshold <- 0.8
threshold_name <- "0.8and1.2thresholds"
seg$NewCNVcall <- "CopyNeutral"
seg[seg$MedianRatio >= upperThreshold, "NewCNVcall"] <- "Gain"
seg[seg$MedianRatio <= lowerThreshold, "NewCNVcall"] <- "Loss"

# Create new column with log2(MedianRatio)
seg$Log2MedianRatio <- log2(seg$MedianRatio)
seg[seg$Log2MedianRatio=="-Inf", "Log2MedianRatio"] <- -10
# Note there are a few segments with MedianRatio=0 that become -Inf; change these to -10


##### Write out segmentation data in IGV SEG format (same format used for svpluscnv) 
igv <- seg[,c("Sample", "Chromosome", "Start", "End", "BinCount", "Log2MedianRatio")]
colnames(igv) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
write.table(igv, paste0("results/IMR5-merge_segmentation_", filter_name, ".seg"), quote=F,
            col.names = T, row.names = F, sep="\t")


##### Set up ratio data 

# Read in ratio data (~600 MB)
ratio <- fread("results/IMR5-merge_ratio.txt", stringsAsFactors = F, head=T, sep="\t")

# Remove bins with Ratio or MedianRatio = -1 (this is ControlFREEC's annotation for NA)
ratio <- subset(ratio, (Ratio != -1 & MedianRatio != -1))

# Multiply Ratio and MedianRatio in sex chromosomes by 2
sex_chrom_logical <- ratio$Chromosome %in% c("X", "Y")
ratio[sex_chrom_logical, "Ratio"] <- ratio[sex_chrom_logical, "Ratio"] * 2
ratio[sex_chrom_logical, "MedianRatio"] <- ratio[sex_chrom_logical, "MedianRatio"] * 2

# Annotate bins with gain or loss status based on MedianRatio thresholds (defined above) 
ratio$NewCNVcall <- "CopyNeutral"
ratio[ratio$MedianRatio >= upperThreshold, "NewCNVcall"] <- "Gain"
ratio[ratio$MedianRatio <= lowerThreshold, "NewCNVcall"] <- "Loss"


##### Filter seg and ratio data to remove segments with <5 bins
# Main purpose is so that this plot doesn't show extra regions that I've removed from my other analyses/plots
# Approach: Filter the segment data, then filter ratio data based on the filtered segment data
# (This will also filter blacklist regions out of the ratio data, since these were already filtered from seg)
## Filter seg
seg_5bins <- subset(seg, BinCount>=5)
## Create logical vector indicating overlap of ratio points with segments 
overlap_vector <- apply(ratio, 1, function(x){
  ratio_chr <- x[1]
  ratio_start <- as.integer(x[2])
  keep_logical <- any(ratio_chr==seg_5bins$Chromosome & ratio_start>=seg_5bins$Start & ratio_start<=seg_5bins$End)
  return(keep_logical)
})
## Filter ratio
ratio_5bins <- ratio[overlap_vector,]


##### Plot ratio and segmentation across the genome

# Define plot theme and options
genomePlot_theme <-
  theme_light() +
  theme(strip.background = element_rect(fill = NA),
        strip.text.x = element_text(colour = "black"),
        strip.text.y = element_text(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(colour = "gray80", fill=NA), # overwrite darker lines with lighter lines
        legend.position = "none")

# Define sample labels
sample_labels <- c("Control clone", "R112*", "R150*", "E287fs")
names(sample_labels) <- c("IMR5-WT", "IMR5-R112", "IMR5-R150", "IMR5-E287")

# Reorder chromosomes and samples for plotting (for both data frames)
ratio_5bins$Chromosome <- factor(ratio_5bins$Chromosome, c(1:22, "X", "Y"))
ratio_5bins$Sample <- factor(ratio_5bins$Sample, c("IMR5-WT", "IMR5-R112", "IMR5-R150", "IMR5-E287"))
seg_5bins$Chromosome <- factor(seg_5bins$Chromosome, c(1:22, "X", "Y"))
seg_5bins$Sample <- factor(seg_5bins$Sample, c("IMR5-WT", "IMR5-R112", "IMR5-R150", "IMR5-E287"))

# Cap ratio and seg values at a max of 2 so they aren't cut off in plots
seg_5bins[seg_5bins$MedianRatio > 2, "MedianRatio"] <- 2
ratio_5bins[ratio_5bins$Ratio > 2, "Ratio"] <- 2

# Plot ratio for each bin and segmentation across genome
p <- ggplot(data=ratio_5bins, aes(x=Start, y=Ratio, color=NewCNVcall)) +
  facet_grid(Sample~Chromosome, scales="free_x", space="free_x", labeller=labeller(Sample=sample_labels)) +
  geom_point(alpha=0.01, shape=".") +
  geom_segment(data=seg_5bins, aes(x=Start, y=MedianRatio, xend=End, yend=MedianRatio), color="black", size=0.5) +
  geom_hline(yintercept=upperThreshold, color="red", size=0.1) +
  geom_hline(yintercept=lowerThreshold, color="blue", size=0.1) +
  scale_color_manual(values=c("gray", "red", "blue")) +
  genomePlot_theme +
  ylab("Copy Number Ratio") +
  ylim(c(0,2))
ggsave(paste0("plots/ratio_genomePlot_", threshold_name, "_", filter_name, "_forPublication.png"), p, width=10, height=4)
