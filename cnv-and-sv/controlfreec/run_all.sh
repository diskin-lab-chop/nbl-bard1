#!/bin/sh

# Run all scripts (needs to run from main analysis directory)

python scripts/01_merge_ratio_files.py
python scripts/02_generate_segmentation_from_ratio.py
bash scripts/03_remove_blacklist_regions_from_segmentation.sh
Rscript scripts/04_finalize_and_plot_segmentation_and_ratio.R

