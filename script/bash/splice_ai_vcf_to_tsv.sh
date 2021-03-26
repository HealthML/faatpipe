#!/bin/bash


if [ $# -ne 2 ]; then
    echo 'need 2 input arguments: VCF and chromosome!'
    exit 1
fi

# this script will only work properly for SNVs !
# the first awk one-liner changes the format and column ordering
# the second awk one-liner filters on the variant effect prediction 

tabix $1 "${2##chr}" | tr '|' '\t' | awk -F "\t" 'BEGIN{OFS = "\t"}{id = $1":"$2":"$4":"$5; start=$2-1; end=$2;  print $1, start, end, id, $9, ".", $10, $11, $12, $13, $14, $15, $16, $17 }' | awk -F "\t" 'BEGIN{OFS = "\t"; print "chrom", "start", "end", "name", "gene", "strand", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL" }{if($7 >= 0.1 || $8 >= 0.1 || $9 >= 0.1 || $10 >= 0.1){print $0}}'