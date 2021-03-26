#!/bin/bash

#  filters the output of the ensembl variant effect predictor, keeps only missense variants 

case "$1" in
*.gz ) 
        zcat $1 | awk 'BEGIN{OFS="\t"}{if($7 ~ /missense/){print $1, $2, $3, $4, $5, $7}else{if($1 ~ /##/){print $0}else{if($1 ~ /#/){print $1, $2, $3, $4, $5, $7}}}}' | grep -v "NMD_transcript_variant"
        ;;
*)
        awk 'BEGIN{OFS="\t"}{if($7 ~ /missense/){print $1, $2, $3, $4, $5, $7}else{if($1 ~ /##/){print $0}else{if($1 ~ /#/){print $1, $2, $3, $4, $5, $7}}}}' $1 | grep -v "NMD_transcript_variant"
        ;;
esac