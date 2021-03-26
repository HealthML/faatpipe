#!/bin/bash


case "$1" in
*.gz ) 
        zcat $1 | awk 'BEGIN{OFS="\t"}{if($NF ~ /^IMPACT=HIGH/){print $1, $2, $3, $4, $5, $7}else{if($1 ~ /##/){print $0}else{if($1 ~ /#/){print $1, $2, $3, $4, $5, $7}}}}'
        ;;
*)
        awk 'BEGIN{OFS="\t"}{if($NF ~ /^IMPACT=HIGH/){print $1, $2, $3, $4, $5, $7}else{if($1 ~ /##/){print $0}else{if($1 ~ /#/){print $1, $2, $3, $4, $5, $7}}}}' $1 
        ;;
esac