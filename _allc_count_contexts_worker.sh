#!/bin/bash 

# autosomes and X chromosome 

input=$1
output=$2

echo "$input -> $output"
zcat $input | awk 'BEGIN{FS="\t";OFS="\t"} \
    ($1!="chrY" && $1!="chrM" && $1!="chrL") \
    {a[$4]+=$5; b[$4]+=$6} \
    END{for(i in a){print i,a[i],b[i]}}' - > $output
gzip $output