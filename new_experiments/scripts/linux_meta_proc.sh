#!/bin/bash

## to pull all blood culture samples out of the barcodes file 
grep -w "culture" ../data/bloodCulture_qiime/S78_barcodes.txt > ../data/bloodCulture_qiime/culture_barcodes.txt

## to pull the column names off of the main barcode file
head -n 1 ../data/bloodCulture_qiime/raw_seqs/S78_barcodes.txt > ../data/bloodCulture_qiime/columns.txt

## and to replace 'SampleID' column with '#SampleID' with sed so qiime doesn't get pissed off
## s/ deliminates the substitution
## /SampleID is what it's currently called
## /#SampleID is what you want to change it to
## /g stands for global replacement and it changes all occurences in the file 
sed 's/SampleID/#SampleID/g' ../data/bloodCulture_qiime/columns.txt > ../data/bloodCulture_qiime/updated_cols.txt

## putting the barcodes and updated column names together
cat ../data/bloodCulture_qiime/updated_cols.txt ../data/bloodCulture_qiime/culture_barcodes.txt > ../data/bloodCulture_qiime/bloodCulture_barcodes.txt



