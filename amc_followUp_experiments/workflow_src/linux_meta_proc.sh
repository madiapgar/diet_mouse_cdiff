#!/bin/bash

## to pull all blood culture samples out of the barcodes file 
grep -w "stool" ./data/SEQ096/S96_barcodes.txt > ./data/SEQ096/stool_barcodes.txt

## to pull the column names off of the main barcode file
head -n 1 ./data/SEQ096/S96_barcodes.txt > ./data/SEQ096/columns.txt

## and to replace 'SampleID' column with '#SampleID' with sed so qiime doesn't get pissed off
## s/ deliminates the substitution
## /SampleID is what it's currently called
## /#SampleID is what you want to change it to
## /g stands for global replacement and it changes all occurences in the file 
sed 's/SampleID/#SampleID/g' ./data/SEQ096/columns.txt > ./data/SEQ096/updated_cols.txt

## putting the barcodes and updated column names together
cat ./data/SEQ096/updated_cols.txt ./data/SEQ096/stool_barcodes.txt > ./data/SEQ096/new_stool_barcodes.txt



