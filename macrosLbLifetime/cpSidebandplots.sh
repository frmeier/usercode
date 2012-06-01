#!/bin/bash

file=filelist.txt
dest=/Users/frankmeier/Documents/psi/lambda/anLambdaBlifetime/notes/AN-12-044/trunk/img_anchain_datamcComparison

ls *_sidebandsubtracted*pdf > $file
#sed -i "" 's/doAnalysis01plots.*\_B0/B0/' $file
#sed -i "" 's/doAnalysis01plots.*\_Lb/Lb/' $file
awk -v dest=$dest '{ f=$1; sub(/doAnalysis01plots[0-9]*_/, "", f); print $1, dest"/"f}' $file > filelist2.txt

cat filelist2.txt | while read LINE
do
    cp -v $LINE
done

