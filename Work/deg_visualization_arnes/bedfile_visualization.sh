#!/bin/bash


for file in *_DEG.bed
do
    name=${file%_DEG.bed}
    grep -v 'CHR' $file > ${name}_removed.bed
    new_file=${name}_removed.bed
    
    makeTagDirectory $name $new_file -format bed
    makeUCSCfile $name -o auto
    # makeBigWig.pl $name mm10 -webdir web/ -url web/ -update

done


