#!/bin/bash

listDirectories=("30" "33" "36" "39" "42" "45" "48" "51" "54" "57" "60")
listFiles=(1067 1068 2694)

for dir in "${listDirectories[@]}"; do
    for file in "${listFiles[@]}"; do
        fileMinus=`expr $file - 1`
        cp "data/c_exp_${dir}/sim$(echo $fileMinus).txt" "data/c_exp_${dir}/sim$(echo $file).txt"
        cp "data/c_exp_${dir}/sim$(echo $fileMinus)_data.txt" "data/c_exp_${dir}/sim$(echo $file)_data.txt"
    done
done