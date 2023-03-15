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

listDirectories=("356" "376" "396" "416" "436" "456" "476" "496" "516" "536" "556")
listFiles=(2694)

for dir in "${listDirectories[@]}"; do
    for file in "${listFiles[@]}"; do
        fileMinus=`expr $file - 1`
        cp "data/c_exp_${dir}/sim$(echo $fileMinus).txt" "data/c_exp_${dir}/sim$(echo $file).txt"
        cp "data/c_exp_${dir}/sim$(echo $fileMinus)_data.txt" "data/c_exp_${dir}/sim$(echo $file)_data.txt"
    done
done

listDirectories=("1" "5" "9" "13" "17" "456")
listFiles=(695 696 955 956 1067 1068 1258 2023 2050 2084 2085 2652 2653 2693 2694 2713 2714 3413 3414 3522 3523 4273 4274 4388 4389)

for dir in "${listDirectories[@]}"; do
    for file in "${listFiles[@]}"; do
        fileMinus=`expr $file - 1`
        cp "data/c_exp_${dir}/sim$(echo $fileMinus).txt" "data/c_exp_${dir}/sim$(echo $file).txt"
        cp "data/c_exp_${dir}/sim$(echo $fileMinus)_data.txt" "data/c_exp_${dir}/sim$(echo $file)_data.txt"
    done
done

listDirectories=("101" "105" "109" "113" "117" "121")
listFiles=(1257 1258 1649 2652 2653)

for dir in "${listDirectories[@]}"; do
    for file in "${listFiles[@]}"; do
        fileMinus=`expr $file - 1`
        cp "data/c_exp_${dir}/sim$(echo $fileMinus).txt" "data/c_exp_${dir}/sim$(echo $file).txt"
        cp "data/c_exp_${dir}/sim$(echo $fileMinus)_data.txt" "data/c_exp_${dir}/sim$(echo $file)_data.txt"
    done
done