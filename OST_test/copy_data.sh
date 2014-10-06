#!/bin/bash

for i in 0 1 2 3 4
do
    echo $i
    cp /scratch3/scratchdirs/houhun/data/vpic-medium/original/eparticle_T22860_1.5_filter.h5p  ./stripe4/test$i.bin
done


