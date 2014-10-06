#!/bin/bash
let CORE=$1*24

mkdir ./$1P
cd ./$1P

cp ../template/batch_submit.sh ./batch_submit$1.sh
sed -i "s!TSHPATH!\/global\/homes\/h\/houhun\/edison\/OST_test\/$1P!g" ./batch_submit$1.sh
sed -i "s!TNCORE!${CORE}!g" ./batch_submit$1.sh
sed -i "s!TCONTACTOST!$1!g" ./batch_submit$1.sh

for i in 64 128 256 512 1024
do
    cp ../template/OSTtest.sh ./OSTtest_${i}M_$1P.sh
    sed -i "s!TDATASIZE!$i!g" ./OSTtest_${i}M_$1P.sh

    if [ $i == 512 ]; then
        REPEAT="3"
    elif [ $i == 1024 ]; then
        REPEAT="3"
    else
        REPEAT="9"
    fi
    
    sed -i "s!TREPEAT!${REPEAT}!g" ./OSTtest_${i}M_$1P.sh
    sed -i "s!TOSTPERCORE!$1!g" ./OSTtest_${i}M_$1P.sh


done
