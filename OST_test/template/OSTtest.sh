#!/bin/bash

EXEPATH="/global/homes/h/houhun/edison/OST_test/bin/new_OST_multitest"
DATAPATH="/scratch3/scratchdirs/houhun/data/OST"

DATASIZE="TDATASIZE"
REPEAT="TREPEAT"
OSTPERCORE="TOSTPERCORE"

echo "aprun -n ${OSTPERCORE} -N 1 ${EXEPATH} ${DATAPATH} ${DATASIZE} ${REPEAT} ${OSTPERCORE}"
aprun -n ${OSTPERCORE} -N 1 ${EXEPATH} ${DATAPATH} ${DATASIZE} ${REPEAT} ${OSTPERCORE}



