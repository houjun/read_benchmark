#!/bin/bash
SHPATH="TSHPATH"
NCORE="TNCORE"
CONTACTOST="TCONTACTOST"

DATASIZE="64"
echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P.csv


DATASIZE=128
echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P.csv


DATASIZE="256"
echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P.csv


DATASIZE="512"
echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_0.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_0.csv

echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_1.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_1.csv

echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_2.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_2.csv

DATASIZE="1024"
echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_0.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_0.csv

echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_1.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_1.csv

echo "qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_2.csv"
qsub -I -V -q debug -lmppwidth=$NCORE -x ${SHPATH}/OSTtest_${DATASIZE}M_${CONTACTOST}P.sh>> ./OST_test_${DATASIZE}M_${CONTACTOST}P_2.csv
