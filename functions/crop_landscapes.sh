#!/bin/bash


DATE="14DEC2022"

CURR_DIR=`pwd`
INF_STAT=("non_inf" "inf")

for inf in ${INF_STAT[*]}; do

cd $CURR_DIR
cd ./figures/${DATE}/${inf}/landscapes/wo_2dose_b+o_adj/optimization_1

PNG_LIST=`ls -lart *.png | awk '{print $9}'`

echo "Cropping landscapes"

for f in $PNG_LIST; do
 #   convert ${f} -crop 2200x1200+500+600 ${f}
    convert ${f} -gravity center -crop 2000x1100-130+150 ${f}
done


done
