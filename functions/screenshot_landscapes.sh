#!/bin/bash

DATE="14DEC2022"

CURR_DIR=`pwd`
INF_STAT=("non_inf" "inf")

for inf in ${INF_STAT[*]}; do

echo $inf
cd $CURR_DIR
cd ./figures/${DATE}/${inf}/landscapes/wo_2dose_b+o_adj/optimization_1

HTML_LIST=`ls -lart *.html | awk '{print $9}'`

echo "Capturing screenshots of html landscapes"

for f in $HTML_LIST; do
    new_name=${f%.*}.png
    echo $new_name
    open $f
    sleep 2
   # screencapture -x -w $new_name
    screencapture -x -m $new_name
done

done