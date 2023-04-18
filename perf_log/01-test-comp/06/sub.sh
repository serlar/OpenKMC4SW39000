#!/bin/sh
PWD=`pwd`
echo $PWD
OUTFILE=origin_out

if [ ! -f "./$OUTFILE" ];then
    echo " output doesn't exit!"
else
    rm -f $OUTFILE
fi

bsub -b -J xulei_kmc -o $OUTFILE -m 1 -q q_sw_expr -n 1 -cgsp 64 -share_size 15000 -host_stack 256  ../../src/spk_swg++ -in input
