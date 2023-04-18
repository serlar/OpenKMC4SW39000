#!/bin/sh
PWD=`pwd`
echo $PWD
OUTFILE=origin_out

if [ ! -f "./$OUTFILE" ];then
    echo " output doesn't exit!"
else
    rm -f $OUTFILE
fi

bsub -b -J xulei_kmc -o $OUTFILE -m 1 -q q_swhfnl -n 184320 -cgsp 64  -share_size 15000 -host_stack 256 -cross_size 32 -ldm_share_mode 4 -ldm_share_size 4 ../../../src/spk_swg++ -in input

