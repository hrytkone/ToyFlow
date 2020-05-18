#!/bin/bash
###############################################
PROG=`basename $0`
if [ $# -ne 5 ]
then
    echo "Usage: ${PROG} events ptdep multiScale extraConvPart bDecay"
    exit;
fi

export evts=$1
export ptdep=$2
export multiScale=$3
export extraConvPart=$4
export bDecay=$5
./run.sh ${evts} ${ptdep} 1.00 ${multiScale} ${extraConvPart} ${bDecay}
./run.sh ${evts} ${ptdep} 0.80 ${multiScale} ${extraConvPart} ${bDecay}
./run.sh ${evts} ${ptdep} 0.65 ${multiScale} ${extraConvPart} ${bDecay}
