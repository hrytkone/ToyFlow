#!/bin/bash
###############################################

export evts=100000
export ptdep=0

#./run.sh events ptdep vnScale multiScale extraConvPart bDecay bptCut"
./run.sh ${evts} ${ptdep} 1.00 0.3        0.0           0.0    0
./run.sh ${evts} ${ptdep} 1.00 1.0        0.8           0.0    0
./run.sh ${evts} ${ptdep} 0.80 1.0        0.0           0.0    0
