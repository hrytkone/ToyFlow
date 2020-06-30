#!/bin/bash
###############################################

export evts=10000
export ptdep=0

#./run.sh events ptdep vnScale multiScale extraConvPart bDecay"
./run.sh ${evts} ${ptdep} 0.80 0.8        0.2           0.0
./run.sh ${evts} ${ptdep} 1.00 1.0        0.0           0.1
./run.sh ${evts} ${ptdep} 1.00 1.0        0.0           0.2
