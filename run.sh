#!/bin/bash
###############################################
PROG=`basename $0`
if [ $# -ne 8 ]
then
    echo "Usage: ${PROG} events ptdep vnScale multiScale extraConvPart bDecay bptCut comment"
    exit;
fi

# If you run multiple parallel runs, remember to change seed number.

export evts=$1
export ptdep=$2
export vnScale=$3
export multiScale=$4
export extraConvPart=$5
export bDecay=$6
export bptCut=$7
export today=20200714
export comment=$8

export fileg0=toyFlow_${today}_${multiScale}Multi_${extraConvPart}extraConvPart_${bDecay}bDecay_PtDep${ptdep}_Gran0_Scale${vnScale}_ptCut${bptCut}_${comment}
export fileg1=toyFlow_${today}_${multiScale}Multi_${extraConvPart}extraConvPart_${bDecay}bDecay_PtDep${ptdep}_Gran1_Scale${vnScale}_ptCut${bptCut}_${comment}

mkdir output/${fileg0}
mkdir output/${fileg1}

# ./toyFlow filename.root nEvents bUsePtDep bUseGran scale multiScale extraConvPart decays seedNum bSaveAsTrees
./toyFlow output/${fileg0}/${fileg0}.root ${evts} ${ptdep} 0 ${vnScale} ${multiScale} ${extraConvPart} ${bDecay} 1000 0 ${bptCut}
./toyFlow output/${fileg1}/${fileg1}.root ${evts} ${ptdep} 1 ${vnScale} ${multiScale} ${extraConvPart} ${bDecay} 1000 0 ${bptCut}

./runFinalize.sh ${fileg0} ${fileg1}
