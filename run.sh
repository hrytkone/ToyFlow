#!/bin/bash
###############################################
PROG=`basename $0`
if [ $# -ne 9 ]
then
    echo "Usage: ${PROG} events ptdep vnScale multiScale extraConvPart bDecay bptCut effTPC comment"
    exit;
fi

# If you run multiple parallel runs, remember to change seed number.

export today=`date +"%Y%m%d"`

export evts=$1
export ptdep=$2
export vnScale=$3
export multiScale=$4
export extraConvPart=$5
export bDecay=$6
export bptCut=$7
export effTPC=$8
export comment=$9

export fileg0=toyFlow_${today}_${multiScale}Multi_${extraConvPart}extraConvPart_${bDecay}bDecay_PtDep${ptdep}_Gran0_Scale${vnScale}_ptCut${bptCut}_effTPC${effTPC}_${comment}
export fileg1=toyFlow_${today}_${multiScale}Multi_${extraConvPart}extraConvPart_${bDecay}bDecay_PtDep${ptdep}_Gran1_Scale${vnScale}_ptCut${bptCut}_effTPC${effTPC}_${comment}

mkdir output/${fileg0}
mkdir output/${fileg1}

#./toyFlow filename.root nEvents bUsePtDep bUseGran scale multiScale extraConvPart decays seedNum bSaveAsTrees bptCuts effTPC
./toyFlow output/${fileg0}/${fileg0}.root ${evts} ${ptdep} 0 ${vnScale} ${multiScale} ${extraConvPart} ${bDecay} 1000 0 ${bptCut} ${effTPC}
./toyFlow output/${fileg1}/${fileg1}.root ${evts} ${ptdep} 1 ${vnScale} ${multiScale} ${extraConvPart} ${bDecay} 1000 0 ${bptCut} ${effTPC}

./runFinalize.sh ${fileg0} ${fileg1}
