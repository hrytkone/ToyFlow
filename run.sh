#!/bin/bash
###############################################
PROG=`basename $0`
if [ $# -ne 6 ]
then
    echo "Usage: ${PROG} events ptdep vnScale multiScale extraConvPart bDecay"
    exit;
fi

export evts=$1
export ptdep=$2
export vnScale=$3
export multiScale=$4
export extraConvPart=$5
export bDecay=$6
export today=20200630

export fileg0=toyFlow_${today}_${multiScale}Multi_${extraConvPart}extraConvPart_${bDecay}bDecay_PtDep${ptdep}_Gran0_Scale${vnScale}
export fileg1=toyFlow_${today}_${multiScale}Multi_${extraConvPart}extraConvPart_${bDecay}bDecay_PtDep${ptdep}_Gran1_Scale${vnScale}

mkdir output/${fileg0}
mkdir output/${fileg1}

# ./toyFlow filename.root nEvents bUsePtDep bUseGran scale multiScale extraConvPart decays seedNum bSaveAsTrees
./toyFlow output/${fileg0}/${fileg0}.root ${evts} ${ptdep} 0 ${vnScale} ${multiScale} ${extraConvPart} ${bDecay} 1000 0
./toyFlow output/${fileg1}/${fileg1}.root ${evts} ${ptdep} 1 ${vnScale} ${multiScale} ${extraConvPart} ${bDecay} 1000 0

root -b -x -l -q "MakeGraphs.C(\"output/${fileg0}/${fileg0}.root\",\"output/${fileg0}/${fileg0}-plot-output-normal.root\")"
root -b -x -l -q "MakeGraphs.C(\"output/${fileg1}/${fileg1}.root\",\"output/${fileg1}/${fileg1}-plot-output-normal.root\")"

root -b -x -l -q "MakeCentralityGraphs.C(\"output/${fileg0}/${fileg0}.root\",\"output/${fileg0}/${fileg0}-plot-output-centrality.root\",6)"
root -b -x -l -q "MakeCentralityGraphs.C(\"output/${fileg1}/${fileg1}.root\",\"output/${fileg1}/${fileg1}-plot-output-centrality.root\",6)"

root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${fileg0}/${fileg0}.root\",\"output/${fileg0}/${fileg0}-plot-output-alice_comp.root\")"
root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${fileg1}/${fileg1}.root\",\"output/${fileg1}/${fileg1}-plot-output-alice_comp.root\")"

root -b -x -l -q "PlotCentralityData3sub_alice_ref.C(\"${fileg0}\",\"${fileg1}\",2,6,${multiScale},${vnScale})"
