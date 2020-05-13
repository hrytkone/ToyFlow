#!/bin/bash
###############################################
PROG=`basename $0`
if [ $# -ne 3 ]
then
    echo "Usage: ${PROG} events ptdep multiScale"
    exit;
fi

export evts=$1
export ptdep=$2
#export gran=$3
export multiScale=$3
export today=20200430_${multiScale}Multi

export scale1=1.00
export scale2=0.80
export scale3=0.65

export file1g0=toyFlow_${today}_PtDep${ptdep}_Gran0_Scale${scale1}
export file1g1=toyFlow_${today}_PtDep${ptdep}_Gran1_Scale${scale1}
export file2g0=toyFlow_${today}_PtDep${ptdep}_Gran0_Scale${scale2}
export file2g1=toyFlow_${today}_PtDep${ptdep}_Gran1_Scale${scale2}
export file3g0=toyFlow_${today}_PtDep${ptdep}_Gran0_Scale${scale3}
export file3g1=toyFlow_${today}_PtDep${ptdep}_Gran1_Scale${scale3}

mkdir output/${file1g0}
mkdir output/${file1g1}
mkdir output/${file2g0}
mkdir output/${file2g1}
mkdir output/${file3g0}
mkdir output/${file3g1}

# ./toyFlow filename.root nEvents bUsePtDep bUseGran scale multiScale seedNum bSaveAsTrees
#./toyFlow output/${file1g0}/${file1g0}.root ${evts} ${ptdep} 0 ${scale1} ${multiScale} 1000 0
#./toyFlow output/${file1g1}/${file1g1}.root ${evts} ${ptdep} 1 ${scale1} ${multiScale} 1000 0
#./toyFlow output/${file2g0}/${file2g0}.root ${evts} ${ptdep} 0 ${scale2} ${multiScale} 1000 0
#./toyFlow output/${file2g1}/${file2g1}.root ${evts} ${ptdep} 1 ${scale2} ${multiScale} 1000 0
#./toyFlow output/${file3g0}/${file3g0}.root ${evts} ${ptdep} 0 ${scale3} ${multiScale} 1000 0
#./toyFlow output/${file3g1}/${file3g1}.root ${evts} ${ptdep} 1 ${scale3} ${multiScale} 1000 0
#
#root -b -x -l -q "MakeGraphs.C(\"output/${file1g0}/${file1g0}.root\",\"output/${file1g0}/${file1g0}-plot-output-normal.root\")"
#root -b -x -l -q "MakeGraphs.C(\"output/${file1g1}/${file1g1}.root\",\"output/${file1g1}/${file1g1}-plot-output-normal.root\")"
#root -b -x -l -q "MakeGraphs.C(\"output/${file2g0}/${file2g0}.root\",\"output/${file2g0}/${file2g0}-plot-output-normal.root\")"
#root -b -x -l -q "MakeGraphs.C(\"output/${file2g1}/${file2g1}.root\",\"output/${file2g1}/${file2g1}-plot-output-normal.root\")"
#root -b -x -l -q "MakeGraphs.C(\"output/${file3g0}/${file3g0}.root\",\"output/${file3g0}/${file3g0}-plot-output-normal.root\")"
#root -b -x -l -q "MakeGraphs.C(\"output/${file3g1}/${file3g1}.root\",\"output/${file3g1}/${file3g1}-plot-output-normal.root\")"
#
#root -b -x -l -q "MakeCentralityGraphs.C(\"output/${file1g0}/${file1g0}.root\",\"output/${file1g0}/${file1g0}-plot-output-centrality.root\",6)"
#root -b -x -l -q "MakeCentralityGraphs.C(\"output/${file1g1}/${file1g1}.root\",\"output/${file1g1}/${file1g1}-plot-output-centrality.root\",6)"
#root -b -x -l -q "MakeCentralityGraphs.C(\"output/${file2g0}/${file2g0}.root\",\"output/${file2g0}/${file2g0}-plot-output-centrality.root\",6)"
#root -b -x -l -q "MakeCentralityGraphs.C(\"output/${file2g1}/${file2g1}.root\",\"output/${file2g1}/${file2g1}-plot-output-centrality.root\",6)"
#root -b -x -l -q "MakeCentralityGraphs.C(\"output/${file3g0}/${file3g0}.root\",\"output/${file3g0}/${file3g0}-plot-output-centrality.root\",6)"
#root -b -x -l -q "MakeCentralityGraphs.C(\"output/${file3g1}/${file3g1}.root\",\"output/${file3g1}/${file3g1}-plot-output-centrality.root\",6)"
#
#root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${file1g0}/${file1g0}.root\",\"output/${file1g0}/${file1g0}-plot-output-alice_comp.root\")"
#root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${file1g1}/${file1g1}.root\",\"output/${file1g1}/${file1g1}-plot-output-alice_comp.root\")"
#root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${file2g0}/${file2g0}.root\",\"output/${file2g0}/${file2g0}-plot-output-alice_comp.root\")"
#root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${file2g1}/${file2g1}.root\",\"output/${file2g1}/${file2g1}-plot-output-alice_comp.root\")"
#root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${file3g0}/${file3g0}.root\",\"output/${file3g0}/${file3g0}-plot-output-alice_comp.root\")"
#root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${file3g1}/${file3g1}.root\",\"output/${file3g1}/${file3g1}-plot-output-alice_comp.root\")"


root -b -x -l -q "PlotCentralityData3sub_alice_ref.C(\"output/${file1g0}/\",2,6,${multiScale},${scale1})"
root -b -x -l -q "PlotCentralityData3sub_alice_ref.C(\"output/${file1g1}/\",2,6,${multiScale},${scale1})"
root -b -x -l -q "PlotCentralityData3sub_alice_ref.C(\"output/${file2g0}/\",2,6,${multiScale},${scale2})"
root -b -x -l -q "PlotCentralityData3sub_alice_ref.C(\"output/${file2g1}/\",2,6,${multiScale},${scale2})"
root -b -x -l -q "PlotCentralityData3sub_alice_ref.C(\"output/${file3g0}/\",2,6,${multiScale},${scale3})"
root -b -x -l -q "PlotCentralityData3sub_alice_ref.C(\"output/${file3g1}/\",2,6,${multiScale},${scale3})"
