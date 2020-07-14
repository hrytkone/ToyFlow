#!/bin/bash
###############################################
PROG=`basename $0`
if [ $# -ne 2 ]
then
    echo "Usage: ${PROG} fileNoGran fileGran"
    exit;
fi

export fileg0=$1
export fileg1=$2

root -b -x -l -q "MakeGraphs.C(\"output/${fileg0}/${fileg0}.root\",\"output/${fileg0}/${fileg0}-plot-output-normal.root\")"
root -b -x -l -q "MakeGraphs.C(\"output/${fileg1}/${fileg1}.root\",\"output/${fileg1}/${fileg1}-plot-output-normal.root\")"

root -b -x -l -q "MakeCentralityGraphs.C(\"output/${fileg0}/${fileg0}.root\",\"output/${fileg0}/${fileg0}-plot-output-centrality.root\",6)"
root -b -x -l -q "MakeCentralityGraphs.C(\"output/${fileg1}/${fileg1}.root\",\"output/${fileg1}/${fileg1}-plot-output-centrality.root\",6)"

root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${fileg0}/${fileg0}.root\",\"output/${fileg0}/${fileg0}-plot-output-alice_comp.root\")"
root -b -x -l -q "MakeCentralityGraphs_alice_ref.C(\"output/${fileg1}/${fileg1}.root\",\"output/${fileg1}/${fileg1}-plot-output-alice_comp.root\")"

root -b -x -l -q "PlotCentralityData3sub_alice_ref.C(\"${fileg0}\",\"${fileg1}\",2,6,${multiScale},${vnScale})"
