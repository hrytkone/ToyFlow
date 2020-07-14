#!/bin/bash
###############################################
PROG=`basename $0`
if [ $# -ne 0 ]
then
    echo "Usage: ${PROG} "
    exit;
fi

export commong0=toyFlow_20200708_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran1_Scale1.00_ptCut0
export commong1=toyFlow_20200708_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran1_Scale1.00_ptCut0

export file01=output/${commong0}_run01/${commong0}_run01.root
export file02=output/${commong0}_run02/${commong0}_run02.root
export file03=output/${commong0}_run03/${commong0}_run03.root
export file04=output/${commong0}_run04/${commong0}_run04.root
export file05=output/${commong0}_run05/${commong0}_run05.root
export file06=output/${commong0}_run06/${commong0}_run06.root
export file07=output/${commong0}_run07/${commong0}_run07.root
export file08=output/${commong0}_run08/${commong0}_run08.root
export file09=output/${commong0}_run09/${commong0}_run09.root
export file10=output/${commong0}_run10/${commong0}_run10.root

mkdir output/${commong0}

hadd -f output/${commong0}/${commong0}.root ${file01} ${file02} ${file03} ${file04} ${file05} ${file06} ${file07} ${file08} ${file09} ${file10} 

export file01=output/${commong1}_run01/${commong1}_run01.root
export file02=output/${commong1}_run02/${commong1}_run02.root
export file03=output/${commong1}_run03/${commong1}_run03.root
export file04=output/${commong1}_run04/${commong1}_run04.root
export file05=output/${commong1}_run05/${commong1}_run05.root
export file06=output/${commong1}_run06/${commong1}_run06.root
export file07=output/${commong1}_run07/${commong1}_run07.root
export file08=output/${commong1}_run08/${commong1}_run08.root
export file09=output/${commong1}_run09/${commong1}_run09.root
export file10=output/${commong1}_run10/${commong1}_run10.root

mkdir output/${commong1}

hadd -f output/${commong1}/${commong1}.root ${file01} ${file02} ${file03} ${file04} ${file05} ${file06} ${file07} ${file08} ${file09} ${file10} 

./runFinalize.sh ${commong0} ${commong1}
