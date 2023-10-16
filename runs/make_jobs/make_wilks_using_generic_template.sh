#!/bin/bash

marray=($(seq -0.6 0.2 1.0))
#umarray=($(seq -1.25 0.2 -0.15)) #ideal -1.3 - -0.15
umarray=0

RUN_FOLDER=/cluster/home/jmical01/uboone/analysis/whipping_star/
XML=uboone_numi_wirecell_apponly_split_overlay_intrinsic.xml
OUTDIR=uboone_apponly_split_CNP
TAG=UBOONE
DETECTORS=M #M, I, S or combine like MI
[ ! -d ${RUN_FOLDER}/runs/scripts_${OUTDIR} ] && mkdir ${RUN_FOLDER}/runs/scripts_${OUTDIR}
[ ! -d ${RUN_FOLDER}/runs/joblogs_${OUTDIR} ] && mkdir ${RUN_FOLDER}/runs/joblogs_${OUTDIR}

m_max=${#marray[@]}
((m_max=m_max-1))
um_max=${#umarray[@]}
echo $m_max $um_max


COUNT=0
for ((m=0; m<$m_max; m++));
do
    for ((u=0; u<$um_max; u++));
    do
        # Calculate the absolute value of u
        begin=${marray[$m]}
        end_id=$(($m+1))
        end=${marray[$end_id]}
        um4=${umarray[$u]}
        ue4=-10

        echo $begin $end $um4 $ue4
        #name=`basename $file`
        name="m41_${begin}_${end}_um4_${um4}_ue4_${ue4}"
        sed -e "s|@@begin@@|${begin}|g" \
            -e "s|@@end@@|${end}|g" \
            -e "s|@@um4@@|${um4}|g" \
            -e "s|@@ue4@@|${ue4}|g" \
            -e "s|@@dir@@|${RUN_FOLDER}|g" \
            -e "s|@@outdir@@|${OUTDIR}|g" \
            -e "s|@@xml@@|${XML}|g" \
            -e "s|@@tag@@|${TAG}|g" \
            -e "s|@@det@@|${DETECTORS}|g" \
            < submit_wilks_generic_template.sh > $RUN_FOLDER/runs/scripts_${OUTDIR}/${name}.sh
        let COUNT=$COUNT+1
   done
done
echo $COUNT
echo ${RUN_FOLDER}/runs/scripts_${OUTDIR}
