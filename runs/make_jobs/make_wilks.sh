#!/bin/bash

#marray=($(seq -2.0 3.0 1.0))
marray=($(seq -0.5 1.5 1.0))
uearray=($(seq -1.3 0.05 -0.15)) #ideal -1.3 - -0.15
#umarray=0

#Given sin^2 (th14) * sin^2 (th24) = 4 |Umu4|^2 |Ue4|^2
# Solve for Umu4, such that sin^2 th24 remains constant
# Set sin^2 (th24) = 0.0011, 0.0032, 0.0050, 0.0079, 0.0093

sin2th24=0.0050


RUN_FOLDER=/cluster/home/jmical01/uboone/analysis/whipping_star/
XML=uboone_numi_wirecell_split_overlay_intrinsic.xml
OUTDIR=uboone_3+1_split_intrinsic_sin${sin2th24}
TAG=UBOONE
DETECTORS=M #M, I, S or combine like MI
[ ! -d ${RUN_FOLDER}/runs/scripts_${OUTDIR} ] && mkdir ${RUN_FOLDER}/runs/scripts_${OUTDIR}
[ ! -d ${RUN_FOLDER}/runs/joblogs_${OUTDIR} ] && mkdir ${RUN_FOLDER}/runs/joblogs_${OUTDIR}
cp run_all.sh ${RUN_FOLDER}/runs/scripts_${OUTDIR}

m_max=${#marray[@]}
((m_max=m_max-1))
ue_max=${#uearray[@]}
echo $m_max $ue_max


COUNT=0
for ((m=0; m<$m_max; m++));
do
    for ((u=0; u<$ue_max; u++));
    do
        # Calculate the absolute value of u
        begin=${marray[$m]}
        end_id=$(($m+1))
        end=${marray[$end_id]}
        ue4=${uearray[$u]}
        logUm4=$(echo "scale=10; sqrt($sin2th24 * (1 - e(2*${ue4}* l(10))))" | bc -l)
        um4=$(echo "scale=4; l($logUm4)/l(10)" | bc -l)

        echo $begin $end $um4 $ue4
        #name=`basename $file`
        name="m41_${begin}_${end}_um4_${um4}_ue4_${ue4}"
        if [ "$u" -eq 0 ]
        then
            sed -e "s|@@begin@@|${begin}|g" \
                -e "s|@@end@@|${end}|g" \
                -e "s|@@um4@@|${um4}|g" \
                -e "s|@@ue4@@|${ue4}|g" \
                -e "s|@@dir@@|${RUN_FOLDER}|g" \
                -e "s|@@outdir@@|${OUTDIR}|g" \
                -e "s|@@xml@@|${XML}|g" \
                -e "s|@@tag@@|${TAG}|g" \
                -e "s|@@det@@|${DETECTORS}|g" \
                < submit_gen_template.sh > $RUN_FOLDER/runs/scripts_${OUTDIR}/generate_${name}.sh
        fi
        sed -e "s|@@begin@@|${begin}|g" \
            -e "s|@@end@@|${end}|g" \
            -e "s|@@um4@@|${um4}|g" \
            -e "s|@@ue4@@|${ue4}|g" \
            -e "s|@@dir@@|${RUN_FOLDER}|g" \
            -e "s|@@outdir@@|${OUTDIR}|g" \
            -e "s|@@xml@@|${XML}|g" \
            -e "s|@@tag@@|${TAG}|g" \
            -e "s|@@det@@|${DETECTORS}|g" \
            < bash_test_template.sh > $RUN_FOLDER/runs/scripts_${OUTDIR}/test_${name}.sh
        let COUNT=$COUNT+1
   done
done
echo $COUNT
echo ${RUN_FOLDER}/runs/scripts_${OUTDIR}
