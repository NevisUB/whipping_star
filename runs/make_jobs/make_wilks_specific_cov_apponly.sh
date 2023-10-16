#!/bin/bash

#marray=($(seq -2.0 3.0 1.0))
marray=($(seq -0.5 1.5 1.0))
uearray=($(seq -2.3 0.1 -0.3)) #ideal -1.3 - -0.15

RUN_FOLDER=/cluster/home/jmical01/uboone/analysis/whipping_star/
XML=icarus_uboone_numi_wirecell_apponly_split_overlay_intrinsic.xml
OUTDIR=icarus_uboone_apponly_split_intrinsic_NOoffdiag
COV=sbnfit_2detectors_file_systematics_split_int_over_nodirt.root
TAG=UBOONE_ICARUS
DETECTORS=MI #M, I, S or combine like MI
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
        um4=0

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
		-e "s|@@cov@@|${COV}|g" \
                -e "s|@@tag@@|${TAG}|g" \
                -e "s|@@det@@|${DETECTORS}|g" \
                < submit_gen_template_specific_cov.sh > $RUN_FOLDER/runs/scripts_${OUTDIR}/generate_${name}.sh
        fi
        sed -e "s|@@begin@@|${begin}|g" \
            -e "s|@@end@@|${end}|g" \
            -e "s|@@um4@@|${um4}|g" \
            -e "s|@@ue4@@|${ue4}|g" \
            -e "s|@@dir@@|${RUN_FOLDER}|g" \
            -e "s|@@outdir@@|${OUTDIR}|g" \
            -e "s|@@xml@@|${XML}|g" \
	    -e "s|@@cov@@|${COV}|g" \
            -e "s|@@tag@@|${TAG}|g" \
            -e "s|@@det@@|${DETECTORS}|g" \
            < submit_test_template_specific_cov.sh > $RUN_FOLDER/runs/scripts_${OUTDIR}/test_${name}.sh
        let COUNT=$COUNT+1
   done
done
echo $COUNT
echo ${RUN_FOLDER}/runs/scripts_${OUTDIR}
