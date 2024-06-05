#!/bin/bash

DIRS="/cluster/home/jmical01/uboone/analysis/whipping_star//runs/scripts_uboone_3+1_split_intrinsic_sin0.00?[1,6]"
FILES="test*.sh"

for d in $DIRS
do
    cd $d
    echo $d
    for f in $FILES
	do
    	sbatch $f
    done
done

