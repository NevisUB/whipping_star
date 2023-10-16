#!/bin/bash

FILES="test*.sh"

for f in $FILES
do
    sbatch $f
    sleep 3
done

