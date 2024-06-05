#!/bin/bash

FILES="test*.sh"

for f in $FILES
do
    sbatch $f
done

