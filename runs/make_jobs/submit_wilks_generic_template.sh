#!/bin/bash 

########### SBATCH Lines for Resource Request ##########

#SBATCH --time=02:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                 # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --mem=12G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name gen_@@tag@@      # you can give your job a name for easier identification (same as -J)
#SBATCH --partition=batch,wongjiradlab
#SBATCH --exclude=i2cmp006,s1cmp001,s1cmp002,s1cmp003,p1cmp005,p1cmp041
#SBATCH --output @@dir@@/runs/joblogs_@@outdir@@/output_m41_@@begin@@_@@end@@_um4_@@um4@@_@@tag@@.log

########### Command Lines to Run ##################


module load singularity/3.5.3

SINGULARITY_IMAGE="/cluster/tufts/wongjiradlabnu//larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_pyspark.sif"
WILKS="@@dir@@/build/bin/sbnfit_wilks"
ARGS='-x @@dir@@/xml/@@xml@@ -b @@begin@@ -e @@end@@ -u @@um4@@ -a @@ue4@@ -l @@det@@ -t @@tag@@'

singularity exec -B /cluster:/cluster $SINGULARITY_IMAGE /bin/bash -c "source /usr/local/root/bin/thisroot.sh && ${WILKS} ${ARGS} -m gen && ${WILKS} $ARGS -m test"

