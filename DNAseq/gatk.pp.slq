#!/bin/bash

# -----------------------[SLURM : START]---------------------------
# -- run name --
#SBATCH --job-name=GATK_PP

#  -- email preferences --
#SBATCH --mail-user=rodrigo.gularte@ulg.ac.be
#SBATCH --mail-type=END,FAIL
#     (ALL = BEGIN, END, FAIL, REQUEUE)


# -- output prefex
#SBATCH --output="pp-%j.out"
#     %j     Job allocation number.
#     %N     Main node name.  


# -- time requierements --
#SBATCH --time=15-00:00:00
# "days-hours:minutes:seconds"
# ** Note that the lower the requested run-time, the higher the
#    chances to get scheduled to 'fill in the gaps' between other
#    jobs. 

# ---- Resources ----
#SBATCH --ntasks-per-node=20 --mem-per-cpu=5G

# ---- Array Control ---- 
# -- format: 0-7:1 0 through 7 by one
#SBATCH --array=0-7

# -- email messages --
#SBATCH --mail-type=END

# -----------------------[SLURM : END]---------------------------

source ~/.bashrc
source ../envar.sh

## module add
module add Java/1.7.0_21
module add R/3.0.2

input=(
$(cat bam.files)
)


# You will need to merge the individual BAM files together before processing with GATK.                                                                  

## Run w/ 5 cores & 5G RAM
java -Xmx5000M \
    -jar ${GATK}/GenomeAnalysisTK.jar  \
    -T RealignerTargetCreator \
    -R $MM10 \
    -I ${input[$SLURM_ARRAY_TASK_ID]} \
    -nt $SLURM_JOB_CPUS_PER_NODE\
    -o ${input[$SLURM_ARRAY_TASK_ID]}.intervals
##    --known $INDEL

## Run w/ 1 core & 5G RAM
java -Xmx4000M -verbose\
    -jar ${GATK}/GenomeAnalysisTK.jar  \
    -R $MM10\
    -T IndelRealigner \
    -targetIntervals ${input[$SLURM_ARRAY_TASK_ID]}.intervals\
    -I ${input[$SLURM_ARRAY_TASK_ID]} \
#    -known $INDEL\
    -o ${input[$SLURM_ARRAY_TASK_ID]}.ir.bam


## Run w/ 10 cores & 5G RAM
java -Xmx5000M\
    -jar ${GATK}/GenomeAnalysisTK.jar  \
    -R $MM10\
    -T BaseRecalibrator\
    -I ${input[$SLURM_ARRAY_TASK_ID]}.ir.bam\
    -knownSites $SNP\
    -o ${input[$SLURM_ARRAY_TASK_ID]}.grp\
    -nct $SLURM_JOB_CPUS_PER_NODE
 
java -Xmx5000M -Djava.io.tmpdir=$GLOBALSCRATCH/tmp \
    -jar ${GATK}/GenomeAnalysisTK.jar  \
    -T PrintReads\
    -R $MM10\
    -I ${input[$SLURM_ARRAY_TASK_ID]}.ir.bam \
    -BQSR ${input[$SLURM_ARRAY_TASK_ID]}.grp\
    -nct $SLURM_JOB_CPUS_PER_NODE\
    -compress 6\
    -o ${input[$SLURM_ARRAY_TASK_ID]}.ir.bqsr.bam



#########################
#####   DO NOT RUN   ####
#########################
## Run w/ 1 core and 5G RAM
#java -Xmx5000M\
# -Djava.io.tmpdir=$GLOBALSCRATCH/tmp\
#    -jar ${GATK}/GenomeAnalysisTK.jar \
#    -R $MM10\
#    -T ReduceReads\
#    -I ${input[$SLURM_ARRAY_TASK_ID]}.ir.bqsr.bam\
#    -o ${input[$SLURM_ARRAY_TASK_ID]}.ir.bqsr.reduced.bam 
