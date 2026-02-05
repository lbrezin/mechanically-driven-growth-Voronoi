#!/bin/bash -l

# Set SCC project
#$ -P qedksk

# -l buyin
#$ -l h_rt=8:00:00
#$ -l gpus=1
#$ -l gpu_c=3.5

#$ -o /projectnb/qedksk/lbrezin/logs
#$ -j y
# -m b

module purge
module load cuda/8.0
module load netcdf/4.6.1
module load boost/1.69.0
module load eigen/3.3.5
module load cgal/4.9.1
module load netcdf-cxx/4.2

#Create directory and move the program to it
cd /projectnb/qedksk/lbrezin/results/
mkdir -p $JOB_ID
cd $JOB_ID/

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

../../Original-cellGPU/cellDivision.out -n 500 -v 0.05 -h 0 -k 1 -i 2 -u 4 -y 4 -a 1 -p 1 -g 1 -x 1 -f 0 -q -1.00 -s 0 -b 1 -d 1 -j 0 -w $SGE_TASK_ID -l 0 -t 200000 
