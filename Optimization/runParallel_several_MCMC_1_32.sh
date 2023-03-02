#!/bin/bash
#
# Use 1 node, i.e. 1x32 cores, for 72 hours
#batch1: 1:32
#SBATCH -N 1
#SBATCH -t 156:00:00
#SBATCH --exclusive

# Update the project if possible. 
#git pull


# Setup parameters for the scripts
modelName='avatarHEALTH'
dataFolder='Data'
#range in the format: 5,6,7  or 1:5 (empty = run all patients)
patRange='1:32'
doParamEst=1
patNums=''
doPV=1
continueMCMC=0
now=$(date '+%Y%m%d-%H%M%S')
method='MCMC'
cdprojpath="cd('/proj/add your path');"

#Use hidden matlab r2017 module to compile AMICII models
module load MATLAB/.R2017b-nsc1
MATLAB_0='matlab -nodesktop -nodisplay -singleCompThread'
module load buildenv-gcc/2018a-eb
${MATLAB_0} -r "GenerateModels; exit"
# run setup script to create results folders
${MATLAB_0} -r "Setup('${patNums}',$doParamEst,'${now}','${patRange}','${dataFolder}',$doPV,'${method}',$continueMCMC); exit"

# Load the Matlab module
module load MATLAB/R2020a-nsc1 
module load parallel/20181122-nsc1
#module add gcc/8.2.0 
MATLAB='matlab -nodesktop -nodisplay -singleCompThread'

# The name of the Matlab script (without .m)
job=EstimateParametersMCMC_HEALTH

# Run main task. Note the explicit "exit" to exit Matlab.
seq 1 32| parallel --ssh=jobsh -S $(hostlist -e -s',' -d  -p "$SLURM_CPUS_ON_NODE/" $SLURM_JOB_NODELIST) ${MATLAB} -r "\"$cdprojpath  ${job}({},'${now}','${patNums}','${patRange}', '${dataFolder}',$continueMCMC);exit\""

# End of script
