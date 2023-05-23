#!/bin/bash
#
# Use 1 node, i.e. 1x32 cores, for 22 hours.
#SBATCH -N 1
#SBATCH -t 22:00:00
#SBATCH --exclusive

# Update the project if possible. 
#git pull


# Setup parameters for the scripts
modelName='avatarHEALTH'
dataFolder='Data'
#range in the format: 5,6,7  or 1:5 (empty = run all patients)
patRange=''
doParamEst=1
patNums='5'
doPV=1
continueMCMC=0
now=$(date '+%Y%m%d-%H%M%S')
method='HEALTH_BP_sensitivity'
cdprojpath="cd('/proj/cardiovascular_modeling/users/x_kajtu/cardiovascularavatar_T2D_HT/Optimization');"

#Use hidden matlab r2017 module to compile AMICII models
module load MATLAB/.R2017b-nsc1
MATLAB_0='matlab -nodesktop -nodisplay -singleCompThread'
module load buildenv-gcc/2018a-eb
cd Modelfiles
${MATLAB_0} -r "GenerateModels; exit"
# run setup script to create results folders
cd ../Simulation
${MATLAB_0} -r "Setup('${patNums}',$doParamEst,'${now}','${patRange}','${dataFolder}',$doPV,'${method}',$continueMCMC); exit"
cd ..

# Load the Matlab module
module load MATLAB/R2020a-nsc1 
module load parallel/20181122-nsc1
#module add gcc/8.2.0 
MATLAB='matlab -nodesktop -nodisplay -singleCompThread'

# The name of the Matlab script (without .m)
job=EstimateParametersESS_HEALTH_changeBP

# ---------SBP max---------
# Run main task. Note the explicit "exit" to exit Matlab.
# random start guess - 40 runs * 1 pat
randomstart=1
SBPchange=16
DBPchange=0
seq 1 40| parallel --ssh=jobsh -S $(hostlist -e -s',' -d  -p "$SLURM_CPUS_ON_NODE/" $SLURM_JOB_NODELIST) ${MATLAB} -r "\"$cdprojpath  ${job}({},'${now}','${patNums}', '${modelName}','${patRange}', '${dataFolder}',$randomstart,$SBPchange,$DBPchange);exit\""

# load start guess - 10 runs * 1 pats
randomstart=0
seq 1 10| parallel --ssh=jobsh -S $(hostlist -e -s',' -d  -p "$SLURM_CPUS_ON_NODE/" $SLURM_JOB_NODELIST) ${MATLAB} -r "\"$cdprojpath  ${job}({},'${now}','${patNums}', '${modelName}','${patRange}', '${dataFolder}',$randomstart,$SBPchange,$DBPchange);exit\""

# End of script
