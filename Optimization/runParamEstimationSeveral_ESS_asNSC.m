%runParamEstimationSeveral_ESS_asNSC
clear
close all

rng('shuffle')

addpath(genpath(['..' filesep]))

%% Set variables
modelName = 'Avatar_HEALTH';
inputpathbase = split(pwd,'cardiovascularavatar_T2D_HT');
dataFolder = [inputpathbase{1},filesep 'cardiovascularavatar_T2D_HT' filesep 'Data'];
patNums = '';%eg: '16', if only one subject should be optimized
patrange = '1:10';%eg: '1:10', if you want to optimize several subjects
patientNums = findPatientRange([],patrange,dataFolder);
numpats = length(patientNums);
doPV=1; % include data of flow in pulmonary vein or not
method = 'HEALTH';
continueMCMC = 0;
doParamEst=1;
datename = datestr(now,'yymmdd-HHMMSS');
Setup(patNums,doParamEst,datename,patrange,dataFolder,doPV,method,continueMCMC)

restoredefaultpath

%% Run optimizations
tic
randomstart = 1; 
disp('***Starting optimization run 1 - random starts***')
parfor i = 1:(40*numpats) %40 per subject
    EstimateParametersESS_HEALTH(i,datename,patNums,modelName,patrange,dataFolder,randomstart) 
end
randomstart = 0; 
disp('***Starting optimization run 2 - from previous results***')
parfor i = 1:(10*numpats) %10 per subject, startguess from previous results
    EstimateParametersESS_HEALTH(i,datename,patNums,modelName,patrange,dataFolder,randomstart) 
end
disp('***Starting final optimization run***')
parfor i = 1:numpats %1 per subject, last run
    EstimateParametersESS_HEALTH_lastrun(i,datename,patNums,modelName,patrange,dataFolder) 
end

%% Information about the time of optimization
timetotalOpt = toc
timetotalOpt_minutes = timetotalOpt/60
timetotalOpt_hours = timetotalOpt_minutes/60


