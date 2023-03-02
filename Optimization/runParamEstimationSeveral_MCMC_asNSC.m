%runParamEstimationSeveral_MCMC_asNSC
clear
close all

rng('shuffle')
addpath(genpath(['..' filesep]))

%% Set variables
modelName = 'Avatar_HEALTH';
patNums = '';%eg: '16'
patrange = '1:32';%eg: 1:5
numberOfSubjects = 32;
inputpathbase = split(pwd,'cardiovascularavatar_T2D_HT');
dataFolder = [inputpathbase{1},filesep 'cardiovascularavatar_T2D_HT' filesep 'Data'];
doPV=1; % include data of flow in pulmonary vein or not
method = 'MCMC';
continueMCMC = 0; %1=load old MCMC optimization file
doParamEst=1;
datename = datestr(now,'yymmdd-HHMMSS');
Setup(patNums,doParamEst,datename,patrange,dataFolder,doPV,method,continueMCMC)

restoredefaultpath

%% Run optimizations
tic
parfor i = 1:numberOfSubjects
    EstimateParametersMCMC_HEALTH(i,datename,patNums,patrange,dataFolder,continueMCMC)
end

%% Information about the time of optimization
timetotalOpt = toc
timetotalOpt_minutes = timetotalOpt/60
timetotalOpt_hours = timetotalOpt_minutes/60


