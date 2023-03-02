function [] =  EstimateParametersMCMC_HEALTH(loopNumber,dateStart,patientNums,patrange,dataFolder,continueMCMC)
%% Setup
doMultiStart = 0;
doSensitivites = 0;

s=rng('shuffle'); %do not note remove this line
d = datestr(now,'YYmmDD-HHMMSSFFF'); %creating different seed for each start
rngSettings=rng(s.Seed+loopNumber*sum(d)); %do not note remove this line

pathbase = split(pwd,'cardiovascularavatar_T2D_HT');
pathbase = [pathbase{1},filesep, 'cardiovascularavatar_T2D_HT'];

addpath(genpath(pathbase))
addpath([pathbase filesep 'Parameters'])
addpath([pathbase filesep 'Modelfiles'])

addpath([pathbase filesep 'Requirements' filesep 'PESTO-1.1.0'])
addpath([pathbase filesep 'Requirements' filesep 'AMICI-0.10.11_SS_eventFix' filesep 'matlab'])
addpath([pathbase filesep 'Requirements' filesep 'MEIGO'])

run([pathbase filesep 'Requirements' filesep 'AMICI-0.10.11_SS_eventFix' filesep 'matlab' filesep 'installAMICI.m'])
run([pathbase filesep 'Requirements' filesep 'MEIGO' filesep 'install_MEIGO.m'])

%% Choose which subject to simulate based on the range of subjects and which loop number you are in
patientNums = findPatientRange(patientNums,patrange,dataFolder);
num = mod(loopNumber,length(patientNums))+1;
patientNum = patientNums{num};
disp(['----------------------- P' patientNum ' ---------------------------------'])

%% Load parameters and data
resultsFolder = [pathbase 'Parameters' filesep 'ESS'];
[~,data,~,pvaluessetup,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues] = setup_simulations_HEALTH({patientNum},resultsFolder,1);
% settings for the cost function:
numHeartBeats = 20;
doPlot = 0;
dispErrors = 0;
%% Set bounds (separate for each subject, since paramdata is individual)
[lbOrig,ubOrig] = loadParamBounds_HEALTH(ind,data.params);
lbOrig = lbOrig';%mcmc wants format n x 1
ubOrig = ubOrig';
lb = log10(lbOrig);
ub = log10(ubOrig);
parameters.min = lb;
parameters.max = ub;
parameters.name = paramNames';%mcmc wants format n x 1
% the toolbox PESTO uses this parameters structure, list all parameters names for automatic plottning
parameters.number = length(parameters.name); % number of parameters

%% Find save folder
savefolders = dir([pathbase filesep 'Parameters' filesep 'MCMC' filesep 'P' patientNum '_*']);
[~,latestfolderInd] = max([savefolders.datenum]);
foldernames = {savefolders.name};
folderpaths = {savefolders.folder};
savefolder = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd})

%% MCMC options
 % Options
optionsPesto = PestoOptions();% loads optimization options and options for everything basically
optionsPesto.obj_type = 'negative log-posterior';  % Works with least squares cost which we typically been using, i.e the goal is the minimize. so keep it as is
optionsPesto.n_starts = 50; % number of multi-starts optimization performed
optionsPesto.mode = 'visual'; %visual it will plot the results itteratively.    * 'text': optimization results for multi-start are printed on screen* 'silent': no output during the multi-start local optimization
optionsPesto.comp_type = 'sequential'; % no parallel computing

% Setup for how awesome you are, if using SBtoolbox, set =1, if using AMICI
% and can return dJ/dp set =2, if superduper-awesome cool and can calculate
% hessian, set to =3 (J=cost)
% i.e, this is the number of output arguments you can have from your
% objective-function, 1=cost only, 2=cost & dJ/dp, 3=cost &dJ/dp & hessian
% (in that order).
if doSensitivites
    optionsPesto.objOutNumber=2;
else
    optionsPesto.objOutNumber=1;
end
    
%% Markov Chain Monte Carlo sampling -- Parameters
% Values for the parameters are sampled by using an Parallel Tempering (PT)
% algorithm. This way, the underlying probability density of the parameter 
% distribution can be captured. Since only one temperature is used, this is
% effectively an adapted Metropolis algorithm single-chain algorithm.

% Building a struct covering all sampling options:
optionsPesto.MCMC = PestoSamplingOptions();
optionsPesto.MCMC.nIterations = 1e5;% This is the most important setting, covering how many iterations to be performed, here i think the minimum is 1e5 and above, especially now, as we only want the really good (i.e significant) parameter vectors
optionsPesto.MCMC.mode = optionsPesto.mode;

optionsPesto.MCMC.saveEach = 1;
if continueMCMC
    continuenames = dir(sprintf('%s%sMCMCrun_P%s_*',savefolder,filesep,patientNum));
    continuesavefile=fullfile(continuenames(1).folder,continuenames(1).name);
    % OBS performRAMPART adds the .mat, if included here it can't find it
    continuesavefile = strrep(continuesavefile,'.mat','');
    fprintf('Continuing mcmc from file %s \n',continuesavefile)
    optionsPesto.MCMC.saveFileName = continuesavefile;
else
    memory_address = system_dependent('getpid'); % create an unuqie name for the save file
    optionsPesto.MCMC.saveFileName = sprintf('%s%sMCMCrun_P%s_%diter-%i%i',savefolder,filesep,patientNum,optionsPesto.MCMC.nIterations,loopNumber,memory_address);
end
optionsPesto.MCMC.PT.nTemps = parameters.number;

%% Using RAMPART  Most of this should be left as is, however look in the file called RAMPARTOPTIONS in @RAMPARTOPTIONS folder to see what each setting does
optionsPesto.MCMC.samplingAlgorithm     = 'RAMPART';
optionsPesto.MCMC.RAMPART.nTemps           =  parameters.number;
optionsPesto.MCMC.RAMPART.exponentT        = 1000;
optionsPesto.MCMC.RAMPART.maxT             = 2000;
optionsPesto.MCMC.RAMPART.alpha            = 0.51;
optionsPesto.MCMC.RAMPART.temperatureNu    = 1e3;
optionsPesto.MCMC.RAMPART.memoryLength     = 1;
optionsPesto.MCMC.RAMPART.regFactor        = 1e-8;
optionsPesto.MCMC.RAMPART.temperatureEta   = 10;

optionsPesto.MCMC.RAMPART.trainPhaseFrac   = 0.1;
optionsPesto.MCMC.RAMPART.nTrainReplicates = 5;

optionsPesto.MCMC.RAMPART.RPOpt.rng                  = 1;
optionsPesto.MCMC.RAMPART.RPOpt.nSample              = floor(optionsPesto.MCMC.nIterations*optionsPesto.MCMC.RAMPART.trainPhaseFrac)-1;
optionsPesto.MCMC.RAMPART.RPOpt.crossValFraction     = 0.2;
optionsPesto.MCMC.RAMPART.RPOpt.modeNumberCandidates = 1:20;
optionsPesto.MCMC.RAMPART.RPOpt.displayMode          = 'silent';%'text';
optionsPesto.MCMC.RAMPART.RPOpt.maxEMiterations      = 100;
%       optionsPesto.MCMC.RAMPART.RPOpt.nDim                 = parameters.number;
optionsPesto.MCMC.RAMPART.RPOpt.nSubsetSize          = 1000;
%       optionsPesto.MCMC.RAMPART.RPOpt.lowerBound           = parameters.min;
%       optionsPesto.MCMC.RAMPART.RPOpt.upperBound           = parameters.max;
% tolMu & tolSigma --> Break if terminiation condition was reached before i == nAlg 'Terminated because movement tolerances were reached.'
%       optionsPesto.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (parameters.max(1)-parameters.min(1));
%       optionsPesto.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (parameters.max(1)-parameters.min(1));
optionsPesto.MCMC.RAMPART.RPOpt.dimensionsToPlot     = [1,2];
%       optionsPesto.MCMC.RAMPART.RPOpt.isInformative        = [1,1,ones(1,optionsPesto.MCMC.RAMPART.RPOpt.nDim-2)];

%% define rest of RAMPART settings
optionsPesto.MCMC.RAMPART.RPOpt.isInformative        = ones(1,length(parameters.min)); %[1,1,ones(1,optionsPesto.MCMC.RAMPART.RPOpt.nDim-2)];
optionsPesto.MCMC.RAMPART.RPOpt.nDim                 = parameters.number;
optionsPesto.MCMC.RAMPART.RPOpt.lowerBound           = parameters.min;
optionsPesto.MCMC.RAMPART.RPOpt.upperBound           = parameters.max;
optionsPesto.MCMC.RAMPART.RPOpt.tolMu                = 1e-4 * (parameters.max(1)-parameters.min(1));
optionsPesto.MCMC.RAMPART.RPOpt.tolSigma             = 1e-2 * (parameters.max(1)-parameters.min(1));

%%
warning('off','all') % AMICI prints error-messages for all simulations (if they fail) so this will fix that annoying orange text form appearing
format long
format compact


if doMultiStart
    %% Performs the multi-start optimization (skip this if the results aren't great (i.e the same cost as you have in your optimal parameter set)
    allparams = pvaluessetup;
    ind.estParams = 1:length(allparams);
    objectiveFunction = @(p) costFunction_HEALTH(p,constants,allparams,simOptions,numHeartBeats,data,ind,dispErrors,doPlot);
    tic
    parameters = getMultiStarts(parameters, objectiveFunction, optionsPesto);
    timemultistart = toc
    %% Here important part, these two lines are only valid if you have done the multi-start opt
    optionsPesto.MCMC.theta0 = parameters.MS.par(:,1);
    optionsPesto.MCMC.sigma0 = 0.5 * inv(squeeze(parameters.MS.hessian(:,:,1)));
else
    % Find startguess if not multistart
    %take latest result from previous ess optimization
    loadfolders = dir(fullfile(pathbase,'Parameters','MCMC',['P' patientNum '_*']));
    [~,latestfolderInd] = max([loadfolders.datenum]);
    foldernames = {loadfolders.name};
    folderpaths = {loadfolders.folder};
    loadfolderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    loadsummaryfile=0;
    [optParam,loadedoptcost,~,loadedconstants] = findBestParams(loadfolderName,loadsummaryfile,length(paramNames));
    fprintf('Estimateparameters: Loaded paramvalues with cost %0.2f from %s\n',loadedoptcost,loadfolderName)
    optParam = max(optParam,lbOrig);%make sure the startguess is within the bounds
    optParam = min(optParam,ubOrig);%make sure the startguess is within the bounds
    allparams = optParam;
    ind.estParams = 1:length(allparams);
    optParam = log10(optParam);
    optParam = max(optParam,lb);%make sure the startguess is within the bounds
    optParam = min(optParam,ub);%make sure the startguess is within the bounds
    
    % if multi-start does not work use best param
    optionsPesto.MCMC.theta0 = optParam;
    % if no sigma given by user, it is set to deafault: this.sigma0 = 1e4 * diag(ones(1,par.number));
    optionsPesto.MCMC.sigma0 = 1e5*eye(length(optParam));
    objectiveFunction = @(p) costFunction_HEALTH(p,constants,allparams,simOptions,numHeartBeats,data,ind,dispErrors,doPlot);
end
save(sprintf('%s/testsave-%i-%s.mat',savefolder,loopNumber,dateStart) ,'dateStart')

%% Sample with MCMC
%Please make sure opt.theta0, the par.number and opt.PT.nTemps are consistent.
tic
% try
    parameters = getParameterSamples(parameters, objectiveFunction, optionsPesto);
% catch ME
%     disp(['ERROR - pat ' patientNum ' failed getParameterSamples:'])
%     disp(ME.identifier)
% end
timeMCMC = toc;
hoursMCMC = toc/3600;

%% Save results & settings to folder
bounds.lb = 10.^lb;
bounds.ub = 10.^ub;
dateEnd = datestr(now,'yymmdd-HHMMSS');
save(sprintf('%s%sparameters-%i.mat',savefolder,filesep,loopNumber),'parameters','bounds')
save(sprintf('%s%soptionsPesto-%i.mat',savefolder,filesep,loopNumber),'optionsPesto')
save(sprintf('%s%srngSettings-%i.mat',savefolder,filesep,loopNumber),'rngSettings')
save(sprintf('%s%stimeforMCMC-%i.mat',savefolder,filesep,loopNumber),'timeMCMC','hoursMCMC')
save(sprintf('%s%sextrainfoMCMC-%i.mat',savefolder,filesep,loopNumber),'constants','constantsNames','ind','paramNames','simOptions','dateStart','dateEnd')
w = warning ('on','all'); % yay error messages again

%% Save all open plots
disp('SAVING PLOTS...')
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for fig = 1:length(FigList)
    currFig = FigList(fig);
    figName = currFig.Name;
    saveas(currFig, fullfile(savefolder,[figName,'.png']))
%     saveas(currFig, fullfile(savefolder,[figName,'.fig'])) %fig takes too
%     much memory
end
    
end
