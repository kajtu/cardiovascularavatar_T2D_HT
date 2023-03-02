function [] =  EstimateParametersESS_HEALTH_lastrun(loopNumber,dateStart,patientNums,modelName,patrange,dataFolder)
s=rng('shuffle'); %do not note remove this line
d = datestr(now,'YYmmDD-HHMMSSFFF'); %creating different seed for each start
s=rng(s.Seed+loopNumber*sum(d)); %do not note remove this line

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


patientNums = findPatientRange(patientNums,patrange,dataFolder);
num = mod(loopNumber,length(patientNums))+1;
patientNum = patientNums{num};
disp(['-----------------------LAST RUN P' patientNum ' ---------------------------------'])


%% ESS options
% MEIGO OPTIONS I (COMMON TO ALL SOLVERS):
opts.ndiverse     = 'auto'; %100; %500; %5; %
opts.maxtime      = 3600; % MAX-Time of optmization, i.e how long the optimization will last %3600=1h
opts.maxeval      = 1e8; % max number of evals, i.e cost function calls
opts.log_var      = [];  %skip this

opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; % local solver, fmincon works in my experience best
opts.local.finish = opts.local.solver; %uses the local solver to check the best p-vector
opts.local.bestx = 0; % read the documentation, think it's best to leave at zero for now.
problem.f   = 'costFunction_HEALTH'; % %name of cost-function
opts.iterprint = 0;

% MEIGO OPTIONS II (FOR ESS AND MULTISTART):
opts.local.iterprint = 0; % prints what going on during optimization

% MEIGO OPTIONS III (FOR ESS ONLY):
opts.dim_refset   = 'auto'; % leave to auto for now

% OPTIONS AUTOMATICALLY SET AS A RESULT OF PREVIOUS OPTIONS:
if(strcmp(opts.local.solver,'fmincon'))
    opts.local.use_gradient_for_finish = 1; %DW: provide gradient to fmincon
else
    opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
end
opts.local.check_gradient_for_finish = 0; %DW: gradient checker


%% Load parameters and data
resultsfolder = [pathbase filesep 'Parameters' filesep 'ESS'];
[~,data,~,~,constants,paramNames,constantsNames,ynames,xnames,simOptions,ind,origParamvalues] = setup_simulations_HEALTH({patientNum},resultsfolder,1);
% settings for the cost function:
numHeartBeats = 20;
doPlot = 0;
dispErrors = 0;
%% Set bounds (separate for each subject, since paramdata is individual)
[lbOrig,ubOrig] = loadParamBounds_HEALTH(ind,data.params);
lb = log10(lbOrig);
ub = log10(ubOrig);
problem.x_L       = lb; % essOPT uses a problem structure where crucial information is specified
problem.x_U       = ub;

%% Find startguess
loadsummaryfile=0;
folders = dir(fullfile(pathbase,'Parameters','ESS',['P' patientNum '_*']));
[~,latestfolderInd] = max([folders.datenum]);
foldernames = {folders.name};
folderpaths = {folders.folder};
currentResultsFolder = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd})
[optParam,loadedoptcost,~,loadedconstants] = findBestParams(currentResultsFolder,loadsummaryfile,length(paramNames)); %load from current results folder
fprintf('Estimateparameters: Loaded startguess from current results folder with cost %0.2f\n',loadedoptcost)
if sum(isnan(optParam))>0 % load from a previous optimization if current result is empty
    disp('Current results folder empty')
    [~,latestfolderInd] = max([folders.datenum]);
    latestfolderInd = latestfolderInd+1; %not the latest = this opt, but the one before
    if latestfolderInd > length(folders)
        latestfolderInd = latestfolderInd-2;
    end
    foldernames = {folders.name};
    folderpaths = {folders.folder};
    loadfolderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    [optParam,loadedoptcost,~,loadedconstants] = findBestParams(loadfolderName,loadsummaryfile,length(paramNames));
    fprintf('Estimateparameters: Loaded startguess with cost %0.2f from %s\n',loadedoptcost,loadfolderName)
end

optParam = max(optParam,lbOrig);%make sure the startguess is within the bounds
optParam = min(optParam,ubOrig);%make sure the startguess is within the bounds
allparams = optParam;
ind.estParams = 1:length(allparams);

optParam = log10(optParam);
optParam = max(optParam,lb);%make sure the startguess is within the bounds
optParam = min(optParam,ub);%make sure the startguess is within the bounds
problem.x_0=optParam;

%% Run optimization
format long
format compact
warning('off','all') % AMICI prints error-messages for all simulations (if they fail) so this will fix that annoying orange text form appearing
save(sprintf('%s/testsave-%i%s-%s.mat',currentResultsFolder,loopNumber,'lastrun',dateStart) ,'dateStart')

%% Solve with ESS
optim_algorithm = 'ess'; % 'multistart'; %  'cess'; % ESS IS BEST EVAH
Results = MEIGO(problem,opts,optim_algorithm,constants,allparams,simOptions,numHeartBeats,data,ind,dispErrors,doPlot); % Run the optimization
%% Save results
fitting_speed     = Results.time(end);
bestcost          = Results.fbest;
optParam          = Results.xbest;

% Save the results.
optParam = 10.^(optParam); % scale bestparam back to it's original size
bounds.lb = 10.^lb;
bounds.ub = 10.^ub;
dateEnd = datestr(now,'yymmdd-HHMMSS');
    save(sprintf('%s%sopt-ESS(%.3f)-%i%s-r%i.mat',currentResultsFolder,filesep,bestcost,loopNumber,randomstring,r) ,'optParam','bestcost','bounds','constants','ind','paramNames','simOptions','dateStart','dateEnd','fitting_speed','modelName','s','Results','problem')

warning ('on','all'); % yay error messages again


end
