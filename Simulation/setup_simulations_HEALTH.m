function [patNums,data,extradata,paramValues,constants,...
    paramNames,constantsNames,ynames,xnames,simulationoptions,inds,...
    origParamvalues,units,loadedcosts,meanParams,medianParams,paramuncertainty] =setup_simulations_HEALTH(patNumsOrig,resultsFolder,doOptimization)

loadMCMC=1;
dataFolder = 'Data';
if ~doOptimization 
    Setup('',0)
end
%names in the model
ynames = {'P_Aortic','Pperipheral','pressGrad MV','P Dmv','mv open','P Dav','av open','pressgrad AV','Ela','Elv','Qcaa','Qpc','P pulmvein','Qpvc','qLA','qLV','pLA','pLV','aaCorr','avCorr','mvCorr','pvCorr','Vla','Vlv','P_Brachial'};
xnames = {'Ppvc','Qpv','Vla','Qmv','Vlv','Qav','Paa','Qaa','Ppc','mv_open','av_open'};

%% select patients (if not specified)
if isempty(patNumsOrig)
    folders = dir([dataFolder '/*_P*']);
    names = {folders.name};
    patNumsOrig = cell(1,length(names));
    for n = 1:length(names)
        patNumsOrig{n} = names{n}(end-5:end-4);
        if strcmp(patNumsOrig{n}(1),'P')
            patNumsOrig{n} = patNumsOrig{n}(2);
        end
    end
end
numPatients = length(patNumsOrig);

%% Load data and parameters
origParamvalues= cell(1,length(patNumsOrig));
data = cell(1,length(patNumsOrig));
extradata = cell(1,length(patNumsOrig));
paramdata = cell(1,length(patNumsOrig));
nparams = 28;
nconstants = 17;
constants = nan(nconstants,length(patNumsOrig));
paramValues = nan(nparams,length(patNumsOrig));
meanParams = nan(nparams,length(patNumsOrig));
medianParams = nan(nparams,length(patNumsOrig));

loadedcosts =zeros(length(patNumsOrig),1);

patNumsNoSim = [];
patNumsNoSimNums = [];
patNums = patNumsOrig;
inds = cell(length(patNumsOrig),1);
for p = 1:length(patNumsOrig)
    patientNum = patNums{p};

    % Format subject number if needed
    if strcmp(patientNum(end),'_')
        patNums{p} = patientNum(1);
    end
    
    % Load data
    [data{p},extradata{p},paramdata{p}] = loadData_HEALTH(patientNum);
    [origParamvalues{p},constants(:,p),paramNames,constantsNames,units,inds{p},data{p}.params] = loadParameters_HEALTH(paramdata{p});

    % Load ESS opt parameters
    if doOptimization || ~loadMCMC
        % find the latest patientfolder with optimizations
        loadresultsfolder = fullfile(resultsFolder, ['P' patientNum '_*']);
        folders = dir(loadresultsfolder);
        if isempty(folders)
            disp(['SETUP: OBS couldnt load opt parameters for p' patientNum, ' (no prev folder), setting to litterature+data values (prev folder: ' loadresultsfolder ')'])
            patNumsNoSim = [patNumsNoSim,p];
            patNumsNoSimNums = [patNumsNoSimNums,{patientNum}];
            paramValues(:,p) = origParamvalues{p};
            loadedcosts(p) = 1e10;
        else
            [~,latestfolderInd] = max([folders.datenum]);
            foldernames = {folders.name};
            folderpaths = {folders.folder};
            
            folderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
            if isempty(dir([folderName,filesep, 'opt-*'])) && doOptimization && length(foldernames) > 1%new empty folder created, so take the next latest one if it exists
                folderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd-1});
            end
            
            if doOptimization
                loadsummaryfile = 0;
            else
                loadsummaryfile = 1;
            end
            loadmean=1;
            [param,loadedcosts(p),~,loadedconstants,meanparam,medianparam] = findBestParams(folderName,loadsummaryfile,length(paramNames),loadmean);
            
            if isnan(param)
                patNumsNoSim = [patNumsNoSim,p];
                patNumsNoSimNums = [patNumsNoSimNums,{patientNum}];
                disp(['SETUP: OBS couldnt load opt parameters for p' patientNum ' (empty folder), setting to litterature+data values'])
                paramValues(:,p) = origParamvalues{p};
            else
                paramValues(:,p) = param;
                meanParams(:,p) = meanparam;
                medianParams(:,p) = medianparam;
            end
        end
         fprintf('SETUP: Loaded parameter values with cost %0.2f for p%s\n',loadedcosts(p),patientNum)
    end
end

if loadMCMC && ~doOptimization
    fprintf('SETUP: Loading MCMC and ESS parameters...\n')
    [patNums,paramValues,loadedcosts,meanParams,paramuncertainty.allokParams,...
    paramuncertainty.allokCosts,paramuncertainty.minValuesparams,...
    paramuncertainty.maxValuesparams,paramuncertainty.numincluded,...
    data,inds,constants,medianParams] = findAllParams(resultsFolder,patNumsOrig,patNums,paramNames,data,inds,constants);
    fprintf('SETUP: Loaded paramvalues from both MCMC and ESS optimizations for %d subjects\n',length(patNums))
else
    paramuncertainty = NaN;
end


%% Simulation settings
% sensi: sensitivity order - 0 = no calculation of sensitivities (default)
% maxsteps: maximum number of integration steps
simulationoptions = amioption('sensi',0,'maxsteps',1e4);
% set simulation tolerances
simulationoptions.atol = 1e-16;
simulationoptions.rtol = 1e-8;

% Calculate initial conditions for each subject based on the parameter values and data
ICs = zeros(length(patNums),length(xnames));
for p = 1:length(patNums)
    [ICs(p,:),~] = iccalcs(paramValues(:,p),constants(:,p),paramNames,constantsNames,data{p},inds{p});
    data{p}.IC = ICs(p,:)';
end

%% Format the function output
if numPatients == 1 && doOptimization
    patNums = patNums{1};
    data = data{1};
    extradata= extradata{1};
    inds= inds{1};
    origParamvalues = origParamvalues{1};
else
    %exclude patients for whom parameters could not be loaded with findBestParams
    patNums(patNumsNoSim)= [];
    data(patNumsNoSim)= [];
    extradata(patNumsNoSim)= [];
    constants(:,patNumsNoSim)= [];
    paramValues(:,patNumsNoSim) = [];
    origParamvalues(patNumsNoSim) = [];
end

end