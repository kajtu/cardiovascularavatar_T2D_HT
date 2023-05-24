%% Main script 
% avatar model applied on data from HEALTH study to evalute hemodynamics in T2D and hypertension
% Script to simulate several individuals and plot the results

close all
clear

disp('*** Re-creating results for model training, validation, and group comparisons. This might take several minutes. ***')


%% Set up: Load data, optimization results, and settings
% Add dependencies to matlab path
addpath(['.' filesep 'Requirements' filesep 'matplotlib']) 
addpath(['.' filesep 'Tools'])
addpath(['.' filesep 'Simulation'])

% Load data and optimizaton results for all subjects
pathbase = split(pwd,'cardiovascularavatar_T2D_HT');
pathbase = [pathbase{1},filesep, 'cardiovascularavatar_T2D_HT'];
resultsFolder = fullfile(pathbase,'Parameters','ESS');
[patNums, data,extradata,bestparamValues,constants,paramNames,...
    constantsNames,ynames,xnames,options,inds,origParamvalues,...
    units,loadedcosts,meanParams,medianParams,paramuncertainty] =setup_simulations_HEALTH([],resultsFolder,0);

% Create a folder to save the resulting figures in
thispath = split(pwd,'cardiovascularavatar_T2D_HT');
plotFolderName = fullfile(thispath{1},'cardiovascularavatar_T2D_HT',['Resultfigures_',datestr(now,'yymmdd_HHMM')]);
mkdir(plotFolderName)

%% Simulate all subjects
datanames = {'MV','AV','AC','PV'};
numflows = length(datanames);
maxSVdiffData = zeros(size(patNums));
simLast = cell(1,length(patNums));
simCostFunc = cell(1,length(patNums));
rmse = zeros(size(patNums));rmseFlow = zeros(size(patNums));
rmsePressure=zeros(size(patNums));rmseParam=zeros(size(patNums));
SVdiffstot = zeros(size(patNums));
cost = zeros(size(patNums));
fittedcost = zeros(size(patNums));
percenterrorFlow = zeros(1,length(patNums));
percenterrorParamsansBP = zeros(1,length(patNums));
percenterrorTotal=zeros(1,length(patNums));
SVdiffs = zeros(length(patNums),numflows);
warning('off','all');
numHeartBeats = 20; %maximum number of heartbeats to simulate
for p = 1:length(patNums)
    fprintf('Simulating subject %d/%d...\n',p,length(patNums))
    % Calculate cost and fit to data
    inds{p}.estParams = 1:length(bestparamValues(:,p));
    [cost(p),simCostFunc{p},~,SVdiffs(p,:),~,rmses,~,percenterror,MAPE(p)] = costFunction_HEALTH_allcalculations(log10(bestparamValues(:,p)),constants(:,p),bestparamValues(:,p),options,numHeartBeats,data{p},inds{p},datanames);
    fittedcost(p) = costFunction_HEALTH(log10(bestparamValues(:,p)),constants(:,p),bestparamValues(:,p),options,numHeartBeats,data{p},inds{p},1,0);
    SVdiffstot(p) = sum(abs(SVdiffs(p,:)));
    rmse(p) = sum(rmses);
    rmseFlow(p) = rmses(numflows+1);
    rmsePressure(p) = rmses(numflows+2);
    rmseParam(p) = rmses(numflows+3);
    percenterrorFlow(p) = percenterror.flow;
    percenterrorParamsansBP(p) = percenterror.params;
    percenterrorTotal(p)=percenterror.sumall;
    
    %calculate difference in stroke volume (SV) in the data
    dataSV = zeros(size(datanames));
    for d = 1:length(datanames)
        dataSV(d) = trapz(data{p}.time,data{p}.(datanames{d})(:,2));
    end
    maxSVdiffData(p) = max(dataSV) - min(dataSV);
    
    % Simulate
    T = constants(inds{p}.T,p);
    step = 0.001;
    simtime = sort([data{p}.time,0:step:T]);
    simtime = unique(simtime);
    options.x0 = data{p}.IC;
    [~,simLast{p}] = simulate_avatarHEALTH(bestparamValues(:,p),constants(:,p),options,numHeartBeats,inds{p},simtime);
end
warning('on','all');


%% Create a table of costs and fit to data (Table 3)
costs = {fittedcost',rmseFlow,rmsePressure,rmseParam,SVdiffstot,maxSVdiffData,percenterrorTotal,MAPE}';
means = zeros(size(costs));mins = zeros(size(costs));maxs = zeros(size(costs));medians = zeros(size(costs));
for c = 1:length(costs)
    means(c) = round(mean(costs{c}),1);
    medians(c) = round(median(costs{c}),1);
    mins(c) = round(min(costs{c}),1);
    maxs(c) = round(max(costs{c}),1);
end
costTableMean = table(means,medians,mins,maxs,'VariableNames',...
    {'Mean','Median','Min','Max'},'RowNames',...
    {'Cost','RMSE blood flow (ml/s)','RMSE pressure (mmHg)','RMSE parameters (-)'...
    'Sum of |SV simulation - SV data| (ml)','Max difference within SV data (ml)',...
    'Total error (%)','MAPE (%)'});

%Save cost table
writetable(costTableMean,fullfile(plotFolderName,'costtable.xlsx'),"WriteRowNames",1)

%% Plot results
% Load parameter bounds
[bounds.lb,bounds.ub] = loadParamBounds_HEALTH(inds{1},[]);
lbs = zeros(length(patNums),length(paramNames));
ubs=lbs;
for p = 1:length(patNums)
[lbs(p,:),ubs(p,:)] = loadParamBounds_HEALTH(inds{p},data{p}.params);
end

% Create table with the parameter bounds (for supplementary)
varnames = {'min lb','max lb','mean lb','min ub','max ub','mean ub'};
lbubtable = table(min(lbs)',max(lbs)',mean(lbs)',min(ubs)',max(ubs)',mean(ubs)','Rownames',paramNames,'Variablenames',varnames);
writetable(lbubtable,fullfile(plotFolderName,'parameter_bounds.xlsx'),"WriteRowNames",1)

% Create a table with constants
constantTable = table(min(constants(5:end-2,:),[],2),max(constants(5:end-2,:),[],2),'RowNames',constantsNames(5:end-2),'Variablenames',{'Min value','Max value'});
writetable(constantTable,fullfile(plotFolderName,'constants.xlsx'),"WriteRowNames",1)

% Load (or simulate)the simulation uncertainty for all subjects
cd Simulation
loadResults=1;
[minmaxSim] = findSimulationUncertainty(patNums,paramuncertainty.allokParams,constants,inds,data,options,bestparamValues,1,loadResults);
cd ..

% Plot the figures
cd Plotting
plotmedian=1;

disp('Plotting fit to data...')
plot_fitToData_uncertainty(minmaxSim,patNums,bestparamValues,paramuncertainty,paramNames,data,inds,fittedcost,plotFolderName);

disp('Plotting fit to data: validation...')
plot_validation_uncertainty(minmaxSim,simLast,data,extradata,patNums,xnames,plotFolderName)

disp('Plotting box plot of parameter differences...')
[testTableParameters,meantableParameters,mediantableParameters,...
    statisticsTableParams,numRejectedHypParameters] = plot_parameterDifferences(paramNames,medianParams,patNums,bounds,inds{p},units,plotFolderName);

disp('Plotting spider plots of parameter differences...')
plot_spiderparams(paramNames,medianParams,patNums,bounds,inds{p},plotFolderName,statisticsTableParams,numRejectedHypParameters,mediantableParameters);

disp('Plotting clusters and PCA results...')
[clusterindex,groupedinbothFinal] = clusterCalculation(medianParams,constants,2,paramNames,constantsNames,patNums,plotFolderName);

disp('Plotting box plot of prediction differences...')
[testTablePredictions] = plot_predictionDifferences(simLast,patNums,ynames,plotFolderName);

disp('Calculating sensitivity to blood pressure measurement...')
load(fullfile(pathbase,'Parameters','bestpatients.mat'),'pbestT2Dh')
pat = strcmp(patNums, pbestT2Dh);
doPlot = 0;
plot_BPsensitivity(data{pat},simLast{pat},medianParams(:,pat),constants(:,pat),inds{pat},ynames,xnames,options,plotFolderName,paramNames,constantsNames,doPlot)

cd ..

%% Calculate subject specific standard deviations vs group standard deviations/interquartile ranges (supplementary)
allParamSDs = zeros(length(paramuncertainty.allokParams),length(paramNames));
allParamIQRs = zeros(length(paramuncertainty.allokParams),length(paramNames));
for p = 1:length(paramuncertainty.allokParams)
    allParamSDs(p,:)=std(paramuncertainty.allokParams{p});
    allParamIQRs(p,:)=iqr(paramuncertainty.allokParams{p});
end

[groups] = loadGroupIndexes(patNums);
groupnames = {'Controls HBP','T2D HBP','Hypertensive HBP','Hypertensive T2D HBP'};
groupinds = {groups.C_NT_home,groups.T2D_NT_home, groups.C_HT_home,groups.T2D_HT_home};
groupSD = zeros(1,length(paramNames));
pooledSubjectSD = zeros(1,length(paramNames));
groupIQR = zeros(1,length(paramNames));
medianIQRvsgroup = zeros(1,length(paramNames));
for param = 1:length(paramNames)
    for g = 1:length(groupinds)
        groupSD(param,g) = meantableParameters.(['sd' paramNames{param}])(strcmp(meantableParameters.Group,groupnames{g}));
        pooledSubjectSD(param,g) = sqrt(sum(allParamSDs(groupinds{g},param).^2)/length(paramNames));
        
        groupIQR(param,g) = mediantableParameters.(['iqr' paramNames{param}])(strcmp(meantableParameters.Group,groupnames{g}));
        medianIQRvsgroup(param,g) = median(100.*allParamIQRs(groupinds{g},param)./groupIQR(param,g));
    end
end

individual_vs_groupSD=table(100*pooledSubjectSD./groupSD,'Rownames',paramNames','Variablenames',{'Pooled SD % of group SD'});
individual_vs_groupIQR=table(medianIQRvsgroup(:,1),medianIQRvsgroup(:,2),medianIQRvsgroup(:,3),medianIQRvsgroup(:,4),'Rownames',paramNames','Variablenames',{'Controls','T2D','Hypertensive','Hypertensive T2D'});

%save the tables
writetable(individual_vs_groupIQR,fullfile(plotFolderName,'individual_vs_groupIQR.xlsx'),"WriteRowNames",1)
writetable(individual_vs_groupSD,fullfile(plotFolderName,'individual_vs_groupSD.xlsx'),"WriteRowNames",1)

%% Save other results
disp('Saving tables with group differences...')
%create excel files with statistics results
comparisons4hbb= {'HBP: control vs hypertensive', ...
    'HBP: T2D vs hypertensive T2D',...
    'HBP: Hypertensive vs hypertensive & T2D',...
    'HBP: control vs T2D',...
    'HBP: control vs hypertensive T2D',...
    'HBP: Hypertensive vs T2D'};
groups4hbp = {'Controls HBP','T2D HBP','Hypertensive HBP','Hypertensive T2D HBP'};
parameters_4groups_hbp = createResultsTable(testTableParameters, groups4hbp, comparisons4hbb,1,fullfile(plotFolderName,'parameters_4groups_hbp'));
predictions_4groups_hbp = createResultsTable(testTablePredictions, groups4hbp, comparisons4hbb,1,fullfile(plotFolderName,'predictions_4groups_hbp'));


fprintf('***\nDONE. The figures and tables are saved in %s.\nTo re-create figures with subject-specific predicitons, run createFigures_hypertensionT2D_predictions.\n***\n',plotFolderName)
