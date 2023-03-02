function [] = createFigures_hypertensionT2D_predictions()
% Simulate and plot subject-specific predictions with uncertainty
disp('*** Re-creating results for subject-specific predictions with uncertainty. This might take several minutes. ***')

%% Setup
close all
clear

% Create a folder to save the resulting figures in
thispath = split(pwd,'cardiovascularavatar_T2D_HT');
plotFolderName = fullfile(thispath{1},'cardiovascularavatar_T2D_HT',['Resultfigures_predictions_',datestr(now,'yymmdd_HHMM')]);
mkdir(plotFolderName)

% Add dependencies to matlab path
addpath '.\Requirements\matplotlib'
addpath '.\Tools'
addpath '.\Simulation'

resultfolder = 'Parameters/ESS';

% Load which subjects to simulate
load('./Parameters/bestpatients.mat','pbestT2Dh','pbestT2Dn','pbestCn','pbestCh')
pC= pbestCn;
pT2D=pbestT2Dn;
pHT= pbestCh;
pHTT2D = pbestT2Dh;

%% Load paramvalues and simulate for one person from each group
%This takes several minutes
disp('*** Simulating uncertainty for the control subject... *** ')
[~,data_CNT,~,bestparams_CNT,constants_CNT,paramNames,constantsNames,ynames,xnames,options,inds_CNT,~,~,~,~,~,paramuncertainty_CNT] = setup_simulations_HEALTH({pC},resultfolder,0);
[~,minmaxSimCNT] = simulatePredictionUncertainty({pC},[],[],paramuncertainty_CNT.allokParams,bestparams_CNT,constants_CNT,inds_CNT,data_CNT,options);

disp('*** Simulating uncertainty for the T2D subject... *** ')
[~,data_T2DNT,~,bestparams_T2DNT,constants_T2DNT,~,~,~,~,~,inds_T2DNT,~,~,~,~,~,paramuncertainty_T2DNT] = setup_simulations_HEALTH({pT2D},resultfolder,0);
[~,minmaxSimT2DNT] = simulatePredictionUncertainty({pT2D},[],[],paramuncertainty_T2DNT.allokParams,bestparams_T2DNT,constants_T2DNT,inds_T2DNT,data_T2DNT,options);

disp('*** Simulating uncertainty for the HT subject... *** ')
[~,data_CHT,~,bestparams_CHT,constants_CHT,~,~,~,~,~,inds_CHT,~,~,~,~,~,paramuncertainty_CHT] = setup_simulations_HEALTH({pHT},resultfolder,0);
[~,minmaxSimCHT] = simulatePredictionUncertainty({pHT},[],[],paramuncertainty_CHT.allokParams,bestparams_CHT,constants_CHT,inds_CHT,data_CHT,options);

disp('*** Simulating uncertainty for the HT+T2D subject... *** ')
[~,data_T2DHT,~,bestparams_T2DHT,constants_T2DHT,~,~,~,~,~,inds_T2DHT,~,~,~,~,~,paramuncertainty_T2DHT] = setup_simulations_HEALTH({pHTT2D},resultfolder,0);
[~,minmaxSimT2DHT] = simulatePredictionUncertainty({pHTT2D},[],[],paramuncertainty_T2DHT.allokParams,bestparams_T2DHT,constants_T2DHT,inds_T2DHT,data_T2DHT,options);

%% Plot the simulations of the subjects from the 4 groups
magmacols = flip(magma(4));
colors = num2cell(magmacols,2);
simulationNames = {'Control','T2D','Hypertensive',sprintf('T2D &\nhypertensive')}';
simulations = {minmaxSimCNT,minmaxSimT2DNT,minmaxSimCHT,minmaxSimT2DHT};
plotPredictionUncertaintyNice(simulations,simulationNames,xnames,ynames,'Subjects from 4 groups',colors)

%% Simulate prediction: control to reduced m2LV in HT+T2D
% Set parameter values for the prediction
allParams = paramuncertainty_CNT.allokParams;%parameters for the control subject
allParams2 = allParams;
bestparams2 = bestparams_CNT;

% set which parameters to change
paramsToChange = {'m2_LV'};
paramsToChange_inds = [find(strcmp(paramNames,'m2_LV'))];

%set the change in m2LV so that the value in the best parameter set is changed to 20
fractionChange = 20/bestparams_CNT(paramsToChange_inds(1),1);

% change of m2LV in all parameter sets based on the same fraction
for pat = 1:size(allParams)
    for param = 1:length(paramsToChange)
        allParams2{pat}(:,paramsToChange_inds(param)) = allParams{pat}(:,paramsToChange_inds(param))*fractionChange(param);
        bestparams2(paramsToChange_inds(param),pat) = bestparams_CNT(paramsToChange_inds(param),pat)*fractionChange(param);
    end
end

% Save information about the changed parameters in the prediciton
hr_range_control = [bestparams_CNT(paramsToChange_inds(1),1);min(allParams{1}(:,paramsToChange_inds(1))); max(allParams{1}(:,paramsToChange_inds(1)))];
hr_range_prediction = [bestparams2(paramsToChange_inds(1),1);min(allParams2{1}(:,paramsToChange_inds(1))); max(allParams2{1}(:,paramsToChange_inds(1)))];
HRpredictiontable = table(hr_range_control,hr_range_prediction,'Variablenames',{'m2_L_V Control','m2_L_V Prediction'},'Rownames',{'Best value','Minimum value','Maximum value'});
writetable(HRpredictiontable,fullfile(plotFolderName,'m2predictiontable.xlsx'),"WriteRowNames",1)

% Simulate the predicted diastolic function. This takes several minutes.
disp('*** Simulating the prediction of reducted diastolic function... *** ')
constants2 = constants_CNT;
[~,minmax_HTT2D_m2] = simulatePredictionUncertainty({pC},[],[],allParams2,bestparams2,constants2,inds_CNT,data_CNT,options);

% Plot results
magmacols = flip(magma(4));
colors = num2cell(magmacols([1,4],:),2);
simulationNamesShort = {'Healthy control',sprintf('Predicted\ndiastolic dysfunction')}';
simulations = {minmaxSimCNT,minmax_HTT2D_m2};
plotPredictionUncertaintyNice(simulations,simulationNamesShort,xnames,ynames,'Predictions healthy to t2d+ht - m2change',colors)

%% Save figures
saveAllFigures(plotFolderName)

fprintf('***\nDONE. The figures and tables are saved in %s.\nTo re-create the rest of the figures and tables, run createFigures_hypertensionT2D.\n***\n',plotFolderName)

end
