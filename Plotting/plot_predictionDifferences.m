function [ttestTablePredictions] = plot_predictionDifferences(simulations,patNums,ynames,plotFolderName)
% Compares model predictions of the maximum and minimum pressures in the LV
% and the mean pressure in the LA between the four groups.

%% Setup colors and font size
magmacols = flip(magma(4));
colors4groups=num2cell(magmacols,2);
global fzlarge 
global fzsmall 
global fzsupersmall 
fzlarge = 14;
fzsmall = 10;
fzsupersmall = 8;

%% Prepare for statistical test - calculate prediction variables
% Maximum pressures
units = {'mmHg','mmHg','mmHg'};
predictionNames= {'maxpLV','minpLV','meanpLA'};
predictionNamesNice = {'Maximum LV pressure','Minimum LV pressure','Mean LA pressure'};
ypredNames = {'pLV','pLV','pLA'};
allVariableValues = zeros(length(patNums),length(predictionNames));
for pred = 1:length(predictionNames)
    predictionStruct.(predictionNames{pred}) = zeros(size(patNums));
    for pat = 1:length(patNums)
        predind = strcmp(ynames,ypredNames{pred});
        if strcmp(predictionNames{pred}(1:3),'max')
            predictionStruct.(predictionNames{pred})(pat) = max(simulations{pat}.y(:,predind));
        elseif strcmp(predictionNames{pred}(1:3),'min')
            predictionStruct.(predictionNames{pred})(pat) = min(simulations{pat}.y(:,predind));
        else
            predictionStruct.(predictionNames{pred})(pat) = mean(simulations{pat}.y(:,predind));
        end
        allVariableValues(pat,pred) = predictionStruct.(predictionNames{pred})(pat);
    end
end

%% Run statistical test to compare groups
% Setup groups
[groups] = loadGroupIndexes(patNums);

% Run test
[ttestTablePredictions,statisticsTable,predMean,...
    ~,~,numRejectedHyp,~,~,...
    ~,statdocumenttable] = runStatisticalComparison(...
    'numeric',patNums,predictionStruct,predictionNames,predictionNames);

% Save a table for the statistical summary document
writetable(statdocumenttable,fullfile(plotFolderName,'statisticssummary_predictions.xlsx'))

%% Box plots
plotMedian=1;
%Compare the 4 groups
comptitles = {'HBP: control vs hypertensive', ...
    'HBP: T2D vs hypertensive T2D',...
    'HBP: Hypertensive vs hypertensive & T2D',...
    'HBP: control vs T2D',...
    'HBP: Hypertensive vs T2D',...
    'HBP: control vs hypertensive T2D'};
rejects = cell(size(comptitles));
for c = 1:length(comptitles)
    rejects{c} = statisticsTable.reject{strcmp(comptitles{c},statisticsTable.groupComparisonNames)};
    rejects{c}(isnan(rejects{c})) = 0;
end
correctedrejection4groups = numRejectedHyp.HBP.rBenjH2>0;
predictionNamesNice(correctedrejection4groups) = strcat(predictionNamesNice(correctedrejection4groups),'*');

groups4 = {groups.C_NT_home,groups.T2D_NT_home,groups.C_HT_home,groups.T2D_HT_home};
meanT2DHBP = table2array(predMean(strcmp(predMean.Group,'T2D HBP'),3:length(predictionNames)+2));%2 first are groupName and N
sdT2DHBP = table2array(predMean(strcmp(predMean.Group,'T2D HBP'),3+length(predictionNames):end-2));
meanT2DhighHBP = table2array(predMean(strcmp(predMean.Group,'Hypertensive T2D HBP'),3:length(predictionNames)+2));%2 first are groupName and N
sdT2DhighHBP = table2array(predMean(strcmp(predMean.Group,'Hypertensive T2D HBP'),3+length(predictionNames):end-2));
meanCHBP = table2array(predMean(strcmp(predMean.Group,'Controls HBP'),3:length(predictionNames)+2));%2 first are groupName and N
sdCHBP = table2array(predMean(strcmp(predMean.Group,'Controls HBP'),3+length(predictionNames):end-2));
meanChighHBP = table2array(predMean(strcmp(predMean.Group,'Hypertensive HBP'),3:length(predictionNames)+2));%2 first are groupName and N
sdChighHBP = table2array(predMean(strcmp(predMean.Group,'Hypertensive HBP'),3+length(predictionNames):end-2));
mean4 = {meanCHBP,meanT2DHBP,meanChighHBP,meanT2DhighHBP};
sd4 = {sdCHBP,sdT2DHBP,sdChighHBP,sdT2DhighHBP};
names.xtick = {'C','T2D','HT','T2D+HT'};
names.xtickangle = 0;
names.legend = {'Control','T2D','Hypertensive','T2D & hypertensive'};
names.figname ='predictions_4groups';
plot4groups_variables(1:length(predictionNamesNice),comptitles,groups4,mean4,sd4,statisticsTable,allVariableValues,predictionNamesNice,units,plotMedian,names,colors4groups)

%% Save figures
saveAllFigures(plotFolderName)

end

