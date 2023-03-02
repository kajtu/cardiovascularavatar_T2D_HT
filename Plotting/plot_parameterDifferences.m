function [ttestTableParameters,meantableParameters,medianTableParameters,...
    statisticsTable,numRejectedHyp] = plot_parameterDifferences(paramNames,allParamValues,patNums,bounds,ind,units,plotFolderName)

%% Set colors
global colnormal
global colhigh
global colT2Dnormal
global colT2Dhigh

magmacols = magma(4);
colors4groupsHBP  = flip(magmacols);
colnormal = colors4groupsHBP(1,:);
colT2Dnormal = colors4groupsHBP(2,:);
colhigh = colors4groupsHBP(3,:);
colT2Dhigh = colors4groupsHBP(4,:);
colors4groupsHBP  = {colnormal,colT2Dnormal,colhigh,colT2Dhigh};

%Set font size in figures
global fzlarge
global fzsmall
fzlarge = 11;
fzsmall = 8;

%% set parameter values to their actual values (as when simulating)
allParamValues(ind.onset_LV,:) = allParamValues(ind.onset_LV,:) - 1.5;
allParamValues(ind.onset_LA,:) = 1 + allParamValues(ind.onset_LV,:) - allParamValues(ind.onset_LA,:);

bounds.lb(ind.onset_LV) = bounds.lb(ind.onset_LV) - 1.5;
bounds.ub(ind.onset_LV) = bounds.ub(ind.onset_LV) - 1.5;

lbLA = bounds.lb(ind.onset_LA) ;
bounds.lb(ind.onset_LA) = 1 + bounds.lb(ind.onset_LV) - bounds.ub(ind.onset_LA);
bounds.ub(ind.onset_LA) = 1 + bounds.ub(ind.onset_LV) - lbLA;

paramNamesOrig = paramNames;
paramNames =  removeUnderscore(paramNames,' ');

%% Load groups
[groups] = loadGroupIndexes(patNums);

%% Run statistical tests
for param = 1:length(paramNamesOrig)
    paramvaluesStruct.(paramNamesOrig{param}) = allParamValues(param,:);
end
[ttestTableParameters,statisticsTable,meantableParameters,~,~,numRejectedHyp,~,medianTableParameters,~,statdocumenttable] = runStatisticalComparison('numeric',patNums,paramvaluesStruct,paramNamesOrig,paramNamesOrig);

%save table for statistical summary document
writetable(statdocumenttable,fullfile(plotFolderName,'statisticssummary_params.xlsx'))


%% Create box plots of the differences between the 4 groups
correctedRejection4groupsHBP = numRejectedHyp.HBP.rBenjH2>0;
paramNamesMarkedHBP = paramNames;
paramNamesMarkedHBP(correctedRejection4groupsHBP) = strcat(paramNamesMarkedHBP(correctedRejection4groupsHBP),'*');
comptitles = {'HBP: control vs hypertensive', ...
    'HBP: T2D vs hypertensive T2D',...
    'HBP: Hypertensive vs hypertensive & T2D',...
    'HBP: control vs T2D',...
    'HBP: Hypertensive vs T2D',...
    'HBP: control vs hypertensive T2D'};
groups4HBP = {groups.C_NT_home,groups.T2D_NT_home,groups.C_HT_home,groups.T2D_HT_home};

names.xtick = {'C','T2D','HT','T2D+HT'};
names.xtickangle = 0;
names.legend = {'Control','T2D','Hypertensive','T2D & hypertensive'};
names.figname = 'Params_Groups_significant_4groups_HBP';
plotMedian = 1;
plot4groups(comptitles,[],groups4HBP,[],[],...
    statisticsTable,allParamValues,paramNamesMarkedHBP,...
    bounds,units,plotMedian,names,colors4groupsHBP)


%% Save figures
saveAllFigures(plotFolderName)

end