function plot_spiderparams(paramNamesPlot,allParamValues,patNums,bounds,ind,plotFolderName,statisticsTable,numRejectedHyp,mediantableParameters)
% Creates two spider plots of parameter values in the 4 groups. 
% If not provided in the input, a statistical test is runned to determine 
% which of the 28 parameters to plot.

%% Set colors
addpath '..\Requirements\spider_plot'

magmacols = magma(4);
colors4groupshbp=flip(magmacols);

%% set parameter values to their actual values used when simulating
allParamValues(ind.onset_LV,:) = allParamValues(ind.onset_LV,:) - 1.5;
allParamValues(ind.onset_LA,:) = 1 + allParamValues(ind.onset_LV,:) - allParamValues(ind.onset_LA,:);

bounds.lb(ind.onset_LV) = bounds.lb(ind.onset_LV) - 1.5;
bounds.ub(ind.onset_LV) = bounds.ub(ind.onset_LV) - 1.5;

lbLA = bounds.lb(ind.onset_LA) ;
bounds.lb(ind.onset_LA) = 1 + bounds.lb(ind.onset_LV) - bounds.ub(ind.onset_LA);
bounds.ub(ind.onset_LA) = 1 + bounds.ub(ind.onset_LV) - lbLA;

paramNamesOrig = paramNamesPlot;
paramNamesPlot = fixUnderscore(paramNamesPlot);

%% physiological order of parameters
PVparams = {'Ppu','Rpu','Cpvc','Lpv','Rpv'};
LAparams = {'Emax_LA','Emin_LA','k_diast_LA','k_syst_LA','m1_LA','m2_LA','onset_LA','V0_LA'};
MVparams = {'Lmv','Rmv'};
LVparams= {'Emax_LV','Emin_LV','k_diast_LV','k_syst_LV','m1_LV','m2_LV','onset_LV','V0_LV'};
AVparams = {'Lav','ElCo'};
Aortaparams = {'Caa','Lao','Rao'};
Perparams = {'Rtot','Ctot','SBPdiff'};
paramNamesPhysOrder = [PVparams,LAparams,MVparams,LVparams,AVparams,Aortaparams,Perparams];

%% functional order of parameters
hrparamsLA = {'Emax_LA','Emin_LA','k_diast_LA','k_syst_LA','m1_LA','m2_LA','onset_LA'};
hrparamsLV = {'Emax_LV','Emin_LV','k_diast_LV','k_syst_LV','m1_LV','m2_LV','onset_LV'};
Rparams = {'Rpu','Rpv','Rmv','Rao','Rtot'};
Cparams = {'Cpvc','Caa','Ctot'};
Lparams = {'Lpv','Lmv','Lav','Lao'};
Vparams = {'V0_LA','V0_LV'};
restparams = {'Ppu','ElCo','SBPdiff'};
paramNamesFunctionOrder = [hrparamsLA,hrparamsLV,Rparams,Cparams,Lparams,Vparams,restparams];


%% Calculate statistics
if nargin < 7
    % Structure data
    for param = 1:length(paramNamesOrig)
        paramvaluesStruct.(paramNamesOrig{param}) = allParamValues(param,:);
    end
    % Run statistical tests
    [~,statisticsTable,meantableParameters,~,~,numRejectedHyp,~,mediantableParameters,~] = runStatisticalComparison('numeric',patNums,paramvaluesStruct,paramNamesOrig,paramNamesOrig);
end

%number of subjects in each group
nControlHypHBP = table2array(mediantableParameters(strcmp(mediantableParameters.Group,'Hypertensive HBP'),2));%meantableParameters
nControlNormHBP = table2array(mediantableParameters(strcmp(mediantableParameters.Group,'Controls HBP'),2));
nT2DHypHBP = table2array(mediantableParameters(strcmp(mediantableParameters.Group,'Hypertensive T2D HBP'),2));
nT2DNormHBP = table2array(mediantableParameters(strcmp(mediantableParameters.Group,'T2D HBP'),2));

% Add * to the predictions that contain significant differences
correctedRejection4groupsHBP = numRejectedHyp.HBP.rBenjH2>0;
paramNamesPlotHBP = paramNamesPlot;
paramNamesPlotHBP(correctedRejection4groupsHBP) = strcat(paramNamesPlot(correctedRejection4groupsHBP),'*');

% extract the median values
medianT2DnormalHBP = table2array(mediantableParameters(strcmp(mediantableParameters.Group,'T2D HBP'),3:length(paramNamesPlot)+2));%2 first are groupName and N
medianT2DhighHBP = table2array(mediantableParameters(strcmp(mediantableParameters.Group,'Hypertensive T2D HBP'),3:length(paramNamesPlot)+2));%2 first are groupName and N
medianCnormalHBP = table2array(mediantableParameters(strcmp(mediantableParameters.Group,'Controls HBP'),3:length(paramNamesPlot)+2));%2 first are groupName and N
medianChighHBP = table2array(mediantableParameters(strcmp(mediantableParameters.Group,'Hypertensive HBP'),3:length(paramNamesPlot)+2));%2 first are groupName and N

%Index of the parameters to plot, with any p<0.05
rejectHTdiabetes_HBP = statisticsTable.reject{strcmp('HBP: T2D vs hypertensive T2D',statisticsTable.groupComparisonNames)};
rejectHTcontrol_HBP = statisticsTable.reject{strcmp('HBP: control vs hypertensive',statisticsTable.groupComparisonNames)};
rejectHTcontrolT2D_HBP = statisticsTable.reject{strcmp('HBP: Hypertensive vs hypertensive & T2D',statisticsTable.groupComparisonNames)};
rejectNTconrolT2D_HBP = statisticsTable.reject{strcmp('HBP: control vs T2D',statisticsTable.groupComparisonNames)};
rejectControlT2DHT_HBP = statisticsTable.reject{strcmp('HBP: control vs hypertensive T2D',statisticsTable.groupComparisonNames)};
rejectHTT2D_HBP = statisticsTable.reject{strcmp('HBP: Hypertensive vs T2D',statisticsTable.groupComparisonNames)};
indSign4hbp = find(or(rejectNTconrolT2D_HBP,or(rejectHTcontrolT2D_HBP,or(rejectHTcontrol_HBP,or(rejectHTdiabetes_HBP,or(rejectControlT2DHT_HBP,rejectHTT2D_HBP))))));

%% PLOT
fzsmall = 10;
P2 = [medianCnormalHBP(indSign4hbp);medianT2DnormalHBP(indSign4hbp);medianChighHBP(indSign4hbp);medianT2DhighHBP(indSign4hbp)];
names2 = paramNamesPlotHBP(indSign4hbp);
% sort in physiological order
[~,nameind] = ismember(paramNamesPhysOrder,paramNamesOrig(indSign4hbp));
nameind(nameind==0) = [];
names2 = names2(nameind);
P2 = P2(:,nameind);
%shift 2 steps
P2 = circshift(P2,[0,-2]);
precision = spider_findprecision(P2);
names2 = circshift(names2,-2);
figure('Name','Spider_4groups','Visible','off') 
set(gcf,'Color','white')
xdim_CM = 25; %A4: 21x29.7
ydim_CM = 14;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiles=tiledlayout(1,2,'Padding','compact','TileSpacing','loose');
nexttile
spider_plot(P2,...
    'AxesLabels', names2,...
    'AxesInterval', 1,...
    'AxesPrecision',precision,...
    'FillOption', {'on', 'on', 'on','on'},...
    'FillTransparency', [0.12, 0.12, 0.12,0.12],...
    'Color', colors4groupshbp,...
    'LineStyle', {'-', '-', '-','-'},...
    'AxesFontSize', fzsmall,...
    'LabelFontSize', fzsmall,...
    'AxesLabelsRotate', 'on',...
    'AxesLabelsEdge', 'none',...
    'AxesLabelsOffset', 0.1,...
    'AxesLimits', [min(P2);max(P2)],...
    'LineWidth',2,'MarkerSize',1);
set(gca,'FontSize',fzsmall)

%  functional
P2 = [medianCnormalHBP(indSign4hbp);medianT2DnormalHBP(indSign4hbp);medianChighHBP(indSign4hbp);medianT2DhighHBP(indSign4hbp)];
names2 = paramNamesPlotHBP(indSign4hbp);
% sort in functional order
[~,nameind] = ismember(paramNamesFunctionOrder,paramNamesOrig(indSign4hbp));
nameind(nameind==0) = [];
names2 = names2(nameind);
P2 = P2(:,nameind);
%shift 2 steps
P2 = circshift(P2,[0,-2]);
precision = spider_findprecision(P2);
names2 = circshift(names2,-2);
nexttile
spider_plot(P2,...
    'AxesLabels', names2,...
    'AxesInterval', 1,...
    'AxesPrecision',precision,...
    'FillOption', {'on', 'on', 'on','on'},...
    'FillTransparency', [0.12, 0.12, 0.12,0.12],...
    'Color', colors4groupshbp,...
    'LineStyle', {'-', '-', '-','-'},...
    'AxesFontSize', fzsmall,...
    'LabelFontSize', fzsmall,...
    'AxesLabelsRotate', 'on',...
    'AxesLabelsEdge', 'none',...
    'AxesLabelsOffset', 0.1,...
    'AxesLimits', [min(P2);max(P2)],...
    'LineWidth',2,'MarkerSize',1);
set(gca,'FontSize',fzsmall) 
legend({sprintf('Control (N = %d)',nControlNormHBP),sprintf('T2D (N = %d)',nT2DNormHBP),sprintf('Hypertensive (N = %d)',nControlHypHBP),sprintf('T2D and hypertensive (N = %d)',nT2DHypHBP)},...
    'FontSize',fzsmall,'Position',[0.4135478003012 0.0984478719265485 0.275818639798489 0.117941176470588]);
legend('boxoff')

%% Save figures
saveAllFigures(plotFolderName)

end