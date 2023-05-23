function plot_BPsensitivity(dataMiddle,simMiddle,medianParamsMiddle,constants,patindex,ynames,xnames,options,plotFolderName,paramNames,constantsNames,doPlot)
%This function calculates model sensitivity to measured blood pressure based on
%parameter estimations for different data values of SBP and DBP (+-16) for the
%T2D+HT subject with best fit to data.

%% Load best params from optimizations for different presure data
SBP = [dataMiddle.SBP-16 dataMiddle.SBP dataMiddle.SBP+16];
DBP = [dataMiddle.DBP-16 dataMiddle.DBP dataMiddle.DBP+16];

pathbase = split(pwd,'cardiovascularavatar_T2D_HT');
pathbase = [pathbase{1},filesep, 'cardiovascularavatar_T2D_HT'];
loadfolder = fullfile(pathbase,'Parameters','BPsensitivity');
allparams.SBP = nan(3,length(medianParamsMiddle));
allparams.DBP = nan(3,length(medianParamsMiddle));
allSims.SBP = cell(3,1);
allSims.DBP = cell(3,1);

allSims.SBP{2} = simMiddle;
allSims.DBP{2} = simMiddle;
allparams.SBP(2,:) = medianParamsMiddle;
allparams.DBP(2,:) = medianParamsMiddle;

% Simulation settings
step = 0.001;
numHeartBeats=20;
T = constants(patindex.T);
simtime = sort([dataMiddle.time,0:step:T]);
simtime = unique(simtime);

% SBP min
folder=dir(fullfile(loadfolder,'*SBPmin'));
[bestparams,~,~,~,~,allparams.SBP(1,:)] = findBestParams(fullfile(folder.folder,folder.name),0,length(medianParamsMiddle),1);%get parameters
data = dataMiddle;
data.SBP = data.SBP-16;
[IC] = iccalcs(bestparams,constants,paramNames,constantsNames,data,patindex);
options.x0 = IC';
warning('off')
[~,allSims.SBP{1}] = simulate_avatarHEALTH(bestparams,constants,options,numHeartBeats,patindex,simtime);%simulate
warning('on')

% SBP max
folder=dir(fullfile(loadfolder,'*SBPmax'));
[bestparams,~,~,~,~,allparams.SBP(3,:)] = findBestParams(fullfile(folder.folder,folder.name),0,length(medianParamsMiddle),1);%get parameters
data = dataMiddle;
data.SBP = data.SBP+16;
[IC] = iccalcs(bestparams,constants,paramNames,constantsNames,data,patindex);
options.x0 = IC';
warning('off')
[~,allSims.SBP{3}] = simulate_avatarHEALTH(bestparams,constants,options,numHeartBeats,patindex,simtime);%simulate
warning('on')

% DBP min
folder=dir(fullfile(loadfolder,'*DBPmin'));
[bestparams,~,~,~,~,allparams.DBP(1,:)] = findBestParams(fullfile(folder.folder,folder.name),0,length(medianParamsMiddle),1);%get parameters
data = dataMiddle;
data.DBP = data.DBP-16;
[IC] = iccalcs(bestparams,constants,paramNames,constantsNames,data,patindex);
options.x0 = IC';
warning('off')
[~,allSims.DBP{1}] = simulate_avatarHEALTH(bestparams,constants,options,numHeartBeats,patindex,simtime);%simulate
warning('on')

% DBP max
folder=dir(fullfile(loadfolder,'*DBPmax'));
[bestparams,~,~,~,~,allparams.DBP(3,:)] = findBestParams(fullfile(folder.folder,folder.name),0,length(medianParamsMiddle),1);%get parameters
data = dataMiddle;
data.DBP = data.DBP+16;
[IC] = iccalcs(bestparams,constants,paramNames,constantsNames,data,patindex);
options.x0 = IC';
warning('off')
[~,allSims.DBP{3}] = simulate_avatarHEALTH(bestparams,constants,options,numHeartBeats,patindex,simtime);%simulate
warning('on')

%% Plot sensitivity to selected variables
if doPlot
figure('Visible', 'off','Name','SBPsensitivity_simulations');
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 20;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(4,3,'Padding','compact','TileSpacing','Compact');
varsToPlot = {'Qpv','Qmv','Qav','Qaa','Vla','Vlv','Ela','Elv','pLA','pLV','P_Aortic'};
for var = 1:length(varsToPlot)
    obsValueMax = zeros(size(SBP));
    obsValueMin = zeros(size(SBP));
    isstate = sum(strcmp(varsToPlot{var},xnames));
    obsValuesAllS = getVariableValues(allSims.SBP,varsToPlot{var},ynames,xnames,isstate);
    [~,tmax] = max(obsValuesAllS(2,:));
    [~,tmin] = min(obsValuesAllS(2,:));
    for i = 1:3
        obsValueMax(i) = obsValuesAllS(i,tmax);
        obsValueMin(i) = obsValuesAllS(i,tmin);
    end
    dObsMax = gradient(obsValueMax,SBP);
    dObsMin = gradient(obsValueMin,SBP);
    nexttile
    hold on
    ylabel(varsToPlot{var})
    xlabel('SBP')
    plot(SBP,obsValueMax,'r*')
    plot(SBP,obsValueMin,'b*')
    quiver(SBP,obsValueMax,[16 16 16],dObsMax*16)
    quiver(SBP,obsValueMin,[16 16 16],dObsMin*16)
end
end

%% Plot sensitivity to all parameters
if doPlot
figure('Visible', 'off','Name','SBPsensitivity_params');
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 20;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(7,4,'Padding','compact','TileSpacing','Compact');
for param = 1:length(paramNames)
    paramValueS = zeros(size(SBP));
    for i = 1:3
        paramValueS(i) = allparams.SBP(i,param);
    end
    dParamS = gradient(paramValueS,SBP);
    nexttile
    hold on
    ylabel(paramNames{param})
    xlabel('SBP')
    plot(SBP,paramValueS,'r*')
    ylim([min(min(paramValueS*0.9),min(paramValueS+dParamS*16)),max(max(paramValueS*1.1),max(paramValueS+dParamS*16))])
    quiver(SBP,paramValueS,[16 16 16],dParamS*16,'MaxHeadSize',0.1*(paramValueS(2)/10))
end
end

%% Calculate sensitivity to all parameters, also normalized
dparamdSBP = nan(length(paramNames),1);
dparamdSBPN = nan(length(paramNames),1);
dparamdDBP = nan(length(paramNames),1);
dparamdDBPN = nan(length(paramNames),1);
for param = 1:length(paramNames)
    paramValueS = nan(size(SBP));
    paramValueD = nan(size(SBP));
    for i = 1:3
        paramValueS(i) = allparams.SBP(i,param);
        paramValueD(i) = allparams.DBP(i,param);
    end
    % 2 is the middle gradient out of 3
    dParamS = gradient(paramValueS,SBP);
    dparamdSBPN(param) = dParamS(2)/paramValueS(2);
    dparamdSBP(param) = dParamS(2);
    
    dParamD = gradient(paramValueD,DBP);
    dparamdDBPN(param) = dParamD(2)/paramValueD(2);
    dparamdDBP(param) = dParamD(2);
end

expl={'Normalized sensitivity to SBP (% per mmHg)','Normalized sensitivity to DBP (% per mmHg)','Sensitivity to SBP (per mmHg)','Sensitivity to DBP (per mmHg)'};
BPsensitivity_params=table(round(100*dparamdSBPN,3),round(100*dparamdDBPN,3),dparamdSBP,dparamdDBP,'RowNames',paramNames,'VariableNames',expl);
writetable(BPsensitivity_params,fullfile(plotFolderName,'BPsensitivity_params.xlsx'),"WriteRowNames",1)

%% Calculate sensitivity to selected variables, also normalized
varsToPlot = {'Qpv','Qmv','Qav','Qaa','Vla','Vlv','Ela','Elv','pLA','pLV','P_Aortic'};
dvarsdSBPmax = nan(length(varsToPlot),1);
dvarsdSBPmin = nan(length(varsToPlot),1);
dvarsdSBPmaxN = nan(length(varsToPlot),1);
dvarsdSBPminN = nan(length(varsToPlot),1);
dvarsdDBPmax = nan(length(varsToPlot),1);
dvarsdDBPmin = nan(length(varsToPlot),1);
dvarsdDBPmaxN = nan(length(varsToPlot),1);
dvarsdDBPminN = nan(length(varsToPlot),1);
for var = 1:length(varsToPlot)
    obsValueMaxS = nan(size(SBP));
    obsValueMinS = nan(size(SBP));
    obsValueMaxD = nan(size(SBP));
    obsValueMinD = nan(size(SBP));
    isstate = sum(strcmp(varsToPlot{var},xnames));
    obsValuesAllS = getVariableValues(allSims.SBP,varsToPlot{var},ynames,xnames,isstate);
    obsValuesAllD = getVariableValues(allSims.DBP,varsToPlot{var},ynames,xnames,isstate);
    [~,tmax] = max(obsValuesAllS(2,:));
    [~,tmin] = min(obsValuesAllS(2,:));
    for i = 1:3
        obsValueMaxD(i) = obsValuesAllD(i,tmax);
        obsValueMinD(i) = obsValuesAllD(i,tmin);
        obsValueMaxS(i) = obsValuesAllS(i,tmax);
        obsValueMinS(i) = obsValuesAllS(i,tmin);
    end
    dvarmax = gradient(obsValueMaxS,SBP);
    dvarmin = gradient(obsValueMinS,SBP);
    dvarsdSBPmax(var)=dvarmax(2);
    dvarsdSBPmin(var)=dvarmin(2);
    dvarsdSBPmaxN(var)=dvarmax(2)/obsValueMaxS(2);
    dvarsdSBPminN(var)=dvarmin(2)/obsValueMinS(2);
    
    dvarmax = gradient(obsValueMaxD,DBP);
    dvarmin = gradient(obsValueMinD,DBP);
    dvarsdDBPmax(var)=dvarmax(2);
    dvarsdDBPmin(var)=dvarmin(2);
    dvarsdDBPmaxN(var)=dvarmax(2)/obsValueMaxD(2);
    dvarsdDBPminN(var)=dvarmin(2)/obsValueMinD(2);
end
expl={'Normalized sensitivity to SBP (% per mmHg)',...
    'Normalized sensitivity to DBP (% per mmHg)',...
    'Sensitivity to SBP (per mmHg)',...
    'Sensitivity to DBP (per mmHg)'};
BPsensitivity_vars = table(round(100*[dvarsdSBPminN;dvarsdSBPmaxN],3),...
    round(100*[dvarsdDBPminN;dvarsdDBPmaxN],3),...
    [dvarsdSBPmin;dvarsdSBPmax],...
    [dvarsdDBPmin;dvarsdDBPmax],...
    'RowNames',[strcat(varsToPlot,' (min)'),strcat(varsToPlot,' (max)')],'VariableNames',expl);
writetable(BPsensitivity_vars,fullfile(plotFolderName,'BPsensitivity_vars.xlsx'),"WriteRowNames",1)

%% Save figures
if doPlot
    saveAllFigures(plotFolderName)
end

end

%% Get variable
function var = getVariableValues(simulations,variablename,ynames,xnames,isx)
if isx
    varind = strcmp(xnames,variablename);
else
    varind = strcmp(ynames,variablename);
end
var = nan(length(simulations),length(simulations{2}.t));
for i = 1:length(simulations)
    if ~isempty(simulations{i})
        if isx
            var(i,:) = simulations{i}.x(:,varind)';
        else
            var(i,:) = simulations{i}.y(:,varind)';
        end
    end
end

end