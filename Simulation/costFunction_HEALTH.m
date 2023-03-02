function [cost] = costFunction_HEALTH(theta,constants,allparams,simulationOptions,numHeartBeats,data,ind,dispErrors,doPlot)
% Calculates the cost for the Avatar model. 

%convert the parameters back after optimization
allparams(ind.estParams) = 10.^theta;
theta = allparams;

T = constants(ind.T);

%% SIMULATE THE MODEL FOR THE GIVEN PARAM
simulationOptions.x0 = data.IC;
datatime = data.MV(:,1)';
simtime = sort([datatime,0:0.001:T]);
simtime = unique(simtime);
datatimeinds = find(ismember(simtime,datatime));
try
    simLast = simulate_avatarHEALTH_short(theta,constants,simulationOptions,numHeartBeats,ind,simtime);
catch myError
    disp("!!Error when simulating!!")
    disp(myError.message)
    cost = 1e99; %Return an extremly large cost if the simulation fails and there is an error.
    return %This command causes the cost function to stop here, and it does not read the subsequent lines
end 


nansum = sum(isnan(simLast.x(:)));
if (simLast.status<0) || (nansum > 0) %if simulation is NaN (ie the simulation failed)
    cost = 1e10;
    doPlot = 0;
    if dispErrors
        disp(['Simulation failed. Nan sum ' num2str(nansum) ', Status ' num2str(simLast.status)])
    end
else
    %% calculate cost for each of the three blood flow measurements
    indsyst = 1:ind.tdiast; %systole
    inddiast = ind.tdiast:length(datatime); %diastole
    
    % which simulations and parameters to be compared with data
    datanames = {'MV','AV','AC','PV'};
    simnames = datanames;
    dataparams = {'Emax_LV','Caa','Ctot','Rtot','ElCo'};
    
    weightedCosts = zeros(1,length(datanames));
    residuals = cell(1,length(datanames));
    compSims = cell(1,length(datanames));
    compTime = cell(1,length(datanames));
    compDatas = cell(1,length(datanames));
    oscillationPunishment = zeros(1,length(datanames));
    highresSims = cell(1,length(datanames));
    highresTime = cell(1,length(datanames));
    pkLocs = cell(1,length(datanames));

    for i = 1:length(datanames)
        if strcmp('MV',datanames{i})
            dataind = inddiast;
        elseif strcmp('PV',datanames{i}) 
            dataind = 1:length(datatime);
        else
            dataind = indsyst;
            Npeaks = 1;
        end
        simind = datatimeinds(dataind);
        compSims{i} = simLast.x(simind,ind.(simnames{i}))';
        compTime{i} = simLast.t(simind);
        indt = find(ismember(round(simLast.t,2),round(compTime{i},2)));
        highresSims{i} = simLast.x(indt(1):indt(end),ind.(simnames{i}))';
        highresTime{i} = simLast.t(indt(1):indt(end));
        
        if ~strcmp('MV',datanames{i}) &&  ~strcmp('PV',datanames{i})
            [pks,pkLocs{i}] = findpeaks(highresSims{i},'NPeaks',Npeaks+1);
            if length(pks) > Npeaks %if more than expected num of peaks: add cost punishment
                oscillationPunishment(i) = 30;
            end
        end
        compDatas{i} = data.(datanames{i})(dataind,:);
        residuals{i} = (compSims{i} - compDatas{i}(:,2)').^2;
        weightedCosts(i) = sum(residuals{i})/mean(compDatas{i}(:,2))/length(compSims{i}); %normalize number of datapoints
    end
    
    %% Brachial/Aortic pressure
    maxAorticPressure = max(simLast.y(:,ind.brachialPressure));
    minAorticPressure = min(simLast.y(:,ind.brachialPressure));
    
    % diffparam determining difference between brachial and aortic systolic
    % pressure (data from anglo-cardiff trial II
    % http://dx.doi.org/10.1161/HYPERTENSIONAHA.107.105445) 8.6 mean value
    lb = 0.1;ub =18.9;
    diffparam = max(min(data.SBP-maxAorticPressure,ub),lb);%parameter determining the difference between brachial and aortic SBP
    
    estimatedSBP = maxAorticPressure + diffparam;
    estimatedDBP = minAorticPressure;

    costSBP = (data.SBP - estimatedSBP)^2;
    costDBP = (data.DBP - estimatedDBP)^2;
    
    costSBPWeighted = costSBP /data.SBP;
    costDBPWeighted= costDBP/data.DBP;
    weightedCosts = [weightedCosts,costSBPWeighted,costDBPWeighted];
    
    %% Parameters
    paramcosts = zeros(1,length(data.params));
    paramcostsWeighted = zeros(1,length(data.params));
    for p = 1:length(dataparams)
        pval = theta(ind.(dataparams{p}));
        dval = data.params.(dataparams{p});
        paramcosts(p) = ((pval - dval(1))^2);
        paramcostsWeighted(p) = ((pval - dval(1))^2)/dval(1);
    end
    weightedCosts = [weightedCosts,paramcostsWeighted];
    
    %% punishment if Emax < Emin (non physiological)
    emaxEminPunishment = 0;
    if theta(ind.Emax_LA) < theta(ind.Emin_LA)
        emaxEminPunishment = 200 + (theta(ind.Emin_LA)-theta(ind.Emax_LA))*10000;
    end
    if theta(ind.Emax_LV) < theta(ind.Emin_LV)
        emaxEminPunishment = emaxEminPunishment + 200 + (theta(ind.Emin_LV)-theta(ind.Emax_LV))*10000;
    end
    
    %% Final cost
    % (AV + AC/2 to even out 2 measurments around systole from AV & AC compared
    % to 1 measurement in diastiole from MV
    SBPweight = 3;
    DBPweight = 3;
    paramweight = 3;
    
    cost = weightedCosts(4) + weightedCosts(1) + (weightedCosts(2) + weightedCosts(3))/2 ...
        +sum(oscillationPunishment) + emaxEminPunishment ...
        + sum(paramcostsWeighted)*paramweight + costSBPWeighted*SBPweight + costDBPWeighted*DBPweight;
end

%% Plotting
if doPlot
    f = figure('Name','Cost plot');
    f.Units = 'Normalized';
    f.OuterPosition = [0 0 1 1 ];
    sgtitle(sprintf('Cost: %0.2f (pure: %0.2f)', cost,sum(weightedCosts)))
    hold on
    for i = 1:length(compDatas)
        subplot(2,3,i)
        if (i == 2 || i == 3)
            title(sprintf('%s: %0.2f (pure %0.2f)',datanames{i},weightedCosts(i)/2,weightedCosts(i)))
        end
        hold on
        plot(compDatas{i}(:,1),compDatas{i}(:,2),'ko-','LineWidth',2)
         if strcmp(datanames{i},'PV')
             stdv = 0.3;
         else
             stdv = 0.15;
         end
        errorbar(compDatas{i}(:,1),compDatas{i}(:,2),ones(1,length(compDatas{i}))*mean(compDatas{i}(:,2)*stdv),'k','LineWidth',1.5)
        plot(highresTime{i},highresSims{i},'r--','LineWidth',2)
        plot(compTime{i},compSims{i},'r*','LineWidth',2)
        for p = 1:length(pkLocs{i})
            plot(highresTime{i}(pkLocs{i}(p)),highresSims{i}(pkLocs{i}(p)),'bo','LineWidth',2)
        end
        legend({'Data','Simulation'})
    end
    
    subplot(2,3,i+1)
    title(sprintf('Params: %0.2f (pure %0.2f)',sum(paramcostsWeighted)*paramweight,sum(paramcostsWeighted)))
    hold on
    for p = 1:length(dataparams)
        pval = theta(ind.(dataparams{p}));
        dval = data.params.(dataparams{p});
        plot(p,log10(pval),'r*','LineWidth',2)
        plot(p,log10(dval(1)),'ko','LineWidth',2)
    end
    xticks(1:length(dataparams))
   xticklabels(dataparams);
   ylabel('Param value (log10)')
   
   subplot(2,3,i+2)
   title(sprintf('Cost SBP: %0.2f (pure %0.2f), Cost DBP: %0.2f (pure %0.2f)',costSBPWeighted*DBPweight,costSBPWeighted,costDBPWeighted*DBPweight,costDBPWeighted))
   hold on
   plot(simLast.t,simLast.y(:,ind.brachialPressure),'r-','LineWidth',2)
   plot(simLast.t(simLast.y(:,ind.brachialPressure) == maxAorticPressure),estimatedSBP,'r*')
   plot(simLast.t(simLast.y(:,ind.brachialPressure) == minAorticPressure),estimatedDBP,'r*')
   yline(data.SBP,'k-','LineWidth',2);
   yline(data.DBP,'k-','LineWidth',2);
   xlim([0 simLast.t(end)])
end

end