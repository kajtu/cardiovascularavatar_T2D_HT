function [cost,simLast,residuals,SVdiffs,costs,rmse,weightedCosts,percenterror,MAPE] = costFunction_HEALTH_allcalculations(theta,constants,allparams,simulationOptions,numHeartBeats,data,ind,datanames)

% Calculates the cost for the Avatar model. 
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
if (simLast.status<0) || (nansum > 0) %if simulation is NaN
    cost = 1e10;
    residuals = NaN;
    SVdiffs = NaN;
    costs = NaN;
    MAPE = NaN;
    simLast = NaN;
    rmse = NaN;
    weightedCosts = NaN;
    percenterror = NaN;
    disp(['Simulation failed. Nan sum ' num2str(nansum) ', Status ' num2str(simLast.status)])
else
    %% calculate cost for each of the four blood flow measurements
    indsyst = 1:ind.tdiast; %syst
    inddiast = ind.tdiast:length(datatime); %diast
    
    if nargin < 9
        datanames = {'MV','AV','AC','PV'};%{'MV','AV','AC','PV'};
    end
    simnames = datanames;
    dataparams = {'Emax_LV','Caa','Ctot','Rtot','ElCo'};
    
    costs = zeros(1,length(datanames));
    rmse = zeros(1,length(datanames));
    weightedCosts = zeros(1,length(datanames));
    residuals = cell(1,length(datanames)+2+length(dataparams));
    MAPE = zeros(1,length(datanames)+2+length(dataparams));
    compSims = cell(1,length(datanames));
    compTime = cell(1,length(datanames));
    compDatas = cell(1,length(datanames));
    SVdiffs = zeros(1,length(datanames));
    oscillationPunishment = zeros(1,length(datanames));
    highresSims = cell(1,length(datanames));
    highresTime = cell(1,length(datanames));
    %pkLocs = cell(1,length(datanames));
    percenterror.all = zeros(1,length(datanames)+2+length(dataparams));
    percenterror.flow = 0;
    alldatasum=0;
    allresidualssum=0;
    for i = 1:length(datanames)
        if strcmp('MV',datanames{i})
            dataind = inddiast;
        elseif strcmp('PV',datanames{i}) %PV
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
        
%         if ~strcmp('MV',datanames{i}) &&  ~strcmp('PV',datanames{i})
%             [pks,pkLocs{i}] = findpeaks(highresSims{i},'NPeaks',Npeaks+1);
%             if length(pks) > Npeaks %if more than expected num of peaks
%                 %disp('oscillationPunishment, but not included in cost')
%                 %oscillationPunishment(i) = 30;
%             end
%         end
        
        compDatas{i} = data.(datanames{i})(dataind,:);
        residuals{i} = (compSims{i} - compDatas{i}(:,2)');
        weightedCosts(i) = sum(residuals{i}.^2)/mean(compDatas{i}(:,2))/length(compSims{i}); %normalize num of datapoints
        
        costs(i) = sum(residuals{i}.^2);
        rmse(i) = sqrt(sum(residuals{i}.^2)/length(compSims{i}));
        % Look at SV (total blood volume through the AV/MV/AC during the cardiac
        % cycle)
        SVdiffs(i) = trapz(compTime{i},compSims{i}) - trapz(compTime{i},compDatas{i}(:,2)');
        percenterror.all(i) = 100*sum(abs(residuals{i}))/sum(compDatas{i}(:,2));
        allresidualssum = allresidualssum+sum(abs(residuals{i}));
        alldatasum = alldatasum+sum(compDatas{i}(:,2));
        
        MAPE(i) = 100*sum(abs(residuals{i}./compDatas{i}(:,2)'));
    end
    flowcosts = costs;
    percenterror.flow = 100*allresidualssum/alldatasum;
    
    %rmse all flow
    rmseflow=0;n=0;
    for i = 1:length(datanames)
        rmseflow = rmseflow+sum(residuals{i}.^2);
        n = length(compSims{i})+n;
    end
    rmseflow = sqrt(rmseflow/n);
    rmse = [rmse,rmseflow];
    
    numdatapoints = n;
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
 
    residuals(length(datanames)+1:length(datanames)+2) = {data.SBP - estimatedSBP,data.DBP - estimatedDBP};
    percenterror.all(length(datanames)+1:length(datanames)+2) = 100*[abs(data.SBP - estimatedSBP)/data.SBP, abs(data.DBP - estimatedDBP)/data.DBP];
    rmse = [rmse,sqrt((costSBP+costDBP)/2)];
    
    MAPE(length(datanames)+1) = 100*abs((data.SBP - estimatedSBP)/data.SBP);
    MAPE(length(datanames)+2) = 100*abs((data.DBP - estimatedDBP)/data.DBP);
    numdatapoints = numdatapoints+2;
    %% parameters
    paramdatasum=0;paramresidualsum=0;
    paramcosts = zeros(1,length(data.params));
    paramcostsWeighted = zeros(1,length(data.params));
    for p = 1:length(dataparams)
        pval = theta(ind.(dataparams{p}));
        dval = data.params.(dataparams{p});
        paramcosts(p) = ((pval - dval(1))^2);
        paramcostsWeighted(p) = ((pval - dval(1))^2)/dval(1);
        residuals(length(datanames)+2+p) = {pval - dval(1)};
        percenterror.all(length(datanames)+2+p) = 100*(abs(pval - dval(1)))/dval(1);
        paramresidualsum = paramresidualsum+abs(pval - dval(1));
        paramdatasum = paramdatasum+dval(1);
        MAPE(length(datanames)+2+p) = 100*abs((dval(1) - pval)/dval(1));
        numdatapoints = numdatapoints+1;
    end
    weightedCosts = [weightedCosts,paramcostsWeighted];
    
    costs = [costs,costSBP+costDBP,paramcosts];
    rmse = [rmse,sqrt(sum(paramcosts)/length(paramcosts))];
    percenterror.params = 100*(paramresidualsum+abs(data.SBP - estimatedSBP)+abs(data.DBP - estimatedDBP))/(data.SBP+data.DBP+paramdatasum);
    percenterror.sumall = 100*(allresidualssum+paramresidualsum+abs(data.SBP - estimatedSBP)+abs(data.DBP - estimatedDBP))/(data.SBP+data.DBP+paramdatasum+alldatasum);
    MAPE = sum(MAPE)/numdatapoints;
    
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
        + sum(paramcostsWeighted)*paramweight + costSBPWeighted*SBPweight + costDBPWeighted*DBPweight;%costMAPWeighted*BPweight;% + costLA;
end

end 