function [minmaxSim,allSims,allParams,allPressure] = findSimulationUncertainty(patNums,parameters,constants,inds,data,options,bestParams,saveResults,loadResults)
%parameters: cell with one or several parameter sets for each patient
%constants: cell with constant values for each patient

%% Load saved results if possible
loadfiles=dir(sprintf('../Parameters/MCMC/minmaxsims_*'));
if loadResults && ~isempty(loadfiles)
    [~,latestfolderInd] = max([loadfiles.datenum]);
    foldernames = {loadfiles.name};
    folderpaths = {loadfiles.folder};
    loadfilename = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    load(loadfilename,'minmaxSim');
    fprintf('findSimulationUncertainty: Loading results from %s \n',loadfilename)
    
    
else
    %% Simulation settings
    numHeartBeats = 20;
    step = 0.001;
    
    %% Setup variables
    minmaxSim.States.min = cell(size(parameters));
    minmaxSim.States.max = cell(size(parameters));
    minmaxSim.States.best = cell(size(parameters));
    
    minmaxSim.Observables.min = cell(size(parameters));
    minmaxSim.Observables.max = cell(size(parameters));
    minmaxSim.Observables.best = cell(size(parameters));
    
    minmaxSim.Observables.all = cell(size(parameters));
    minmaxSim.States.all = cell(size(parameters));
    
    minmaxSim.pressurenames = {'pLA','pLV','pAo','pPV'};
    minmaxSim.pressureindexes = [17,18,1,13];
    for pr = 1:length(minmaxSim.pressurenames)
        minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).best = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min = cell(size(parameters));
        minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).best = cell(size(parameters));
    end
    
    allParams = cell(size(parameters));
    allPressure.SBP = cell(size(parameters));
    allPressure.DBP  = cell(size(parameters));
    allSims = cell(size(parameters));
    
    %% Run simulations
    w = warning('off','all');
    for p = 1:size(parameters,2)%for each subject
        T = constants(inds{p}.T,p);
        simtime = sort([data{p}.time,0:step:T]);
        simtime = unique(simtime);
        options.x0 = data{p}.IC;
        
        parameters{p} = parameters{p}';
        %simulate best paramset
        if ~isempty(bestParams)
            simLast = simulate_avatarHEALTH_short(bestParams(:,p),constants(:,p),options,numHeartBeats,inds{p},simtime);
            [estimatedSBP,estimatedDBP] = brachialpressure(simLast,inds{p},data{p});
            minmaxSim.States.best{p} = simLast.x;
            minmaxSim.Observables.best{p} = simLast.y;
            minmaxSim.SBP.best{p} = estimatedSBP;
            minmaxSim.DBP.best{p} = estimatedDBP;
            SVsim=trapz(simLast.t,simLast.x(:,inds{p}.AV));
            
            EDVsim = max(simLast.x(:,inds{p}.LV));
            ESVsim = min(simLast.x(:,inds{p}.LV));
            EFsim=100*(EDVsim-ESVsim)/EDVsim;%100 * SV/EDV = 100*(EDV-ESV)/EDV
            minmaxSim.SV.best{p} = SVsim;
            minmaxSim.EF.best{p} = EFsim;
            
            minmaxSim.States.min{p} = simLast.x;
            minmaxSim.States.max{p} = simLast.x;
            minmaxSim.Observables.min{p} = simLast.y;
            minmaxSim.Observables.max{p} = simLast.y;
            minmaxSim.time{p} = simLast.t;
            
            minmaxSim.SBP.max{p} = estimatedSBP;
            minmaxSim.SBP.min{p} = estimatedSBP;
            minmaxSim.DBP.max{p} = estimatedDBP;
            minmaxSim.DBP.min{p} = estimatedDBP;
            minmaxSim.SV.max{p} = SVsim;
            minmaxSim.SV.min{p} = SVsim;
            minmaxSim.EF.max{p} = EFsim;
            minmaxSim.EF.min{p} = EFsim;
            
            for pr = 1:length(minmaxSim.pressurenames)
                thispressure = simLast.y(:,minmaxSim.pressureindexes(pr));
                minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max{p} = mean(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min{p} = mean(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).best{p} = mean(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max{p} = max(thispressure)-min(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min{p} = max(thispressure)-min(thispressure);
                minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).best{p} = max(thispressure)-min(thispressure);
            end
            
            %save the best simulation and parameter set for
            %sensitivity analysis
            minmaxSim.SBP.bestParams{p} = bestParams(:,p);
            minmaxSim.SBP.bestSim{p} = simLast;
            minmaxSim.DBP.bestParams{p} = bestParams(:,p);
            minmaxSim.DBP.bestSim{p} = simLast;
        end
        
        %reduce number of simulated params if they are too many
        numparams=size(parameters{p},2);
        if numparams>2000
            fprintf('The number of simulated parameter sets is reduced from %d to %d\n',numparams,2000)
            ind=round(linspace(1,numparams,2000));
            %always include last paramset
            paramsToSim = parameters{p}(:,[ind,numparams]);
        else
            paramsToSim = parameters{p};
        end
        
        % Simulate all parameter sets
        fprintf('Simulating %d parameter sets...\n',size(paramsToSim,2))
        for s = 1:size(paramsToSim,2)%for each parameter set
%             fprintf('Simulating paramset %d/%d\n',s,size(paramsToSim,2))
            simLast = simulate_avatarHEALTH_short(paramsToSim(:,s),constants(:,p),options,numHeartBeats,inds{p},simtime);
            [estimatedSBP,estimatedDBP] = brachialpressure(simLast,inds{p},data{p});
            SVsim=trapz(simLast.t,simLast.x(:,inds{p}.AV));
            
            EDVsim = max(simLast.x(:,inds{p}.LV));
            ESVsim = min(simLast.x(:,inds{p}.LV));
            EFsim=100*(EDVsim-ESVsim)/EDVsim;%100 * SV/EDV = 100*(EDV-ESV)/EDV
            if isempty(bestParams) && s == 1
                minmaxSim.States.min{p} = simLast.x;
                minmaxSim.States.max{p} = simLast.x;
                minmaxSim.Observables.min{p} = simLast.y;
                minmaxSim.Observables.max{p} = simLast.y;
                minmaxSim.time{p} = simLast.t;
                minmaxSim.SBP.max{p} = estimatedSBP;
                minmaxSim.SBP.min{p} = estimatedSBP;
                minmaxSim.DBP.max{p} = estimatedDBP;
                minmaxSim.DBP.min{p} = estimatedDBP;
                minmaxSim.SV.max{p} = SVsim;
                minmaxSim.SV.min{p} = SVsim;
                minmaxSim.EF.max{p} = EFsim;
                minmaxSim.EF.min{p} = EFsim;
                
                minmaxSim.States.best{p} = NaN(size(simLast.x));
                minmaxSim.Observables.best{p} = NaN(size(simLast.y));
                minmaxSim.SBP.best{p} = NaN;
                minmaxSim.DBP.best{p} = NaN;
                minmaxSim.SV.best{p} = NaN;
                minmaxSim.EF.best{p} = NaN;
                
                for pr = 1:length(minmaxSim.pressurenames)
                    thispressure = simLast.y(:,minmaxSim.pressureindexes(pr));
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max{p} = mean(thispressure);
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min{p} = mean(thispressure);
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).best{p} = NaN;
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max{p} = max(thispressure)-min(thispressure);
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min{p} = max(thispressure)-min(thispressure);
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).best{p} = NaN;
                end
                
                for obs = 1:size(simLast.y,2)
                    minmaxSim.all{obs}{p} =  zeros(size(parameters{p},1),length(simtime));
                end
                
                %save the best simulation and parameter set for
                %sensitivity analysis
                minmaxSim.SBP.bestParams{p} = paramsToSim(:,s);
                minmaxSim.SBP.bestSim{p} = simLast;
                minmaxSim.DBP.bestParams{p} = paramsToSim(:,s);
                minmaxSim.DBP.bestSim{p} = simLast;
            else
                % Save the maximum and minimum values for each simulated variable in a struct
                minmaxSim.States.min{p} = min(simLast.x,minmaxSim.States.min{p});
                minmaxSim.States.max{p} = max(simLast.x,minmaxSim.States.max{p});
                minmaxSim.Observables.min{p} = min(simLast.y,minmaxSim.Observables.min{p});
                minmaxSim.Observables.max{p} = max(simLast.y,minmaxSim.Observables.max{p});
                minmaxSim.SBP.max{p} = max(estimatedSBP,minmaxSim.SBP.max{p});
                minmaxSim.SBP.min{p} = min(estimatedSBP,minmaxSim.SBP.min{p});
                minmaxSim.DBP.max{p} = max(estimatedDBP,minmaxSim.DBP.max{p});
                minmaxSim.DBP.min{p} = min(estimatedDBP,minmaxSim.DBP.min{p});
                minmaxSim.SV.max{p} = max(SVsim,minmaxSim.SV.max{p});
                minmaxSim.SV.min{p} = min(SVsim,minmaxSim.SV.min{p});
                minmaxSim.EF.max{p} = max(EFsim,minmaxSim.EF.max{p});
                minmaxSim.EF.min{p} = min(EFsim,minmaxSim.EF.min{p});
                
                for pr = 1:length(minmaxSim.pressurenames)
                    thispressure = simLast.y(:,minmaxSim.pressureindexes(pr));
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max{p} = max(mean(thispressure),minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).max{p});
                    minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min{p} = min(mean(thispressure),minmaxSim.([minmaxSim.pressurenames{pr} 'mean']).min{p});
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max{p} = max(max(thispressure)-min(thispressure),minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).max{p});
                    minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min{p} = min(max(thispressure)-min(thispressure),minmaxSim.([minmaxSim.pressurenames{pr} 'diff']).min{p});
                end
            end
            
            if max(estimatedSBP,minmaxSim.SBP.max{p}) == estimatedSBP
                %save the maximum simulation and parameter set for
                %sensitivity analysis
                minmaxSim.SBP.maxParams{p} = paramsToSim(:,s);
                minmaxSim.SBP.maxSim{p} = simLast;
            end
            if max(estimatedDBP,minmaxSim.DBP.max{p}) == estimatedDBP
                %save the maximum simulation and parameter set for
                %sensitivity analysis
                minmaxSim.DBP.maxParams{p} = paramsToSim(:,s);
                minmaxSim.DBP.maxSim{p} = simLast;
            end
            if min(estimatedSBP,minmaxSim.SBP.max{p}) == estimatedSBP
                %save the minimum simulation and parameter set for
                %sensitivity analysis
                minmaxSim.SBP.minParams{p} = paramsToSim(:,s);
                minmaxSim.SBP.minSim{p} = simLast;
            end
            if min(estimatedDBP,minmaxSim.DBP.max{p}) == estimatedDBP
                %save the minimum simulation and parameter set for
                %sensitivity analysis
                minmaxSim.DBP.minParams{p} = paramsToSim(:,s);
                minmaxSim.DBP.minSim{p} = simLast;
            end
            
            allSims{p}{s} = simLast;
            allPressure.SBP{p}(s) = estimatedSBP;
            allPressure.DBP{p}(s) = estimatedDBP;
        end
        allParams{p} = paramsToSim;
    end
    w = warning('on','all');
    if saveResults
        save(sprintf('Parameters/MCMC/minmaxsims_%s',datestr(now,'yymmdd-hhMM')),'minmaxSim','patNums');
    end
end

end

function [estimatedSBP,estimatedDBP] = brachialpressure(simLast,ind,data)
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
end
