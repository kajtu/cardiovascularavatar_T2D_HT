function [patNumsSims,bestparams,bestcosts,meanparams,allokParams,...
    allokCosts,minValuesparams,maxValuesparams,numsincluded,data,inds,constants,medianparams] = findAllParams(resultsFolder,patNumsOrig,patNumsPlot,paramNames,data,inds,constants)

%% Load saved results if possible
%load saved results if it exists
loadResults = 1;
loadfiles=dir(sprintf('./Parameters/MCMC/allokparams_*'));
if loadResults && ~isempty(loadfiles)
    % Always load the latest results
    [~,latestfolderInd] = max([loadfiles.datenum]);
    foldernames = {loadfiles.name};
    folderpaths = {loadfiles.folder};
    loadfilename = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
    load(loadfilename,'patNumsSims','bestparams','bestcosts','meanparams','medianparams','allokParams',...
        'allokCosts','minValuesparams','maxValuesparams','numsincluded','data','inds','constants');
    fprintf('findAllParams: Loading existing results from %s \n',loadfilename)
    
    % In patient-specific predictions (eg. when running createFigures_hypertensionT2D_predictions),
    % only the results for one subject at a time is wanted.
    % Thus, we remove all other subjects:
    removenow = ~ismember(patNumsSims,patNumsPlot);
    if sum(removenow) > 0
        disp('findAllParams: removing all subjects not wanted. Wanted subjects:')
        disp(patNumsPlot)
        allokParams(removenow) = [];
        allokCosts(removenow) = [];
        bestparams(:,removenow) = [];
        bestcosts(removenow) = [];
        meanparams(:,removenow) = [];
        medianparams(:,removenow) = [];
        data(removenow) = [];
        constants(:,removenow) = [];
        minValuesparams(:,removenow) = [];
        maxValuesparams(:,removenow) = [];
        inds(removenow) = [];
        numsincluded.MCMC(removenow) = [];
        numsincluded.ESS(removenow) = [];
        patNumsSims(removenow) = [];
    end
        
else
    %pre-define vectors
    patNumsSims = [];
    maxValuesparams = zeros(length(paramNames),length(patNumsOrig));
    minValuesparams= zeros(length(paramNames),length(patNumsOrig));
    allokParams = cell(1,length(patNumsOrig));
    allokCosts = cell(1,length(patNumsOrig));
    meanparams = zeros(length(paramNames),length(patNumsOrig));
    medianparams = zeros(length(paramNames),length(patNumsOrig));
    bestparams= zeros(length(paramNames),length(patNumsOrig));
    bestcosts= zeros(length(patNumsOrig),1);
    numsincluded.MCMC = nan(1,length(patNumsOrig));
    numsincluded.ESS = nan(1,length(patNumsOrig));
    
    %settings for findBestParams
    loadmean=0;
    loadsummaryfile=0;
    
    % Load parameters for each subject p
    for p = 1:size(patNumsOrig,2)
        patientNum = patNumsOrig{p};
        
        %% Load ESS params
        % find latest patientfolder with ESS optimizations
        loadresultsfolder = [resultsFolder '/P' patientNum '_*'];
        folders = dir(loadresultsfolder);
        if isempty(folders)
            disp(['findAllParams: OBS couldnt load ESS parameters for p' patientNum, ' (no prev folder)'])
            foundESS = false;
        else
            foundESS = true;
            [~,latestfolderInd] = max([folders.datenum]);
            foldernames = {folders.name};
            folderpaths = {folders.folder};
            
            folderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
            if isempty(dir([folderName,'/opt-*'])) && doOptimization && length(foldernames) > 1%new empty folder created, so take the next latest one if it exists
                folderName = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd-1});
            end
            [bestessParam,bestessCost,allparamsESS,~,~] = findBestParams(folderName,loadsummaryfile,length(paramNames),loadmean);
        end
        %% Load MCMC params
        % Find patient folder
        folders = dir(['Parameters/MCMC/P' patientNum '_*']);
        if length(folders) < 1
            foundMCMC = false;
            bestMCMCcost = bestessCost;
        else
            %find latest folder
            [~,latestfolderInd] = max([folders.datenum]);
            foldernames = {folders.name};
            folderpaths = {folders.folder};
            folder = fullfile(folderpaths{latestfolderInd},foldernames{latestfolderInd});
            
            % Load min & max for all parameters in this patient folder
            files = dir(fullfile(folder,'parameters-*'));
            if isempty(files)
                foundMCMC = false;
                bestMCMCcost = bestessCost;
            else
                foundMCMC = true;
                load(fullfile(files(1).folder,files(1).name),'parameters');
                allMCMCparams = 10.^(parameters.S.par);
                allMCMCcosts = -parameters.S.logPost;
                [bestMCMCcost,bestMCMCind] = min(allMCMCcosts);
                bestMCMCparams = allMCMCparams(:,bestMCMCind);
                
                %best of mcmc and ess
                if bestMCMCcost < bestessCost
                    %fprintf('MCMC opt lower cost (MCMC %0.2f vs ESS %0.2f)\n',bestMCMCcost,bestessCost)
                    bestparams(:,p) = bestMCMCparams;
                    bestcosts(p) = bestMCMCcost;
                else
                    %fprintf('ESS opt lower cost (MCMC %0.2f vs ESS %0.2f)\n',bestMCMCcost,bestessCost)
                    bestparams(:,p) = bestessParam';
                    bestcosts(p) = bestessCost;
                end
            end
        end
        
        if foundESS
            patNumsSims = [patNumsSims,{patNumsPlot{p}}];
            %% Combine all params
            if foundMCMC
                allparams = [allparamsESS(:,2:end);allMCMCparams'];
                allcosts = [allparamsESS(:,1);allMCMCcosts];
            else
                disp('findAllParams: no MCMC params found')
                allparams = allparamsESS(:,2:end);
                allcosts = allparamsESS(:,1);
            end
            
            %include all params within +10% of best MCMC cost
            okinds = allcosts <= bestMCMCcost*1.10;
            allokParams{p} = allparams(okinds,:);
            allokCosts{p} = allcosts(okinds,:);
            meanparams(:,p) = mean(allokParams{p});
            medianparams(:,p) = median(allokParams{p});
            
            minValuesparams(:,p) = min(allokParams{p},[],1);
            maxValuesparams(:,p) = max(allokParams{p},[],1);
            
            if foundMCMC
                numsincluded.MCMC(p) = sum(allMCMCcosts <= bestMCMCcost*1.1);
            else
                numsincluded.MCMC(p) = NaN;
            end
            numsincluded.ESS(p) = sum(allparamsESS(:,1) <= bestMCMCcost*1.1);
        end
    end
    
    %% remove subjects where no parameters were loaded
    if isempty(patNumsSims)
        removeindex = 1:length(patNumsPlot);
    else
        removeindex = ~ismember(patNumsPlot,patNumsSims);
        allokParams(removeindex) = [];
        allokCosts(removeindex) = [];
        bestparams(:,removeindex) = [];
        bestcosts(removeindex) = [];
        meanparams(:,removeindex) = [];
        medianparams(:,removeindex) = [];
        data(removeindex) = [];
        constants(:,removeindex) = [];
        minValuesparams(:,removeindex) = [];
        maxValuesparams(:,removeindex) = [];
        inds(removeindex) = [];
        numsincluded.MCMC(removeindex) = [];
        numsincluded.ESS(removeindex) = [];
    end
    
    %% Save results
    savefilename = sprintf('./Parameters/MCMC/allokparams_%s',datestr(now,'yymmdd-HHMM'));
    save(savefilename,'patNumsSims','bestparams','bestcosts','meanparams','medianparams','allokParams',...
        'allokCosts','minValuesparams','maxValuesparams','numsincluded','removeindex','data','inds','constants');
    
end



end