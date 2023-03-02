function [bestparams,bestcostAll,allparams,bestconstants,meanparams,medianparams] = findBestParams(folderName,loadsummaryfile,numParams,loadMean)
% Load parameters with lowest cost
if nargin <3
    numParams = 22;
else
    numParams = numParams+1;
end

meanparams = NaN;
medianparams = NaN;
files = dir(fullfile(folderName,'opt-*'));
summaryfile = dir(fullfile(folderName,'bestparam*'));
if isempty(files)
    disp('findBestParams: ERROR - emtpy load folder')
    bestparams  =NaN;
    bestcostAll = NaN;
    allparams  =NaN;
    bestconstants = NaN;
elseif ~loadsummaryfile || (loadsummaryfile && isempty(summaryfile)) || loadMean
    bestcostAll = 1e100;
    allparams = zeros(length(files),numParams);
    for f = 1:length(files)
        load(fullfile(files(f).folder,files(f).name),'optParam','bestcost','constants')
        allparams(f,1) = bestcost;
        if f == 1 && length(optParam) ~= length(allparams(f,2:end)) %if wrong param length
            allparams = zeros(length(files),length(optParam)+1);
            disp('findBestParams: OBS length of parameter vectors does not match!')
        end
        allparams(f,2:end) = optParam;
        if bestcost < bestcostAll
            bestparams = optParam;
            bestcostAll = bestcost;
            bestconstants = constants;
        end
    end
    
    okinds = allparams(:,1) <= bestcostAll*1.10;%1.05
    meanparams = mean(allparams(okinds,2:end));
    medianparams = median(allparams(okinds,2:end));
    
    if loadsummaryfile
        save(sprintf('%s/bestparam(%0.3f).mat',files(f).folder,bestcostAll),'bestparams','bestcostAll','bestconstants','allparams')
    end
else
    %load summaryfile
    fprintf('Findbestparams: Loaded summaryfile %s \n',summaryfile(1).name)
    load(fullfile(summaryfile(1).folder,summaryfile(1).name),'bestparams','bestcostAll','bestconstants','allparams')
end
end