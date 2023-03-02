function [minmaxSim,minmaxSim2] = simulatePredictionUncertainty(patNums,allParams,bestparamsESS,allParams2,bestparamsESS2,constants,inds,data,options)
w = warning('off','all');
% Simulate using all minmax parameter sets
saveResults=0;
loadResults=0;
if isempty(allParams)
    minmaxSim = NaN;
else
    [minmaxSim] = findSimulationUncertainty(patNums,allParams,constants,inds,data,options,bestparamsESS,saveResults,loadResults);
end
[minmaxSim2] = findSimulationUncertainty(patNums,allParams2,constants,inds,data,options,bestparamsESS2,saveResults,loadResults);
w = warning('on','all');
end
