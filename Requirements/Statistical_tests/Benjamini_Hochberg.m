function [reject,resulttable] = Benjamini_Hochberg(pvalues,testnames,falseDiscoveryRate)
% Function written by Kajsa Tunedal 2022.
%Performs benjamini hochberg correction.
%Also implemented in the fdr_BH in multiple testing toolbox on matlab file
%share - which seems to give the same result.
indnan = isnan(pvalues);
if sum(indnan) == length(pvalues)
    resulttable = NaN;
    reject = nan(size(pvalues));
elseif sum(~indnan) == 1
    %only one comparison valid - no correction needed
    resulttable = table;
    [resulttable.pvalueSorted,resulttable.sortindex] = sort(pvalues);
    resulttable.comparison = testnames(resulttable.sortindex)';
    resulttable.rank = [1,nan(1,length(pvalues)-1)]';
    resulttable.criticalValue = [falseDiscoveryRate,nan(1,length(pvalues)-1)]';
    resulttable.reject = nan(size(pvalues));
    resulttable.reject(~indnan) = pvalues(~indnan) <= falseDiscoveryRate;
    reject = pvalues;
    reject(~isnan(reject)) = pvalues(~indnan) <= falseDiscoveryRate;
else
    resulttable = table;
    % Sort p-values and rank them from smallest to largest
    [resulttable.pvalueSorted,resulttable.sortindex] = sort(pvalues);
    resulttable.comparison = testnames(resulttable.sortindex)';
    resulttable.rank = [1:length(pvalues)]';
    
    % Create index to go back to the original order of comparisons
    [~,resulttable.indexback] = ismember(pvalues,resulttable.pvalueSorted);
    
    % Critical value: (k/m)*Q, where Q is the given false discovery rate (number of expected false
    % discoveries / total number of discoveries)
    resulttable.criticalValue = falseDiscoveryRate.*(resulttable.rank./length(testnames(~indnan)));
    
    % Find the largest pvalue that is smaller than the critical value
    breakpoint = resulttable.pvalueSorted < resulttable.criticalValue;
    breakpoint = find(breakpoint,1,'last');
    % All values with rank above and equal to the breakpoint are considered
    % significant (even if their p-value would be larger than the critical value)
    resulttable.reject = zeros(size(pvalues));
    resulttable.reject(1:breakpoint) = 1;
    reject = resulttable.reject;
    
    if sum(indnan)>0
        indnansorted = isnan(resulttable.pvalueSorted);
        resulttable.rank(indnansorted) = NaN;
        resulttable.criticalValue(indnansorted) = NaN;
        resulttable.reject(indnansorted) = NaN;
        reject(indnan) = NaN;
        reject(~indnan)= resulttable.reject(resulttable.indexback(resulttable.indexback~=0));
    else
            reject = reject(resulttable.indexback);
    end
    
end
end