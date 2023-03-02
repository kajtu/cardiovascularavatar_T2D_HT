function [groupTable2,statisticsTable,meanTable,stattestTable] = ttestTable(comparison,data,datanames,tablenames,statisticsTable,meanTable,groupTable2,indexes,toTest,stattestTable)
% calculate unpaired t-test (Two-sample t-test with pooled or unpooled
% variance estimate.) H = ttest2(X,Y). This function performs an unpaired
% two-sample t-test.
% for each variable in the tablenames list
for i = 1:length(toTest)
    t = toTest(i);
    d1 =data.(datanames{t})(indexes{comparison.compinds(1)});
    d2 = data.(datanames{t})(indexes{comparison.compinds(2)});
    if sum(isnan(d1)) == length(d1) || sum(isnan(d2)) == length(d2) %if one of the groups is full of NaN values
        meanTable.(tablenames{t})(comparison.pind) = NaN;
        statisticsTable.reject{comparison.statind}(t) = NaN;
        statisticsTable.statistics{comparison.statind}{t} = NaN;
        groupTable2.(tablenames{t}){comparison.pind} = ' - ';
        statisticsTable.pvalues{comparison.statind}(t) = NaN;
        statisticsTable.rejectBonferroni{comparison.statind}(t) = NaN;
        stattestTable.(tablenames{t}){comparison.pind} = ' - ';
    else
        d1 = d1(~isnan(d1));
        d2 = d2(~isnan(d2));
        [statisticsTable.reject{comparison.statind}(t),...
            meanTable.(tablenames{t})(comparison.pind),...
            statisticsTable.confidenceIntervals{comparison.statind}{t},...
            statisticsTable.statistics{comparison.statind}{t}] = ...
            ttest2(d1,d2,'Alpha',0.05,'Vartype','equal');
        if meanTable.(tablenames{t})(comparison.pind) < 0.0001 %0.0001
            groupTable2.(tablenames{t}){comparison.pind}  = '<0.0001*';
        elseif meanTable.(tablenames{t})(comparison.pind) <= 0.05
            groupTable2.(tablenames{t}){comparison.pind}  = sprintf('%0.4f*',meanTable.(tablenames{t})(comparison.pind));%0.3f
        else
            groupTable2.(tablenames{t}){comparison.pind}  = sprintf('%0.4f',meanTable.(tablenames{t})(comparison.pind));
        end
        statisticsTable.pvalues{comparison.statind}(t) = meanTable.(tablenames{t})(comparison.pind);
    end
    stattestTable.(tablenames{t}){comparison.pind} = 'unpaired ttest';

end

end
