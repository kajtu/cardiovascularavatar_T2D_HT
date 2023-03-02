function [groupTable,statisticsTable,meanTable,stattestTable] = categoricaltestTable(comparison,data,datanames,tablenames,statisticsTable,meanTable,groupTable,indexes,stattestTable)

% calculate categorical test
for t = 1:length(tablenames)
    bothgroups = indexes{comparison.compinds(1)} | indexes{comparison.compinds(2)};
    totalGroupInd = indexes{comparison.compinds(1)}(bothgroups);
    totalVariables = data.(datanames{t})(bothgroups);
    if sum(totalVariables) == length(totalVariables) %groups are the same for this variable. sum(totalVariables) == 0 || 
        statisticsTable.reject{comparison.statind}(t)= NaN;
        meanTable.(tablenames{t})(comparison.pind)= NaN;
        statisticsTable.statistics{comparison.statind}{t} = NaN;
        groupTable.(tablenames{t}){comparison.pind}  = NaN;
        statisticsTable.pvalues{comparison.statind}(t) = NaN;
        stattestTable.(tablenames{t}){comparison.pind} = ' - ';
    else
        %calculate expected values from chi2
        [dochi2,chi2res.chi2,chi2res.p,chi2res.expected,chi2res.calculated] = chi2calc(totalGroupInd,totalVariables);
        %use chi2 if possible (large enough groups), otherwise use fisher
        if ~dochi2
            %use fisher
            [statisticsTable.reject{comparison.statind}(t),...
                meanTable.(tablenames{t})(comparison.pind),...
                statisticsTable.statistics{comparison.statind}{t}] = fishertest(chi2res.calculated);
            stattestTable.(tablenames{t}){comparison.pind} = 'fisher';
        else
            %use chi2
            statisticsTable.reject{comparison.statind}(t) = chi2res.p < 0.005;
            meanTable.(tablenames{t})(comparison.pind) = chi2res.p;
            statisticsTable.statistics{comparison.statind}{t} = chi2res;
            stattestTable.(tablenames{t}){comparison.pind} = 'chi2';
        end
        statisticsTable.statistics{comparison.statind}{t}.chi2 = chi2res.chi2;
        statisticsTable.statistics{comparison.statind}{t}.expected = chi2res.expected;
        statisticsTable.statistics{comparison.statind}{t}.calculated = chi2res.calculated;
        statisticsTable.statistics{comparison.statind}{t}.chi2pval = chi2res.p;
        if meanTable.(tablenames{t})(comparison.pind) < 0.0001
            groupTable.(tablenames{t}){comparison.pind}  = '<0.0001*';
        elseif meanTable.(tablenames{t})(comparison.pind) <= 0.05
            groupTable.(tablenames{t}){comparison.pind}  = sprintf('%0.4f*',meanTable.(tablenames{t})(comparison.pind));
        else
            groupTable.(tablenames{t}){comparison.pind}  = sprintf('%0.4f',meanTable.(tablenames{t})(comparison.pind));
        end
        statisticsTable.pvalues{comparison.statind}(t) = meanTable.(tablenames{t})(comparison.pind);
    end

end


%%%%
    function [dochi2,chi2calc,pcalc,expected,calculated] = chi2calc(rows,columns)
        % gives the same results as:
        %[chi2res.tbl, chi2res.chi2,chi2res.p,chi2res.labels] = crosstab(totalGroupInd,totalVariables);
        tabletotal = length(rows);
        expected(1,1) = (sum(~rows)*sum(~columns))/tabletotal;
        expected(1,2) = (sum(~rows)*sum(columns))/tabletotal;
        expected(2,1) = (sum(rows)*sum(~columns))/tabletotal;
        expected(2,2) = (sum(rows)*sum(columns))/tabletotal;
        
        calculated(1,1) = sum(~rows & ~columns);
        calculated(1,2) = sum(~rows & columns);
        calculated(2,1) = sum(rows & ~columns);
        calculated(2,2) = sum(rows & columns);
        chi2calc = sum(sum(((calculated-expected).^2)./expected));
        df = 1;%(r-1)*(c-1) = (2-1)*(2-1)
        pcalc = 1-chi2cdf(chi2calc,df);
        
        %chi2 is ok if all expected values are >=1 and at least 50% >= 5
        if sum(expected(:)>=5) >= 2 && sum(expected(:)>=1) == 4
            dochi2 = 1;
        else
            dochi2 = 0;
        end
    end

end