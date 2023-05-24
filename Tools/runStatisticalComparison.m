function [groupTable2,statisticsTable,meanTable,checkMean,checkNormality,...
    numRejectedHyp,rejecttables,medianTable,stattestTable,...
    statdocumenttable] = runStatisticalComparison(testName,patNums,data,datanames,tablenames)
% grouptable2: with strings eg "mean +-sd"
% groupTable: with pure results, mean std and pvalues separatly
% statisticsTable: with all statistics saved: statistics, confidence
% intervals, rejection true/false and p-values

pathbase = split(pwd,'cardiovascularavatar_T2D_HT');
addpath(fullfile(pathbase{1},'cardiovascularavatar_T2D_HT/Requirements/Statistical_tests'))

%% load group indexes based on patNums
[groups] = loadGroupIndexes(patNums);

%% setup tables
groupNames = {'Controls HBP','T2D HBP',...
    'Hypertensive HBP','Hypertensive T2D HBP'};

compNames = {
    'HBP: control vs hypertensive', ...
    'HBP: T2D vs hypertensive T2D',...
    'HBP: Hypertensive vs hypertensive & T2D',...
    'HBP: control vs T2D',...
    'HBP: control vs hypertensive T2D',...
    'HBP: Hypertensive vs T2D'};

tablenamesSD = strcat('sd',tablenames);
tablenamesIQR = strcat('iqr',tablenames);

meanTable = table();
indexes = {groups.C_NT_home,groups.T2D_NT_home,...
    groups.C_HT_home,groups.T2D_HT_home}';

meanTable.Group = [groupNames,compNames]';
for g = 1:length(indexes)
    meanTable.N(g) = sum(indexes{g});
end
meanTable.N(g+1:end)  = nan(length(compNames),1);

groupTable2 = table;
groupTable2.Group = meanTable.Group;
groupTable2.N = meanTable.N;
groupTable2.correctedP = cell(size(meanTable.N));

%normality check tables
checkMean = table;
checkMean.Group = meanTable.Group;
checkMeanall = table;
checkMeanall.names = tablenames';
checkMeanall.mean = nan(size(tablenames))';
checkMeanall.std = nan(size(tablenames))';
checkMeanall.perc = nan(size(tablenames))';

checkNormality.isnormal = nan(length(tablenames),length(indexes)+1);
checkNormality.pvalNormality = nan(length(tablenames),length(indexes)+1);
checkNormality.sdPercMean = nan(length(tablenames),length(indexes)+1);
checkNormality.groupNames = [groupNames,{'all'}];
checkNormality.tablenames = tablenames;

%mean and median tables
for t = 1:length(tablenames)
    meanTable.(tablenames{t}) =nan(size(meanTable.Group));
    groupTable2.(tablenames{t}) = cell(size(groupTable2.Group));
    checkMean.(tablenames{t}) =nan(size(meanTable.Group));
end
medianTable = meanTable;
for t = 1:length(tablenamesSD)
    meanTable.(tablenamesSD{t}) =nan(size(meanTable.Group));
    medianTable.(tablenamesIQR{t}) = nan(size(medianTable.Group)); 
end
medianTable.mediansubjects = cell(size(medianTable.Group));
medianTable.groupsubjects = cell(size(medianTable.Group));
for g = 1:length(medianTable.Group)
    medianTable.mediansubjects{g} = cell(1,length(tablenames));
    medianTable.groupsubjects{g} = cell(1,length(tablenames));
end
meanTable.meansubjects=medianTable.mediansubjects;
meanTable.groupsubjects=medianTable.groupsubjects;

%table for which statistical test that was performed
stattestTable = groupTable2(:,[1:2,4:end]);

% create statistics table
groupComparisonNames= compNames;
statistics = cell(size(groupComparisonNames));
confidenceIntervals = cell(size(groupComparisonNames));
pvalues = cell(size(groupComparisonNames));
reject = pvalues;
for n = 1:length(groupComparisonNames)
    statistics{n} =cell(size(tablenames));
    confidenceIntervals{n} =cell(size(tablenames));
    pvalues{n} = nan(size(tablenames));
    reject{n} =zeros(size(tablenames));%logical
end

rejectBenjH2=reject;
statisticsTable = table(groupComparisonNames,statistics,confidenceIntervals,pvalues,reject,rejectBenjH2);



%% fill tables - calculate mean, median and std for each group and each variable in the tablenames
for ind = 1:length(indexes)
    for t = 1:length(tablenames)
        %check for normal distribution using a Shapiro-Wilk parametric hypothesis test of composite normality
        if ~strcmp('categorical',testName)
            if ind == 1
                try
                [checkNormality.isnormal(t,length(indexes)+1), checkNormality.pvalNormality(t,length(indexes)+1), ~] = swtest(data.(datanames{t}), 0.005);
                catch e
                    fprintf('Shapiro-w could not run for %s! The message was:\n%s\n',tablenames{t},e.message);
                    checkNormality.isnormal(t,length(indexes)+1) = 0;
                    checkNormality.pvalNormality(t,length(indexes)+1) = 0;
                end
                checkNormality.sdPercMean(t,length(indexes)+1) = 100* nanstd(data.(datanames{t}))/nanmean(data.(datanames{t}));
            end
            thisvariable=data.(datanames{t})(indexes{ind});
            if length(thisvariable(~isnan(thisvariable))) > 2
                try
                    [checkNormality.isnormal(t,ind), checkNormality.pvalNormality(t,ind), ~] = swtest(thisvariable, 0.005);
                catch e
                    fprintf('Shapiro-w could not run for %s! The message was:\n%s\n',tablenames{t},e.message);
                    checkNormality.isnormal(t,ind)=0;
                    checkNormality.pvalNormality(t,ind)=NaN;
                end
            else
                checkNormality.isnormal(t,ind)=0;
                checkNormality.pvalNormality(t,ind)=NaN;
            end
            checkNormality.sdPercMean(t,ind) = 100* nanstd(thisvariable)/nanmean(thisvariable);
        end
        
        allvalues = data.(datanames{t})(indexes{ind});
        patNumsind = patNums(indexes{ind});
        medianTable.groupsubjects{ind}{t} = patNumsind;
        meanTable.groupsubjects{ind}{t} = patNumsind;
        
        %calculate mean & standard deviation
        meanTable.(tablenames{t})(ind) = mean(allvalues,'omitnan');
        meanTable.(tablenamesSD{t})(ind) = std(allvalues,'omitnan');
        [mindiff,meanPat] = min(abs(allvalues- meanTable.(tablenamesSD{t})(ind)));
        meanTable.meansubjects{ind}{t} = patNumsind{meanPat};%mean subject for the variable t in group ind
        
        %calculate median& interquartile range
        medianTable.(tablenames{t})(ind) = median(allvalues,'omitnan');
        medianTable.(tablenamesIQR{t})(ind) = iqr(allvalues);%or use quantile(allvalues,[0.25 0.75])
        [mindiff,medianPat] = min(abs(allvalues-medianTable.(tablenames{t})(ind)));
        medianTable.mediansubjects{ind}{t} = patNumsind{medianPat};%median subject for the variable t in group ind
    end
end

%% Define comparison indexes for all group comparisons
comparisons = {};
comparisongroups = cell(6,2);

% four groups based on home blood pressure
[comparisongroups{1:6,1}] = deal('fourgroups_hbp');
[comparisongroups{1:6,2}] = deal({'Controls HBP','Hypertensive HBP','Hypertensive T2D HBP','T2D HBP'});
           
comparisons{1}.name = 'HBP: control vs hypertensive';
comparisons{1}.compinds = [find(strcmp(meanTable.Group,'Controls HBP')),find(strcmp(meanTable.Group,'Hypertensive HBP'))];

comparisons{2}.name = 'HBP: T2D vs hypertensive T2D';
comparisons{2}.compinds = [find(strcmp(meanTable.Group,'Hypertensive T2D HBP')),find(strcmp(meanTable.Group,'T2D HBP'))];

comparisons{3}.name = 'HBP: Hypertensive vs hypertensive & T2D';
comparisons{3}.compinds = [find(strcmp(meanTable.Group,'Hypertensive HBP')),find(strcmp(meanTable.Group,'Hypertensive T2D HBP'))];

comparisons{4}.name = 'HBP: control vs T2D';
comparisons{4}.compinds = [find(strcmp(meanTable.Group,'T2D HBP')),find(strcmp(meanTable.Group,'Controls HBP'))];

comparisons{5}.name = 'HBP: control vs hypertensive T2D';
comparisons{5}.compinds = [find(strcmp(meanTable.Group,'Controls HBP')),find(strcmp(meanTable.Group,'Hypertensive T2D HBP'))];

comparisons{6}.name = 'HBP: Hypertensive vs T2D';
comparisons{6}.compinds = [find(strcmp(meanTable.Group,'Hypertensive HBP')),find(strcmp(meanTable.Group,'T2D HBP'))];

    
for c = 1:length(comparisons)
    comparisons{c}.pind = find(strcmp(meanTable.Group,comparisons{c}.name));
    comparisons{c}.statind = find(strcmp(statisticsTable.groupComparisonNames,comparisons{c}.name));
end

checkNormality.isnormalComparisongroups = nan(length(tablenames),length(comparisons));
checkNormality.isnormalComparisongroupNames = cell(1,length(comparisons));
checkNormality.isnormalComparisongroupNamesAll = cell(1,length(comparisons));

%% make the comparisons
lastgroup = '';
for c = 1:length(comparisons)
    comparison = comparisons{c};
    switch testName
        case 'nonparametric'
            indnormal = zeros(size(tablenames));
            [groupTable2,statisticsTable,meanTable,stattestTable] = nonparametrictestTable(comparison,data,datanames,tablenames,statisticsTable,meanTable,groupTable2,indexes,1:length(tablenames),stattestTable);
        case 'ttest'
            indnormal = ones(size(tablenames));
            [groupTable2,statisticsTable,meanTable,stattestTable] = ttestTable(comparison,data,datanames,tablenames,statisticsTable,meanTable,groupTable2,indexes,1:length(tablenames),stattestTable);
        case 'numeric'
            compgroups = ismember(checkNormality.groupNames,[comparisongroups{c,2},{'all'}]);%indexes for the groups in the comparison. 'all' is to test for all subjects
            indnormal = sum(checkNormality.isnormal(:,compgroups),2) == length(comparisongroups{c,2})+1; %check if all groups in the comparison are normally distributed
            [groupTable2,statisticsTable,meanTable,stattestTable] = ttestTable(comparison,data,datanames,tablenames,statisticsTable,meanTable,groupTable2,indexes,find(indnormal),stattestTable);
            [groupTable2,statisticsTable,meanTable,stattestTable] = nonparametrictestTable(comparison,data,datanames,tablenames,statisticsTable,meanTable,groupTable2,indexes,find(~indnormal),stattestTable);
        case 'categorical'
            indnormal = ones(size(tablenames));
            [groupTable2,statisticsTable,meanTable,stattestTable] = categoricaltestTable(comparison,data,datanames,tablenames,statisticsTable,meanTable,groupTable2,indexes,stattestTable);
    end
    checkNormality.isnormalComparisongroupNamesAll{c} = comparisongroups{c,2};
    checkNormality.isnormalComparisongroups(:,c) = indnormal;
    if ~strcmp(comparisongroups{c,1},lastgroup)
        checkNormality.isnormalComparisongroupNames{c} = comparisongroups{c,1};
    end
    lastgroup = comparisongroups{c,1};
    
    %add pvalues to the median table as well
    for t = 1:length(tablenames)
        medianTable.(tablenames{t})(comparison.pind) = meanTable.(tablenames{t})(comparison.pind);
    end
    checkNormality.normalDistributed_4groups = indnormal;
end

%% put mean, median or percentage values into the grouptable
for ind = 1:length(indexes)
    for t = 1:length(tablenames)
        isnormal = checkNormality.normalDistributed_4groups(t);
        if strcmp('categorical',testName)
            groupTable2.(tablenames{t}){ind} = sprintf('%0.2f%%',meanTable.(tablenames{t})(ind)*100);
            stattestTable.(tablenames{t}){ind} = 'categorical';
        elseif isnormal
            % put together a string with mean and std
            if isnan(meanTable.(tablenames{t})(ind))
                groupTable2.(tablenames{t}){ind} = ' - ';
            else
                groupTable2.(tablenames{t}){ind} = sprintf('%0.2f (+-%0.2f)',meanTable.(tablenames{t})(ind),meanTable.(tablenamesSD{t})(ind));
            end
            stattestTable.(tablenames{t}){ind} = 'normal';
        elseif  ~isnormal
            %string with median and interquartile range
            if isnan(medianTable.(tablenames{t})(ind))
                groupTable2.(tablenames{t}){ind} = ' - ';
            else
                groupTable2.(tablenames{t}){ind} = sprintf('%0.2f (%0.2f)',medianTable.(tablenames{t})(ind),medianTable.(tablenamesIQR{t})(ind));
            end
            stattestTable.(tablenames{t}){ind} = 'notnormal';
        end
    end
end

%% correct p-values 
compinds = zeros(6,1);pinds= zeros(6,1);
rownames = {
    'HBP: control vs hypertensive',...
    'HBP: T2D vs hypertensive T2D',...
    'HBP: Hypertensive vs hypertensive & T2D',...
    'HBP: control vs T2D',...
    'HBP: control vs hypertensive T2D',...
    'HBP: Hypertensive vs T2D'};
rejectBenjH2= nan(1,length(rownames));

for r = 1:length(rownames)
    compinds(r) = find(strcmp(statisticsTable.groupComparisonNames,rownames{r}));
    pinds(r) = find(strcmp(meanTable.Group,rownames{r}));
end
rejectany4hbp = zeros(length(tablenames),2);
for t = 1:length(tablenames)
    pvals = zeros(size(compinds));
    reject = zeros(size(compinds));
    for c = 1:length(compinds)
        pvals(c) = statisticsTable.pvalues{compinds(c)}(t);
        reject(c) = statisticsTable.reject{compinds(c)}(t);
    end
    % correct pvalues Benjamini-Hochberg
    [rejectBenjH2(1:6),resulttableBenjH] = Benjamini_Hochberg(pvals(1:6),rownames(1:6),0.05);
        
    for c = 1:length(compinds)
        statisticsTable.rejectBenjH2{compinds(c)}(t) = rejectBenjH2(c);
    end
    
    %sum pvalues from differenttests in a table
        rejecttables.(tablenames{t}) = table(pvals,rejectBenjH2','VariableNames',...
        {'pOrig','rBenjH2'},'Rownames',rownames);
        rejectany4hbp(t,:) = [sum(reject(1:6)),sum(rejectBenjH2(1:6))];

    %replace table with the corrected values for reject or not
    for c = 1:length(compinds)
        if isnan(rejectBenjH2(c))
            groupTable2.(tablenames{t}){pinds(c)}  = ' - ';
        elseif rejectBenjH2(c) && statisticsTable.pvalues{compinds(c)}(t) < 0.0001 %0.001
            groupTable2.(tablenames{t}){pinds(c)}  = '<0.0001*';%0.001
        elseif rejectBenjH2(c)
            groupTable2.(tablenames{t}){pinds(c)}  = sprintf('%0.4f*',statisticsTable.pvalues{compinds(c)}(t));%0.3f
        else
            groupTable2.(tablenames{t}){pinds(c)}  = sprintf('%0.4f',statisticsTable.pvalues{compinds(c)}(t));
        end
        groupTable2.correctedP{pinds(c)} = 'BenjaminiHochberg';
    end

end
numRejectedHyp.HBP = array2table(rejectany4hbp,'Rownames',tablenames,'VariableNames',{'r','rBenjH2'});


%% Create table for statistics document
tosummarize = {comparisons{1:6}};

statdocumenttable = table;
groups4hbp = {'Controls HBP','T2D HBP','Hypertensive HBP','Hypertensive T2D HBP'};
groupnamesshort = {'C','T2D','HT','T2D+HT'};
statdocumenttable.variable = cell(length(tablenames)*length(tosummarize),1);
statdocumenttable.group = cell(length(tablenames)*length(tosummarize),1);
statdocumenttable.mean = cell(length(tablenames)*length(tosummarize),1);
statdocumenttable.sdgroup = cell(length(tablenames)*length(tosummarize),1);
statdocumenttable.sd = cell(length(tablenames)*length(tosummarize),1);
statdocumenttable.ngroup = cell(length(tablenames)*length(tosummarize),1);
statdocumenttable.n = cell(length(tablenames)*length(tosummarize),1);
n=1;
for t =1:length(tablenames)
    isnormal=checkNormality.isnormalComparisongroups(t,1:6);
    allarenormal = sum(isnormal)==length(isnormal);
    statdocumenttable.variable{n} = tablenames{t};
    for g = 1:length(groups4hbp)
        statdocumenttable.group{n} = groupnamesshort{g};
        statdocumenttable.sdgroup{n} = groupnamesshort{g};
        statdocumenttable.ngroup{n} = groupnamesshort{g};
        groupind = strcmp(medianTable.Group,groups4hbp{g});
        %median, SD,n
        statdocumenttable.n{n} = medianTable.N(groupind);
        if strcmp(testName,'categorical')
            statdocumenttable.mean{n} =sprintf('%0.2f%%',meanTable.(tablenames{t})(groupind)*100);
            statdocumenttable.sd{n} = '-';
        elseif allarenormal
            statdocumenttable.mean{n} = meanTable.(tablenames{t})(groupind);
            statdocumenttable.sd{n} = meanTable.(tablenamesSD{t})(groupind);
        else
            statdocumenttable.mean{n} = medianTable.(tablenames{t})(groupind);
            statdocumenttable.sd{n} = medianTable.(tablenamesIQR{t})(groupind);
        end
        n = n+1;
    end
    n = n+2;%add two empty rows
end

%pvals
compnamesshort = {'C vs HT','T2D vs HT+T2D','HT vs HT+T2D','C vs T2D','C vs HT+T2D','HT vs T2D'};
statdocumenttable.comp = cell(length(tablenames)*length(tosummarize),1);
statdocumenttable.p = cell(length(tablenames)*length(tosummarize),1);
statdocumenttable.testtype = cell(length(tablenames)*length(tosummarize),1);
n=1;
for t = 1:length(tablenames)
    isnormal=checkNormality.isnormalComparisongroups(t,1:6);
    %pvals
    for c = 1:length(tosummarize)
        comp = tosummarize{c}.name;
        statdocumenttable.comp{n} = compnamesshort{c};
        statdocumenttable.p{n} = groupTable2.(tablenames{t}){find(strcmp(meanTable.Group,comp))};
        if strcmp(testName,'categorical')
            statdocumenttable.testtype{n} = [stattestTable.(tablenames{t}){tosummarize{c}.pind},', Benjamini Hochberg'];
        elseif isnormal(c)
            statdocumenttable.testtype{n} = 't-test, Benjamini Hochberg';
        else
            statdocumenttable.testtype{n} = 'Wilcoxon ranked sum, Benjamini Hochberg';
        end
        n = n+1;
    end
end

end



