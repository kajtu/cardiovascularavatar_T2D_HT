function [groupTable2,statstable,meanTable,medianTable,checkMean,checkNormality,meanTableCategorical,checkMeanCategorical,allscapisdata,stattestTable] = cohortTable(patNums,paramdata,plotFolderName)

%% Load data
inputpathbase = split(pwd,'cardiovascularavatar_HEALTH_T2D_HT');
loadpath = [inputpathbase{1},'/cardiovascularavatar_HEALTH_T2D_HT/Data'];

load(fullfile(loadpath,'scapisdata_HEALTH.mat'),'scapisdata_HEALTH')
load(fullfile(loadpath,'bpData.mat'),'bpData');

%% sort scapisdata in the same order at patNums (and remove patients not in patnums)
patinds = zeros(size(patNums));
patindsBPdata = zeros(size(patNums));
for p = 1:length(patNums)
    patientNumnoP=removeUnderscore(patNums(p));
    patientNum = ['P' patientNumnoP{1}];
    patinds(p) = find(strcmp(scapisdata_HEALTH.healthIDs,patientNum));
        patindsBPdata(p) = find(strcmp(bpData.patnums,patNums(p)));
end
scapisdata_HEALTH_cohort = scapisdata_HEALTH(patinds,:);
allscapisdata = scapisdata_HEALTH(patinds,:);

% HR used  (during MRI)
scapisdata_HEALTH_cohort.HR = zeros(size(scapisdata_HEALTH_cohort.HR)); %removes the other HR data
for p = 1:length(patNums)
    scapisdata_HEALTH_cohort.HR(p) = round(60/paramdata{p}.T,1);
end

% BP used (before MRI)
scapisdata_HEALTH_cohort.SBP = bpData.SBP(patindsBPdata);
scapisdata_HEALTH_cohort.DBP = bpData.DBP(patindsBPdata);


%remove cells from numerical values
varnames = scapisdata_HEALTH_cohort.Properties.VariableNames;
for i = 1:length(varnames)
    if iscell(scapisdata_HEALTH_cohort.(varnames{i}))
        if isnumeric(scapisdata_HEALTH_cohort.(varnames{i}){1})
            scapisdata_HEALTH_cohort.(varnames{i}) = [scapisdata_HEALTH_cohort.(varnames{i}){:}]';
        end
    end
end

%% Run statistics
%numeric variables
datanames = {'age','diabDuration','weight','height','BMI','HR','SBP','DBP','SBPhome','DBPhome','EAratio','Eeprimratio','CapillaryP_glucose1','hbA1c','HDL','LDL','TG'};
tablenames =  {'Age','DiabetesDuration','Weight','Height','BMI','HR','SBP','DBP','SBPhome','DBPhome','EAratio','Eeprimratio','CapillaryPGlucose','hbA1c','HDL','LDL','TG'};
[groupTable2,statstable,meanTable,checkMean,checkNormality,~,~,medianTable,stattestTable,statdocumenttable] = runStatisticalComparison('numeric',patNums,scapisdata_HEALTH_cohort,datanames,tablenames);%ttest

%categorical variables
categoricalDatanames = {'sex','smoking','insulinUse'};
categoricalTablenames= {'Sex_men','Smoking','InsulinUse'};
for n = 1:length(categoricalDatanames)
        scapisdata_HEALTH_cohort.(categoricalDatanames{n}) = logical(scapisdata_HEALTH_cohort.(categoricalDatanames{n}));
end
[groupTable2Categorical,statstableCategorical,meanTableCategorical,checkMeanCategorical,checkNormalityCategorical,~,~,medianTableCategorical,stattestTableCat,statdocumenttableCategorical] = runStatisticalComparison('categorical',patNums,scapisdata_HEALTH_cohort,categoricalDatanames,categoricalTablenames);

%combine tables
stattestTable = [stattestTable,stattestTableCat(:,3:end)];
groupTable2 = [groupTable2,groupTable2Categorical(:,4:end)];
statstable = [statstable;statstableCategorical];
checkNormality = [checkNormality;checkNormalityCategorical];
statdocumenttable = [statdocumenttable;statdocumenttableCategorical];

%save table for statistical summary document
writetable(statdocumenttable,fullfile(plotFolderName,'statisticssummary_cohort.xlsx'))
end