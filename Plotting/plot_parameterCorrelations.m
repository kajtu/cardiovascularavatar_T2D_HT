function plot_parameterCorrelations(paramValues,constants,paramNames,simLast,ind,data,extradata,patNums,ynames,plotFolderName)

%% Set values
%Pre-define vectors
Caas = zeros(1,length(patNums));
EmaxLVs = zeros(1,length(patNums));
Cpvcs = zeros(1,length(patNums));
Raos = zeros(1,length(patNums));
m2LVs = zeros(1,length(patNums));
m1LVs = zeros(1,length(patNums));
k_syst_LVs=zeros(1,length(patNums));
k_diast_LVs=zeros(1,length(patNums));

HbA1c=zeros(1,length(patNums));HR = zeros(1,length(patNums));
diabduration=zeros(1,length(patNums));
SBPhome=zeros(1,length(patNums));DBPhome=zeros(1,length(patNums));
SBPduringMRI_data=zeros(1,length(patNums));

% Get values for each subject
for p = 1:length(patNums)
    % Parameters
    EmaxLVs(p) = paramValues(ind.Emax_LV,p);
    Raos(p) = paramValues(strcmp(paramNames,'Rao'),p);
    Caas(p) = paramValues(ind.Caa,p);
    m2LVs(p) = paramValues(strcmp(paramNames,'m2_LV'),p);
    m1LVs(p) = paramValues(strcmp(paramNames,'m1_LV'),p);
    k_syst_LVs(p) = paramValues(strcmp(paramNames,'k_syst_LV'),p);
    k_diast_LVs(p) = paramValues(strcmp(paramNames,'k_diast_LV'),p);
    Cpvcs(p) =  paramValues(strcmp(paramNames,'Cpvc'),p);
    
    % Long- and shortterm markers
    HR(p) = 60/constants(ind.T,p);
    HbA1c(p) = extradata{p}.SCAPIS.hbA1c;
    diabduration(p) = extradata{p}.SCAPIS.diabDuration;
    SBPhome(p) = extradata{p}.HomeBP.SBPmean;
    DBPhome(p) = extradata{p}.HomeBP.DBPmean;
    SBPduringMRI_data(p)=data{p}.SBP;
end

%% Calculate and plot correlations
figure('Visible', 'off','Name','Correlations_summary');
set(gcf,'Color','white')
xdim_CM = 20;
ydim_CM = 15;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
t = tiledlayout(3,2,'TileSpacing','Compact','Padding','compact');
corrvarnames = {'m2LV','Cpvc','k_s_y_s_t_L_V','k_d_i_a_s_t_L_V','m1LV','EmaxLV','Rao','Caa'};
corrvars = {m2LVs,Cpvcs, k_syst_LVs,k_diast_LVs,m1LVs,EmaxLVs,Raos,Caas};
p = zeros(6,length(corrvars)); 
r2 = zeros(6,length(corrvars));
rho = zeros(6,length(corrvars));
numobs = zeros(6,length(corrvars));
testtype=cell(6,length(corrvars));
nexttile(t);
[p(1,:),r2(1,:),rho(1,:),numobs(1,:),testtype(1,:)] = createSummaryCorrelationPlot('Correlation with heartrate during MRI',corrvarnames,corrvars,HR);

nexttile(t);
x = SBPduringMRI_data;
[p(2,:),r2(2,:),rho(2,:),numobs(2,:),testtype(2,:)] =createSummaryCorrelationPlot('Correlation with SBP during MRI',corrvarnames,corrvars,x);

nexttile(t);
[p(3,:),r2(3,:),rho(3,:),numobs(3,:),testtype(3,:)] =createSummaryCorrelationPlot('Correlation with HbA1c',corrvarnames,corrvars,HbA1c);

nexttile(t);
[p(4,:),r2(4,:),rho(4,:),numobs(4,:),testtype(4,:)] =createSummaryCorrelationPlot('Correlation with Diabetes duration',corrvarnames,corrvars,diabduration);

nexttile(t);
[p(5,:),r2(5,:),rho(5,:),numobs(5,:),testtype(5,:)] =createSummaryCorrelationPlot('Correlation with home SBP',corrvarnames,corrvars,SBPhome);

nexttile(t);
[p(6,:),r2(6,:),rho(6,:),numobs(6,:),testtype(6,:)] =createSummaryCorrelationPlot('Correlation with home DBP',corrvarnames,corrvars,DBPhome);

%% Save figures
saveAllFigures(plotFolderName)

%% Create summary table of the correlations
correlationTablesummary = table;
corrnames = {'HR','SBP MR','Hba1c','diabdur','SBPhome','DBPhome'};
correlationTablesummary.corr = cell(length(corrnames)*length(corrvarnames),1);
correlationTablesummary.corrvar = cell(length(corrnames)*length(corrvarnames),1);
correlationTablesummary.n = cell(length(corrnames)*length(corrvarnames),1);
correlationTablesummary.rho = cell(length(corrnames)*length(corrvarnames),1);
correlationTablesummary.p = cell(length(corrnames)*length(corrvarnames),1);
correlationTablesummary.r2 = cell(length(corrnames)*length(corrvarnames),1);
correlationTablesummary.testtype = cell(length(corrnames)*length(corrvarnames),1);
n=1;
for c = 1:length(corrnames)
    correlationTablesummary.corr{n} = corrnames{c};
    for v = 1:length(corrvarnames)
        correlationTablesummary.corrvar{n} = corrvarnames{v};
        if p(c,v) < 0.0001
            correlationTablesummary.p{n} = '<0.0001';
        else
            correlationTablesummary.p{n} = sprintf('%0.4f',p(c,v));
        end
        correlationTablesummary.n{n} = numobs(c,v);
        correlationTablesummary.testtype{n} = testtype{c,v};
        correlationTablesummary.r2{n} = sprintf('%0.4f',r2(c,v));
        correlationTablesummary.rho{n} = sprintf('%0.4f',rho(c,v));
        n=n+1;
    end
end

%save table for statistical summary document
writetable(correlationTablesummary,fullfile(plotFolderName,'statisticssummary_correlations.xlsx'))


end
