function [clusterfinaltable,groupedinbothFinal] = clusterCalculation(paramValues,constants,nGroups,paramNames,constantsNames,patNums,plotFolderName)
% PCA analyis and clustering with k-means and hiearchical clustering --> 2
% clusters.
% Based on the method in https://doi.org/10.1113/JP281845 jones et al 2021

%% Load groups
[groups] = loadGroupIndexes(patNums);

m2_lvs = paramValues(strcmp(paramNames,'m2_LV'),:);
Cpvcs = paramValues(strcmp(paramNames,'Cpvc'),:);

%% Set colors
%magma
addpath '..\Requirements\matplotlib'
magmacols = flip(magma(4));
colors4groups = magmacols;
fourgroups = {groups.C_NT_home,groups.T2D_NT_home,groups.C_HT_home,groups.T2D_HT_home};

%% Prepare the parameter matrix: standardize the values
%centering:
% "Before any of the
% clustering methods are applied, each column is centred
% by subtracting the average of each column from each
% element in that column. " jones et al 2021
% "Because our clinical data and
% optimized parameters had different units within their
% respective matrices, we normalized each clinical measure
% or parameter by its standard deviation." jones et al 2021

% standardize the params
standardizedParams = paramValues';
for c = 1:size(standardizedParams,2)
    standardizedParams(:,c) = (standardizedParams(:,c)-mean(standardizedParams(:,c)))./std(standardizedParams(:,c));
end

%% PCA
% PCA centers the data by itself, but we already did that, so we set it to
% false to not do it twice
[coeff,score,latent,tsquared,explained,mu] = pca(standardizedParams,'Centered',false); 
twofirstpcsexplainedpercent = sum(explained(1:2));
firstPCAexplainedpercent = explained(1);
secondPCAexplainedpercent = explained(2);

pca1= score(:,1);
pca2 = score(:,2);

% Plot result
figure('Visible', 'off','Name','clustering_PCA')
xdim_CM = 17;
ydim_CM = 16;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(1,2,'TileSpacing','loose','Padding','compact');
breakpoint = 9;
nexttile
hold on
plot(1:length(explained),latent,'-*')
xline(breakpoint);
xlabel('Num component')
ylabel('eigenvalue')
title(sprintf('%0.2f%% explained by the first %d components',sum(explained(1:breakpoint)),breakpoint))

nexttile
hold on
plot(1:length(explained),explained,'-*')
plot(1:length(explained),cumsum(explained),'--*')
yline(50);
fiftypercexplaind = find(cumsum(explained)>50,1);
xline(fiftypercexplaind);
xlabel('Num component')
ylabel('explained%')
title(sprintf('%0.2f%% explained by the first %d components',sum(explained(1:fiftypercexplaind)),fiftypercexplaind))
legend({'Explained % per component','Cumulative explained %'},'Position',[0.638865999370241 0.188457301332931 0.30171072432079 0.0570247919106287])

% Create table with all coefficients
pcacoefficienttable = table([coeff(:,1);firstPCAexplainedpercent],[coeff(:,2);secondPCAexplainedpercent],[coeff(:,1)+coeff(:,2);twofirstpcsexplainedpercent],'Rownames',[paramNames,{'% explained'}],'Variablenames',{'PCA1','PCA2','Sum'});
pcacoefficienttablefull = array2table(coeff,'Rownames',paramNames);

% Create biplot to show spread of coefficients in the first 2 PCs
figure('Visible', 'off','Name','cluster_PCA_biplot')
xdim_CM = 17; 
ydim_CM = 16;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])

h=biplot(coeff(:,1:2),'Scores',score(:,1:2),'VarLabels',paramNames,'MarkerSize',15,'color',[0.3 0.3 0.3]);
hold on
idx = zeros(size(fourgroups{1}));
idx(fourgroups{1}) = 1;
idx(fourgroups{2}) = 2;
idx(fourgroups{3}) = 3;
idx(fourgroups{4}) = 4;
for p = 1:size(standardizedParams,1)
    h(p + size(standardizedParams,2)*3).MarkerEdgeColor = colors4groups(idx(p),:); % change the color of data 
end

%% hiearchical clustering
% "To do this in MATLAB, the linkage function is used, and the Ward metric (Ward, 1963) is selected to group the two clusters together at each step that minimize the total in-cluster variation. Using the dendrogram, we partitioned our patients into two clusters by cutting the dendrogram halfway between the second-from-last and last linkages."jones et al 2021 
Y = pdist(standardizedParams);
Z = linkage(Y,'ward');
%Column 1 and 2 of Z contain cluster indices linked in pairs to form a binary tree. 
%"Using the dendrogram, we
% partitioned our patients into two clusters by cutting the
% dendrogram halfway between the second-from-last and
% last linkages." jones et al 2021
hiearchicalClusters = cluster(Z,'maxclust',nGroups);

figure('Visible', 'off','Name','Clustering_hiearchical_Dendrogram')
dendrogram(Z)

%make sure the largest cluster is named cluster 1 to be able to more easily
%compare it to the k-means clustering
numones = sum(hiearchicalClusters == 1);
numtwoes = sum(hiearchicalClusters == 2);
if numtwoes>numones
    hiearchicalClusters(hiearchicalClusters == 1) = 20;%will be set to 2
    hiearchicalClusters(hiearchicalClusters == 2) = 1;
    hiearchicalClusters(hiearchicalClusters == 20) = 2;
end


%% k-means clustering
%"This method is dependent on the random initial cluster centroids selected, 
%so we run this process 20 times and select the clustering result that has
%the smallest total cluster variance" jones et al 2021
numrepeats = 100;%20 was not stable enough with this data, increased to 100.
clusters = cell(1,numrepeats);
withinClusterDistancesum = zeros(1,numrepeats);
for i = 1:numrepeats
    [clusters{i},C,sumd] = kmeans(standardizedParams,nGroups);
    withinClusterDistancesum(i) = sum(sumd);
end
[mindistanceKmeans,mindistanceInd] = min(withinClusterDistancesum);
kmeansClusters = clusters{mindistanceInd};

%make sure the largest cluster is named cluster 1 to be able to more easily
%compare it to the hierarchical clustering
numones = sum(kmeansClusters == 1);
numtwoes = sum(kmeansClusters == 2);
if numtwoes>numones
    kmeansClusters(kmeansClusters == 1) = 20;%will be set to 2
    kmeansClusters(kmeansClusters == 2) = 1;
    kmeansClusters(kmeansClusters == 20) = 2;
end

%% Combining hearchical and kmeans clusters
%if clustered into the same group in both hiearchical and kmeans clusters,
%they are in the same group
groupedinboth = nan(length(patNums),1);
for p = 1:length(patNums)
    for pcomp = 1:length(patNums)
        %if not looked at before and not comparing a patient with itself
        if isnan(groupedinboth(p)) &&(pcomp ~= p)
            %if they are in same group in both methods
            if (hiearchicalClusters(p) == hiearchicalClusters(pcomp)) && ...
                    (kmeansClusters(p) == kmeansClusters(pcomp))
                groupedinboth(p) = str2double([num2str(kmeansClusters(p)) num2str(hiearchicalClusters(p))]);
                groupedinboth(pcomp) = str2double([num2str(kmeansClusters(p)) num2str(hiearchicalClusters(p))]);
                break
            end
        end
    end
end
allgroups = unique(groupedinboth(~isnan(groupedinboth)));

%% Using the convex pca hull to check if subjects are HT+T2D-like (creating a  HT+T2D group vs a non-HT+T2D group)
a = pca1(groups.T2D_HT_home);
b = pca2(groups.T2D_HT_home);
pcaHullT2DHT = convhull(a,b); %computes the 2-D convex hull of the points in column vectors x and y.
inHTT2Dhull = inpolygon(pca1,pca2,a(pcaHullT2DHT),b(pcaHullT2DHT));

HTgroups = zeros(size(inHTT2Dhull));
HTgroups(inHTT2Dhull) = 4;%within the hull of hypertension AND diabetes
HTgroups(groups.T2D_HT_home) = 3;%hypertension AND diabetes only

groupedinbothFinal = nan(length(patNums),1);
for p = 1:length(patNums)
    for pcomp = 1:length(patNums)
        %if not looked at before and not comparing a patient with itself
        if isnan(groupedinbothFinal(p)) &&(pcomp ~= p)
            %if they are in same group in both methods
            if (HTgroups(p) == HTgroups(pcomp)) && ...
                    (groupedinboth(p) == groupedinboth(pcomp))
                groupedinbothFinal(p) = str2double([num2str(groupedinboth(p)) num2str(HTgroups(p))]);
                groupedinbothFinal(pcomp) = str2double([num2str(groupedinboth(p)) num2str(HTgroups(p))]);
                break
            end
        end
    end
end

%find clusters with HT+T2D subjects in, and select the largest one as HTT2Dcluster
HTT2Dclusters = [113,223];%we do not want 123 or 213 since they are not consistently clustered.
groupsize = size(HTT2Dclusters);
for i = 1:length(HTT2Dclusters)
    groupsize(i) = sum(groupedinbothFinal == HTT2Dclusters(i));
end
[numInHTT2Dcluster,indmax] = max(groupsize);
HTT2Dcluster = HTT2Dclusters(indmax);

diseaselikecluster = HTT2Dcluster-3+4;%should end with 4 but the first two groupings should be the same as the HT+T2D-like cluster
numIndiseaselikecluster = sum(groupedinbothFinal == diseaselikecluster);

if HTT2Dcluster  == 113
    numHTT2D_notinHtT2Dcluster = sum(groupedinbothFinal == 223);
    nondiseasecluster = 220;
else
    numHTT2D_notinHtT2Dcluster = sum(groupedinbothFinal == 113);
    nondiseasecluster = 110;
end
numInnondiseaselikecluster = sum(groupedinbothFinal == nondiseasecluster);

numHTT2D_NCC = sum(groups.T2D_HT_home) - numHTT2D_notinHtT2Dcluster - numInHTT2Dcluster;

%inconsistently clustered with k-means and hiearchical
NCCdiseaselikecluster = HTT2Dcluster-3;%in disease cluster but not in disease hull (ends with 0)
NCCnondiseaselikecluster = nondiseasecluster+4;%not in disease cluster but in disease hull (ends with 4)
NCCclusters = [120,210,124,214,NCCdiseaselikecluster,NCCnondiseaselikecluster];

%summarize new groups
finalclusters = cell(size(groupedinbothFinal));
for f = 1:length(finalclusters)
    if HTgroups(f) == 3
        finalclusters{f} = 'HTT2D';
    elseif ismember(groupedinbothFinal(f),NCCclusters)
        finalclusters{f} = 'NCC';
    elseif groupedinbothFinal(f) == nondiseasecluster
        finalclusters{f} = 'non-HTT2D-like';
    elseif groupedinbothFinal(f) == diseaselikecluster
        finalclusters{f} = 'HTT2D-like';
    elseif isnan(groupedinbothFinal(f))
        finalclusters{f} = 'NCC';
    else
        disp(groupedinbothFinal(f))
    end
end
numinNCCcluster = sum(ismember(finalclusters,'NCC'));

clusternamesfinal = {'NCC','non-HTT2D-like','HTT2D-like','HTT2D'};
clustercolorsfinal = {[0.5 0.5 0.5],[0.2 0.5 0.2],[0.2 0.2 0.5],colors4groups(4,:)};
clustermarkersfinal = {'+','^','o','*'};
clustermarkers = {'*','+','x','o','^','v','<','>'};

clustercolors = {[0.1 0.4 0.1],[0.2 0.2 0.5],[0.5 0.2 0.2],[0.2 0.2 0.2],[0.2 0.2 0.2]};

% Create result table
clusterfinaltable = table(groups.fourgroupsHBP,groupedinbothFinal,kmeansClusters,hiearchicalClusters,inHTT2Dhull,finalclusters,'Rownames',patNums,'Variablenames',{'Subject group','kmeans, hiearchical,hull','kmeans','hiearchical','PCA hull','final clustered group'});
clusterfinaltablesorted = sortrows(clusterfinaltable,1);

%% Plot final clusters
figure('Visible', 'off','Name','clusters_final')
global fzlarge 
global fzsmall 
global fzsupersmall 
fzlarge = 14;
fzsmall = 10;
fzsupersmall = 8;
set(gcf,'Color','white')
xdim_CM = 17;
ydim_CM = 14;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(2,2,'TileSpacing','loose','Padding','compact');

%PCAs
nexttile
hold on
for i = 1:length(fourgroups)
    a = pca1(fourgroups{i});
    b = pca2(fourgroups{i});
    hullN = convhull(a,b); 
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1);
end
for i = 1:length(allgroups)
    plot(pca1(groupedinboth==allgroups(i)),pca2(groupedinboth==allgroups(i)),clustermarkers{i},'color',[0.5 0.5 0.5],'LineWidth',1,'Markersize',5);
end
set(gca,'FontSize',fzsmall)
xlabel('1st Principal Component','FontSize',fzsmall)
ylabel('2nd Principal Component','FontSize',fzsmall)
hold off

%cpvc
nexttile
hold on
for i = 1:length(fourgroups)
    a = m2_lvs(fourgroups{i});
    b = Cpvcs(fourgroups{i});
    hullN = convhull(a,b); 
    diseasegroups(i) = fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1);
end
for i = 1:length(allgroups) 
    clustergroups(i)=plot(m2_lvs(groupedinboth==allgroups(i)),Cpvcs(groupedinboth==allgroups(i)),clustermarkers{i},'color',[0.5 0.5 0.5],'LineWidth',1,'Markersize',5);
    namestring = num2str(allgroups(i));
    clusternames{i} = sprintf('Cluster %d (%s,%s)',i,namestring(1),namestring(2));
end
set(gca,'FontSize',fzsmall)
xlabel('LV relaxation rate (-)','FontSize',fzsmall)
ylabel(sprintf('Pulmonary venous compliance\n(mmHg/mL)'),'FontSize',fzsmall)
legend([clustergroups,diseasegroups],[clusternames,{'Control hull','T2D hull','HT hull','HT & T2D hull'}],'FontSize',7,...
    'Position',[0.802190942086316 0.833414502411537 0.144634523714312 0.129629626042313])
hold off


nexttile
hold on
for i = 1:length(fourgroups)
    a = pca1(fourgroups{i});
    b = pca2(fourgroups{i});
    hullN = convhull(a,b); 
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1)
end
for i = 1:length(clusternamesfinal)
    a = pca1(ismember(finalclusters,clusternamesfinal{i}));
    b = pca2(ismember(finalclusters,clusternamesfinal{i}));
    plot(a,b,clustermarkersfinal{i},'color',clustercolorsfinal{i},'LineWidth',1,'Markersize',5);
end
set(gca,'FontSize',fzsmall)
xlabel('1st Principal Component','FontSize',fzsmall)
ylabel('2nd Principal Component','FontSize',fzsmall)
hold off

%cpvc
nexttile
hold on
for i = 1:length(fourgroups)
    a = m2_lvs(fourgroups{i});
    b = Cpvcs(fourgroups{i});
    hullN = convhull(a,b); 
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1)
end
for i = 1:length(clusternamesfinal)
    a = m2_lvs(ismember(finalclusters,clusternamesfinal{i}));
    b = Cpvcs(ismember(finalclusters,clusternamesfinal{i}));
    plot(a,b,clustermarkersfinal{i},'color',clustercolorsfinal{i},'LineWidth',1,'Markersize',5);
end
set(gca,'FontSize',fzsmall)
xlabel('LV relaxation rate (-)','FontSize',fzsmall)
ylabel(sprintf('Pulmonary venous compliance\n(mmHg/mL)'),'FontSize',fzsmall)
l = legend([{'Control hull','T2D hull','HT hull','T2D+HT hull'},clusternamesfinal],...
    'FontSize',7,'NumColumns',2,'Position',[0.683478143201293 0.382888873886038 0.24 0.0840000000000001]);%[0.22 0.57 0.24 0.084]


%% PLOT all clusters
global fzlarge 
global fzsmall 
global fzsupersmall 
fzlarge = 14;
fzsmall = 10;
fzsupersmall = 8;

clustercols = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};

figure('Visible', 'off','Name','Clustering')
set(gcf,'Color','white')
xdim_CM = 21;
ydim_CM = 25;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(3,2,'TileSpacing','tight','Padding','compact');
clf
%----------- kmeans and hiearchical plotted against PCAs ---------------
nexttile
title('Grouped on kmeans','FontSize',fzlarge)
hold on
clusternames = cell(1,nGroups);
for i = 1:length(fourgroups)
    a = pca1(fourgroups{i});
    b = pca2(fourgroups{i});
    hullN = convhull(a,b);
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1.5)
end
for i = 1:nGroups
    hold on
    clusternames{i} = num2str(i);
    clusterplot(i) = plot(pca1(kmeansClusters==i),pca2(kmeansClusters==i),'*','LineWidth',2,'color',clustercols{i});
end
set(gca,'FontSize',fzsmall)
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
l = legend(clusterplot,clusternames,'FontSize',fzsupersmall);
l.BoxFace.ColorType='truecoloralpha';
l.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
hold off

nexttile
title('Grouped on hiearchical','FontSize',fzlarge)
hold on
for i = 1:length(fourgroups)
    a = pca1(fourgroups{i});
    b = pca2(fourgroups{i});
    hullN = convhull(a,b); 
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1.5)
end
for i = 1:nGroups
    hold on
    clusterplot(i) = plot(pca1(hiearchicalClusters==i),pca2(hiearchicalClusters==i),'*','LineWidth',2,'color',clustercols{i});
    clusternames{i} = num2str(i);
end
set(gca,'FontSize',fzsmall)
xlabel('1st Principal Component','FontSize',fzsmall)
ylabel('2nd Principal Component','FontSize',fzsmall)
l = legend(clusterplot,clusternames,'FontSize',fzsupersmall);
l.BoxFace.ColorType='truecoloralpha';
l.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
hold off

%--------------- kmeans and hiearcical plotted against selected parameters ---------------
nexttile
hold on
title('Grouped on kmeans','FontSize',fzlarge)
for i = 1:length(fourgroups)
    a = Cpvcs(fourgroups{i});
    b = m2_lvs(fourgroups{i});
    hullN = convhull(a,b); 
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1.5)
end
for i = 1:nGroups
    clusterplot(i) =plot(Cpvcs(kmeansClusters==i),m2_lvs(kmeansClusters==i),'*','LineWidth',2,'color',clustercols{i});
    clusternames{i} = num2str(i);
end
set(gca,'FontSize',fzsmall)
xlabel('Cpvc','FontSize',fzsmall)
ylabel('m2 LV','FontSize',fzsmall)
l = legend(clusterplot,clusternames,'FontSize',fzsupersmall);
l.BoxFace.ColorType='truecoloralpha';
l.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
hold off

nexttile
hold on
title('Grouped on hiearchical','FontSize',fzlarge)
for i = 1:length(fourgroups)
    a = Cpvcs(fourgroups{i});
    b = m2_lvs(fourgroups{i});
    hullN = convhull(a,b); 
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1.5)
end
for i = 1:nGroups
     clusterplot(i) = plot(Cpvcs(hiearchicalClusters==i),m2_lvs(hiearchicalClusters==i),'*','LineWidth',2,'color',clustercols{i});
     clusternames{i} = num2str(i);
end
set(gca,'FontSize',fzsmall)
xlabel('Cpvc','FontSize',fzsmall)
ylabel('m2 LV','FontSize',fzsmall)
l = legend(clusterplot,clusternames,'FontSize',fzsupersmall);
l.BoxFace.ColorType='truecoloralpha';
l.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
hold off

% --------------- Grouped in both methods (final clusters) ---------------
nexttile
hold on
title('Grouped in both methods','FontSize',fzlarge)
for i = 1:length(fourgroups)
    a = Cpvcs(fourgroups{i});
    b = m2_lvs(fourgroups{i});
    hullN = convhull(a,b); 
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.1,'EdgeColor',colors4groups(i,:),'LineWidth',1.5)
end
for i = 1:length(allgroups) 
    clusterplot(i) =plot(Cpvcs(groupedinboth==allgroups(i)),m2_lvs(groupedinboth==allgroups(i)),clustermarkersfinal{i},'color',clustercolors{i},'LineWidth',1,'MarkerSize',5);
    clusternames{i} = num2str(allgroups(i));
end
set(gca,'FontSize',fzsmall)
xlabel('Cpvc','FontSize',fzsmall)
ylabel('m2 LV','FontSize',fzsmall)
l = legend(clusterplot,clusternames,'FontSize',fzsupersmall);
l.BoxFace.ColorType='truecoloralpha';
l.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
hold off


nexttile
title('Grouped in both methods','FontSize',fzlarge)
hold on
for i = 1:length(fourgroups)
    a = pca1(fourgroups{i});
    b = pca2(fourgroups{i});
    hullN = convhull(a,b); %computes the 2-D convex hull of the points in column vectors x and y.
    fill(a(hullN),b(hullN),colors4groups(i,:),'FaceAlpha',0.07,'EdgeColor',colors4groups(i,:),'LineWidth',1)
end

for i = 1:length(allgroups)
    clusterplot(i) = plot(pca1(groupedinboth==allgroups(i)),pca2(groupedinboth==allgroups(i)),clustermarkersfinal{i},...
        'color',clustercolors{i},'LineWidth',1,'Markersize',5);
    clusternames{i} = num2str(allgroups(i));
end
set(gca,'FontSize',fzsmall)
xlabel('1st Principal Component','FontSize',fzsmall)
ylabel('2nd Principal Component','FontSize',fzsmall)
l = legend(clusterplot,clusternames,'FontSize',fzsupersmall);
l.BoxFace.ColorType='truecoloralpha';
l.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');

% Add legend
nums = [numIndiseaselikecluster;numInnondiseaselikecluster;numinNCCcluster;numInHTT2Dcluster;numHTT2D_notinHtT2Dcluster;numHTT2D_NCC];
cluster_counts = table([nums;sum(nums)],...
    'Rownames',{'HT+T2D-like cluster','non-HT+T2D-like cluster','NCC cluster','HT+T2D in HT+T2D cluster','HT+T2D not in HT+T2D cluster','HT+T2D in NCC cluster','Sum'},'Variablenames',{'Number of subjects'});

%% Save figures and tables
saveAllFigures(plotFolderName)

writetable(pcacoefficienttablefull,fullfile(plotFolderName,'cluster_pcacoefficients_all.xlsx'),"WriteRowNames",1)
writetable(pcacoefficienttable,fullfile(plotFolderName,'cluster_pcacoefficients.xlsx'),"WriteRowNames",1)
writetable(clusterfinaltablesorted,fullfile(plotFolderName,'clusterfinaltablesorted.xlsx'),"WriteRowNames",1)
writetable(cluster_counts,fullfile(plotFolderName,'cluster_counts.xlsx'),"WriteRowNames",1)


