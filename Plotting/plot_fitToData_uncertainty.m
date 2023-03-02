function plot_fitToData_uncertainty(minmaxSim,patNums,bestparams,paramuncertainty,paramNames,data,inds,bestcost,plotFolderName)

%% Define colors, setup
paramcol = [21,96,122]./255;%blue
paramcol2 = [255 181 95]./255;%yellow
colorParams = paramcol;

bluescale = {[208,209,230]./255,[200,188,220]./255,[166,189,219]./255,[116,169,207]./255,[54,144,192]./255,[5,112,176]./255,[3,78,123]./255};
namesavatar = {'pulmonary','LA','MV','LV','AV','AA','peripheral'};
colororder = {'MV','AV','AA','pulmonary','LA','LV','peripheral'};
[~,colind] = ismember(colororder,namesavatar);
bluescaleordered = bluescale(colind);
colors = bluescaleordered;

% set font sizes
fz1 = 10;
fzsmall = 7;
letters = 'A':'P';

%load groups
[groups] = loadGroupIndexes(patNums);

%% Find the subjects with worst, median, and best fit to data
comparisoncost = bestcost;
% best
[bestcost,pbest] = min(comparisoncost);
% median
pmedianorder = round(length(comparisoncost)/2);
[sorttotalcost,sortpatstotalcost] = sort(comparisoncost);
mediancost = sorttotalcost(pmedianorder);
pmedian = sortpatstotalcost(pmedianorder);
% worst
[worstcost,pworst] = max(comparisoncost);

%% Find the subjects with best fit to data within each group
% hypertension is defined both during MR and home BP
% best t2d+hypertensive
costT2Dhyp = comparisoncost(groups.T2D_HT_home & groups.hypertensionMR);
[bestcostT2Dh,pbestT2Dh] = min(costT2Dhyp);
patNumsT2DHT = patNums(groups.T2D_HT_home & groups.hypertensionMR);
pbestT2Dh=patNumsT2DHT{pbestT2Dh};

% best t2d
costT2Dnormal = comparisoncost(groups.T2D_NT_home & groups.nothypertensionMR);
[bestcostT2Dn,pbestT2Dn] = min(costT2Dnormal);
patNumsT2D = patNums(groups.T2D_NT_home & groups.nothypertensionMR);
pbestT2Dn=patNumsT2D{pbestT2Dn};

% best control
costCnormal = comparisoncost(groups.C_NT_home & groups.nothypertensionMR);
[bestcostCn,pbestCn] = min(costCnormal);
patNumsC= patNums(groups.C_NT_home & groups.nothypertensionMR);
pbestCn=patNumsC{pbestCn};

% best hypertensive
costChigh = comparisoncost(groups.C_HT_home & groups.hypertensionMR);
[bestcostCh,pbestCh] = min(costChigh);
patNumsHT= patNums(groups.C_HT_home & groups.hypertensionMR);
pbestCh=patNumsHT{pbestCh};


save('../Parameters/bestpatients.mat','pbestT2Dh','pbestT2Dn','pbestCn','pbestCh')


%% Plot the blood flow simulation vs data for 3 example subjects
selectedNums = [pbest,pmedian,pworst];
titles = {sprintf('Best fit to data (V=%0.1f)',bestcost),...
    sprintf('Median fit to data (V=%0.1f)',mediancost),...
    sprintf('Worst fit to data (V=%0.1f)',worstcost)};

figure('Visible', 'off','Name','Flows - data vs simulation 3 examples uncertainty');
set(gcf,'Color','white')
xdim_CM = 21;
ydim_CM = 6;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiles=tiledlayout(1,3,'TileSpacing','loose','Padding','compact');
for pat = 1:length(selectedNums)
    p = selectedNums(pat);
    nexttile(tiles)
    hold on
    simtime = minmaxSim.time{p};
    datatime = data{p}.MV(:,1);
    tend = datatime(end);
    tstart = datatime(1);
    indd = data{p}.indtdiast;
    
    %Plot the data
    cdat(1)=plot(datatime(indd:end,1),data{p}.MV(indd:end,2),'o','color',colors{1},'LineWidth',1.5);
    cdat(2)=plot(datatime(1:indd),data{p}.AV(1:indd,2),'^','color',colors{2},'LineWidth',1.5);
    cdat(3)=plot(datatime(1:indd),data{p}.AC(1:indd,2),'square','color',colors{3},'LineWidth',1.5);
    cdat(4)=plot(datatime,data{p}.PV(:,2),'diamond','color',colors{4},'LineWidth',1.5);
    
    % Plot the simulation uncertainty
    t = [simtime;flipud(minmaxSim.time{p})];
    MV = [minmaxSim.States.min{p}(:,inds{p}.MV);flipud(minmaxSim.States.max{p}(:,inds{p}.MV))];
    fill(t,MV,colors{1},'FaceAlpha',0.6,'EdgeColor','none')
    AA = [minmaxSim.States.min{p}(:,inds{p}.AC);flipud(minmaxSim.States.max{p}(:,inds{p}.AC))];
    fill(t,AA,colors{3},'FaceAlpha',0.8,'EdgeColor','none')
    AV = [minmaxSim.States.min{p}(:,inds{p}.AV);flipud(minmaxSim.States.max{p}(:,inds{p}.AV))];
    fill(t,AV,colors{2},'FaceAlpha',0.6,'EdgeColor','none')
    PV = [minmaxSim.States.min{p}(:,inds{p}.PV);flipud(minmaxSim.States.max{p}(:,inds{p}.PV))];
    fill(t,PV,colors{4},'FaceAlpha',0.6,'EdgeColor','none')
    
    % Plot the simulation that best fit to data
    plot(simtime,minmaxSim.States.best{p}(:,inds{p}.MV),'-','color',0.9.*colors{1},'LineWidth',2)
    plot(simtime,minmaxSim.States.best{p}(:,inds{p}.AV),'-','color',0.9.*colors{2},'LineWidth',2)
    plot(simtime,minmaxSim.States.best{p}(:,inds{p}.AC),'-','color',0.9.*colors{3},'LineWidth',2)
    plot(simtime,minmaxSim.States.best{p}(:,inds{p}.PV),'-','color',0.9.*colors{4},'LineWidth',2)
    
    xlim([tstart tend])
    xticks(0:0.1:tend)
    alllabels = strsplit(num2str(0:0.1:tend));
    alllabels(2:2:end-2) = {''};
    xticklabels(alllabels)
    xtickangle(0)
    ylim([-100 425])
    yticks(-100:100:425)
    yticklabels({'-100','0','','200','','400','','600'})
    
    set(gca,'FontSize',fz1)
    xlabel('Time (s)','FontSize',fz1)
    ylabel('Flow volume (ml/s)','FontSize',fz1)
    title(titles{pat},'FontSize',fz1)
    TilePos = tiles.Children.InnerPosition;
    letter = annotation('textbox',[TilePos(1)-0.045 TilePos(2)+1.15*TilePos(4) .015 .015],'String',letters(pat),'Linestyle','none','FitBoxToText','on','BackgroundColor','none');
end
s=plot(NaN,NaN,'k-','linewidth',2);
legend([s,cdat],{'Simulation','Data MV','Data AV','Data AA','Data PV'},'FontSize',8,'NumColumns',2,...
                'Position',[0.120675725046687 0.719565923782884 0.2 0.17])
legend('boxoff')


%% Plot the fit to data in the one subjects with best overall fit to data (used in the method figure)
% Fit to parameter data
figure('Visible', 'off','Name','Parameters fit to data (best fit)');
xdim_CM = 10; 
ydim_CM = 6;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
dataparams = {'Emax_LV','Caa','Ctot','Rtot','ElCo'};
pat = pbest;
ind = inds{pat};
ymin = 0;
ymax = 0;
hold on
[lb,ub] = loadParamBounds_HEALTH(ind,data{pat}.params);
lb = log10(lb);ub=log10(ub);

maxparpat = paramuncertainty.maxValuesparams(:,pat);
minparpat = paramuncertainty.minValuesparams(:,pat);
plotparpat = bestparams(:,pat);
for p = 1:length(dataparams)
    pvalmin = minparpat(ind.(dataparams{p}));
    pvalmax = maxparpat(ind.(dataparams{p}));
    pvalbest = plotparpat(ind.(dataparams{p}));
    dval = data{pat}.params.(dataparams{p});
    emin = abs(log10(pvalbest)-log10(pvalmin));
    emax = abs(log10(pvalmax)-log10(pvalbest));
    errorbar(p,log10(pvalbest),emin,emax,'.','color',colorParams,'LineWidth',1.5,'MarkerSize',20)
    plot(p,log10(dval(1)),'ko','LineWidth',1)
    plot([p-0.5,p+0.5],[lb(ind.(dataparams{p})),lb(ind.(dataparams{p}))],'k--')
    plot([p-0.5,p+0.5],[ub(ind.(dataparams{p})),ub(ind.(dataparams{p}))],'k--')
    ymin = min(ymin,lb(ind.(dataparams{p})));
    ymax = max(ymax,ub(ind.(dataparams{p})));
end
%sbp
emin = abs(log10(minmaxSim.SBP.best{p})-log10(minmaxSim.SBP.min{p}));
emax = abs(log10(minmaxSim.SBP.max{p})-log10(minmaxSim.SBP.best{p}));
errorbar(length(dataparams)+1,log10(minmaxSim.SBP.best{p}),emin,emax,'.','color',colorParams,'LineWidth',1.5,'MarkerSize',20)
plot(length(dataparams)+1,log10(data{pat}.SBP),'ko','LineWidth',1)

%dbp
emin = log10(minmaxSim.DBP.best{p})-log10(minmaxSim.DBP.min{p});
emax = log10(minmaxSim.DBP.max{p})-log10(minmaxSim.DBP.best{p});
errorbar(length(dataparams)+2,log10(minmaxSim.DBP.best{p}),max(0,emin),max(0,emax),'.','color',colorParams,'LineWidth',1.5,'MarkerSize',20)
plot(length(dataparams)+2,log10(data{pat}.DBP),'ko','LineWidth',1)
ymax = max([ymax,log10(data{pat}.SBP),log10(minmaxSim.SBP.max{p})]);
ymin = min([ymin,log10(data{pat}.DBP),log10(minmaxSim.DBP.min{p})]);
    
ylim([ymin*1.1,ymax*1.1])
xlim([0,length(dataparams)+1+2])
xticks(1:length(dataparams)+2)
dataparamsNames = {'Emax LV','Caa','Ctot','Rtot','ElCo','SBP','DBP'};
xticklabels(dataparamsNames);
ylabel(sprintf('Parameter value\n(log10)'))
set(gca,'FontSize',12)

% Fit to blood flow data
figure('Visible', 'off','Name','Flows fit to data (best fit)');
set(gcf,'Color','white')
xdim_CM = 24; 
ydim_CM = 6;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
p = pbest;
hold on
simtime = minmaxSim.time{p};
datatime = data{p}.time;
tend = datatime(end);
tstart = datatime(1);
indd = data{p}.indtdiast;
% Plot data
plot(datatime(indd:end),data{p}.MV(indd:end,2),'o','color',colors{1},'LineWidth',1.5);
plot(datatime(1:indd),data{p}.AV(1:indd,2),'^','color',colors{2},'LineWidth',1.5);
plot(datatime(1:indd),data{p}.AC(1:indd,2),'square','color',colors{3},'LineWidth',1.5);
plot(datatime,data{p}.PV(:,2),'diamond','color',colors{4},'LineWidth',1.5);

% Plot simulation uncertainty
t = [simtime;flipud(minmaxSim.time{p})];
MV = [minmaxSim.States.min{p}(:,inds{p}.MV);flipud(minmaxSim.States.max{p}(:,inds{p}.MV))];
fill(t,MV,colors{1},'FaceAlpha',0.6,'EdgeColor','none')
AA = [minmaxSim.States.min{p}(:,inds{p}.AC);flipud(minmaxSim.States.max{p}(:,inds{p}.AC))];
fill(t,AA,colors{3},'FaceAlpha',0.8,'EdgeColor','none')
AV = [minmaxSim.States.min{p}(:,inds{p}.AV);flipud(minmaxSim.States.max{p}(:,inds{p}.AV))];
fill(t,AV,colors{2},'FaceAlpha',0.6,'EdgeColor','none')
PV = [minmaxSim.States.min{p}(:,inds{p}.PV);flipud(minmaxSim.States.max{p}(:,inds{p}.PV))];
fill(t,PV,colors{4},'FaceAlpha',0.6,'EdgeColor','none')

% Plot simulation with best fit to data
plot(simtime,minmaxSim.States.best{p}(:,inds{p}.MV),'-','color',0.9.*colors{1},'LineWidth',2)
plot(simtime,minmaxSim.States.best{p}(:,inds{p}.AV),'-','color',0.9.*colors{2},'LineWidth',2)
plot(simtime,minmaxSim.States.best{p}(:,inds{p}.AC),'-','color',0.9.*colors{3},'LineWidth',2)
plot(simtime,minmaxSim.States.best{p}(:,inds{p}.PV),'-','color',0.9.*colors{4},'LineWidth',2)

xlim([tstart round(tend,1)])
ymax=  max([minmaxSim.States.max{p}(:,inds{p}.AV);minmaxSim.States.max{p}(:,inds{p}.AC);data{p}.AV(:,2);data{p}.AC(:,2)]);
ymin = min([minmaxSim.States.min{p}(:,inds{p}.AC);minmaxSim.States.min{p}(:,inds{p}.AV);minmaxSim.States.min{p}(:,inds{p}.PV);data{p}.PV(:,2)]);
ylim([ymin,ymax])
if round(ymin,-1)<0
    yticks([round(ymin,-1),[0:100:round(ymax,-2)]])
else
    yticks([0:100:round(ymax,-2)])
end
xticks([tstart,tend/2,round(tend,1)]);
xticklabels({num2str(tstart),'',num2str(round(tend,1))})
xlabel('Time (s)')
ylabel('Flow (ml/s)')
set(gca,'FontSize',12)

%% Create plots of fit to parameter data for all subjects
% Collect the parameter data from all sujbects
dvalsDataparams = cell(size(dataparams));
pvalsDataparams = cell(size(dataparams));

for param = 1:length(dataparams)
    pvals = zeros(size(patNums));
    dvals = zeros(size(patNums));
    for pat = 1:length(patNums)
        theta = bestparams(:,pat);
        pvals(pat) = theta(ind.(dataparams{param}));
        dval = data{pat}.params.(dataparams{param});
        dvals(pat) = dval(1);
    end
    dvalsDataparams{param} = dvals;
    pvalsDataparams{param} = pvals;
end

DBP = zeros(size(patNums));
SBP = DBP;dataDBP=DBP;dataSBP=DBP;
for p = 1:length(patNums)
    DBP(p) = minmaxSim.DBP.best{p};
    SBP(p) = minmaxSim.SBP.best{p};
    dataDBP(p) = data{p}.DBP;
    dataSBP(p) = data{p}.SBP;
end

% Create bland altman plot of parameters fit to data 
figure('Visible', 'off','Name','Parameters fit to data compact bland altman uncertainty');
LW=1;
set(gcf,'Color','white')
xdim_CM = 21; 
ydim_CM = 15;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
dataparams = {'Emax_LV','Caa','Ctot','Rtot','ElCo'};
dataparamsTitles = {'Emax_L_V','Caa','Ctot','Rtot','ElCo'};
dataparamsUnits = {'mmHg/mL','mL/mmHg','mL/mmHg','mmHg·s/mL','cm^2'};
tiles=tiledlayout(2,3,'TileSpacing','loose','Padding','compact');
t1=nexttile;
TilePos = tiles.Children.InnerPosition;
annotation('textbox',[TilePos(1)-0.035 TilePos(2)+1.05*TilePos(4) .015 .015],'String',letters(1),'Linestyle','none','FitBoxToText','on','BackgroundColor','none');
hold on

%SBP
averageS = (dataSBP+SBP)./2;
differenceS = (dataSBP-SBP);
averageD = (dataDBP+DBP)./2;
differenceD = (dataDBP-DBP);
maxx = max(averageS)*1.1;
minx = min(averageD)*0.9;
sbpplot=plot(averageS,differenceS,'.','color',paramcol2,'LineWidth',1,'Markersize',10);
yline(mean(differenceS),'-','color',paramcol2,'LineWidth',LW*2);
high = mean(differenceS)+std(differenceS)*1.96;
low = mean(differenceS)-std(differenceS)*1.96;
fill([minx,minx,maxx,maxx],[low,high,high,low],paramcol2,'FaceAlpha',0.3,'LineStyle','None')

%DBP
dbpplot=plot(averageD,differenceD,'.','color',paramcol,'LineWidth',1,'Markersize',10);
yline(mean(differenceD),'-','color',paramcol,'LineWidth',LW*2);
high = mean(differenceD)+std(differenceD)*1.96;
low = mean(differenceD)-std(differenceD)*1.96;
fill([minx,minx,maxx,maxx],[low,high,high,low],paramcol,'FaceAlpha',0.15,'LineStyle','None')
for p = 1:length(patNums)
    diff1 = dataSBP(p)-minmaxSim.SBP.min{p};
    diff2 = dataSBP(p)-minmaxSim.SBP.max{p};
    diffmax = max([diff1,diff2,differenceS(p)]);
    diffmin = min([diff1,diff2,differenceS(p)]);
    es=errorbar(averageS(p),differenceS(p),differenceS(p)-diffmin,diffmax-differenceS(p),'.','color',paramcol2,'LineWidth',0.9);
    set([es.Bar, es.Line], 'ColorType', 'truecoloralpha', 'ColorData', [es.Line.ColorData(1:3); 255*0.3]);
    diff1 = dataDBP(p)-minmaxSim.DBP.min{p};
    diff2 = dataDBP(p)-minmaxSim.DBP.max{p};
    diffmax = max([diff1,diff2,differenceD(p)]);
    diffmin = min([diff1,diff2,differenceD(p)]);
    ed=errorbar(averageD(p),differenceD(p),differenceD(p)-diffmin,diffmax-differenceD(p),'.','color',paramcol,'LineWidth',0.9);
    set([ed.Bar, ed.Line], 'ColorType', 'truecoloralpha', 'ColorData', [ed.Line.ColorData(1:3); 255*0.3]);
end 
xlabel('Average (mmHg)','FontSize',fz1)
ylabel(sprintf('Brachial pressure (mmHg)\nData - Simulation'),'FontSize',fz1)
xlim([minx maxx])
set(gca,'FontSize',fz1)

% The data-based parameters
for param = 1:length(dataparams)
    t=nexttile;
    paramind = strcmp(paramNames,dataparams{param});
    createBlandAltmanPlot_uncertainty(t,tiles,dvalsDataparams{param},pvalsDataparams{param},paramcol,[0 0 0],LW,['Average (' dataparamsUnits{param} ')'],sprintf('%s (%s)\nData - Parameter value', dataparamsTitles{param},dataparamsUnits{param}),letters(param+1),paramuncertainty.allokParams,paramind)
    set(gca,'FontSize',fz1)
end

legend(t1,[ed,sbpplot,dbpplot],{'Uncertainty','SBP','DBP'},'FontSize',7,...
    'Position',[0.196000225159964 0.889744944420004 0.12 0.07]);

%% Plot the fit to blood flow data in ALL subjects (supplementary figure)
realPatNumsTitles = 0;
LWsim = 1.5;
LWdata = 1;
figure('Visible', 'off','Name','Flows - data vs simulation ALL uncertainty');
set(gcf,'Color','white')
xdim_CM = 20;
ydim_CM = (20/21)*29.7;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
c=6;
r = ceil(length(patNums)/c);
tiles = tiledlayout(r,c);
tiles.Padding = 'none';
tiles.TileSpacing = 'tight';
for pnum = 1:length(sortpatstotalcost)
    p = sortpatstotalcost(pnum);
     nexttile
     hold on
    simtime = minmaxSim.time{p};
    datatime = data{p}.MV(:,1);
    tend = datatime(end);
    indd = data{p}.indtdiast;
    
    maxflow = max([minmaxSim.States.max{p}(:,inds{p}.MV);...
        minmaxSim.States.max{p}(:,inds{p}.AV);...
        minmaxSim.States.max{p}(:,inds{p}.AC);...
        minmaxSim.States.max{p}(:,inds{p}.PV);...
        data{p}.MV(indd:end,2);data{p}.AV(1:indd,2);data{p}.AC(1:indd,2);data{p}.PV(:,2)]);
    
    % Plot simulation uncertainty
    t = [100*simtime./tend;flipud(100*simtime./tend)];
    MV = [minmaxSim.States.min{p}(:,inds{p}.MV);flipud(minmaxSim.States.max{p}(:,inds{p}.MV))];
    fill(t,100*MV./maxflow,colors{1},'FaceAlpha',0.6,'EdgeColor','none')
    AA = [minmaxSim.States.min{p}(:,inds{p}.AC);flipud(minmaxSim.States.max{p}(:,inds{p}.AC))];
    fill(t,100*AA./maxflow,colors{3},'FaceAlpha',0.8,'EdgeColor','none')
    AV = [minmaxSim.States.min{p}(:,inds{p}.AV);flipud(minmaxSim.States.max{p}(:,inds{p}.AV))];
    fill(t,100*AV./maxflow,colors{2},'FaceAlpha',0.6,'EdgeColor','none')
    PV = [minmaxSim.States.min{p}(:,inds{p}.PV);flipud(minmaxSim.States.max{p}(:,inds{p}.PV))];
    fill(t,100*PV./maxflow,colors{4},'FaceAlpha',0.6,'EdgeColor','none')
    
    % Plot simulation with best fit to data
    plot(100*simtime./tend,minmaxSim.States.best{p}(:,inds{p}.MV)./maxflow,'-','color',0.9.*colors{1},'LineWidth',LWsim)
    plot(100*simtime./tend,minmaxSim.States.best{p}(:,inds{p}.AV)./maxflow,'-','color',0.9.*colors{2},'LineWidth',LWsim)
    plot(100*simtime./tend,minmaxSim.States.best{p}(:,inds{p}.AC)./maxflow,'-','color',0.9.*colors{3},'LineWidth',LWsim)
    plot(100*simtime./tend,minmaxSim.States.best{p}(:,inds{p}.PV)./maxflow,'-','color',0.9.*colors{4},'LineWidth',LWsim)
    
    % Plot data
    d(1)=plot(100*datatime(indd:end,1)./tend,100*data{p}.MV(indd:end,2)./maxflow,'o','color',colors{1},'LineWidth',LWdata,'MarkerSize',2);
    d(2)=plot(100*datatime(1:indd)./tend,100*data{p}.AV(1:indd,2)./maxflow,'^','color',colors{2},'LineWidth',LWdata,'MarkerSize',2);
    d(3)=plot(100*datatime(1:indd)./tend,100*data{p}.AC(1:indd,2)./maxflow,'square','color',colors{3},'LineWidth',LWdata,'MarkerSize',2);
    d(4)=plot(100*datatime./tend,100*data{p}.PV(:,2)./maxflow,'diamond','color',colors{4},'LineWidth',LWdata,'MarkerSize',2);
    
    xlim([0 100])
    ylim([-30 100])
    set(gca,'FontSize',fzsmall)
    if pnum > c*(r-1)
        xlabel('% of cardiac cycle','FontSize',fzsmall)
        xticks([0 50 100])
    else
        xticks([])
    end
    if mod(pnum,c) == 1
        yticks([-30 0 50 100])
        ylabel('% of max flow','FontSize',fzsmall)
    else
        yticks([])
    end
    titlename = groupTitle(groups,p,patNums,realPatNumsTitles,num2str(pnum)); 
    title(titlename,'FontSize',fzsmall)
    hold off
end
tiles.Padding = 'none';
tiles.TileSpacing = 'tight';
hold on
sim = plot(NaN,NaN,'k-','Linewidth',1);
unc = fill(NaN,NaN,[0.5 0.5 0.5]);
legend([sim,unc,d],{'Simulation: best fit to data',...
    'Simulation uncertainty','Data MV','Data AV','Data AA','Data PV'},'FontSize',8,'Position',[0.415961204816109 0.00345895415058883 0.23015872530994 0.0846585570594995])

%% Save all figures
saveAllFigures(plotFolderName)

clear minmaxSim

end