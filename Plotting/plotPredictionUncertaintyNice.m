function plotPredictionUncertaintyNice(simulations,simulationNames,xnames,ynames,figureName,colors)
orange = [255 123 83]./255;
blue = [0 176 240]./255;
orange2 = [246 202 184]./255;
blue2 = [173 210 231]./255;
liuLilac = [137 129 211]./255;
liuOrange = [255 100 65]./255;
liuBlue = [0 185 231]./255;
if nargin <6
    colors = {[0 0 0],liuOrange,liuBlue,liuLilac,orange2,blue2};
end
fz1 = 8.5;

r = 3;c=2;
figure('Name',figureName,'Visible','off')
set(gcf,'Color','white')
set(gcf, 'InvertHardCopy', 'off'); % setting 'grid color reset' off
xdim_CM = 17;
ydim_CM = 12.75;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
tiledlayout(r,c,'TileSpacing','loose','Padding','compact');

nexttile
yInd = strcmp('mvCorr',ynames);
hold on
for s = 1:length(simulations)
    for p = 1:length(simulations{s}.time)
        tend = simulations{s}.time{p}(end);
        t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
        yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
        fill(t,yvals,colors{s},'FaceAlpha',0.4,'EdgeColor','none');
        plot(100.*simulations{s}.time{p}./tend,simulations{s}.Observables.best{p}(:,yInd),'-','color',colors{s},'LineWidth',2.2)
    end
end
xlabel('% of cardiac cycle','FontSize',fz1)
ylabel('Flow in mitral valve (ml/s)','FontSize',fz1)
set(gca,'FontSize',fz1)

nexttile
hold on
yInd = strcmp('avCorr',ynames);
for s = 1:length(simulations)
    for p = 1:length(simulations{s}.time)
        tend = simulations{s}.time{p}(end);
        t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
        yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
        sims(s)=fill(t,yvals,colors{s},'FaceAlpha',0.4,'EdgeColor','none');
        plot(100.*simulations{s}.time{p}./tend,simulations{s}.Observables.best{p}(:,yInd),'-','color',colors{s},'LineWidth',2.2)
    end
end
xlabel('% of cardiac cycle','FontSize',fz1)
ylabel('Flow in aortic valve (ml/s)','FontSize',fz1)
set(gca,'FontSize',fz1)

t1=nexttile;
colororder(t1,{'k','k'})%y axis colors
yInd = strcmp('Ela',ynames);
hold on
yyaxis right
ymax = 0;
for s = 1:length(simulations)
    for p = 1:length(simulations{s}.time)
        tend = simulations{s}.time{p}(end);
        t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
        yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
        fill(t,yvals,colors{s},'FaceAlpha',0.4,'EdgeColor','none')
        plot(100.*simulations{s}.time{p}./tend,simulations{s}.Observables.best{p}(:,yInd),'-','color',colors{s},'LineWidth',2.2)
        ymax = max([ymax,max(simulations{s}.Observables.max{p}(:,yInd))]);
    end
end
ylim([-ymax ymax])

yyaxis left
yInd = strcmp('Elv',ynames);
ymax=0;
for s = 1:length(simulations)
    for p = 1:length(simulations{s}.time)
        tend = simulations{s}.time{p}(end);
        t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
        yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
        fill(t,yvals,colors{s},'FaceAlpha',0.4,'EdgeColor','none')
        plot(100.*simulations{s}.time{p}./tend,simulations{s}.Observables.best{p}(:,yInd),'-','color',colors{s},'LineWidth',2.2)
        ymax = max([ymax,max(simulations{s}.Observables.max{p}(:,yInd))]);
    end
end
ylim([0 ymax*2])
xlabel('% of cardiac cycle','FontSize',fz1)
yyaxis left
ylabel(sprintf('Time-varying elastance\nin LV (mmHg/ml)'),'FontSize',fz1)
yyaxis right
ylabel(sprintf('Time-varying elastance\nin LA (mmHg/ml)'),'FontSize',fz1)
set(gca,'FontSize',fz1)
xlim([0 100])

yparamsToPlot = {'P_Aortic','pLA','pLV'};
ylabels = {'Aortic pressure (mmHg)','Pressure in LA (mmHg)','Pressure in LV (mmHg)'};
for y = 1:length(yparamsToPlot)
    nexttile
    yInd = strcmp(yparamsToPlot{y},ynames);
    hold on
    for s = 1:length(simulations)
        for p = 1:length(simulations{s}.time)
            tend = simulations{s}.time{p}(end);
            t = 100.*[simulations{s}.time{p}./tend;flipud(simulations{s}.time{p}./tend)];
            yvals = [simulations{s}.Observables.min{p}(:,yInd);flipud(simulations{s}.Observables.max{p}(:,yInd))];
            fill(t,yvals,colors{s},'FaceAlpha',0.4,'EdgeColor','none')
            plot(100.*simulations{s}.time{p}./tend,simulations{s}.Observables.best{p}(:,yInd),'-','color',colors{s},'LineWidth',2.2)
        end
    end
    xlabel('% of cardiac cycle','FontSize',fz1)
    ylabel(ylabels{y},'FontSize',fz1)
    set(gca,'FontSize',fz1)
    xlim([0 100])
end
l=legend(sims,simulationNames,'Location','NorthEast','FontSize',8);
set(l,'position',[0.765074192528902 0.802419119653766 0.189153435526702 0.154320983396847]);

end