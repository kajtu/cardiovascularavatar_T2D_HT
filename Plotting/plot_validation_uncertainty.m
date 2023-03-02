function []=plot_validation_uncertainty(minmaxSim,paramuncertainty,bestparams,constants,simLast,indexes,options,data,extradata,patNums,xnames,plotFolderName)
%% Setup colors
colornames = {'MV','AV','AA','pulmonary','LA','LV','peripheral'};
namesavatar = {'pulmonary','LA','MV','LV','AV','AA','peripheral'};
bluescale = {[208,209,230]./255,[200,188,220]./255,[166,189,219]./255,[116,169,207]./255,[54,144,192]./255,[5,112,176]./255,[3,78,123]./255};
[~,colind] = ismember(colornames,namesavatar);
bluescaleordered = bluescale(colind);
colors = bluescaleordered;
colorLV = colors{strcmp('LV',colornames)};


%% Extract the simulated and measured stroke volumes for each subject
SVmvDat = zeros(size(patNums));
SVacDat = zeros(size(patNums));
SVavDat = zeros(size(patNums));
SVmriDat = zeros(size(patNums));
SVsim = zeros(size(patNums));

for p = 1:length(patNums)
    SVmriDat(p)=data{p}.EDV-data{p}.ESV;
    SVmvDat(p)=extradata{p}.SV_MV;
    SVavDat(p)=extradata{p}.SV_AV;
    SVacDat(p)=extradata{p}.SV_AC;
    SVsim(p)=trapz(simLast{p}.t,simLast{p}.x(:,strcmp(xnames,'Qav')));
end

%% Create valiation plot with prediction vs data
figure('Visible', 'off','Name','validation_predictionVSdata_uncertainty');
set(gcf,'Color','white')
xdim_CM = 21;
ydim_CM = 20;
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
LW = 1.5;
tiles=tiledlayout(4,2,'Padding','compact','TileSpacing','Compact');
t=nexttile;
createValidationPlot_uncertainty(t,tiles,SVmvDat,SVsim,colors{1},[0.3 0.3 0.3],LW,'Volume (ml)','A',minmaxSim.SV)
title('SV mitral valve (blood flow data)')

t=nexttile;
createValidationPlot_uncertainty(t,tiles,SVavDat,SVsim,colors{2},[0.3 0.3 0.3],LW,'Volume (ml)','B',minmaxSim.SV)
title('SV aortic valve (blood flow data)')

t=nexttile;
createValidationPlot_uncertainty(t,tiles,SVacDat,SVsim,colors{3},[0.3 0.3 0.3],LW,'Volume (ml)','C',minmaxSim.SV)
title('SV ascending aorta (blood flow data)')

t=nexttile;
createValidationPlot_uncertainty(t,tiles,SVmriDat,SVsim,colorLV,[0.3 0.3 0.3],LW,'Volume (ml)','D',minmaxSim.SV)
title('Validation: SV left ventricle (3D cine MRI data)')

legend([t.Children(1),t.Children(2)],'Data','Simulation uncertainty','Position',[0.554274030341607 0.673662299047534 0.198336696744927 0.0419312177890192])

%% Save figures
saveAllFigures(plotFolderName)

end


%%%%%%%%%%%%%%%%%%%%%%%%
function createValidationPlot_uncertainty(t,tiles,S1,S2,dotcolor,linecolor,linewidth,ylabelname,lettertowrite,uncertainty)

average = (S1+S2)./2;
nonans = ~isnan(average) & ~isnan(S1) & ~isnan(S2);
average = average(nonans);
S1 = S1(nonans);
S2 = S2(nonans);

[~,patorder] = sort(average);
S2sorted = S2(patorder);
unc_min_sorted = uncertainty.min(patorder);
unc_max_sorted = uncertainty.max(patorder);

hold on
% Plot simulation
for p = 1:length(S2)
    minp = unc_min_sorted{p};
    maxp = unc_max_sorted{p};
    errorbar(t,p,S2sorted(p),S2sorted(p)-minp,maxp-S2sorted(p),'.','color',dotcolor,'LineWidth',1,'Markersize',linewidth*10)
end
% Plot data
plot(t,1:length(S1),S1(patorder),'o','color',linecolor,'LineWidth',1,'Markersize',linewidth*4);

ylabel(t,ylabelname)
xlim(t,[0 length(patorder)+1])
xticks([1 round(length(patorder)/2) length(patorder)])
xticklabels({'1','Subject number',num2str(length(average))})

TilePos = tiles.Children.InnerPosition;
chnum = size(tiles.Children);
if chnum(1)>1
    h=TilePos(2)+1.15*TilePos(4);
else
    h=TilePos(2)+1.05*TilePos(4);
end
if h>1
    h=0.98;
end
letter = annotation('textbox',[TilePos(1)-0.035 h .015 .015],'String',lettertowrite,'Linestyle','none','FitBoxToText','on','BackgroundColor','none');


end
