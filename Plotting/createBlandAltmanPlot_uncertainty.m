function createBlandAltmanPlot_uncertainty(t,tiles,data,simulation,dotcolor,linecolor,linewidth,xlabelname,ylabelname,lettertowrite,uncertainty,paramind)

isvariable = isstruct(uncertainty);
%uncertainty is saved in allokparams if parameter, or in minmaxsims.variable if
%variable

average = (data+simulation)./2;
difference = (data-simulation);
average = average(~isnan(average));
difference = difference(~isnan(difference));

maxx = max(average)*1.1;
minx = min(average)*0.9;
high = mean(difference)+std(difference)*1.96;
low = mean(difference)-std(difference)*1.96;

plot(t,average,difference,'.','color',dotcolor,'LineWidth',1,'Markersize',10);
hold on
%uncertainty
for p = 1:length(data)
    if isvariable
        minpdiff = data(p)-uncertainty.min{p};
        maxpdiff = data(p)-uncertainty.max{p};
    else
        minpdiff = data(p)-min(uncertainty{p}(:,paramind));
        maxpdiff = data(p)-max(uncertainty{p}(:,paramind));
    end
    diffmax = max([minpdiff,maxpdiff,difference(p)]);
    diffmin = min([minpdiff,maxpdiff,difference(p)]);
    e=errorbar(average(p),difference(p),abs(difference(p)-diffmin),abs(diffmax-difference(p)),'color',dotcolor,'LineWidth',0.9);
    set([e.Bar, e.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e.Line.ColorData(1:3); 255*0.3]);
end
yline(t,mean(difference),'-','color',linecolor,'LineWidth',linewidth);
yline(t,low,'--','color',linecolor,'LineWidth',linewidth);
yline(t,high,'--','color',linecolor,'LineWidth',linewidth);

ylabel(t,ylabelname)
xlabel(t,xlabelname)
xlim(t,[minx maxx])

% Print letter
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
annotation('textbox',[TilePos(1)-0.035 h .015 .015],'String',lettertowrite,'Linestyle','none','FitBoxToText','on','BackgroundColor','none');

end