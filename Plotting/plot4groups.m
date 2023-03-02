function plot4groups(comptitles,indSignComp,groups4,mean4,sd4,statisticsTable,allParamValues,paramNamesMarked,bounds,units,plotMedian,names,colors)
% plots box plots or errorbar plots of 4 groups.

% Font sizes
global fzlarge 
global fzsmall

% Find significantly different parameters
allParamValues = allParamValues';
rejects = cell(1,6);rejectsfind=rejects;
rejectsCorrected = cell(1,6);
for c = 1:6
    if ~isempty(comptitles{c})
        rejects{c} = statisticsTable.reject{strcmp(comptitles{c},statisticsTable.groupComparisonNames)};
        rejectsfind{c}=rejects{c};
        rejectsCorrected{c} = statisticsTable.rejectBenjH2{strcmp(comptitles{c},statisticsTable.groupComparisonNames)};
    else
        rejects{c} = nan(1,length(paramNamesMarked));
        rejectsfind{c} = zeros(1,length(paramNamesMarked));
        rejectsCorrected{c} = nan(1,length(paramNamesMarked));
    end
end
indSignThis = find(or(rejectsfind{1},or(rejectsfind{2},or(rejectsfind{3},or(rejectsfind{4},or(rejectsfind{5},rejectsfind{6}))))));
if isempty(indSignComp)
    differentInds =  [];
    indSign = indSignThis;
else
    indSign = unique([indSignComp,indSignThis],'stable');
    differentInds = indSignThis(~ismember(indSignThis,indSignComp));
end


% Plot
if plotMedian
    figure('Visible', 'off','Name',[names.figname '_median'])
else
    figure('Visible', 'off','Name',[names.figname '_mean'])
end
set(gcf,'Color','white')
if length(indSign) > 16
    tiles = tiledlayout(7,4,'TileSpacing','tight','Padding','compact');
    xdim_CM = 17;
    ydim_CM = 29;
    rowboundary = 6*4;
    legendposition=[0.0809834738580603 0.106126702239051 0.798908480268682 0.0620588235294119];
else
    tiles = tiledlayout(5,4,'TileSpacing','tight','Padding','compact');
    xdim_CM = 17.2;
    ydim_CM = 12.5+(12.5/4);
    rowboundary = 4*4;
    legendposition=[0.0436898230644095 0.143818554090903 0.798908480268682 0.0620588235294118];
end
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
letters = ['C':'Z' 'a':'z'];

% Create a plot for each significant parameter
for param = 1:length(indSign)
    t=nexttile;
    hold on
    % Plot values for all of the subjects
    % control normotensive
    swarmchart(0.1.*ones(size(allParamValues(groups4{1},indSign(param)))),allParamValues(groups4{1},indSign(param)),2,colors{1},'*','XJitterWidth',0.08)
    % t2d normotensive
    swarmchart(0.2.*ones(size(allParamValues(groups4{2},indSign(param)))),allParamValues(groups4{2},indSign(param)),2,colors{2},'*','XJitterWidth',0.08)
    % control hypertensive
    swarmchart(0.3.*ones(size(allParamValues(groups4{3},indSign(param)))),allParamValues(groups4{3},indSign(param)),2,colors{3},'*','XJitterWidth',0.08)
    % t2d hypertensive
    swarmchart(0.4.*ones(size(allParamValues(groups4{4},indSign(param)))),allParamValues(groups4{4},indSign(param)),2,colors{4},'*','XJitterWidth',0.08)
    
    if plotMedian
         % box plot of median values
        cnormal = boxchart(0.1.*ones(size(allParamValues(groups4{1},indSign(param)))),allParamValues(groups4{1},indSign(param)),'BoxWidth',0.08,'BoxFaceColor',colors{1},'MarkerStyle','+','MarkerColor',[0 0 0],'Notch','off');
        t2dnormal = boxchart(0.2.*ones(size(allParamValues(groups4{2},indSign(param)))),allParamValues(groups4{2},indSign(param)),'BoxWidth',0.08,'BoxFaceColor',colors{2},'MarkerStyle','+','MarkerColor',[0 0 0],'Notch','off');
        chigh = boxchart(0.3.*ones(size(allParamValues(groups4{3},indSign(param)))),allParamValues(groups4{3},indSign(param)),'BoxWidth',0.08,'BoxFaceColor',colors{3},'MarkerStyle','+','MarkerColor',[0 0 0],'Notch','off');
        t2dhigh = boxchart(0.4.*ones(size(allParamValues(groups4{4},indSign(param)))),allParamValues(groups4{4},indSign(param)),'BoxWidth',0.08,'BoxFaceColor',colors{4},'MarkerStyle','+','MarkerColor',[0 0 0],'Notch','off');
    else
        %error bars mean+-sd
        cnormal = errorbar(0.1,mean4{1}(indSign(param)),sd4{1}(indSign(param)),'*','color',colors{1},'LineWidth',2);
        t2dnormal=errorbar(0.2,mean4{2}(indSign(param)),sd4{2}(indSign(param)),'*','color',colors{2},'LineWidth',2);
        chigh=errorbar(0.3,mean4{3}(indSign(param)),sd4{3}(indSign(param)),'*','color',colors{3},'LineWidth',2);
        t2dhigh=errorbar(0.4,mean4{4}(indSign(param)),sd4{4}(indSign(param)),'*','color',colors{4},'LineWidth',2);
    end
    xlim([0 0.5])
    
    %find minimim and maximum y values
    if plotMedian
        ymax = max([max(allParamValues(groups4{1},indSign(param))),max(allParamValues(groups4{2},indSign(param))),max(allParamValues(groups4{3},indSign(param))),max(allParamValues(groups4{4},indSign(param)))]);
        ymin = min([min(allParamValues(groups4{1},indSign(param))),min(allParamValues(groups4{2},indSign(param))),min(allParamValues(groups4{3},indSign(param))),min(allParamValues(groups4{4},indSign(param)))]);
    else
        max1 = 1.05*max(mean4{1}(indSign(param))+sd4{1}(indSign(param)),mean4{3}(indSign(param))+sd4{3}(indSign(param)));
        min1 = min(mean4{1}(indSign(param))-sd4{1}(indSign(param)),mean4{3}(indSign(param))-sd4{3}(indSign(param)));
        max2 = 1.05*max(mean4{2}(indSign(param))+sd4{2}(indSign(param)),mean4{4}(indSign(param))+sd4{4}(indSign(param)));
        min2 = min(mean4{2}(indSign(param))-sd4{2}(indSign(param)),mean4{4}(indSign(param))-sd4{4}(indSign(param)));
        ymax = max([allParamValues(:,indSign(param));max1;max2]);
        ymin = min([allParamValues(:,indSign(param));min1;min2]);
    end
    
    % calculate which y bounds to use and where to plot the p-values
    if bounds.ub(indSign(param)) < ymax*1.2 %bounds close enough to values
        plotboundsU = 1;
        ymax = max([ymax*1.1,bounds.ub(indSign(param))*1.04]);
    else
        ymax = ymax*1.1;
        plotboundsU = 0;
    end
    
    if bounds.lb(indSign(param)) > ymin*0.9 %bounds close enough to values
        plotboundsL = 1;
        ymin = min([ymin,bounds.lb(indSign(param))*0.95]);
    else
        ymin = ymin*0.9;
        plotboundsL = 0;
    end
    
    y1 = 1;
    ytext = 0.12;
    y2 = (y1+ytext+0.15);
    y3 = (y2+ytext+ 0.15);
    y4 = (y3+ytext+ 0.15);
    pvalsize = 8;
    
    ylim([ymin,ymax*(y4+ytext+0.01)])
    

    linecoordinates = {[0.1 0.28],[0.22 0.4],[0.32 0.4],[0.1 0.18],[0.22 0.3],[0.1 0.4]};
    yvals = [y3,y2,y3,y2,y1,y4];
    xvals = [0.06,0.24,0.32,0.06,0.24,0.06];
    
    % get p-values and print them
    pvals = zeros(1,6);
    for i = 1:6
        if ~isempty(comptitles{i})
            pvals(i) = statisticsTable.pvalues{strcmp(comptitles{i},statisticsTable.groupComparisonNames)}(indSign(param));
            print_pval(t,xvals(i),yvals(i),linecoordinates{i},pvals(i),indSign(param),rejects{i},rejectsCorrected{i},pvalsize,ymax,ytext)
        else
            pvals(i) = NaN;
        end
    end
    
    % plot parameter bounds
    if plotboundsU
        bound=yline(bounds.ub(indSign(param)),'--');
    else
        bound=plot(NaN,NaN,'k--');
    end
    if plotboundsL
        yline(bounds.lb(indSign(param)),'--');
    end
    
    sign = plot(NaN,NaN,'w*','LineWidth',1);
    sign2 = plot(NaN,NaN,'w*','LineWidth',1);
    indvalue = plot(NaN,NaN,'k.','LineWidth',2);
    outlier = plot(NaN,NaN,'k+','LineWidth',1);
    xticklabels({});
    xticks([0.1 0.2 0.3 0.41])
    lastrow = (floor(length(indSign)/4))*4;
    if mod(length(indSign),4) == 0
        lastrow = (floor(length(indSign)/4)-1)*4;
    end
    if param > lastrow
        xticklabels(names.xtick)
        xtickangle(names.xtickangle)
    else
        xticklabels({});
    end
    ylabel(units.param(indSign(param)),'FontName','Calibri')
    ax = gca;
    ax.FontSize = fzsmall;
    ax.FontName = 'Calibri';
    if ismember(indSign(param),differentInds)
        title(paramNamesMarked(indSign(param)),'FontSize',fzlarge,'color',[0 0 0.7],'FontName','Calibri')
    elseif  sum(pvals(~isnan(pvals))<0.05) == 0
        title(paramNamesMarked(indSign(param)),'FontSize',fzlarge,'color',[1 0 0],'FontName','Calibri')
    else
        title(paramNamesMarked(indSign(param)),'FontSize',fzlarge,'FontName','Calibri')
    end
    TilePos = tiles.Children.InnerPosition;
    letter = annotation('textbox',[TilePos(1)-0.035 TilePos(2)+1.15*TilePos(4) .015 .015],'String',letters(param),'Linestyle','none','FitBoxToText','on','BackgroundColor','none','FontSize',fzlarge*1.1,'FontName','Calibri');
end

% Plot legend
if plotMedian && (floor(length(indSign)/4))*4 <= rowboundary
    legend([cnormal,t2dnormal,chigh,t2dhigh,bound,indvalue,outlier,sign],{...
        sprintf('%s (N=%d)',names.legend{1},sum(groups4{1})),...
        sprintf('%s (N=%d)',names.legend{2},sum(groups4{2})),...
        sprintf('%s (N=%d)',names.legend{3},sum(groups4{3})),...
        sprintf('%s (N=%d)',names.legend{4},sum(groups4{4})),...
        'Litterature-based parameter bounds','Individual subject value','Outlier','*Significant with correction'},...
        'FontSize',fzsmall,'Numcolumns',3,'Position',legendposition,'FontName','Calibri')
elseif (floor(length(indSign)/4))*4 <= rowboundary
    legend([cnormal,t2dnormal,chigh,t2dhigh,bound,indvalue],{...
        sprintf('%s (N=%d)',names.legend{1},sum(groups4{1})),...
        sprintf('%s (N=%d)',names.legend{2},sum(groups4{2})),...
        sprintf('%s (N=%d)',names.legend{3},sum(groups4{3})),...
        sprintf('%s (N=%d)',names.legend{4},sum(groups4{4})),...
        'Litterature-based parameter bounds','Individual subject value'},...
        'FontSize',fzsmall,'Numcolumns',3,'Position',legendposition,'FontName','Calibri')
end

end