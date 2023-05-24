function plot4groups_variables(indSignComp,comptitles,groups4,mean4,sd4,statisticsTable,allValues,varnamesMarked,units,plotMedian,names,colors)
% Define font sizes
fzlarge = 10;
fzsmall = 8;

% Find the significantly different variables
rejects = cell(1,6);
rejectsCorrected = cell(1,6);
for c = 1:6
    if ~isempty(comptitles{c})
        rejects{c} = statisticsTable.reject{strcmp(comptitles{c},statisticsTable.groupComparisonNames)};
        try
            rejectsCorrected{c} = statisticsTable.rejectBenjH2{strcmp(comptitles{c},statisticsTable.groupComparisonNames)};
        catch
            rejectsCorrected{c} = nan(1,length(varnamesMarked));
        end
    else
        rejects{c} = nan(1,length(varnamesMarked));
        rejectsCorrected{c} = nan(1,length(varnamesMarked));
    end
end
rejectsfind = rejects;
for c = 1:6
    rejectsfind{c}(isnan(rejects{c})) = 0;
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
    tiledlayout(7,4,'TileSpacing','tight','Padding','compact');
    xdim_CM = 21;
    ydim_CM = 29;
    rowboundary = 6*4;
    legendposition=[0.0809834738580603 0.106126702239051 0.798908480268682 0.0620588235294119];
elseif length(indSign) > 4
    tiles = tiledlayout(5,4,'TileSpacing','tight','Padding','compact');
    xdim_CM = 21;
    ydim_CM = 15+(15/4);
    rowboundary = 4*4;
    legendposition=[0.0436898230644095 0.143818554090903 0.798908480268682 0.0620588235294118];
else
    tiles = tiledlayout(2,2,'TileSpacing','tight','Padding','compact');
    xdim_CM = 9.5;
    ydim_CM = 7;
    rowboundary = 1*2;
    legendposition=[0.0436898230644095 0.143818554090903 0.798908480268682 0.0620588235294118];
end
set(gcf,'Units','centimeters','Position',[0 0 xdim_CM ydim_CM])
set(gcf,'PaperUnits', 'centimeters', 'PaperSize', [xdim_CM, ydim_CM])
letters = ['A':'Z' 'a':'z'];
% create a subplot for each variable
for var = 1:length(indSign)
    t=nexttile;
    hold on
    % Plot variable values for all subjects
    % control normotensive
    swarmchart(0.1.*ones(size(allValues(groups4{1},indSign(var)))),allValues(groups4{1},indSign(var)),2,colors{1},'*','XJitterWidth',0.08)
    % t2d normotensive
    swarmchart(0.2.*ones(size(allValues(groups4{2},indSign(var)))),allValues(groups4{2},indSign(var)),2,colors{2},'*','XJitterWidth',0.08)
    % control hypertensive
    swarmchart(0.3.*ones(size(allValues(groups4{3},indSign(var)))),allValues(groups4{3},indSign(var)),2,colors{3},'*','XJitterWidth',0.08)
    % t2d hypertensive
    swarmchart(0.4.*ones(size(allValues(groups4{4},indSign(var)))),allValues(groups4{4},indSign(var)),2,colors{4},'*','XJitterWidth',0.08)
    
    if plotMedian
        % box plot of median values
        cnormal = boxchart(0.1.*ones(size(allValues(groups4{1},indSign(var)))),allValues(groups4{1},indSign(var)),'BoxWidth',0.08,'BoxFaceColor',colors{1},'MarkerStyle','+','MarkerColor',[0 0 0],'Notch','off');
        t2dnormal = boxchart(0.2.*ones(size(allValues(groups4{2},indSign(var)))),allValues(groups4{2},indSign(var)),'BoxWidth',0.08,'BoxFaceColor',colors{2},'MarkerStyle','+','MarkerColor',[0 0 0],'Notch','off');
        chigh = boxchart(0.3.*ones(size(allValues(groups4{3},indSign(var)))),allValues(groups4{3},indSign(var)),'BoxWidth',0.08,'BoxFaceColor',colors{3},'MarkerStyle','+','MarkerColor',[0 0 0],'Notch','off');
        t2dhigh = boxchart(0.4.*ones(size(allValues(groups4{4},indSign(var)))),allValues(groups4{4},indSign(var)),'BoxWidth',0.08,'BoxFaceColor',colors{4},'MarkerStyle','+','MarkerColor',[0 0 0],'Notch','off');
    else
        %error bars mean+-sd
        cnormal = errorbar(0.1,mean4{1}(indSign(var)),sd4{1}(indSign(var)),'*','color',colors{1},'LineWidth',2);
        t2dnormal=errorbar(0.2,mean4{2}(indSign(var)),sd4{2}(indSign(var)),'*','color',colors{2},'LineWidth',2);
        chigh=errorbar(0.3,mean4{3}(indSign(var)),sd4{3}(indSign(var)),'*','color',colors{3},'LineWidth',2);
        t2dhigh=errorbar(0.4,mean4{4}(indSign(var)),sd4{4}(indSign(var)),'*','color',colors{4},'LineWidth',2);
    end
    xlim([0 0.5])
    
    %find minimim and maximum y values
    max1 = 1.05*max(mean4{1}(indSign(var))+sd4{1}(indSign(var)),mean4{3}(indSign(var))+sd4{3}(indSign(var)));
    min1 = min(mean4{1}(indSign(var))-sd4{1}(indSign(var)),mean4{3}(indSign(var))-sd4{3}(indSign(var)));
    max2 = 1.05*max(mean4{2}(indSign(var))+sd4{2}(indSign(var)),mean4{4}(indSign(var))+sd4{4}(indSign(var)));
    min2 = min(mean4{2}(indSign(var))-sd4{2}(indSign(var)),mean4{4}(indSign(var))-sd4{4}(indSign(var)));
    ymax = max([allValues(:,indSign(var));max1;max2]);
    ymin = min([allValues(:,indSign(var));min1;min2]);
    ymin = ymin*0.9;
    ymax = ymax*1.05;
    
    
    % calculate which y bounds to use and where to plot the p-values
    y1 = 1;
    ytext = 0.09;
    y2 = (y1+ytext+0.1);
    y3 = (y2+ytext+ 0.1);
    y4 = (y3+ytext+ 0.1);
    pvalsize = 8;
    
    ylim([ymin,ymax*(y4+ytext+0.01)])
    
    pvals = zeros(1,6);

    linecoordinates = {[0.1 0.28],[0.22 0.4],[0.32 0.4],[0.1 0.18],[0.22 0.3],[0.1 0.4]};
    yvals = [y3,y2,y3,y2,y1,y4];
    xvals = [0.06,0.24,0.32,0.06,0.24,0.06];
    
     % get p-values and print them
    for i = 1:6
        if ~isempty(comptitles{i})
            pvals(i) = statisticsTable.pvalues{strcmp(comptitles{i},statisticsTable.groupComparisonNames)}(indSign(var));
            print_pval(t,xvals(i),yvals(i),linecoordinates{i},pvals(i),indSign(var),rejects{i},rejectsCorrected{i},pvalsize,ymax,ytext)
        else
            pvals(i) = NaN;
        end
    end
    
    sign = plot(NaN,NaN,'w*','LineWidth',1);
    indvalue = plot(NaN,NaN,'k.','LineWidth',2);
    outlier = plot(NaN,NaN,'k+','LineWidth',1);
    xticklabels({});
    xticks([0.1 0.2 0.3 0.41])
    lastrow = (floor(length(indSign)/4))*4;%4*4
    if mod(length(indSign),4) == 0
        lastrow = (floor(length(indSign)/4)-1)*4;
    end
    if var > lastrow
        xticklabels(names.xtick)
        xtickangle(names.xtickangle)
    else
        xticklabels({});
    end
    ylabel(sprintf('%s\n%s',varnamesMarked{indSign(var)},units{indSign(var)}))
    
    ax = gca;
    ax.FontSize = fzsmall;
    ax.FontName = 'Calibri';
end

% Plot legend
if plotMedian && (floor(length(indSign)/4))*4 <= rowboundary
        legend([cnormal,t2dnormal,chigh,t2dhigh,indvalue,outlier,sign],{...
        sprintf('%s (N=%d)',names.legend{1},sum(groups4{1})),...
        sprintf('%s (N=%d)',names.legend{2},sum(groups4{2})),...
        sprintf('%s (N=%d)',names.legend{3},sum(groups4{3})),...
        sprintf('%s (N=%d)',names.legend{4},sum(groups4{4})),...
        'Individual subject value','Outlier','*Significant with correction'},...
        'FontSize',fzsmall,'Position',[0.5 0.07 0.3 0.397])
elseif (floor(length(indSign)/4))*4 <= rowboundary
    legend([cnormal,t2dnormal,chigh,t2dhigh,indvalue],{...
        sprintf('%s (N=%d)',names.legend{1},sum(groups4{1})),...
        sprintf('%s (N=%d)',names.legend{2},sum(groups4{2})),...
        sprintf('%s (N=%d)',names.legend{3},sum(groups4{3})),...
        sprintf('%s (N=%d)',names.legend{4},sum(groups4{4})),...
        'Individual subject value'},...
        'FontSize',fzsmall,'Numcolumns',3,'Position',legendposition)
end

end