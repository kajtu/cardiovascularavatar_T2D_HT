function [pValues,r2Values,rhos,numObservations,testtype] = createSummaryCorrelationPlot(titlename,corrvarnames,corrvars,xorig)
paramcolor = [21,96,122]./255;%blue
paramcolor2 = [255 181 95]./255;%yellow
colororder([paramcolor;paramcolor2.*0.8])
pValues = zeros(size(corrvars));
r2Values = zeros(size(corrvars));
rhos = zeros(size(corrvars));
numObservations = zeros(size(corrvars));
testtype = cell(size(corrvars));

title(titlename)
hold on
for corrvar = 1:length(corrvars)
    x = xorig;
    y = corrvars{corrvar};
    
    %remove NaN values
    CheckNaN_X = isnan(x);
    CheckNaN_Y = isnan(y);
    y(CheckNaN_X)= NaN;
    x(CheckNaN_Y)= NaN;
    x(isnan(x)) = [];
    y(isnan(y)) = [];
    
    % linear correlation to get r2
    mdl = fitlm(x,y);
    numObservations(corrvar) = length(x);
    r2 = mdl.Rsquared.Ordinary;
    
    %check normal distribution
    try
        [isnormalx, pvalNormality, ~] = swtest(x, 0.005);
        [isnormaly, pvalNormality, ~] = swtest(y, 0.005);
    catch
        disp('sw normality test failed')
        isnormalx = 0;isnormaly = 0;
    end
    
    % Calculate p-value for the correlation
    if isnormalx && isnormaly
        %pearson correlation for normally distributed parameters
        [rhos(corrvar),p] = corr(x',y','type','Pearson');
        testtype{corrvar} = 'Pearson';
    else
        %spearman correlation for unparmaetric variables
        [rhos(corrvar),p] = corr(x',y','type','Spearman');
        testtype{corrvar} = 'Spearman';
    end

    yyaxis left
    bar(corrvar-0.01,1-p,'FaceColor',paramcolor,'EdgeColor','none','FaceAlpha',0.8)
    
    yyaxis right
    bar(corrvar+0.01,r2,'FaceColor',paramcolor2,'EdgeColor','none','FaceAlpha',0.5)
    plot(corrvar,r2,'.','color',paramcolor2,'Markersize',10)
    
    yyaxis left
    plot(corrvar,1-p,'.','color',paramcolor,'Markersize',10)
    
    pValues(corrvar) = p;
    r2Values(corrvar) = r2;
end
xticks(1:length(corrvars))
xticklabels(corrvarnames)
xtickangle(-20)
yyaxis left
ylabel('p-value')
ylim([0.8 1])
% ylim([0 0.2])

yyaxis right
ylim([0 1])
ylabel('r^2')

ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

end