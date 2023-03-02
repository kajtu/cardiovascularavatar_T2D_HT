function saveAllFigures(plotFolderName,saveFigs)

if nargin < 2 || (nargin == 2 && saveFigs)
    disp('Saving figures...')    
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for fig = 1:length(FigList)
        currFig = FigList(fig);
        figName = currFig.Name;
        exportgraphics(currFig,[fullfile(plotFolderName,figName) '.pdf'],'ContentType','vector')
        close(currFig)
        %fprintf('%0.1f \n',100*fig/length(FigList))
    end
else
    disp('Not saving any figures.')
end

end