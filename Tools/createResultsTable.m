function nicetable = createResultsTable(origintable, groups, comparisons,savefile,filename)

grouprows = ismember(origintable.Group,groups);
comprows = ismember(origintable.Group,comparisons);
nicetable = rows2vars(origintable(grouprows | comprows,:));

if savefile
    writetable(nicetable,[filename,'.xlsx'])
end

end


