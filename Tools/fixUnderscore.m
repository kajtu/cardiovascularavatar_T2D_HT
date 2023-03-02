function fixedStrings =  fixUnderscore(stringList)
%takes a cell array with strings and finds any underscores in the strings
%outputs the cell array fixedStrings with the same strings with the extra
%underscores. Useful when creating texts to plot.

fixedStrings = stringList;
for s = 1:length(stringList)
    if ismember('_',stringList{s})
        ind = strfind(stringList{s},'_');
        part1 = stringList{s}(1:ind(1)-1);
        part2 = stringList{s}(ind(1):end);
        part2 = join(split(part2(2:end),''),'_');
        part2 = part2{1};
        word = [part1 part2(1:end-1)];
        doubleind = strfind(word,'___');
        if ~isempty(doubleind)
            word(doubleind) = [];
            word(doubleind) = ' ';
        end
        fixedStrings{s} = word;
    end
end

end


