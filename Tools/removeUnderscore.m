function fixedStrings =  removeUnderscore(stringList,replacement)
%takes a cell array with strings and removes any underscores in the strings
%outputs the cell array fixedStrings with the same strings without the
%underscores.

fixedStrings = stringList;
for s = 1:length(stringList)
    if ismember('_',stringList{s})
        ind = strfind(stringList{s},'_');
        if nargin < 2
            fixedStrings{s}(ind) = [];
        else
            fixedStrings{s}(ind) = replacement;
        end
    end
end

end


