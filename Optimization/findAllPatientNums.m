function patNums = findAllPatientNums(dataFolder,range)
if nargin <1
    dataFolder = 'Data';
end

folders = dir([dataFolder '/*_P*']);
names = {folders.name};
patNums = cell(1,length(names));
for n = 1:length(names)
    patNums{n} = names{n}(end-5:end-4);
    if strcmp(patNums{n}(1),'P')
        patNums{n} = patNums{n}(2);
    end
end

if nargin > 1
    patNums = patNums(range);
end

end