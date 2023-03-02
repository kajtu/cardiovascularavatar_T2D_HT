function patientNums = findPatientRange(patientNums,patRange,dataFolder)

if isempty(patientNums)
    disp('patientNums is empty')
    if isempty(patRange)
        patientNums = findAllPatientNums(dataFolder);
    elseif ismember(':',patRange)
        nums = split(patRange,':');
        patRange = str2double(nums{1}):str2double(nums{2});
        patientNums = findAllPatientNums(dataFolder,patRange);
    elseif ismember(',',patRange)
        nums = split(patRange,',');
        patRange = str2double(nums);
        patientNums = findAllPatientNums(dataFolder,patRange);
    else
        patientNums = findAllPatientNums(dataFolder,patRange);
    end
elseif ischar(patientNums) %if only one pat, eg '47' adn not {'47','48'}
    patientNums = {patientNums};
end

end