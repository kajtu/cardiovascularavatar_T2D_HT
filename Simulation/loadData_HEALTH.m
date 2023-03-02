function [data,extradata,paramdata] = loadData_HEALTH(patientNum)

inputpathbase = split(pwd,'cardiovascularavatar_T2D_HT');
inputpath = [inputpathbase{1},filesep 'cardiovascularavatar_T2D_HT' filesep 'Data'];
try
    load(fullfile(inputpath,['inputdata_avatarmodel_P',patientNum,'.mat']),'data','extradata','paramdata')
catch
    load(fullfile(inputpath,['inputdata_avatarmodel_',patientNum,'.mat']),'data','extradata','paramdata')
end

end
