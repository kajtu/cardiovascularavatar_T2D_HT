function [groups] = loadGroupIndexes(usedPatNums)
% load groups
inputpathbase = split(pwd,'cardiovascularavatar_T2D_HT');
datapath = [inputpathbase{1},'cardiovascularavatar_T2D_HT',filesep,'Data'];
load(fullfile(datapath,'scapisdata_HEALTH.mat'),'SCAPISdata_HEALTH')

usedPatNums=removeUnderscore(usedPatNums);

%% sort and group on scapis home pressure and T2D data
% healthIDs are anonymized IDs connected to the input data and optimized parameters for each subject.
patinds = zeros(size(usedPatNums));
usedPatNumsP = cell(size(usedPatNums));
for p = 1:length(usedPatNums)
    patientNum=removeUnderscore(usedPatNums(p));
    patientNum = ['P' patientNum{1}];
    usedPatNumsP{p} = patientNum;
    patinds(p) = find(strcmp(SCAPISdata_HEALTH.healthIDs,patientNum));
end
patnumcheck = strcmp(SCAPISdata_HEALTH.healthIDs(patinds),usedPatNumsP');
if sum(patnumcheck) ~= length(usedPatNums)
    disp('OBS! check the group indexes')
end

groups.hypertensionHome = logical(SCAPISdata_HEALTH.HTgroup(patinds));
groups.nothypertensionHome = ~groups.hypertensionHome;
groups.T2D = SCAPISdata_HEALTH.T2Dgroup(patinds);
groups.Control = ~groups.T2D;

groups.hypertensionMR = logical(SCAPISdata_HEALTH.HTgroupMR(patinds));
groups.nothypertensionMR = ~groups.hypertensionMR;

% four groups but highBP defined as home bp (SBP >=135 or DBP>=85)
groups.T2D_HT_home = groups.T2D & groups.hypertensionHome;
groups.C_HT_home = groups.Control & groups.hypertensionHome;
groups.T2D_NT_home = groups.T2D & groups.nothypertensionHome;
groups.C_NT_home = groups.Control & groups.nothypertensionHome;

% create indexes for the four groups 
groups.fourgroupsHBP = cell(size(groups.T2D_HT_home));
groups.fourgroupsHBP(groups.C_NT_home) = {'HBP C'};
groups.fourgroupsHBP(groups.T2D_NT_home) = {'HBP T2D'};
groups.fourgroupsHBP(groups.C_HT_home) = {'HBP HT'};
groups.fourgroupsHBP(groups.T2D_HT_home) = {'HBP T2D & HT'};

end