function []=Setup(patientNums,doParamEst,datename,patRange,dataFolder,doPV,method,continueMCMC)

pathbase = split(pwd,'cardiovascularavatar_T2D_HT');
pathbase = [pathbase{1},filesep, 'cardiovascularavatar_T2D_HT'];

addpath(genpath(pathbase))
addpath([pathbase filesep 'Parameters'])
addpath([pathbase filesep 'Modelfiles'])

addpath([pathbase filesep 'Requirements' filesep 'PESTO-1.1.0'])
addpath([pathbase filesep 'Requirements' filesep 'AMICI-0.10.11_SS_eventFix' filesep 'matlab'])
addpath([pathbase filesep 'Requirements' filesep 'MEIGO'])

run([pathbase filesep 'Requirements' filesep 'AMICI-0.10.11_SS_eventFix' filesep 'matlab' filesep 'installAMICI.m'])
run([pathbase filesep 'Requirements' filesep 'MEIGO' filesep 'install_MEIGO.m'])

%METHOD:
% inverse profile likelihood: 'invPL'
% inverse prediction profile likelihood: 'invPredPL'
% profile likelihood: 'PL'
% "normal" parameter optimization: 'opt'


if doParamEst
    %% if no patientnums given, find all
    patientNums = findPatientRange(patientNums,patRange,dataFolder);
    
    %% create folders to save results in
    % Find names of results folder depending on what type of param
    % estimation method
    if nargin > 5
        if strcmp(method,'invPL') || strcmp(method,'PL')
            if doPV
                resultsfolder = 'Parameters/ParameterUncertaintyPV';
            else
                resultsfolder = 'Parameters/ParameterUncertainty';
            end
        elseif strcmp(method,'invPredPL')
            if doPV
                resultsfolder = 'Parameters/UncertaintyPV';
            else
                resultsfolder = 'Parameters/Uncertainty';
            end
        elseif strcmp(method,'MCMC')
                resultsfolder = 'Parameters/MCMC';
        elseif strcmp(method,'HEALTH')
            resultsfolder = 'Parameters/ESS';
        elseif doPV
            resultsfolder = 'Parameters/Fitted to PV';
        else
            resultsfolder = 'Parameters';
        end
    else
        resultsfolder = 'Parameters';
    end
    resultsfolder = fullfile(pathbase,resultsfolder);
    
    disp(['num pats: ' num2str(length(patientNums))])
    if continueMCMC
        disp('continuing on old mcmc run')
    else
        for p = 1:length(patientNums)
            if strcmp(method,'invPL')
                folder = fullfile(resultsfolder, ['P' patientNums{p}]);
            else
                folder = fullfile(resultsfolder, ['P' patientNums{p} '_' datename]);
            end
            if isempty(dir(folder))
                mkdir(folder)
                disp(['folder created: ' folder])
            else
                disp(['using existing folder: ' folder])
            end
        end
    end
    clear folder
end

disp('Path setup done.')
end
