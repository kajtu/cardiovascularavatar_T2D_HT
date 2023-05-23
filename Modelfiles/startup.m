%% Adds paths to toolboxes needed for simulations

pathbase = split(pwd,'cardiovascularavatar_T2D_HT');
pathbase = [pathbase{1},filesep, 'cardiovascularavatar_T2D_HT'];

addpath(genpath(pathbase))
addpath([pathbase filesep 'Modelfiles'])

addpath([pathbase filesep 'Requirements' filesep 'PESTO-1.1.0'])
addpath([pathbase filesep 'Requirements' filesep 'AMICI-0.10.11_SS_eventFix' filesep 'matlab'])
addpath([pathbase filesep 'Requirements' filesep 'MEIGO'])

run([pathbase filesep 'Requirements' filesep 'AMICI-0.10.11_SS_eventFix' filesep 'matlab' filesep 'installAMICI.m'])
run([pathbase filesep 'Requirements' filesep 'MEIGO' filesep 'install_MEIGO.m'])