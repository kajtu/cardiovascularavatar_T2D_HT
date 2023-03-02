% Compile model files into mexfiles
% This automatically creates the files: amai_avatar_HEALTH_fast.mex.. and simulate_avatar_HEALTH_fast.m
clear 
clear mex
startup

path = pwd;

amiwrap('avatar_HEALTH_fast','avatar_HEALTH_syms_fast',path); 

disp('model generated')


