function [lb,ub] = loadParamBounds_HEALTH(ind,paramUncertainty)
%difference from casas et al: 1e-3-1 instead of -0.5 - 0.5
%(convertion: onset_LV = (onset_LV_old*T + T)/T
%onset_LV: 1:2 instead of -0.5:0.5 --->take -1.5 when simulating
%onset LV: reduced the interval -> -0.1 - 0.1 (1.4-1.6) instead of 1-2
%onset LA: adapted based on LV before simulation: onsetLA = 1 + onsetLV - onsetLA
% -> new bounds: from 0.425:1 to 0.12:0.25
% See specific table in the manuscript.

% Litterature values
littnames = {'Rtot' 'Ctot' 'ElCo' 'Caa' 'Emax_LA' 'Emax_LV' 'Emin_LA' 'Emin_LV' 'Lao' 'Lav' 'Lmv' 'Ppu' 'Rao' 'Rmv' 'k_diast_LA' 'k_diast_LV' 'k_syst_LA' 'k_syst_LV' 'm1_LA' 'm1_LV' 'm2_LA' 'm2_LV' 'onset_LA' 'onset_LV','V0_LA','V0_LV'};
littVals = [0.94,1.48,4.7945,0.1, 0.170,3,0.080,0.080,0.00050,0.0004,0.00020,7.40,0.040,0.0037510,0.18,0.45200,0.110,0.2690,1.320,1.32,13.1,27.4,0.15,0.001]; %, 3, 10

lb = zeros(size(littVals));
ub = zeros(size(littVals));
ub([1,13,14]) = littVals([1,13,14]).*5; %resistances (casas et al 2)  %factor3 was slightly too low ub for Rao, so increased it
lb([1,13,14]) = littVals([1,13,14])./5;%factor 3 too low for rmv -> broadened the lower bound too
ub(9:11) = littVals(9:11).*6; % inertance (casas et al 5)
lb(9:11) = littVals(9:11)./6;
ub([2,4]) = littVals([2,4]).*7; % compliance (casas et al 6)
lb([2,4]) = littVals([2,4])./7;

ub(5:8) = littVals(5:8).*3; %Emax Emin LA & LV (casas et al 2) (Emax_LV set from data)
lb([6,8]) = littVals([6,8])./3; %Emax and Emin LV lb same as before (casas et al 2) (but Emax_LV set from data)
lb([5,7]) = littVals([5,7])./6; %Emax and Emin LA lb lowered to allow for no LA contraction in diseased patients

ub(15:16) = littVals(15:16).*3; % kdiast LA & LV (casas et al 2)
lb(15:16) = littVals(15:16)./3;
ub(16) = 0.9; % kdiast LV ub set to 0.9 (was unneccesarily wide)
ub(17:18) = littVals(17:18).*5; % ksyst LA & LV (casas et al 2). 3 before
lb(17:18) = littVals(17:18)./5; %

ub(19:22) = littVals(19:22).*3; % m1 and m2 LA & LV (casas et al 2)
lb(19:22) = littVals(19:22)./3; % 

ub(12) = littVals(12).*3; % Ppu prev: 12
lb(12) = littVals(12)./3; % Ppu prev: 6

ub(23:24) = [0.25,1.6]; % onset La & Lv, set so that they are physiologically reasonable close
lb(23:24) = [0.12,1.4]; % onset La & Lv

% Parameters for pulmonary veins % casas et al: kept constant since no data in pulmonary veins
pvconstants = [4,0.01,0.002,0.00050];%{'Cpvc','Rpu','Rpv','Lpv'}
pvlb  = pvconstants./[7,5,5,6]; %c,r,r,l 
pvub  = pvconstants.*[7,5,5,6]; %c,r,r,l
lb = [pvlb,lb]; 
ub = [pvub,ub];

% Caa, Aao, EOA_ao, Emax_LV, Ctot and Rtot are calculated from data and 
% get their bounds from there, but still above 0 for the lower bound  
if nargin >1 && ~isempty(paramUncertainty)
    multfactor = 0.75; %+-75% of the data value
    lb(ind.Caa)= max(0.01,paramUncertainty.Caa-paramUncertainty.Caa*multfactor);
    ub(ind.Caa)= paramUncertainty.Caa+paramUncertainty.Caa*multfactor;

    lb(ind.ElCo)=max(0.1,paramUncertainty.ElCo-paramUncertainty.ElCo*multfactor);
    ub(ind.ElCo)=paramUncertainty.ElCo+paramUncertainty.ElCo*multfactor;
    
    lb(ind.Emax_LV)=max(0.1,paramUncertainty.Emax_LV-paramUncertainty.Emax_LV*multfactor);
    ub(ind.Emax_LV)=paramUncertainty.Emax_LV+paramUncertainty.Emax_LV*multfactor;
    
    lb(ind.Ctot)=max(0.1,paramUncertainty.Ctot-paramUncertainty.Ctot*multfactor);
    ub(ind.Ctot)=paramUncertainty.Ctot+paramUncertainty.Ctot*multfactor;
    
    lb(ind.Rtot)=max(0.1,paramUncertainty.Rtot-paramUncertainty.Rtot*multfactor);
    ub(ind.Rtot)=paramUncertainty.Rtot+paramUncertainty.Rtot*multfactor;
end

%make sure no bounds are below or equal to 0
lb(lb<=0) = 1e-7;

end