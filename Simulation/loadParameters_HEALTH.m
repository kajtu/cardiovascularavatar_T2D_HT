function [paramValues,constants,paramNames,constantsNames,units,ind,dataParams] = loadParameters_HEALTH(paramdata)
% Set parameter values

%% Initial parameter values from literature
%Density
rho_blood=1.06; %g/ml

%Pulmonary vessels (from Sun, also in Tanné)
Ppu=7.4;
Rpu=0.01;
Rpvc=0.01;
Cpvc=4;
Lpv=0.0005;
Rpv=0.002;

%Varying elastance function
m1_LA=1.32;    %()         Mynard et al.
m2_LA=13.1;    %()         Mynard et al.
Emin_LA=0.08;  %(mmHg/mL)  Mynard et al.
Emax_LA=0.17;  %(mmHg/mL)  Mynard et al.
V0_LA=3;       %(mL)       Mynard et al.
Ks_LA=10e-9;   %(s/mL)     Mynard et al.
onset_LA=0.15; %()         modified value - 0.85  Mynard et al.
k_syst_LA=0.110;%(s) Systolic time constant (Mynard et al.)
k_diast_LA=0.180;%(s)Diastolic time constant (Mynard et al.)

%Viscous resistance
RLAvisc=0.0001;

%Mitral valve parameters
Rmv=3.751e-3;
Lmv=0.0002;

%Left ventricle parameters

%Varying elastance function
m1_LV=1.32;            %()         Mynard et al.
m2_LV=27.4;            %()         Mynard et al.
Emin_LV=0.08;          %(mmHg/mL)  Mynard et al.
V0_LV=10;              %(mL)       Mynard et al.
Ks_LV=4e-9;            %(s/mL)     Mynard et al.
onset_LV=0.001+1.5;    %()Mynard et al.    %log10(0) = -inf --> set to 1e-3 (+1.5 for opt)
k_syst_LV=0.269;
k_diast_LV=0.452;

%Calculate the normalizing factor (Mynard et al.)
norm_factor_LV = calc_norm_factor(paramdata.T,k_syst_LV,k_diast_LV,m1_LV,m2_LV);
norm_factor_LA= calc_norm_factor(paramdata.T,k_syst_LA,k_diast_LA,m1_LA,m2_LA);

%Viscous resistance
RLVvisc=0.0001;         %Sun et al.

%Aortic valve parameters
Lav=0.0004;             %mmHg*s^2/mL

%Ascending aorta parameters
Raa=0.01;
Rao=0.04;
Lao=5e-4;

%Peripheral vessels parameters
Rpc=0.01;

%% Secondary parameters computed from patient data
HR=60*1/paramdata.T;          %Heart rate (beats/min)
PP=paramdata.SBP-paramdata.DBP;         %Peak pressure
MAP=paramdata.DBP+(1/3+HR*0.0012)*(paramdata.SBP-paramdata.DBP); % Beléns calculation of MAP, takes HR in consideration azminia M, Trivedi A, Molnar J, Elbzour M, Guerrero M, Salem Y, Ahmed A, Khosla S, Lubell DL . Validation of a new formula for mean arterial pressure calculation: the new formula is superior to the standard formula. Catheter Cardiovasc Interv 2004; 63: 419–425.

%% Calculate parameters on a subject-specific basis
%Calcualted automatically based on the input measurements:
B_av=rho_blood/(2*paramdata.ElCo^2);
maxQav = max(paramdata.AV_Flow);
deltaPconst = 0.06/133.322;% tagit bort 0.9 här.
deltaP_av = B_av*(maxQav^2)*deltaPconst; %Pressure gradient accross the aortic valve
k=0.9; %Constant to relate systolic blood pressure and end systolic blood pressure: Pes=k*Ps
Emax_LV=k*(paramdata.SBP+deltaP_av)/(paramdata.ESV_seg-V0_LV); %Maximal elastance of the LV


meanQav = mean(paramdata.AV_Flow);
Ctot=meanQav/PP;  %Total compliance of the system

Rtot=MAP/meanQav;     %Total resistance of the system


%% Data-based parameters
% Caa, Aao, EOA_ao, Emax_LV, Ctot and Rtot are calculated from data 
dataParams.Emax_LV = Emax_LV;
dataParams.Ctot = Ctot;
dataParams.Rtot = Rtot;
dataParams.ElCo = paramdata.ElCo;
dataParams.Caa = paramdata.Caa;%pwv fr full aorta, volume fr ascao

%% Set parameters to be optimized and constant parameter values
paramNames = {'Cpvc' 'Rpu' 'Rpv' 'Lpv' 'Rtot' 'Ctot' 'ElCo' 'Caa' 'Emax_LA' 'Emax_LV' 'Emin_LA' 'Emin_LV' 'Lao' 'Lav'...
    'Lmv' 'Ppu' 'Rao' 'Rmv' 'k_diast_LA' 'k_diast_LV' 'k_syst_LA'...
    'k_syst_LV' 'm1_LA' 'm1_LV' 'm2_LA' 'm2_LV' 'onset_LA' 'onset_LV'};
constantsNames = {'aaCorr' 'avCorr' 'mvCorr' 'tdiast' 'Ks_LA' 'Ks_LV'...
    'V0_LA' 'V0_LV' 'RLAvisc' 'RLVvisc' 'Raa' 'Rpc' 'Rpvc' 'T' 'rho_blood' 'norm_factor_LA' 'norm_factor_LV'};
constants = [0 0 0 paramdata.tdiast Ks_LA Ks_LV V0_LA V0_LV RLAvisc RLVvisc Raa Rpc Rpvc paramdata.T rho_blood norm_factor_LA norm_factor_LV];

%R: 'mmHg*s/ml'
%C:'ml/mmHg
%L: 'mmHg*s^2/ml'
units.param = {'ml/mmHg','mmHg*s/ml','mmHg*s/ml','mmHg*s^2/ml','mmHg*s/ml',...
    'ml/mmHg','cm^2','ml/mmHg','ml/mmHg','ml/mmHg','ml/mmHg','ml/mmHg',...
    'mmHg*s^2/ml','mmHg*s^2/ml','mmHg*s^2/ml','mmHg','mmHg*s/ml','mmHg*s/ml',...
    '-','-','-','-','-','-','-','-','s','s'};
units.constant = {'ml','ml','ml','s','s/ml','s/ml','ml','ml','mmHg*s/ml','mmHg*s/ml','mmHg*s/ml','mmHg*s/ml','s','g/ml','-','-'};

Caalitterature=0.1;
Emax_LVlitterature = 3;
paramValues = [Cpvc Rpu Rpv Lpv Rtot Ctot paramdata.ElCo Caalitterature Emax_LA Emax_LVlitterature Emin_LA Emin_LV Lao Lav Lmv Ppu Rao Rmv k_diast_LA k_diast_LV k_syst_LA k_syst_LV m1_LA m1_LV m2_LA m2_LV onset_LA onset_LV];

%% indices for easy access to model parameters, states and observables
ind.MV = 4;
ind.AV = 6;
ind.AC = 8;
ind.LV = 5;
ind.PV = 2;
ind.aaCorrVol = 19;
ind.avCorrVol = 20;
ind.mvCorrVol = 21;
ind.pvCorrVol = 22;
ind.Vla = 23;
ind.aorticPressure = 1;
ind.brachialPressure = 25;

ind.k_syst_LV = strcmp(paramNames,'k_syst_LV');
ind.k_diast_LV = strcmp(paramNames,'k_diast_LV');
ind.m1_LV = strcmp(paramNames,'m1_LV');
ind.m2_LV= strcmp(paramNames,'m2_LV');
ind.k_syst_LA = strcmp(paramNames,'k_syst_LA');
ind.k_diast_LA = strcmp(paramNames,'k_diast_LA');
ind.m1_LA = strcmp(paramNames,'m1_LA');
ind.m2_LA = strcmp(paramNames,'m2_LA');

ind.Emax_LA = strcmp(paramNames,'Emax_LA');
ind.Emin_LA = strcmp(paramNames,'Emin_LA');
ind.Emax_LV = strcmp(paramNames,'Emax_LV');
ind.Emin_LV = strcmp(paramNames,'Emin_LV');

ind.onset_LV = strcmp(paramNames,'onset_LV');
ind.onset_LA = strcmp(paramNames,'onset_LA');

ind.Caa = strcmp(paramNames,'Caa');
ind.ElCo = strcmp(paramNames,'ElCo');
ind.Ctot = strcmp(paramNames,'Ctot');
ind.Rtot = strcmp(paramNames,'Rtot');

ind.avCorr = strcmp(constantsNames,'avCorr');
ind.mvCorr = strcmp(constantsNames,'mvCorr');
ind.aaCorr = strcmp(constantsNames,'aaCorr');

ind.tdiast = paramdata.indtdiast;

ind.T = strcmp(constantsNames,'T');



end