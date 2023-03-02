function [ICs,dICs] = iccalcs(paramValues,constants,paramNames,constantsNames,data,ind,simulation)
% Calculate the initial conditions for a specific subject.

%% Pre-calculations otherwise made in the simulate function
T = constants(strcmp('T',constantsNames));

% Calculate normalizing factors for elastance function based on the parameters
constants(end) = calc_norm_factor(T,paramValues(ind.k_syst_LV),paramValues(ind.k_diast_LV),paramValues(ind.m1_LV),paramValues(ind.m2_LV));
constants(end-1) = calc_norm_factor(T,paramValues(ind.k_syst_LA),paramValues(ind.k_diast_LA),paramValues(ind.m1_LA),paramValues(ind.m2_LA));

%onset_LV: range set to 1-2 instead of -0.5 to 0.5 --> take onset_LV-1.5
paramValues(ind.onset_LV) = paramValues(ind.onset_LV) - 1.5;
%onset LA: adapted after onset LV to be close enough
paramValues(ind.onset_LA) = 1 + paramValues(ind.onset_LV) - paramValues(ind.onset_LA);

%% Set values 
%Params
Cpvc = paramValues(strcmp('Cpvc',paramNames));
Rpu = paramValues(strcmp('Rpu',paramNames));
Rpv = paramValues(strcmp('Rpv',paramNames));
Lpv  = paramValues(strcmp('Lpv',paramNames));
Rtot  = paramValues(strcmp('Rtot',paramNames));
Ctot  = paramValues(strcmp('Ctot',paramNames));
ElCo  = paramValues(strcmp('ElCo',paramNames));
Caa  = paramValues(strcmp('Caa',paramNames));
Emax_LA  = paramValues(strcmp('Emax_LA',paramNames));
Emax_LV  = paramValues(strcmp('Emax_LV',paramNames));
Emin_LA  = paramValues(strcmp('Emin_LA',paramNames));
Emin_LV  = paramValues(strcmp('Emin_LV',paramNames));
Lao  = paramValues(strcmp('Lao',paramNames));
Lav  = paramValues(strcmp('Lav',paramNames));
Lmv  = paramValues(strcmp('Lmv',paramNames));
Ppu  = paramValues(strcmp('Ppu',paramNames));
Rao  = paramValues(strcmp('Rao',paramNames));
Rmv  = paramValues(strcmp('Rmv',paramNames));
k_diast_LA  = paramValues(strcmp('k_diast_LA',paramNames));
k_diast_LV  = paramValues(strcmp('k_diast_LV',paramNames));
k_syst_LA  = paramValues(strcmp('k_syst_LA',paramNames));
k_syst_LV  = paramValues(strcmp('k_syst_LV',paramNames));
m1_LA  = paramValues(strcmp('m1_LA',paramNames));
m1_LV  = paramValues(strcmp('m1_LV',paramNames));
m2_LA  = paramValues(strcmp('m2_LA',paramNames));
m2_LV  = paramValues(strcmp('m2_LV',paramNames));
onset_LA  = paramValues(strcmp('onset_LA',paramNames));
onset_LV  = paramValues(strcmp('onset_LV',paramNames));
% V0_LA  = paramValues(strcmp('V0_LA',paramNames));
% V0_LV = paramValues(strcmp('V0_LV',paramNames));

%Constants
V0_LA  = constants(strcmp('V0_LA',constantsNames));
V0_LV = constants(strcmp('V0_LV',constantsNames));
tdiast  = constants(strcmp('tdiast',constantsNames));
Ks_LA = constants(strcmp('Ks_LA',constantsNames));
Ks_LV = constants(strcmp('Ks_LV',constantsNames));
RLAvisc = constants(strcmp('RLAvisc',constantsNames));
RLVvisc = constants(strcmp('RLVvisc',constantsNames));
Raa = constants(strcmp('Raa',constantsNames));
Rpc = constants(strcmp('Rpc',constantsNames));
Rpvc = constants(strcmp('Rpvc',constantsNames));
rho_blood = constants(strcmp('rho_blood',constantsNames));
norm_factor_LA = constants(strcmp('norm_factor_LA',constantsNames));
norm_factor_LV= constants(strcmp('norm_factor_LV',constantsNames));



Cpc = Ctot-Caa;
Rpr = Rtot-Rao;

%% From data:
Qpv = data.PV(1,2);
Qmv = 0;
Qav = max(data.AV(1,2),0);
Qaa = max(data.AC(1,2),0);
Paa = data.DBP - Raa*(Qav - Qaa); %P_aortic= Raa*(Qav - Qaa)  + Paa
Vlv = data.EDV;
if isnan(data.LaESV(2))
    %disp('No LA ESV value in data - setting to 65 in initial conditions calculation.')
    Vla = 65*2; %Mean laesv: 67.35, median 62, min 29, max 163. Setting missing value --> 65.
else
    Vla = data.LaESV(2)*2; %OBS: Simple estimation
end

mv_open = 0;
if Qav > 0
    av_open = 1;
else
    av_open = 0;
end

%% Derived based on data-vars
qLA = Qpv - Qmv;
qLV = Qmv - Qav;
Qcaa = Qav - Qaa;
dPaa = Qcaa/Caa;

%% The rest
% By assuming values to Qin/Qpvc/dPpvc and to Qpr/Qpc/dPpc, the rest of the initial conditions can be set from that and from data 

% Qpvc = dPpvc*Cpvc;
% Qin = Qpv + Qpvc;
Qin = 1.05*Qpv; % assumption
Qpvc = Qin - Qpv;
dPpvc = Qpvc/Cpvc;

% Qpc = dPpc*Cpc;
% Qpr = Qaa - Qpc;
Qpr = Qaa*0.95; % assumption
Qpc = Qaa - Qpr;
dPpc = Qpc/Cpc;

% From the assumptions we can now get:
Ppc = Qpr*Rpr - Rpc*Qpc;
Ppvc = Ppu - Rpu*Qin - Rpvc*Qpvc;
%dQaa = (Paa + Raa*Qcaa - Rao*Qaa - Rpc*Qpc - Ppc)/Lao; %set below
%dQpv = (Ppvc + Rpvc*Qpvc - Rpv*Qpv - P_LA)/Lpv; % set below

%%
Ts_LA = k_syst_LA*T;
Td_LA = k_diast_LA*T;
ula = T - onset_LA*T;

Ts_LV = k_syst_LV*T;
Td_LV = k_diast_LV*T;
ulv = T - onset_LV*T;

if ulv < T
    nLV = 0;
else
    nLV = 1;
end
% nLV = am_stepfun(ulv,T,1,2*T,2);
tlv = ulv - nLV*T;
g1_LV = (tlv/Ts_LV)^m1_LV;
g2_LV = (tlv/Td_LV)^m2_LV;
Elv = (1/norm_factor_LV)*(Emax_LV-Emin_LV)*(g1_LV/(1+ g1_LV)) *(1./(1 + g2_LV)) + Emin_LV;

if ula < T
    nLA = 0;
else
    nLA = 1;
end
% nLA= am_stepfun(ula,T,1,2*T,2);
tla = ula - nLA*T;
g1_LA = (tla/Ts_LA)^m1_LA;
g2_LA = (tla/Td_LA)^m2_LA;
Ela = (1/norm_factor_LA)*(Emax_LA-Emin_LA)*(g1_LA/(1+ g1_LA)) *(1/(1 + g2_LA)) + Emin_LA;

P_LA = Ela*(Vla - V0_LA)*(1-Ks_LA*qLA);
P_LV = Elv*(Vlv - V0_LV)*(1-Ks_LV*qLV);

Vf =0.00001;
Ron = 0.00001; %1.0000e-05
Goff = 1e-08;
Pvalve = Vf*(1-Ron*Goff);

P_Dmv = mv_open*(Qmv*Ron + Pvalve)  + (1-mv_open)*(Qmv*(1/Goff));
Rav = Qav * rho_blood/(2*ElCo^2)* (0.06/133.322);
P_Dav = av_open*(Qav*Ron + Pvalve)  + (1-av_open)*(Qav*(1/Goff));
Qpc_comp = (Qaa*Rpr - Ppc)/(Rpr+Rpc);


Qpr_comp = (Ppc + Rpc*dPpc*Cpc)/Rpr;
Qpvc_comp = (Ppu-Rpu*Qpv - Ppvc)/(Rpu+Rpvc);
Qin_comp = (Ppu-Rpvc*dPpvc*Cpvc - Ppvc)/Rpu;

% Set derivatives at t=0
dPpvc_comp = (Ppu - Rpu*Qpv - Ppvc)/(Cpvc*(Rpu + Rpvc)); % Ppvc
dQpv = (Ppvc + Rpvc*Qpvc - Rpv*Qpv - P_LA)/Lpv; %Qpv
dVla = qLA;%Qpv - Qmv; %Vla
dQmv = (P_LA - (RLAvisc+Rmv)*Qmv - P_Dmv - P_LV)/Lmv;%Qmv
dVlv = qLV;%Qmv-Qav ; % Vlv
dQav = (P_LV - (RLVvisc+Rav)*Qav - P_Dav - Raa*Qcaa - Paa)/Lav; %Qav
% dPaa = (Qav-Qaa)/Caa; %Paa 
dQaa = (Paa + Raa*Qcaa - Rao*Qaa - Rpc*Qpc - Ppc)/Lao; %Qaa
dPpc_comp = (Qaa*Rpr - Ppc)/(Cpc*(Rpr+Rpc));

%% Sum up
dICs = [dPpvc,dQpv,dVla,dQmv,dVlv,dQav,dPaa,dQaa,dPpc,0,0];
ICs = [Ppvc,Qpv,Vla,Qmv,Vlv,Qav,Paa,Qaa,Ppc,mv_open,av_open];

if sum(isnan(ICs)) > 0
    disp('!!! Initial conditions contains NaN values. Check the data. !!!')
end

%% Check that all conditions are fullfilled
values = [Qpr,Qpc,dPpvc,dPpc,Qpvc,Qin]';
compValues = [Qpr_comp,Qpc_comp,dPpvc_comp,dPpc_comp,Qpvc_comp,Qin_comp]';
compareTable = table(values, compValues,(values-compValues)./values);%diff = 0

kcl1 = Qin - Qpvc-Qpv;
kcl2 = Qpv - qLA - Qmv;
kcl3 = Qmv - qLV - Qav;
kcl4 = Qav - Qcaa - Qaa;
kcl5 = Qaa - Qpc - Qpr;
kvl1 = Ppu-Rpu*Qin - Rpvc*Qpvc - Ppvc;
kvl2 = Ppvc + Rpvc*Qpvc - Lpv*dQpv - Rpv*Qpv - P_LA;
kvl3 = P_LA - RLAvisc*Qmv - P_Dmv - Lmv*dQmv - Rmv*Qmv - P_LV;
kvl4 = P_LV - RLVvisc*Qav - P_Dav - Lav*dQav - Rav*Qav - Raa*Qcaa - Paa;
kvl5 = Paa + Raa*Qcaa - Rao*Qaa - Lao*dQaa - Rpc*Qpc - Ppc;
kvl6 = Ppc + Rpc*Qpc - Qpr*Rpr;
kcl_c1 = Qin - dPpvc*Cpvc - Qpv;
kcl_c2 = Qav - dPaa*Caa - Qaa;
kcl_c3 = Qaa - dPpc*Cpc - Qpr;

lawTable = table(kcl1,kcl2,kcl3,kcl4,kcl5,kvl1,kvl2,kvl3,kvl4,kvl5,kvl6,kcl_c1,kcl_c2,kcl_c3);

if nargin > 6
    simVScalc = table(simulation.x0,simulation.x(1,:)',ICs',ICs'-simulation.x(1,:)')
end

end