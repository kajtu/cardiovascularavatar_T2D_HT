function [model] = avatar_HEALTH_syms_fast()

% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'lin';

%%
% STATES
% create state syms 
syms Ppvc Qpv Vla Qmv Vlv Qav Paa Qaa Ppc mv_open av_open
% create state vector
model.sym.x = [Ppvc Qpv Vla Qmv Vlv Qav Paa Qaa Ppc mv_open av_open];
%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms        
syms Cpvc Rpu Rpv Lpv Rtot Ctot ElCo Caa Emax_LA Emax_LV Emin_LA Emin_LV Lao Lav Lmv Ppu Rao Rmv k_diast_LA k_diast_LV k_syst_LA k_syst_LV m1_LA m1_LV m2_LA m2_LV onset_LA onset_LV
% create parameter vector 
model.sym.p = [Cpvc Rpu Rpv Lpv Rtot Ctot ElCo Caa Emax_LA Emax_LV Emin_LA Emin_LV Lao Lav Lmv Ppu Rao Rmv k_diast_LA k_diast_LV k_syst_LA k_syst_LV m1_LA m1_LV m2_LA m2_LV onset_LA onset_LV];
%%  
% CONSTANTS ( for these no sensitivities will be computed )

% create constant parameter syms
syms  aaCorr avCorr mvCorr tdiast Ks_LA Ks_LV  V0_LA V0_LV RLAvisc RLVvisc Raa Rpc Rpvc T rho_blood norm_factor_LA norm_factor_LV

% create constant parameter vector  
model.sym.k = [ aaCorr avCorr mvCorr tdiast Ks_LA Ks_LV  V0_LA V0_LV RLAvisc RLVvisc Raa Rpc Rpvc T rho_blood norm_factor_LA norm_factor_LV];

%%
% SYSTEM EQUATIONS
% create symbolic variable for time
syms t 
model.sym.xdot = sym(zeros(size(model.sym.x)));

%% Elastance equations
syms P_LA P_LV

Ts_LA = k_syst_LA*T;
Td_LA = k_diast_LA*T;
ula = t - onset_LA*T;

Ts_LV = k_syst_LV*T;
Td_LV = k_diast_LV*T;
ulv = t - onset_LV*T;

nLV = am_if(am_gt(ulv,T),1,0) + am_if(am_gt(ulv,T*2),1,0);
tlv = ulv - nLV*T;
g1_LV = (tlv/Ts_LV)^m1_LV;
g2_LV = (tlv/Td_LV)^m2_LV;
Elv = (1/norm_factor_LV)*(Emax_LV-Emin_LV)*(g1_LV/(1+ g1_LV)) *(1./(1 + g2_LV)) + Emin_LV;

nLA = am_if(am_gt(ula,T),1,0) + am_if(am_gt(ula,T*2),1,0);
tla = ula - nLA*T;
g1_LA = (tla/Ts_LA)^m1_LA;
g2_LA = (tla/Td_LA)^m2_LA;
Ela = (1/norm_factor_LA)*(Emax_LA-Emin_LA)*(g1_LA/(1+ g1_LA)) *(1/(1 + g2_LA)) + Emin_LA;

qLA = Qpv - Qmv;
P_LA = Ela*(Vla - V0_LA)*(1-Ks_LA*qLA);

qLV = Qmv - Qav;
P_LV = Elv*(Vlv - V0_LV)*(1-Ks_LV*qLV);

%% Pulmonary venous system
model.sym.xdot(1) = (Ppu - Rpu*Qpv - Ppvc)/(Cpvc*(Rpu + Rpvc)); % Ppvc pulm pressure change
Qpvc = ((Ppu - Rpu*Qpv - Ppvc)/(Rpu + Rpvc));
model.sym.xdot(2) = (Ppvc + Rpvc*Qpvc - Rpv*Qpv - P_LA)/Lpv; %Qpv pulm flow change

P_Pulmvein =P_LA + Qpv*Rpv;

%% Left atrium
model.sym.xdot(3) = Qpv - Qmv; %Vla

%% Mitral valve
Vf =0.00001;
Ron = 0.00001; %1.0000e-05
Goff = 1e-08;
Pvalve = Vf*(1-Ron*Goff);%1.0000e-05

P_Dmv = mv_open*(Qmv*Ron + Pvalve)  + (1-mv_open)*(Qmv*(1/Goff));
model.sym.xdot(4) = (P_LA - (RLAvisc+Rmv)*Qmv - P_Dmv - P_LV)/Lmv;%Qmv

%% Left ventricle
model.sym.xdot(5) = Qmv-Qav ; % Vlv

%% Aortic valve
%ElCo = (EOA_av*A_ao)/(A_ao-EOA_av);
Rav = Qav * rho_blood/(2*ElCo^2)* (0.06/133.322);
Qcaa = Qav - Qaa;

P_Dav = av_open*(Qav*Ron + Pvalve)  + (1-av_open)*(Qav*(1/Goff));
model.sym.xdot(6) = (P_LV - (RLVvisc+Rav)*Qav - P_Dav - Raa*Qcaa - Paa)/Lav; %Qav


%% Ascending aorta
model.sym.xdot(7) = (Qav-Qaa)/Caa; %Paa 
P_Aortic = Raa*(Qav - Qaa)  + Paa; %Raa constant
Cpc = Ctot-Caa;
Rpr = Rtot-Rao;
Qpc = (Qaa*Rpr - Ppc)/(Rpr+Rpc);
model.sym.xdot(8) = (Paa + Raa*Qcaa - Rao*Qaa - Rpc*Qpc - Ppc)/Lao; %Qaa


P_Brachial = P_Aortic;%conversion of SBP using SBPdiff is done outside the model

%% Peripheral arteries/Vascular bed
model.sym.xdot(9) = (Qaa*Rpr - Ppc)/(Cpc*(Rpr+Rpc));   %Ppc
Pperipheral = Ppc + Rpc*Qpc;

model.sym.xdot(10) = 0;%mv_open
model.sym.xdot(11) = 0;%av_open

%%
% INITIAL CONDITIONS
% Ppvc Qpv Vla Qmv Vlv Qav Paa Qaa Ppc 
model.sym.x0(1) = 10.6845;
model.sym.x0(2) = 11.3164;
model.sym.x0(3) = 89.4388; %LA
model.sym.x0(4) = 11.3164;
model.sym.x0(5) = 17.5520; % LV
model.sym.x0(6) = 11.3164;
model.sym.x0(7) = 10.6218;
model.sym.x0(8) = 11.3164;
model.sym.x0(9) = 9.7165; %Ppc
model.sym.x0(10) = 1;
model.sym.x0(11) = 1;

model.sym.dx0 = sym(zeros(size(model.sym.x)));


% OBSERVALES
%'P ascAo','P peripheral','Pressgrad MV calc','Dmv','mv_open','Dav','av open','pressgrad AV calc','Ela','Elv'
model.sym.y = sym(zeros(25,1)); 
model.sym.y(1) = P_Aortic;
model.sym.y(2) = Pperipheral;
model.sym.y(3) = P_LA - RLAvisc*Qmv - P_LV; %pressGrad_MV
model.sym.y(4) = P_Dmv;
model.sym.y(5) = mv_open;
model.sym.y(6) = P_Dav;
model.sym.y(7) = av_open;
model.sym.y(8) = P_LV + RLVvisc*Qav + Raa*Caa*((Qav-Qaa)/Caa);%pressGrad_AV;
model.sym.y(9) = Ela;
model.sym.y(10) = Elv;
model.sym.y(11) = Qcaa;
model.sym.y(12) = Qpc;
model.sym.y(13) = P_Pulmvein;
model.sym.y(14) = Qpvc;
model.sym.y(15) = qLA;
model.sym.y(16) = qLV;
model.sym.y(17) = P_LA;
model.sym.y(18) = P_LV;
model.sym.y(19) = Qaa + am_max(0,am_min(1,Qaa-5))*am_max(0,(tdiast-tlv)/abs(tdiast-tlv))*aaCorr;
model.sym.y(20) = Qav + am_max(0,am_min(1,Qav-5))*am_max(0,(tdiast-tlv)/abs(tdiast-tlv))*avCorr;
model.sym.y(21) = Qmv + am_max(0,am_min(1,Qmv-5))*am_max(0,(tlv-tdiast)/abs(tlv-tdiast))*mvCorr;
model.sym.y(22) = Qpv; 
model.sym.y(23) = Vla; 
model.sym.y(24) = Vlv; 
model.sym.y(25) = P_Brachial;

%%
%-- EVENTS
%syms t

% events fire when there is a -zero crossing of the root function
% events are used to automatically open and close the valves
%Ppvc Qpv Vla Qmv Vlv Qav Paa Qaa Ppc mv_open av_open

model.event(1) = amievent(P_Dmv - Vf,[0 0 0 0 0 0 0 0 0 1 0],[]);
model.event(2) = amievent(-P_Dmv + Vf,[0 0 0 0 0 0 0 0 0 -1 0],[]);
model.event(3) = amievent(P_Dav - Vf,[0 0 0 0 0 0 0 0 0 0 1],[]);
model.event(4) = amievent(-P_Dav + Vf,[0 0 0 0 0 0 0 0 0 0 -1],[]);

end