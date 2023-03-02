function norm_factor = calc_norm_factor(T,k_syst,k_diast,m1,m2)
N=200;      %Length of the time axis vector
T_axis=linspace(0,T,N);

Ts = k_syst*T;
Td = k_diast*T;
g1=(T_axis/Ts).^m1;
g2=(T_axis/Td).^m2;
prod=(g1./(1+g1)).*(1./(1+g2));
norm_factor = max(prod);

end