function [simAmi,sol] = simulate_avatarHEALTH(theta,constants,options,numHeartBeats,ind,simtime)
% Do a simulation of the avatar model. First, set parameter values, then
% simulate a number of heartbeats one at a time.
%step = 0.001; %good resolution when plotting simulation
%step = T/39; %comparing to data (since there are 40 timeframes)
%ind is a struct containing indexes for the parameters
%options contains simulation options
T = constants(ind.T);

% Calculate normalizing factors for elastance function based on the parameters
norm_factor_LV = calc_norm_factor(T,theta(ind.k_syst_LV),theta(ind.k_diast_LV),theta(ind.m1_LV),theta(ind.m2_LV));
norm_factor_LA= calc_norm_factor(T,theta(ind.k_syst_LA),theta(ind.k_diast_LA),theta(ind.m1_LA),theta(ind.m2_LA));
constants(end) = norm_factor_LV;
constants(end-1) = norm_factor_LA;

%onset_LV: range set to 1-2 instead of -0.5 to 0.5 --> take onset_LV-1.5
theta(ind.onset_LV) = theta(ind.onset_LV) - 1.5;
%onset LA: adapted after onset LV to be close enough
theta(ind.onset_LA) = 1 + theta(ind.onset_LV) - theta(ind.onset_LA);

% First simulation of one heartbeat - to get the sizes of x and y
simtime = simtime+T;
options.tstart = T;
sol = simulate_avatar_HEALTH_fast(simtime,theta, constants, [], options);

options.x0 = sol.x(end,:)';
simlen = length(sol.t)-1;
simAmi.x = zeros(simlen*numHeartBeats,size(sol.x,2));
simAmi.y = zeros(simlen*numHeartBeats,size(sol.y,2));
simAmi.t = zeros(simlen*numHeartBeats,size(sol.t,2));

simAmi.x(1:length(sol.t),:) = sol.x;
simAmi.y(1:length(sol.t),:) = sol.y;
simAmi.t(1:length(sol.t)) = sol.t;

% The rest of the simulations
for i = 1:numHeartBeats-1
    sol = simulate_avatar_HEALTH_fast(simtime,theta, constants, [], options);

    %start next sim with end values from this sim
    options.x0 = sol.x(end,:)';
    
    %save this sim together with the rest
    simAmi.x(simlen*i+2:simlen*(i+1)+1,:) = sol.x(2:end,:);
    simAmi.y(simlen*i+2:simlen*(i+1)+1,:) = sol.y(2:end,:);
    simAmi.t(simlen*i+2:simlen*(i+1)+1)   = sol.t(2:end)+T*i;
end

 % To compare with data we use the last heartbeat when a "steady state" is established.
 % We need to round the time to 3 digits, and start at t=0 to be able to
 % compare with data.
sol.t = round(sol.t-T,3);

end