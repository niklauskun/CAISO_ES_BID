clc
clear all %#ok<CLALL> 
close all
tic
%% setting
load('WALNUT.mat')
Ys = 2016; % year start
Ye = 2016; % year end
Is = (Ys-2016)*365+1;
Ie = (Ye-2015)*365;
T = 1440; %totl periods
H = 24;
Ts = H/T; % time interval between two periods (hour)

E = 1;
Pr = 0.25*E; % normalized power rating wrt energy rating
P = Pr*Ts; % actual power rating taking time step size into account: now 10MW/40MWh
ed = .001; % SoC sample granularity
e0 = .0; % initial SoC
ef = .0; % final SoC target level, use 0 if none
Ne = floor(1/ed)+1; % number of SOC samples
seg_num = 5;% value function segment number
sigma = 0;

% % constant parameters
% eta = .9; % efficiency
% c = 25; % marginal discharge cost - degradation

% nonlinear parameters
Elist = [0.2 0.2 0.2 0.2 0.2];
etalist = [0.8 0.85 0.9 0.85 0.8];
clist = [25 20 15 20 25];
eta = zeros(Ne,1);
c = zeros(Ne,1);
eta(1) = etalist(1); % set value for 0 SoC
c(1) = clist(1); % set value for 0 SoC
ss = 1;
for i = 1:numel(Elist)
    se = ss + (Ne-1)*Elist(i);
    eta(ss+1:se) = etalist(i);
    c(ss+1:se) = clist(i);
    ss = ss + (Ne-1)*Elist(i);
end

%% bid calculation
vpAvg = zeros(seg_num,(Ie-Is+1)*H);
vbAvg = zeros(seg_num,(Ie-Is+1)*H);
for i = Is:Ie
    [vpAvg(:,(i-1)*H+1:i*H),vbAvg(:,(i-1)*H+1:i*H)] = value_fcn_calcu(Q(:,i),seg_num, Ne, T, H, c, P, eta, ed, ef, sigma);
    fprintf('Day %i... \n',i);
end
writematrix(vpAvg,'WALNUT_2016_Discharge.csv');
writematrix(vbAvg,'WALNUT_2016_Charge.csv');
timeElapse = toc;