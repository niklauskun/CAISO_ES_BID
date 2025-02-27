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
Pmax = Pr*Ts; % actual max power rating taking time step size into account: now 0.25MW/1MWh
ed = .001; % SoC sample granularity
e0 = .0; % initial SoC
ef = .0; % final SoC target level, use 0 if none
Ne = floor(1/ed)+1; % number of SOC samples
seg_num = 5;% value function segment number
rng('default')
sigma = 0;
noise = sigma*randn(size(Q,1), size(Q,2));
Q = Q+noise;
% % constant parameters
% eta = .9; % efficiency
% c = 20; % marginal discharge cost - degradation
% P1 = Pmax;
% P2 = Pmax;

% nonlinear parameters
Elist = [0.2 0.2 0.2 0.2 0.2];
etalist = [0.8 0.85 0.9 0.85 0.8];
clist = [25 20 15 20 25];
P1list = [0.7 1.0 1.0 0.9 0.5]*Pmax;
P2list = [0.5 0.8 1.0 1.0 1.0]*Pmax;
eta = zeros(Ne,1);
c = zeros(Ne,1);
P1 = zeros(Ne,1);
P2 = zeros(Ne,1);
eta(1) = etalist(1); % set value for 0 SoC
c(1) = clist(1); % set value for 0 SoC
P1(1) = P1list(1); % set value for 0 SoC
P2(1) = P2list(1); % set value for 0 SoC
ss = 1;
for i = 1:numel(Elist)
    se = ss + (Ne-1)*Elist(i);
    eta(ss+1:se) = etalist(i);
    c(ss+1:se) = clist(i);
    P1(ss+1:se) = P1list(i);
    P2(ss+1:se) = P2list(i);
    ss = ss + (Ne-1)*Elist(i);
end

%% bid calculation
vpAvg = zeros(seg_num,(Ie-Is+1)*H);
vbAvg = zeros(seg_num,(Ie-Is+1)*H);
for i = Is:Ie
    [vpAvg(:,(i-1)*H+1:i*H),vbAvg(:,(i-1)*H+1:i*H)] = value_fcn_calcu(Q(:,i),seg_num, Ne, T, H, c, P1, P2, eta, ed, ef);
    fprintf('Day %i... \n',i);
end
writematrix(vpAvg,'WALNUT_2016_Discharge.csv');
writematrix(vbAvg,'WALNUT_2016_Charge.csv');
timeElapse = toc;