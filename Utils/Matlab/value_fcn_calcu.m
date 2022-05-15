function [vpAvg vbAvg] = value_fcn_calcu(lambda,seg_num,Ne, T, H, c, P1, P2, eta, ed, ef)

vEnd = zeros(Ne,1);  % generate value function samples
vEnd(1:floor(ef*Ne)) = 1e3; % use 100 as the penalty for final discharge level
v = zeros(Ne, T+1); % initialize the value function series
% v(1,1) is the marginal value of 0% SoC at the beginning of day 1
% V(Ne, T) is the maringal value of 100% SoC at the beginning of the last operating day
v(:,end) = vEnd; % update final value function
% process index
es = (0:ed:1)';
Ne = numel(es);
% calculate soc after charge vC = (v_t(e+P*eta))
eC = es + P2.*eta;
% round to the nearest sample
iC = ceil(eC/ed)+1;
iC(iC > (Ne+1)) = Ne + 2;
iC(iC < 2) = 1;
% calculate soc after discharge vC = (v_t(e-P/eta))
eD = es - P1./eta;
% round to the nearest sample
iD = floor(eD/ed)+1;
iD(iD > (Ne+1)) = Ne + 2;
iD(iD < 2) = 1;

price = repelem(lambda,T/numel(lambda));
for t = T:-1:1 %     % Recalculate ES value fcn
    vi = v(:,t+1);
    vo = CalcValueNoUnc(price(t), c, eta, vi, ed, iC, iD);
    v(:,t) = vo;
end
vp = v./eta + c;
vb = v.*eta;
%% modify here if need non-uniform segments
vpAvgt = zeros(seg_num,T+1);
vbAvgt = zeros(seg_num,T+1);
NN = (Ne-1)/seg_num;
for i = 1:seg_num
    vpAvgt(i,:) = mean(vp((i-1)*NN + (1:(NN+1)),:));
    vbAvgt(i,:) = mean(vb((i-1)*NN + (1:(NN+1)),:));
end

vpAvg = zeros(seg_num,H);
vbAvg = zeros(seg_num,H);
NH = T/H;
for j = 1:H
    vpAvg(:,j) = mean(vpAvgt(:,(j-1)*NH + (1:(NH+1))),2);
    vbAvg(:,j) = mean(vbAvgt(:,(j-1)*NH + (1:(NH+1))),2);
end
end