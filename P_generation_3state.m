%% Setting paths
addpath(genpath(pwd));

%%
bar.e            = 50;
bar.i            = 50;
param.ne         = 75;
param.ni         = 25;
param.p_ee       = 0.15;
param.p_ie       = 0.5;
param.p_ei       = 0.5;
param.p_ii       = 0.4;
param.s_ee       = 20;
param.s_ie       = 8;
param.s_ei       = 20;
param.s_ii       = 20;
param.M          = 100;
param.Mr         = 66;
param.lambda_e   = 1/7;
param.lambda_i   = 1/7;
%%
V=load('C++/outputs/model_full/V_small.txt');

res.V_e=V(:,1:param.ne);
res.V_i=V(:,param.ne+1:param.ne+param.ni);
%%
P = P_generation_3state_statistics(res,param,bar);