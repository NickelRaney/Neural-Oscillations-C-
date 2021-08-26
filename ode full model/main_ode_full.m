%% Setting paths
addpath('ode full model');
%% model_L
alpha=1;
beta=1;
param.ne=300;
param.ni=100;
param.M        = 100*beta;
param.Mr       = 66*beta;
% param.p_ee     = 0.15;
% param.p_ie     = 0.5;
% param.p_ei     = 0.42;
% param.p_ii     = 0.4;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
% param.s_ee     = 5/alpha*beta;
% param.s_ie     = 2/alpha*beta;
% param.s_ei     = 4.91/alpha*beta;
% param.s_ii     = 4.91/alpha*beta;
param.s_ee     = 5*0.15/alpha*beta;
param.s_ie     = 2*0.5/alpha*beta;
param.s_ei     = 4.91*0.415/alpha*beta;
param.s_ii     = 4.91*0.4/alpha*beta;

param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.tau_re=0;
param.tau_ri=0;
param.duration = 1;
param.delta_time = 0.1;

param.factor_Leak=0;
param.LeakE = 20;

param.LeakI = 16.7;
param.factor_Leak = inf;
param.ns_ee=alpha;
param.ns_ie=alpha;
param.ns_ei=alpha;
param.ns_ii=alpha;

param2=param;
param2.gridsize=0.1;

%%
res_mif=model_L(param);
rasterplot(res_mif,param);

%%
tic;
res_ode=ode_full(param2);
toc;
rasterplot(res_ode,param2);
%%
set(gcf,'Position',[10,10,1000,300]);
name='pei=0.415-s=n-r=0';    
title(name);
saveas(gcf,['ode full model/output/Raster-',name,'.fig']);