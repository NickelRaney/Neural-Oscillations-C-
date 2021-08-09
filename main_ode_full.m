%% Setting paths
addpath('module');
%% model_ODE
param.ne       = 300;
param.ni       = 100;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.415;
param.s_ii     = 4.91*0.4;

param.M        = 100;
param.Mr       = 66;
param.lambda_e = 7;
param.lambda_i = 7;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.duration = 1;
param.delta_time = 0.1;
param.dt = 0.01;
%%
tic;
res=model_ode2(param);
toc;

%% model_L
alpha=1;
beta=1;
param.ne=300;
param.ni=100;
param.M        = 100*beta;
param.Mr       = 66*beta;
% param.p_ee     = 0.15;
% param.p_ie     = 0.5;
% param.p_ei     = 0.415;
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
tic;
res_ode=ode_full(param2);
toc;

%%
s_path='./';
compare_video(param,res,res_ode,s_path)
%%
rasterplot(res_ode,param2);
%%
set(gcf,'Position',[10,10,1000,300]);
name='pei=0.415-s=n-r=0';    
title(name);
saveas(gcf,['ode full model/output/Raster-',name,'.fig']);