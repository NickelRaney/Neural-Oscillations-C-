%% Path Initialization

addpath('module');

%% Model Parameters Initilization

param.ne       = 300;
param.ni       = 100;
param.M        = 100;
param.Mr       = 66;
param.p_ee     = 0.15;
param.p_ie     = 0.5;
param.p_ei     = 0.5;
param.p_ii     = 0.4;
% param.p_ee     = 1;
% param.p_ie     = 1;
% param.p_ei     = 1;
% param.p_ii     = 1;
param.s_ee     = 5;
param.s_ie     = 2;
param.s_ei     = 4.91;
param.s_ii     = 4.91;
% param.s_ee     = 5*0.15;
% param.s_ie     = 2*0.5;
% param.s_ei     = 4.91*0.42;
% param.s_ii     = 4.91*0.4;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.5;
param.s_ii     = 4.91*0.4;
param.lambda_e = 7000;
param.lambda_i = 7000;
param.tau_ee   = 1.3;
param.tau_ie   = 0.95;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.tau_re    = 0.5;
param.tau_ri    = 2.5;
param.duration = 3;
param.delta_time = 0.1;
param.gridsize=0.1;

%%

ve = zeros(1,300);
vi = zeros(1,100);
he = zeros(2,300);
hi = zeros(2,100);
tic;
res_lif=model_LIF3(param,ve,vi,he,hi);
toc;

%%
alpha          = 20;
beta           = 20;
param.ne       = 300;
param.ni       = 100;
param.p_ee     = 0.15;
param.p_ie     = 0.5;
param.p_ei     = 0.5;
param.p_ii     = 0.4;
param.s_ee     = 5/alpha*beta;
param.s_ie     = 2/alpha*beta;
param.s_ei     = 4.91/alpha*beta;
param.s_ii     = 4.91/alpha*beta;
param.ns_ee    = alpha;
param.ns_ie    = alpha;
param.ns_ei    = alpha;
param.ns_ii    = alpha;
param.tau_ri   = 2.5;
param.tau_re   = 2.5;
param.M        = 100*beta;
param.Mr       = 66*beta;
param.lambda_e = 7000*beta;
param.lambda_i = 7000*beta;
param.tau_ee   = 1.3;
param.tau_ie   = 0.95;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.duration = 3;
param.LeakE = 20;
param.LeakI = 16.7;
param.factor_Leak = inf;
param.delta_time = 0.1;
param.sdbin = 2.5;
param.spectrogram_timewindow = 200;
param.w = 1;
param.frequency_range = [5,100];
param.duration = 3;

%%

tic;
res_mif=model_L(param);
toc;

%% 

figure;
subplot(2,1,1);
rasterplot(res_lif, param);
title('LIF');
hold on;
subplot(2,1,2);
rasterplot(res_mif, param);
title('MIF');
