%ode3
param.ne       = 300;
param.ni       = 100;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.4;
param.s_ii     = 4.91*0.4;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;

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

tic;
res=model_ode3(param);
toc;
