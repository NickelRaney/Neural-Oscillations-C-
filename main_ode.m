%% Setting paths
addpath('module');
%% ode2
param.ne       = 300;
param.ni       = 100;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.435;
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
param.delta_time = 0.5;
param.dt = 0.01;
%%
tic;
res=model_ode2(param);
toc;

%% ode3
param.ne       = 300;
param.ni       = 100;
param.s_ee     = 5;
param.s_ie     = 2;
param.s_ei     = 4.91;
param.s_ii     = 4.91;
param.p_ee     = 0.15;
param.p_ie     = 0.5;
param.p_ei     = 0.5;
param.p_ii     = 0.4;

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
res=model_ode3(param);
toc;

%%
s_path='./';
ode_video(param,res,s_path)

%%
param.lambda_e = 7000;
param.lambda_i = 7000;
param.tau_re=0;
param.tau_ri=0;
param.ns_ee=1;
param.ns_ie=1;
param.ns_ei=1;
param.ns_ii=1;
param.factor_Leak=0;
param.LeakE = 20;
param.LeakI = 16.7;
param.factor_Leak = inf;
param.delta_time = 0.1;

tic;
res_mif=model_L(param);
toc;

%%
s_path='./';
compare_video(param,res,res_mif,s_path)
%%
range=[1:300];
subplot(2,2,1);
plot(sum(res_mif.HE(range,1:300),2))
hold on;
plot(res.h(range,1)*300);
xlabel('HEE');
legend('MIF','ode');

subplot(2,2,2);
plot(sum(res_mif.HE(range,301:400),2))
hold on;
plot(res.h(range,2)*100);
xlabel('HIE');

subplot(2,2,3);
plot(sum(res_mif.HI(range,1:300),2))
hold on;
plot(res.h(range,3)*300);
xlabel('HEI');

subplot(2,2,4);
plot(sum(res_mif.HI(range,301:400),2))
hold on;
plot(res.h(range,4)*100);
xlabel('HII');