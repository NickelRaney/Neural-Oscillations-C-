%% Setting paths
addpath(genpath(pwd));

%%
param.ne       = 75;
param.ni       = 25;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.5;
param.s_ii     = 4.91*0.4;
param.tau_ri   = 0;
param.tau_re   = 0;
param.tau_ee   = 1.3;
param.tau_ie   = 0.95;
param.tau_ei   = 4.5;
param.tau_ii   = 4.5;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
param.M        = 100;
param.Mr       = 66;
param.lambda_e = 7000;
param.lambda_i = 7000;
param.duration = 50;
param.LeakE    = 20;
param.LeakI    = 16.7;
param.factor_Leak = inf;
param.delta_time  = 0.1;
param.sdbin       = 2.5;
param.spectrogram_timewindow = 200;
param.w               = 1;
param.frequency_range = [5,100];
param.bar             = 45;

%%
tic;
res=model_L(param);
toc;


%%
param1 = param;
param1.s_ei     = 4.4*0.5;
param1.duration = 3;
tic;
res1=model_L(param1);
toc;
figure;
rasterplot(res1,param1);
set(gcf,'Position',[50,50,1200,200])
%%
res1.tHE = sum(res1.HE, 2);
res1.tHI = sum(res1.HI, 2);
res1.tGE = sum(res1.VE>50,2);
res1.tGI = sum(res1.VI>50,2);
trajectory_plot(res1,2);

%%
bar.low_e            = 20;
bar.high_e           = 70;
bar.low_i            = 20;
bar.high_i           = 70;
P = P_generation_3state_lowbase_statistics(res,param,bar);

%%
param2 = param;
param2.duration = 3;
tic
res2 = model_LRNL(param2,P);
toc

figure;
rasterplot(res2, param2);
set(gcf,'Position',[50,50,1200,200]);

%%
res2.tHE = sum(res2.HE, 2);
res2.tHI = sum(res2.HI, 2);
res2.tGE = sum(res2.VE>-1,2);
res2.tGI = sum(res2.VI>-1,2);

trajectory_plot(res2,1);
