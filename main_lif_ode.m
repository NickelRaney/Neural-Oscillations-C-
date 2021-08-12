%% Setting path

addpath('module');

%% parameters for Model_ODE

param.ne       = 300;
param.ni       = 100;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.5;
param.s_ii     = 4.91*0.4;
param.M        = 100;
param.Mr       = 66;
param.lambda_e = 7;
param.lambda_i = 7;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei   = 4.5;
param.tau_ii   = 4.5;
param.duration = 1;
param.delta_time = 5;
param.dt = 0.01;

%% Model_ode

tic;
res_ode=model_ode5(param);
toc;

%% Parameters for Model_lif

alpha = 1;
beta = 1;
param.ne = 300;
param.ni = 100;
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
param.s_ei     = 4.91*0.5/alpha*beta;
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


%% Model_lif

tic;
res_lif=model_LIF(param2);
toc;

%% Compare distributions versus time

s_path='./';
compare_video(param,res_ode,res_lif,s_path);

%% Rasterplot of Model_lif

rasterplot(res_lif,param2);

%% Save Rasterplot

set(gcf,'Position',[10,10,1000,300]);
name='pei=0.415-s=n-r=0';    
title(name);
saveas(gcf,['ode full model/output/Raster-',name,'.fig']);

%% SD plots
delta_time = param.delta_time/1000;
size = floor(param2.duration/delta_time);
sd_lif = zeros(size,2);
spike_E = res_lif.spike(2:end,1:300);
spike_I = res_lif.spike(2:end,301:400);

for i = 1:size
    sd_lif(i,1) = sum(sum(spike_E< i*delta_time & spike_E>= (i-1)*delta_time))/delta_time;
    sd_lif(i,2) = sum(sum(spike_I< i*delta_time & spike_I>= (i-1)*delta_time))/delta_time;
end

t1 = param.delta_time:param.delta_time:param.duration*1000;
t2 = param.delta_time:param.delta_time: size*param.delta_time;
max_fr1 = max(max(res_ode.sd));
max_fr2 = max(max(sd_lif(2:end,2)));
max_fr = max([max_fr1 max_fr2]) *1.2;

figure;
subplot(1,2,1);
plot(t1, res_ode.sd(:,1));
hold on;
plot(t2(2:end), sd_lif(2:end,1));
ylim([0, max_fr]);
xlabel('Time (ms)');
ylabel('Firing Rate');
subplot(1,2,2);
plot(t1, res_ode.sd(:,2));
hold on;
plot(t2(2:end), sd_lif(2:end,2));
ylim([0, max_fr]);
xlabel('Time (ms)');
legend('ODE', 'LIF');

sgtitle('Spike Density');
set(gcf,'Position',[10,10,1000,400]);

%% LIF experiment
for i =1:101
    tau = 0.41 + (i-1)*0.001;
    param2.s_ei     = 4.91*tau/alpha*beta;
    tic;
    res_lif=model_LIF(param2);
    toc;
    figure;
    rasterplot(res_lif,param2);
    set(gcf,'Position',[10,10,1000,300]);
    name=['p_ei=',num2str(tau)];    
    title(name);
    saveas(gcf,['figure/LIF/Raster-',name,'.jpg']);
    saveas(gcf,['figure/LIF/Raster-',name,'.fig']);
    close;
end