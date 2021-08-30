%% Test_delta_mean
% This is a scrpit to measure the effect of MFE 
% on the differences of means based on the 
% Leakage Integrated Firing (LIF) model. 

%% Path Initialization
addpath('module');

%% Model Parameters Initilization
param.ne=300;
param.ni=100;
param.M        = 100;
param.Mr       = 66;
% param.p_ee     = 0.15;
% param.p_ie     = 0.5;
% param.p_ei     = 0.42;
% param.p_ii     = 0.4;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
% param.s_ee     = 5;
% param.s_ie     = 2;
% param.s_ei     = 4.91;
% param.s_ii     = 4.91;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.42;
param.s_ii     = 4.91*0.4;
param.lambda_e = 7000;
param.lambda_i = 7000;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.tau_re    = 0;
param.tau_ri    = 0;
param.duration = 1;
param.delta_time = 0.1;
param.gridsize=0.1;

%% Single Experiment to check the width of MFE
tic;
init.ve = zeros(1,300);
init.vi = zeros(1,100);
init.he = zeros(2,300);
init.hi = zeros(2,100);
res_lif=model_LIF2(param,init);
toc;
rasterplot(res_lif,param);

%% Test Parameters Initilization
single_parameter_test_time = 200;
lowerbound = -10;
upperbound = 25;
gridsize = 0.5;

%% Tests For Multiple Parameters
size_params = (upperbound - lowerbound)/gridsize +1;
params = lowerbound:gridsize:upperbound;
delta_m = zeros(size_params, single_parameter_test_time);
%res_m = cell(size_params, single_parameter_test_time);

for iter1 = 1:size_params
    diff = params(iter1);
    disp(['iteration ',num2str(iter1), ': Difference = ', num2str(diff)]);
    param.duration=0.05;
    t_start=tic;
    parfor trial=1:single_parameter_test_time
        he=zeros(2,300);
        hi=zeros(2,100);
        ve=sqrt(75)*randn(1,300)+diff;
        vi=sqrt(80)*randn(1,100);
        m=max(max(ve),max(vi))-100;
        ve=ve-m-1;
        vi=vi-m-1;
        res_lif=model_LIF2(param,ve,vi,he,hi);
        mfe_init=[];
        i=1;
        while res_lif.spikecount_e(i)==0
            i=i+1;
        end
        i=i+80;
        while res_lif.spikecount_e(i)==0
            i=i+1;
        end
        j=1;
        while res_lif.spikecount_i(j)==0
            j=j+1;
        end
        j=j+80;
        while res_lif.spikecount_i(j)==0
            j=j+1;
        end
        ind=min(i,j);
        delta_m(iter1,trial)=mean(res_lif.VE(ind,:))-mean(res_lif.VI(ind,:));
        %res_m{iter1,trial}=res_lif;
    end   
    t_end = toc(t_start);
    disp(['CPU time: ', num2str(t_end),'s']);
end

%% Plots
p_ei = 0.415;
path = 'figure/LIF/Delta_mean/delta_m-';
name = ['sd-p_ei=', num2str(p_ei)];
mean_m=mean(delta_m,2);
se_m=sqrt(var(delta_m')/single_parameter_test_time);
figure;
errorbar(params,mean_m,se_m);
xlabel('m_{VE}-m_{VI} before MFE');
ylabel('m_{VE}-m_{VI} after MFE');
set(gcf,'Position',[10,10,800, 600]);
title(['LIF model:', newline,' P_{EI}=', num2str(p_ei)]);

%%
delta_m1 = zeros(size_params, single_parameter_test_time);
for iter1 = 1:size_params
    for trial = 1:single_parameter_test_time
        delta_m1(iter1, trial) =  params(iter1);
    end
end
a = reshape(delta_m1, size_params*single_parameter_test_time,1);
b = reshape(delta_m, size_params*single_parameter_test_time,1);
scatter(a, b,'.');
xlabel('m_{VE}-m_{VI} before MFE');
ylabel('m_{VE}-m_{VI} after MFE');
title(['LIF model:', newline,' P_{EI}=', num2str(p_ei)]);

%%
saveas(gcf,[path,name,'.jpg']);
saveas(gcf,[path,name,'.fig']);