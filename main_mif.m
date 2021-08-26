%% Setting paths
addpath('module');

%%
param.ne       = 300;
param.ni       = 100;
param.p_ee     = 0.15;
param.p_ie     = 0.5;
param.p_ei     = 0.5;
param.p_ii     = 0.4;
param.s_ee     = 5;
param.s_ie     = 2;
param.s_ei     = 4.91;
param.s_ii     = 4.91;
param.ns_ee    = 1;
param.ns_ie    = 1;
param.ns_ei    = 1;
param.ns_ii    = 1;
param.tau_ri   = 2.5;
param.tau_re   = 2.5;
param.M        = 100;
param.Mr       = 66;
param.lambda_e = 7000;
param.lambda_i = 7000;
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

%% Example
% Setting Parameters
save_bool = false;
param1 = param;
param1.factor_Leak = inf;
beta = 10;
alpha = 10;
param1.M = 100*beta;
param1.Mr = 66*beta;
extra_name = 'p=1-ref=0-M=2000-ns=20-tee=1.4-tie=1.2-p_ei=4.2';
% param1.p_ee     = 0.15;
% param1.p_ie     = 0.5;
% param1.p_ei     = 0.4;
% param1.p_ii     = 0.4;
param1.p_ee     = 1;
param1.p_ie     = 1;
param1.p_ei     = 1;
param1.p_ii     = 1;
param1.tau_ee   = 1.4;
param1.tau_ie   = 1.2;
param1.tau_ei   = 4.5;
param1.tau_ii   = 4.5;
param1.s_ee     = 5*0.15/alpha * beta;
param1.s_ie     = 2*0.5/alpha * beta;
param1.s_ei     = 2.65/alpha* beta;
param1.s_ii     = 4.91*0.4/alpha* beta;
% param1.s_ee     = 5/alpha * beta;
% param1.s_ie     = 2/alpha * beta;
% param1.s_ei     = 4.91/alpha* beta;
% param1.s_ii     = 4.91/alpha* beta;
param1.ns_ee    = alpha;
param1.ns_ie    = alpha;
param1.ns_ei    = alpha;
param1.ns_ii    = alpha;
param1.tau_ri   = 2.5;
param1.tau_re   = 2.5;
param1.delta_time  = 0.1;
param1.duration  = 3;
param1.lambda_e = 7000*beta;
param1.lambda_i = 7000*beta;
bar = 50*beta;
tic;
res = model_LCI(param1);
toc;

%%
if param1.factor_Leak == inf
    model ='B';
else
    model = 'L';
end
save_name = ['M=',model,'-n=', num2str(param1.ne+ param1.ni),...
    '-t=', num2str(param1.duration)]; 
model = ['model_', model];
if isempty(extra_name) == 0
    save_name = [save_name,'-', extra_name];
end

save_path = ['outputs//', model, '//'];
plots_save_path = [save_path, save_name, '//'];

%
res.tHI = sum(res.HI,2);
res.tHE = sum(res.HE,2);
res.NGE = sum(res.VE>bar, 2);
res.NGI = sum(res.VI>bar, 2);
res.NRE = sum(res.VE==-param1.Mr-1,2);
res.NRI = sum(res.VI==-param1.Mr-1,2);

if save_bool
if exist(plots_save_path,'dir') == 0
    mkdir(plots_save_path);
end
save([plots_save_path, 'data.mat'], 'res', 'param1');
end

%%
fr = firing_rate(res, param1);
SSI = spike_synchrony_index(res, param1);

% Raster
subplot(4,1,1);
rasterplot(res, param1);
title({['Raster-', save_name], ['FrE=', num2str(fr.e), ' FrI=', num2str(fr.i), ' SSI=', num2str(SSI)]});
box on;
x1 = param1.duration*1000-3000;
x2 = param1.duration*1000;
if param1.duration >= 3
    xlim([x1, x2]);
end

% HE, HI traj
subplot(4,1,2);
box on;
l1=plot(res.time*1000, res.tHE,'r');
hold on;
l2=plot(res.time*1000, res.tHI,'b');
ylabel('Count');
legend('HE','HI');
if param1.duration >= 3
    xlim([x1, x2]);
end
title('H-Trajectory');

% NG
subplot(4,1,3);
box on;
l2=plot(res.time*1000, res.NGE/param1.ne,'r');
hold on;
l2=plot(res.time*1000, res.NGI/param1.ni,'b');
if param1.duration >= 3
    xlim([x1, x2]);
end
legend('N_{GE}',  'N_{GI}');
title('NG-Trajectory');
ylabel('Percentage');

%NR
subplot(4,1,4);
box on;
l3=plot(res.time*1000, res.NRE/param1.ne,'r');
hold on;
l4=plot(res.time*1000, res.NRI/param1.ni,'b');
if param1.duration >= 3
    xlim([x1, x2]);
end
legend('N_{RE}', 'N_{RI}');
title('NR-Trajectory');
ylabel('Percentage');
set(gcf,'Position',[10,10,2000,1200]);
if save_bool
    saveas(gcf,[plots_save_path, 'Raster-', save_name,'.png']);
    saveas(gcf,[plots_save_path, 'Raster-', save_name,'.fig']);
end

%% Spectrogram
sd = spikedensity(res, param1);
subplot(2,1,1);
spectrogram(sd.e, param1);
title('Spec-E');
subplot(2,1,2);
spectrogram(sd.i, param1);
title('Spec-I');
set(gcf,'Position',[10,10,2000,600]);
sgtitle(['Spec-', save_name]);
if save_bool
    saveas(gcf,[plots_save_path, 'Spec-', save_name,'.png']);
    saveas(gcf,[plots_save_path, 'Spec-', save_name,'.fig']);
end
%%
% Trajectory
figure;
subplot(1,2,1);
a=plot3(res.NGE, res.NGI, res.tHI,'b');
a.Color(4)=0.06;
xlabel('N_{GE}');
ylabel('N_{GI}');
zlabel('H_I');
view([-110,  30]);
grid on;
set(gca,'fontsize',15,'fontname','Arial');
subplot(1,2,2);
a=plot3(res.NGE, res.NGI, res.tHE,'b');
a.Color(4)=0.06;
xlabel('N_{GE}');
ylabel('N_{GI}');
zlabel('H_E');
grid on;
view([-110,  30]);
set(gcf,'Position',[10,10,1400,600]);
%sgtitle(['Trajectory-', save_name]);
set(gca,'fontsize',15,'fontname','Arial');
if save_bool
    saveas(gcf,[plots_save_path, 'Tr-', save_name,'.png']);
    saveas(gcf,[plots_save_path, 'Tr-', save_name,'.fig']);
end
%%
res.tHE = sum(res.HE, 2);
res.tHI = sum(res.HI, 2);
res.tGE = sum(res.VE>bar,2);
res.tGI = sum(res.VI>bar,2);
trajectory_plot(res,param1,plots_save_path, 0.3);