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

%% Example
% Setting Parameters
save_bool = true;
param1 = param;
param1.factor_Leak = inf;
param1.M = 1000;
param1.Mr = 660;
extra_name = 'p=1-M=1000';
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
param1.s_ee     = 7.5;
param1.s_ie     = 10;
param.lambda_e = 70000;
param.lambda_i = 70000;
param1.s_ei   = 25.0;
param1.s_ii   = 19.6;
bar = 500;
tic;
res = model_L(param1);
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
title({['Raster-'], ['FrE=', num2str(fr.e), ' FrI=', num2str(fr.i), ' SSI=', num2str(SSI)]});
box on;
if param1.duration >= 3000
    xlim([param1.duration-3000, param1.duration]);
end

% HE, HI traj
subplot(4,1,2);
box on;
l1=plot(res.time, res.tHE,'r');
hold on;
l2=plot(res.time, res.tHI,'b');
ylabel('Count');
legend('HE','HI');
title('H-Trajectory');

% NG
subplot(4,1,3);
box on;
l2=plot(res.time, res.NGE/param1.ne,'r');
hold on;
l2=plot(res.time, res.NGI/param1.ni,'b');
legend('N_{GE}',  'N_{GI}');
title('NG-Trajectory');
ylabel('Percentage');

%NR
subplot(4,1,4);
box on;
l3=plot(res.time, res.NRE/param1.ne,'r');
hold on;
l4=plot(res.time, res.NRI/param1.ni,'b');
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
    saveas(gcf,[save_path, 'Tr-', save_name,'.png']);
    saveas(gcf,[save_path, 'Tr-', save_name,'.fig']);
end

%%
for i=2:8
save_bool = true;
param1 = param;

param1.tau_re = 0.5*i;
% if param1.tau_ri ==0
%     param1.tau_ri=0.1;
% end
extra_name = ['tre=',num2str(param1.tau_re)];
bar = 50;
tic;
res = model_L(param1);
toc;

if param.factor_Leak == inf
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
else
    rmdir(plots_save_path,'s');
    mkdir(plots_save_path);
end
save([plots_save_path, 'data.mat'], 'res', 'param1');
end

%
fr = firing_rate(res, param1);
SSI = spike_synchrony_index(res, param1);
disp(['FrE=',num2str(fr.e),' FrI=',num2str(fr.i)]);
disp(['SSI=', num2str(SSI)]);

% Raster
f=figure;
subplot(4,1,1);
rasterplot(res, param1);
title({['Raster-', save_name], ['FrE=', num2str(fr.e), ' FrI=', num2str(fr.i), ' SSI=', num2str(SSI)]});
box on;

% HE, HI traj
subplot(4,1,2);
box on;
l1=plot(res.time, res.tHE);
l1.MarkerEdgeColor='r';
hold on;
l2=plot(res.time, res.tHI);
xlim([0,1000*param1.duration]);
l2.MarkerEdgeColor = 'b';
ylabel('Count');
legend('HE','HI');
title('H-Trajectory');

% NG
subplot(4,1,3);
box on;
l2=plot(res.time, res.NGE/param1.ne);
l2.MarkerEdgeColor='b';
hold on;
l2=plot(res.time, res.NGI/param1.ni);
l2.MarkerEdgeColor = 'r';
xlim([0,1000*param1.duration]);
legend('N_{GE}',  'N_{GI}');
title('NG-Trajectory');
ylabel('Percentage');

%NR
subplot(4,1,4);
box on;
l3=plot(res.time, res.NRE/param1.ne);
l3.MarkerEdgeColor='b';
hold on;
l4=plot(res.time, res.NRI/param1.ni);
l4.MarkerEdgeColor = 'r';
xlim([0,1000*param1.duration]);
legend('N_{RE}', 'N_{RI}');
title('NR-Trajectory');
ylabel('Percentage');
set(gcf,'Position',[10,10,2000,1200]);
if save_bool
    saveas(gcf,[plots_save_path, 'Raster-', save_name,'.png']);
    saveas(gcf,[plots_save_path, 'Raster-', save_name,'.fig']);
end
close(figure(gcf)) 

% Spectrogram
f=figure;
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
close(figure(gcf)) 
% Trajectory
f=figure;
subplot(1,2,1);
a=plot3(res.NGE, res.NGI, res.tHI,'b');
a.Color(4)=0.06;
xlabel('N_{GE}');
ylabel('N_{GI}');
zlabel('H_I');
view([-110,  30]);
grid on;
subplot(1,2,2);
a=plot3(res.NGE, res.NGI, res.tHE,'b');
a.Color(4)=0.06;
xlabel('N_{GE}');
ylabel('N_{GI}');
zlabel('H_E');
grid on;
view([-110,  30]);
set(gcf,'Position',[10,10,1400,600]);
sgtitle(['Trajectory-', save_name]);
if save_bool
    saveas(gcf,[plots_save_path, 'Tr-', save_name,'.png']);
    saveas(gcf,[plots_save_path, 'Tr-', save_name,'.fig']);
end
close(figure(gcf)) 
end