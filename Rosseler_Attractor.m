%% Full Model
tic;
system('.\\models\\NOModel_full.exe');
toc;

%% L Model
tic;
system('.\\models\\NOModel_L_full.exe');
toc;

%% Parameters
ne = 300;
ni = 100;
duration = 10000;
model = 'L';
extra_name = 'no-in-noise';
save_bool = true;
bar = 50;

%% Data Loading

if isempty(model)
save_name = ['M=','Full','-n=', num2str(ne+ni),'-t=', num2str(duration/1000)];
model = 'model_full';
else
  save_name = ['M=',model,'-n=', num2str(ne+ni),'-t=', num2str(duration/1000)]; 
  model = ['model_', model, '_full'];
end

if isempty(extra_name) == 0
    save_name = [save_name,'-', extra_name];
end

spike_path = ['outputs//', model, '//spike.txt'];
H_path = ['outputs//', model, '//H.txt'];
V_path = ['outputs//', model, '//V.txt'];
save_path = ['outputs//', model, '//'];
param_file = [model,'_params.txt'];

save_path = [save_path, save_name, '//'];
save_path1 = [save_path, 'data//'];

% Generate res
spike_som = load(spike_path);
spike = zeros(10000,ne+ni);
for i=1:size(spike_som,1)
    neuron_i = spike_som(i,2);
    spike(1,neuron_i+1) = spike(1,neuron_i+1) + 1;
    spike(spike(1,neuron_i+1)+1,neuron_i+1) = spike_som(i,1)*1000;
end
res.spike = spike;

% Generate param
param.ne = ne;
param.ni = ni;
param.duration = duration;
param.sdbin = 2.5;
param.spectrogram_timewindow = 200;
param.w = 1;
param.frequency_range = [5,100];

% Load H, V
H = load(H_path);
V = load(V_path);
H_E = H(:,2) + H(:,3);
H_I = H(:,4) + H(:,5);
V_E = V(:, 1:param.ne);
V_I = V(:, 1+param.ne: param.ni + param.ne);
N_GE = sum(V_E > bar, 2);
N_GI = sum(V_I > bar, 2);

% Discard first 100
H_I = H_I(100:end);
H_E = H_E(100:end);
N_GE = N_GE(100:end);
N_GI = N_GI(100:end);

res.H_I = H_I;
res.H_E = H_E;
res.N_GE = N_GE;
res.N_GI = N_GI;

%% Data Saving
if save_bool
if exist(save_path,'dir') == 0
    mkdir(save_path);
end
if exist(save_path1,'dir') == 0
    mkdir(save_path1);
end
    copyfile(param_file,save_path1);
    copyfile(spike_path, save_path1);
    copyfile(H_path, save_path1);
    copyfile(V_path, save_path1);
end

if save_bool
save([save_path1, save_name,'.mat'], 'res', 'param');
end

%% Rasterplot
rasterplot(res, param);
fr = firing_rate(res, param);
SSI = spike_synchrony_index(res, param);
title({['Raster-', save_name], ['FrE=', num2str(fr.e), ' FrI=', num2str(fr.i), ' SSI=', num2str(SSI)]});
box on;
if param.duration >= 3000
    xlim([param.duration-3000, param.duration]);
end
set(gcf,'Position',[10,10,2000,300]);
if save_bool
    saveas(gcf,[save_path, 'Raster-', save_name,'.png']);
    saveas(gcf,[save_path, 'Raster-', save_name,'.fig']);
end

%% Spectrogram
sd = spikedensity(res, param);
subplot(2,1,1);
spectrogram(sd.e, param);
title('Spec-E');
subplot(2,1,2);
spectrogram(sd.i, param);
title('Spec-I');
set(gcf,'Position',[10,10,2000,600]);
sgtitle(['Spec-', save_name]);
if save_bool
    saveas(gcf,[save_path, 'Spec-', save_name,'.png']);
    saveas(gcf,[save_path, 'Spec-', save_name,'.fig']);
end

%% Generaete GIF for Poincare_map

tr =[N_GE, N_GI, H_I];
plane = [1,0,0; 0,1,0];
mark = 1;
gif_name = [save_path,'PM-NGE-NGI-HI-',save_name,'.gif'];
for i= 1:500
    height = min(H_I) + i* (max(H_I)-min(H_I))/500;
    pm = poincare_map(tr,plane,height);
    title(['H_I= ', num2str(height)]);
    a=patch(pm(:,1), pm(:,2),pm(:,2)/max(pm(:,2)),'edgecolor','flat','facecolor','none');
    xlim([0, max(N_GE)]);
    ylim([0, max(N_GI)]);
    xlabel('N_{GE}');
    ylabel('N_{GI}');
    F = getframe(gcf);
    im = frame2im(F);
    [I, map] =rgb2ind(im, 256);
    if mark == 1
        imwrite(I,map, gif_name, 'GIF', 'Loopcount',inf,'DelayTime',0.01);
        mark = mark + 1;
    else
        imwrite(I,map, gif_name, 'WriteMode','append','DelayTime',0.01);
    end
    colorbar
    caxis([0 1]);
    pause(0.005);
    delete(a);
end

%% Generate GIF for V_E V_I

times = 500;
mark = 1;
set(gcf,'Position',[10,10,1000,500]);
gif_name = [save_path,'Vd-',save_name,'.gif'];
for i = 1:times
    t = size(V_E,1)-500 + i;
    V_E_temp = V_E(t, :);
    V_I_temp = V_I(t, :);
    subplot(1,2,1);
    h1 = histogram(V_E_temp, 'Normalization','probability');
    h1.FaceColor = 'b';
    h1.BinEdges = [-70:5:100];
    xlim([-70,100]);
    xlabel('V_E');
    ylim([0, 0.6]);
    subplot(1,2,2);
    h2 = histogram(V_I_temp, 'Normalization','probability');
    h2.FaceColor = 'r';
    h2.BinEdges = [-70:5:100];
    xlim([-70,100]);
    ylim([0, 0.6]);
    xlabel('V_I');
    sgtitle({save_name,['t=', num2str(t*0.1)]});
    F = getframe(gcf);
    im = frame2im(F);
    [I, map] =rgb2ind(im, 256);
    if mark == 1
        imwrite(I,map, gif_name,'GIF', 'Loopcount',inf,'DelayTime',0.01);
        mark = mark + 1;
    else
        imwrite(I,map, gif_name,'WriteMode','append','DelayTime',0.01);
    end
    pause(0.005);
    delete(h1);
    delete(h2);
end

%% 3D trajectory Anime
figure('Name','TrajIllus')
set(gcf,'Position',[10,10,500,500]);
a=plot3(N_GE, N_GI, H_E,'b');
a.Color(4)=0.03;
zlabel('H^E');
ylabel('N_{GI}');
xlabel('N_{GE}');
grid on;
view([-80, 30]);
set(gca,'fontsize',15,'fontname','Arial');
hold on
ShowSize = 30;
for tInd = 100:length(N_GE) - ShowSize
    Win = tInd:tInd+ShowSize-1;
    a1 = plot3(N_GE(Win), N_GI(Win), H_E(Win), 'r','LineWidth',3);
    pause(0.02);
    delete(a1);
end

%%
figure('Name','TrajIllus')
set(gcf,'Position',[10,10,1000,1000]);
a=plot3(N_GE, N_GI, H_I,'b');
a.Color(4)=0.03;
zlabel('H^I');
ylabel('N_{GI}');
xlabel('N_{GE}');
grid on;
view([-40, 30]);
set(gca,'fontsize',15,'fontname','Arial');

hold on
ShowSize = 30;

for tInd = 100:length(N_GE) - ShowSize
    Win = tInd:tInd+ShowSize-1;
    a1 = plot3(N_GE(Win), N_GI(Win), H_I(Win), 'r','LineWidth',3);
    pause(0.02)
    delete(a1)
end

