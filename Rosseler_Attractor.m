%% Parameters

ne = 300;
ni = 100;
duration = 10000;
model = 'B';
extra_name = '';
save_bool = true;
bar = 50;

%% Run Model

cmd_commend = ['.\\models\\NOModel_', model, '.exe'];
tic;
system(cmd_commend);
toc;

%% Data Loading

save_name = ['M=',model,'-n=', num2str(ne+ni),'-t=', num2str(duration/1000)]; 
model = ['model_', model];
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
V_E = V_E(100:end, :);
V_I = V_I(100:end, :);
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

figure;
tr =[N_GE, N_GI, H_I];
plane = [1,0,0; 0,1,0];
mark = 1;

gif_name = [save_path,'PM-NGE-NGI-HI-',save_name,'.gif'];

for i= 1:500
    height = min(H_I) + i* (max(H_I)-min(H_I))/500;
    [pm, pt] = poincare_map(tr,plane,height);
    a=scatter(pm(:,1), pm(:,2),25, pt/max(pt),'filled');
    xlim([0, max(N_GE)]);
    ylim([0, max(N_GI)]);
    xlabel('N_{GE}');
    ylabel('N_{GI}');
    colorbar
    caxis([0 1]);
    title(['H_I= ', num2str(height)]);
    F = getframe(gcf);
    im = frame2im(F);
    [I, map] =rgb2ind(im, 256);
    if mark == 1
        imwrite(I,map, gif_name, 'GIF', 'Loopcount',inf,'DelayTime',0.01);
        mark = mark + 1;
    else
        imwrite(I,map, gif_name, 'WriteMode','append','DelayTime',0.01);
    end
    pause(0.005);
    delete(a);
end
%% Generate another version of Poincare_map

figure;
tr =[N_GE, N_GI, H_I];
plane = [1,0,0; 0,1,0];
mark = 1;

gif_name = [save_path,'PM2-NGE-NGI-HI-',save_name,'.gif'];

for i= 30:500
    height = min(H_I) + i* (max(H_I)-min(H_I))/500;
    pm2 = poincare_map2(tr,plane,height);
    s=scatter(pm2(:,1), pm2(:,2), '.');
    a=cell(1);
    for j = 1:10:size(pm2,1)-1
    a{j}=PlotLineArrow(gca, [pm2(j,1), pm2(j+1,1)], [pm2(j,2), pm2(j+1,2)], 'b', 'r');
    end
    xlabel('N_{GE}');
    ylabel('N_{GI}');
    title(['H_I= ', num2str(height)]);
    F = getframe(gcf);
    im = frame2im(F);
    [I, map] =rgb2ind(im, 256);
    if mark == 1
        imwrite(I,map, gif_name, 'GIF', 'Loopcount',inf,'DelayTime',0.01);
        mark = mark + 1;
    else
        imwrite(I,map, gif_name, 'WriteMode','append','DelayTime',0.01);
    end
    pause(0.005);
    delete(s);  
    for j=1:size(a,2)
        delete(a{j});
    end
end

%% Generate mean Poincare_map
figure;
tr =[N_GE, N_GI, H_I];
plane = [1,0,0; 0,1,0];
mark = 1;

gif_name = [save_path,'PM2-mean-NGE-NGI-HI-',save_name,'.gif'];
pm3 = [];
for i= 30:500
    height = min(H_I) + i* (max(H_I)-min(H_I))/500;
    pm2 = poincare_map2(tr,plane,height);
    X=pm2(1:size(pm2,1)-1,:)';
    K=9;
    [D,N] = size(X);
    y=zeros(1,N);
    X2 = sum(X.^2,1);
    distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
    size(distance)
    [sorted,index] = sort(distance);
    neighborhood = index(2:(1+K),:);
    neighborhood = [1:1:size(pm2,1)-1;neighborhood];
    for j=1:size(pm2,1)-1
        pm3(j,1)=mean(pm2(neighborhood(:,j),1));
        pm3(j,2)=mean(pm2(neighborhood(:,j),2));
        pm3(j,3)=mean(pm2(neighborhood(:,j)+1,1));
        pm3(j,4)=mean(pm2(neighborhood(:,j)+1,2));
    end
    s=scatter(pm2(:,1), pm2(:,2), '.');
    a=cell(1);
    for j = 1:5:size(pm3,1)-1
        a{j}=PlotLineArrow(gca, [pm3(j,1), pm3(j,3)], [pm3(j,2), pm3(j,4)], 'b', 'r');
    end
    xlabel('N_{GE}');
    ylabel('N_{GI}');
    title(['H_I= ', num2str(height)]);
    F = getframe(gcf);
    im = frame2im(F);
    [I, map] =rgb2ind(im, 256);
    if mark == 1
        imwrite(I,map, gif_name, 'GIF', 'Loopcount',inf,'DelayTime',0.01);
        mark = mark + 1;
    else
        imwrite(I,map, gif_name, 'WriteMode','append','DelayTime',0.01);
    end
    pause(0.005);
    delete(s);  
    for j=1:size(a,2)
        delete(a{j});
    end
end

%% Generate GIF for V_E V_I

figure;
times = 500;
mark = 1;
set(gcf,'Position',[10,10,1000,600]);
gif_name = [save_path,'Vd-',save_name,'.gif'];
subplot(1,3,3);
a=plot3(N_GE, N_GI, H_I,'b');
a.Color(4)=0.03;
grid on;
hold on;
ShowSize = 30;
view([-60,  30]);
t= 1000;
Win = t:t+ShowSize-1;
a1 = plot3(N_GE(Win), N_GI(Win), H_I(Win), 'r','LineWidth',3);
for i = 1:times
    t = size(V_E,1)-5*(times +10) + 5*i;
    V_E_temp = V_E(t, :);
    V_I_temp = V_I(t, :);
    subplot(2,3,1);
    h1 = histogram(V_E_temp, 'Normalization','probability');
    h1.FaceColor = 'b';
    h1.BinEdges = [-70:5:100];
    xlim([-70,100]);
    xlabel('V_E');
    ylim([0, 0.6]);
    subplot(2,3,4);
    h2 = histogram(V_I_temp, 'Normalization','probability');
    h2.FaceColor = 'r';
    h2.BinEdges = [-70:5:100];
    xlim([-70,100]);
    ylim([0, 0.6]);
    xlabel('V_I');
    kkk= subplot(2,3,[2,3,5,6]);
    kkk.Position = kkk.Position + [0.05 0 0 0.0];
    Win = t:t+ShowSize-1;
    a=plot3(N_GE, N_GI, H_I,'b');
    a.Color(4)=0.03;
    xlabel('N_{GE}');
    ylabel('N_{GI}');
    zlabel('H_I');
    grid on;
    hold on;
    view([-110,  30]);
    delete(a1);
    a1 = plot3(N_GE(Win), N_GI(Win), H_I(Win), 'r','LineWidth',3);
    hold on;
    pause(0.02)
    sgtitle({save_name,['t=', num2str(t*0.1)]});
    F = getframe(gcf);
    im = frame2im(F);
    [I, map] =rgb2ind(im, 256);
    if mark == 1
        imwrite(I,map, gif_name,'GIF', 'Loopcount', inf,'DelayTime', 0.01);
        mark = mark + 1;
    else
        imwrite(I,map, gif_name,'WriteMode','append','DelayTime', 0.01);
    end
    pause(0.005);
    delete(h1);
    delete(h2);
end

