%% Parameters

ne = 300;
ni = 100;
<<<<<<< Updated upstream
duration = 10000;
model = 'B';
extra_name = '';
=======
duration = 3000;
model = 'B';
extra_name = 'HitEI=2-SEI=2.45-PEI=1';
>>>>>>> Stashed changes
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
HE_path = ['outputs//', model, '//HE.txt'];
HI_path = ['outputs//', model, '//HI.txt'];
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
HE = load(HE_path);
HI = load(HI_path);
V = load(V_path);
time_rec = V(:,1);
V = V(:,2:end);
HE = HE(:,2:end);
HI = HI(:,2:end);
t_HE = sum(HE,2);
t_HI = sum(HI,2);
V_E = V(:, 1:param.ne);
V_I = V(:, 1+param.ne: param.ni + param.ne);
N_GE = sum(V_E > bar, 2);
N_GI = sum(V_I > bar, 2);
N_RE = sum(V_E==-67, 2);
N_RI = sum(V_I==-67, 2);

% Discard first 100
time_rec = time_rec*1000;
res.t_HI = t_HI;
res.t_HE = t_HE;
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
    copyfile(HE_path, save_path1);
    copyfile(HI_path, save_path1);
    copyfile(V_path, save_path1);
end

if save_bool
save([save_path1, save_name,'.mat'], 'res', 'param');
end

%% Rasterplot


fr = firing_rate(res, param);
SSI = spike_synchrony_index(res, param);

% Raster
subplot(4,1,1);
rasterplot(res, param);
title({['Raster-', save_name], ['FrE=', num2str(fr.e), ' FrI=', num2str(fr.i), ' SSI=', num2str(SSI)]});
box on;
if param.duration >= 3000
    xlim([param.duration-3000, param.duration]);
end

% HE, HI traj
subplot(4,1,2);
box on;
l1=plot(time_rec, t_HE);
l1.MarkerEdgeColor='r';
hold on;
l2=plot(time_rec, t_HI);
l2.MarkerEdgeColor = 'b';
ylabel('Count');
legend('HE','HI');
title('H-Trajectory');

% NG
subplot(4,1,3);
box on;
l1=plot(time_rec, N_GE/ne);
l1.MarkerEdgeColor='b';
hold on;
l2=plot(time_rec, N_GI/ni);
l2.MarkerEdgeColor = 'r';
legend('N_{GE}',  'N_{GI}');
title('NG-Trajectory');
ylabel('Percentage');

%NR
subplot(4,1,4);
box on;
l3=plot(time_rec, N_RE/ne);
l3.MarkerEdgeColor='b';
hold on;
l4=plot(time_rec, N_RI/ni);
l4.MarkerEdgeColor = 'r';
legend('N_{RE}', 'N_{RI}');
title('NR-Trajectory');
ylabel('Percentage');
set(gcf,'Position',[10,10,2000,1200]);
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

<<<<<<< Updated upstream
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
=======
% figure;
% tr =[N_GE, N_GI, H_I];
% plane = [1,0,0; 0,1,0];
% mark = 1;
% 
% gif_name = [save_path,'PM-NGE-NGI-HI-',save_name,'.gif'];
% 
% for i= 1:500
%     height = min(H_I) + i* (max(H_I)-min(H_I))/500;
%     [pm, pt] = poincare_map(tr,plane,height);
%     a=scatter(pm(:,1), pm(:,2),25, pt/max(pt),'filled');
%     xlim([0, max(N_GE)]);
%     ylim([0, max(N_GI)]);
%     xlabel('N_{GE}');
%     ylabel('N_{GI}');
%     colorbar
%     caxis([0 1]);
%     title(['H_I= ', num2str(height)]);
%     F = getframe(gcf);
%     im = frame2im(F);
%     [I, map] =rgb2ind(im, 256);
%     if mark == 1
%         imwrite(I,map, gif_name, 'GIF', 'Loopcount',inf,'DelayTime',0.01);
%         mark = mark + 1;
%     else
%         imwrite(I,map, gif_name, 'WriteMode','append','DelayTime',0.01);
%     end
%     pause(0.005);
%     delete(a);
% end
>>>>>>> Stashed changes

%% Generate GIF for V_E V_I

% figure;
% times = 500;
% mark = 1;
% set(gcf,'Position',[10,10,1600,600]);
% gif_name = [save_path,'Vd-',save_name,'.gif'];
% subplot(2,10,[4,5,9,10]);
% a=plot3(N_GE, N_GI, t_HI,'b');
% a.Color(4)=0.06;
% grid on;
% hold on;
% ShowSize = 30;
% view([-60,  30]);
% t= 1000;
% Win = t:t+ShowSize-1;
% a1 = plot3(N_GE(Win), N_GI(Win), t_HI(Win), 'r','LineWidth',3);
% max_HE = max(max(HE))+2;
% max_HI = max(max(HI))+2;
% for i = 1:times
%     t = size(V_E,1)-5*(times +10) + 5*i;
%     V_E_temp = V_E(t, :);
%     V_I_temp = V_I(t, :);
%     H_E_temp = HE(t,:);
%     H_I_temp = HI(t,:);
%     
%     %V_E distribution
%     subplot(2,5,1);
%     h1 = histogram(V_E_temp, 'Normalization','probability');
%     h1.FaceColor = 'b';
%     h1.BinEdges = [-70:5:100];
%     xlim([-70,100]);
%     xlabel('V_E');
%     ylim([0, 0.6]);
%     
%     %H_EE distribution
%     subplot(2,5,2);
%     h2 = histogram(H_E_temp(1:ne), 'Normalization','probability');
%     h2.FaceColor = 'b';
%     h2.BinEdges = [0:2:max_HE];
%     ylim([0, 1]);
%     xlabel('H_{EE}');
%     
%     %H_IE distribution
%     subplot(2,5,3);
%     h3 = histogram(H_E_temp(ne+1:ne+ni), 'Normalization','probability');
%     h3.FaceColor = 'b';
%     h3.BinEdges = [0:2:max_HE];
%     ylim([0, 1]);
%     xlabel('H_{IE}');
%     
%     %V_I distribution
%     subplot(2,5,6);
%     h4 = histogram(V_I_temp, 'Normalization','probability');
%     h4.FaceColor = 'r';
%     h4.BinEdges = [-70:5:100];
%     xlim([-70,100]);
%     ylim([0, 0.6]);
%     xlabel('V_I');
%     
%     %H_EI distribution
%     subplot(2,5,7);
%     h5 = histogram(H_I_temp(1:ne), 'Normalization','probability');
%     h5.FaceColor = 'r';
%     h5.BinEdges = [0:2:max_HI];
%     ylim([0, 1]);
%     xlabel('H_{EI}');
%     
%     %H_II distribution
%     subplot(2,5,8);
%     h6 = histogram(H_I_temp(ne+1:ne+ni), 'Normalization','probability');
%     h6.FaceColor = 'r';
%     h6.BinEdges = [0:2:max_HI];
%     ylim([0, 1]);
%     xlabel('H_{II}');
%     
%     
%     kkk= subplot(2,5,[4,5,9,10]);
%     kkk.Position = kkk.Position + [0.05 0 0 0.0];
%     Win = t:t+ShowSize-1;
%     a=plot3(N_GE, N_GI, t_HI,'b');
%     a.Color(4)=0.06;
%     xlabel('N_{GE}');
%     ylabel('N_{GI}');
%     zlabel('H_I');
%     grid on;
%     hold on;
%     view([-110,  30]);
%     delete(a1);
%     a1 = plot3(N_GE(Win), N_GI(Win), t_HI(Win), 'r','LineWidth',3);
%     hold on;
%     pause(0.01)
%     sgtitle({save_name,['t=', num2str(t*0.1)]});
%     F = getframe(gcf);
%     im = frame2im(F);
%     [I, map] =rgb2ind(im, 256);
%     if mark == 1
%         imwrite(I,map, gif_name,'GIF', 'Loopcount', inf,'DelayTime', 0.01);
%         mark = mark + 1;
%     else
%         imwrite(I,map, gif_name,'WriteMode','append','DelayTime', 0.01);
%     end
%     pause(0.005);
%     delete(h1);
%     delete(h2);
% end

%% Trajectory
subplot(1,2,1);
a=plot3(N_GE, N_GI, t_HI,'b');
a.Color(4)=0.06;
xlabel('N_{GE}');
ylabel('N_{GI}');
zlabel('H_I');
view([-110,  30]);
grid on;
subplot(1,2,2);
a=plot3(N_GE, N_GI, t_HE,'b');
a.Color(4)=0.06;
xlabel('N_{GE}');
ylabel('N_{GI}');
zlabel('H_E');
grid on;
view([-110,  30]);
set(gcf,'Position',[10,10,1400,600]);
sgtitle(['Trajectory-', save_name]);
if save_bool
    saveas(gcf,[save_path, 'Tr-', save_name,'.png']);
    saveas(gcf,[save_path, 'Tr-', save_name,'.fig']);
end