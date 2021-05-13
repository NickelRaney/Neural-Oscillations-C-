%% Load N_GE, N_GI, H_E, H_I

ne = 300;
ni = 100;
duration = 10000;
bar = 50;
model = 'B';
extra_name = '';
load_name = ['M=',model,'-n=', num2str(ne+ni),'-t=', num2str(duration/1000)]; 
model = ['model_', model];
if isempty(extra_name) == 0
    load_name = [load_name,'-', extra_name];
end
H_path = ['outputs//', model, '//', load_name,'//data//H.txt'];
V_path = ['outputs//', model, '//', load_name,'//data//V.txt'];

% Load H, V
H = load(H_path);
V = load(V_path);
H_E = H(:,2) + H(:,3);
H_I = H(:,4) + H(:,5);
V_E = V(:, 1:ne);
V_I = V(:, 1+ne: ni + ne);
N_GE = sum(V_E > bar, 2);
N_GI = sum(V_I > bar, 2);

% Discard first 100
V_E = V_E(100:end, :);
V_I = V_I(100:end, :);
H_I = H_I(100:end);
H_E = H_E(100:end);
N_GE = N_GE(100:end);
N_GI = N_GI(100:end);


%% 3D trajectory Anime
figure('Name','TrajIllus')
set(gcf,'Position',[10,10,1000,1000]);
a=plot3(N_GE, N_GI, H_E,'b');
a.Color(4)=0.03;
zlabel('H^E');
ylabel('N_{GI}');
xlabel('N_{GE}');
grid on;
view([-80, 20]);
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
a.Color(4)=0.08;
zlabel('H^I');
ylabel('N_{GI}');
xlabel('N_{GE}');
% xlim([150,300]);
% ylim([30,80]);
% zlim([2000,15000]);
grid on;
view([-80, 20]);
set(gca,'fontsize',15,'fontname','Arial');

hold on
ShowSize = 30;

for tInd = 100:length(N_GE) - ShowSize
    Win = tInd:tInd+ShowSize-1;
    a1 = plot3(N_GE(Win), N_GI(Win), H_I(Win), 'r','LineWidth',3);
    pause(0.02)
    delete(a1)
end

%%
figure;
times = 500;
mark = 1;
set(gcf,'Position',[10,10,1000,600]);
subplot(1,3,3);
a=plot3(N_GE, N_GI, H_I,'b');
a.Color(4)=0.03;
grid on;
hold on;
ShowSize = 30;
view([-110,  30]);
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
    subplot(2,3,[2,3,5,6]);
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
    sgtitle(['t=', num2str(t*0.1)]);

    pause(0.005);
    delete(h1);
    delete(h2);
end
%%
H_I_small = H_I(1:10:end);
i=0;
[m,p]=max(H_I_small(1:20));
pos=[p];
m_h=[m];
c_p=p;
while c_p+30<length(H_I_small)
    [m,p]=max(H_I_small(c_p+10:c_p+30));
    m_h=[m_h,m];
    c_p=c_p+10+p-1;
    pos=[pos,c_p];
end
%%
scatter(m_h(1:end-1),m_h(2:end))
%%
H_E_small = H_E(1:10:end);
i=0;
[m,p]=max(H_E_small(1:20));
pos=[p];
m_h=[m];
c_p=p;
while c_p+30<length(H_E_small)
    [m,p]=max(H_E_small(c_p+10:c_p+30));
    m_h=[m_h,m];
    c_p=c_p+10+p-1;
    pos=[pos,c_p];
end
%%
scatter(m_h(1:end-1),m_h(2:end))