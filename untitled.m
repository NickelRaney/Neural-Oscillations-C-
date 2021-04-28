
%% Generate GIF for V_E V_I

times = 500;
mark = 1;
set(gcf,'Position',[10,10,1000,500]);
gif_name = [save_path,'Vd-',save_name,'.gif'];
subplot(2,3,[2,3,5,6]);
view([-80,30]);
a=plot3(N_GE, N_GI, H_E,'b');
a.Color(4)=0.03;
grid on;
hold on;
ShowSize = 30;
t = 1000;
Win = t:t+ShowSize-1;
a1 = plot3(N_GE(Win), N_GI(Win), H_I(Win), 'r','LineWidth',3);
for i = 1:times
    t = size(V_E,1)-1000 + i;
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
    a=plot3(N_GE, N_GI, H_E,'b');
    a.Color(4)=0.03;
    hold on;
    grid on;
    Win = t:t+ShowSize-1;
    delete(a1);
    a1 = plot3(N_GE(Win), N_GI(Win), H_I(Win), 'r','LineWidth',3);
    hold on;
    pause(0.02)
    sgtitle({save_name,['t=', num2str(t*0.1)]});
%     F = getframe(gcf);
%     im = frame2im(F);
%     [I, map] =rgb2ind(im, 256);
%     if mark == 1
%         imwrite(I,map, gif_name,'GIF', 'Loopcount',inf,'DelayTime',0.01);
%         mark = mark + 1;
%     else
%         imwrite(I,map, gif_name,'WriteMode','append','DelayTime',0.01);
%     end
    pause(0.005);
    delete(h1);
    delete(h2);
end