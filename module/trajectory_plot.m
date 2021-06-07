function [] = trajectory_plot(res,param, save_path,start)
N_GI = res.tGI;
N_GE = res.tGE;
tH_I = res.tHI;
V_E = res.VE;
V_I = res.VI;
H_E = res.HE;
H_I = res.HI;

mVE = min(min(V_E));
MVE = max(max(V_E));
mVI = min(min(V_I));
MVI = max(max(V_I));
mHE = min(min(H_E));
MHE = max(max(H_E));
mHI = min(min(H_I));
MHI = max(max(H_I));
gap_VE = max(floor((MVE-mVE)/30),0.5);
gap_VI = max(floor((MVI-mVI)/30),0.5);
gap_HE = max(floor((MHE-mHE)/30),0.5);
gap_HI = max(floor((MHI-mHI)/30),0.5);


hhh=figure;
times = 500;
set(gcf,'Position',[10,10,1200,1200]);
grid on;
hold on;
ShowSize = 30;
ind = find(res.time>=start);
start_t = ind(1);


ax1 = subplot(4,4,1);
ax2 = subplot(4,4,2);
ax3 = subplot(4,4,5);
ax4 = subplot(4,4,6);
ax5 = subplot(4,4,[3,4,7,8]);
    a=plot3(ax5,res.tGE, res.tGI, res.tHI,'b');
    a.Color(4)=0.1;
    xlabel('N_{GE}');
    ylabel('N_{GI}');
    zlabel('H_I');
    grid on;
    hold(ax5,'on');
Win = start_t:start_t+ShowSize-1;
a1 = plot3(ax5, N_GE(Win), N_GI(Win), H_I(Win), 'r','LineWidth',3);
hold(ax5,'on');
view(ax5, [-110,30]);

ax6=subplot(4,4,[9,10]);
rasterplot(res, param);
box(ax6,'on');
hold(ax6,'on');

ax7=subplot(4,4,[11,12]);
box on;
l1=plot(res.time*1000, res.tHE,'r');
hold(ax7,'on');
l2=plot(res.time*1000, res.tHI,'b');
ylabel('Count');
legend([l1,l2],'H_E','H_I');
hold(ax7,'on');

ax8=subplot(4,4,[13,14]);
box on;
l3=plot(res.time*1000, res.NGE/param.ne,'r');
hold on;
l4=plot(res.time*1000, res.NGI/param.ni,'b');
legend([l3,l4],'N_{GE}',  'N_{GI}');
hold(ax8,'on');

ax9=subplot(4,4,[15,16]);
box on;
l5=plot(res.time*1000, res.NRE/param.ne,'r');
hold on;
l6=plot(res.time*1000, res.NRI/param.ni,'b');
legend([l5,l6],'N_{RE}', 'N_{RI}');
hold(ax9,'on');

    t1 = 0;
    xlim(ax6, [t1-100,t1+100]);
    a2 = line(ax6, [t1,t1],[0,param.ne+param.ni]);
    set(get(get(a2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    hold(ax6,'on');
    xlim(ax7, [t1-100,t1+100]);
    a3 = line(ax7, [t1,t1],[0,4*10^6]);
    set(get(get(a3, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    hold(ax7,'on');
    xlim(ax8, [t1-100,t1+100]);
    a4 = line(ax8, [t1,t1],[0,1]);
set(get(get(a4, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    hold(ax8,'on');
    xlim(ax9, [t1-100,t1+100]);
    a5 = line(ax9, [t1,t1],[0,1]);
    set(get(get(a5, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    hold(ax9,'on'); 
writeObj = VideoWriter([save_path,'tr.avi']);
open(writeObj);
frame = getframe(hhh);
writeVideo(writeObj,frame.cdata);
for i = 1:times
    t = start_t+5*i;
    V_E_temp = V_E(t, :);
    V_I_temp = V_I(t, :);
    H_I_temp = H_I(t, :);
    H_E_temp = H_E(t, :);
    
    
    h1 = histogram(ax1,V_E_temp, 'Normalization','probability');
    h1.FaceColor = 'b';
    h1.BinEdges = [mVE:gap_VE:MVE];
    xlabel(ax1,'V_E');
    ylim(ax1,[0, 1]);
    
%     subplot(4,4,2);
    h1 = histogram(ax2,H_E_temp, 'Normalization','probability');
    h1.FaceColor = 'b';
    h1.BinEdges =  [mHE:gap_HE:MHE];
    xlabel(ax2,'H_E');
    ylim(ax2,[0, 1]);
    
%     subplot(4,4,5);
    h2 = histogram(ax3, V_I_temp, 'Normalization','probability');
    h2.FaceColor = 'r';
    h2.BinEdges =  [mVI:gap_VI:MVI];
    ylim(ax3,[0, 1]);
    xlabel(ax3,'V_I');
    
%     subplot(4,4,6);
    h2 = histogram(ax4, H_I_temp, 'Normalization','probability');
    h2.FaceColor = 'r';
    h2.BinEdges =  [mHI:gap_HI:MHI];
    ylim(ax4,[0, 1]);
    xlabel(ax4,'H_I');
    
%     ax = subplot(4,4,[3,4,7,8]);
% ax5 = subplot(4,4,[3,4,7,8]);
%     a=plot3(ax5,res.tGE, res.tGI, res.tHI,'b');
%     a.Color(4)=0.05;
%     xlabel('N_{GE}');
%     ylabel('N_{GI}');
%     zlabel('H_I');
%     grid on;
%     hold(ax5,'on');
    Win = max((t-ShowSize/2),1):(t+ShowSize/2);
    delete(a1);
    a1 = plot3(ax5, N_GE(Win), N_GI(Win), tH_I(Win), 'r','LineWidth',3);
    hold(ax5,'on');
    view(ax5,[-110,30]);
    
    delete(a2);
    delete(a3);
    delete(a4);
    delete(a5);
    t1 = res.time(t)*1000;
    xlim(ax6, [t1-100,t1+100]);
    a2 = line(ax6, [t1,t1],[0,param.ne+param.ni]);
    set(get(get(a2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    box(ax6,'on');
    hold(ax6,'on');
    
    xlim(ax7, [t1-100,t1+100]);
    a3 = line(ax7, [t1,t1],[0,4*10^6]);
    set(get(get(a3, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    hold(ax7,'on');
    xlim(ax8, [t1-100,t1+100]);
    a4 = line(ax8, [t1,t1],[0,1]);
set(get(get(a4, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    hold(ax8,'on');
    xlim(ax9, [t1-100,t1+100]);
    a5 = line(ax9, [t1,t1],[0,1]);
    set(get(get(a5, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    hold(ax9,'on');
    
    pause(0.01)
    sgtitle(['t= ', num2str(res.time(t)*1000),' ms']);
    pause(0.01);
    
    frame = getframe(hhh);
writeVideo(writeObj,frame.cdata);
    delete(h1);
    delete(h2);
end
close(writeObj);
end

