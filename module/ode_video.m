function []= ode_video(param, res, save_path)
hhh=figure;
npis = res.npi;
npes = res.npe;
peak_es = res.peak_e;
peak_is = res.peak_i;
set(gcf,'Position',[10,10,1200,1200]);
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);
writeObj = VideoWriter([save_path,'ode.avi']);
open(writeObj);
for i=1:size(npis)
    npe = npes(i);
    npi = npis(i);
    peak_e = peak_es(i,:);
    peak_i = peak_is(i,:);
    peak_e = reshape(peak_e,[3,10]);
    peak_i = reshape(peak_i,[3,10]);
    x_e = phinv(peak_e(2,1), peak_e(3,1));
    x_e_l = phinv(peak_e(2,1), 0.0013);
    xe1 = x_e_l:0.001:x_e;
    ye1 = exp(-xe1.^2./(2*peak_e(2,1)))./(sqrt(2*pi*peak_e(2,1)));
    xe1 = xe1 + peak_e(1,1);
    plot(ax1,xe1,ye1,'Color','r');
    if npe > 1
        for j=2:npe
            hold(ax1, 'on');
            m = peak_e(1,j);
            v = peak_e(2,j);
            r = peak_e(3,j);
            xe = -3*sqrt(v):0.01:3*sqrt(v);
            ye = exp(-xe.^2./(2*v))./(sqrt(2*pi*v));
            xe = xe + m;
            ye = ye * r;
            plot(ax1,xe,ye,'Color','r');
        end
    end
    r2 = peak_e(3,npe+1);
    if r2>0
        hold(ax1,'on');
        x_e_u = phinv(peak_e(2,1),peak_e(3,1)+r2);
        xe2 = x_e:0.01:x_e_u;
        ye2 = exp(-xe2.^2./(2*peak_e(2,1)))./(sqrt(2*pi*peak_e(2,1)));
        xe2 =xe2-x_e;
        plot(ax1,xe2,ye2,'Color','r');
        hold(ax1,'on');
        a1 = ye2(1);
        a2 = ye2(end);
        a1 = 0:0.0001:a1;
        a2 = 0:0.0001:a2;
        b1 = ones(1,size(a1,2));
        b2 = ones(1,size(a2,2));
        plot(ax1,0*b1,a1,'Color','r');
        hold(ax1,'on');
        plot(ax1,(x_e_u-x_e)*b2,a2,'Color','r');
    end
    xlim(ax1,[-66,100]);
    ylim(ax1,[0,0.2]);
    
    x_i = phinv(peak_i(2,1), peak_i(3,1));
    x_i_l = phinv(peak_i(2,1), 0.0013);
    xi1 = x_i_l:0.001:x_i;
    yi1 = exp(-xi1.^2./(2*peak_i(2,1)))./(sqrt(2*pi*peak_i(2,1)));
    xi1 = xi1 + peak_i(1,1);
    plot(ax2,xi1,yi1,'Color','b');
    
    if npi > 1
        for j=2:npi
            hold(ax2, 'on');
            m = peak_i(1,j);
            v = peak_i(2,j);
            r = peak_i(3,j);
            xi = -3*sqrt(v):0.01:3*sqrt(v);
            yi = exp(-xi.^2./(2*v))./(sqrt(2*pi*v));
            xi = xi + m;
            yi = yi * r;
            plot(ax2,xi,yi,'Color','b');
        end
    end
    r2 = peak_i(3,npi+1);
    if r2>0
        hold(ax2,'on');
        x_i_u = phinv(peak_i(2,1),peak_i(3,1)+r2);
        xi2 = x_i:0.01:x_i_u;
        yi2 = exp(-xi2.^2./(2*peak_i(2,1)))./(sqrt(2*pi*peak_i(2,1)));
        xi2 =xi2-x_i;
        plot(ax2,xi2,yi2,'Color','b');
        hold(ax1,'on');
        a1 = yi2(1);
        a2 = yi2(end);
        a1 = 0:0.0001:a1;
        a2 = 0:0.0001:a2;
        b1 = ones(1,size(a1,2));
        b2 = ones(1,size(a2,2));
        plot(ax2,0*b1,a1,'Color','b');
        hold(ax1,'on');
        plot(ax2,(x_i_u-x_i)*b2,a2,'Color','b');
    end
    
    xlim(ax2,[-66,100]);
    ylim(ax2,[0,0.2]);
    sgtitle(['t = ', num2str(i*param.delta_time),' ms']);
    xlabel(ax1,'V_E');
xlabel(ax2, 'V_I');
ylabel(ax1,'Percentage');
ylabel(ax2,'Percentage');
    frame = getframe(hhh);
    writeVideo(writeObj,frame.cdata);
    hold(ax2,'off');
    hold(ax1,'off');
end
close(writeObj);

    function x = phinv(v,h)
        % note: assume original distribution has zero mean so that there is
        % no m variable.
        h=min(h,0.9987);
        h=max(h,0.0013);
        x = erfinv((2*h-1))*sqrt(2)*sqrt(v);
    end
end

