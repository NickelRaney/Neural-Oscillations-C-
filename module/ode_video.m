function []= ode_video(param,  res, save_path)
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
    xe1 = x_e:0.001:x_e_l;
    ye1 = exp(xe1.^2./(2*peak_e(2,1)))./(sqrt(2*pi*peak_e(2,1)));
    xe1 = xe1 + peak_e(1,1);
    plot(ax1,xe1,ye1);
    xlim(ax1,[-66,100]);
    x_i = phinv(peak_i(2,1), peak_i(3,1));
    x_i_l = phinv(peak_i(2,1), 0.0013);
    xi1 = x_i:0.001:x_i_l;
    yi1 = exp(xi1.^2./(2*peak_i(2,1)))./(sqrt(2*pi*peak_i(2,1)));
    xi1 = xi1 + peak_i(1,1);
    plot(ax2,xi1,yi1);
    xlim(ax2,[-66,100]);
    sgtitle(['t=', num2str(i*param.delta_time)]);
    frame = getframe(hhh);
    writeVideo(writeObj,frame.cdata);
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

