function [spectrogram]=spectrogram(sd,param)

%unit of timewindow and bin is milisecond while unit of frequency is times/second.
tw=param.spectrogram_timewindow;
bin=param.sdbin;
fr=param.frequency_range;
duration = param.duration;
initial_index = ceil((2.5)/ bin);
N=tw/bin;
kernel=(exp(-2*pi*i*bin/1000)*ones(1,N)).^(1:N)';
num_spec= length(sd) - N+1 - initial_index+1;
spectrogram=zeros(fr(2)-fr(1)+1,num_spec);
for j=1:num_spec
    for k=1:fr(2)-fr(1)+1
        spectrogram(k,j)=abs(bin/1000 * sd(initial_index+j -1: initial_index+j+N-2)*(kernel.^(k+fr(1)-1))...
            /sqrt(tw/1000))^2;
    end
end

%smoothing
norm(spectrogram,'fro')
spectrogram1=conv2(spectrogram,ones(1,8)/8, 'valid');
frs   = linspace(fr(1), fr(2), fr(2)-fr(1)+1);
times = 0: num_spec;
times = times *bin + 0;
figure
imagesc(times, frs, spectrogram1);
set(gca,'YDir','normal');
xlabel('Time(ms)','fontsize',11);
ylabel('Freq(Hz)','fontsize',11);
hcb=colorbar;
%set(hcb, 'YTick',[0,100,200,300],'Position',[0.89,0.41, 0.03, 0.24]);
ylabel(hcb,'(spikes/sec)^2/Hz','fontsize',11);
set(gca,'fontsize',11);
%set(gca,'xtick',[]);
%colormap turbo;
%caxis([0 300]);
%title(name);
set(gcf,'Position',[10,10,1500,300]);
end




         
    
    