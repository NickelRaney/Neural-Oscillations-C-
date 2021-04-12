tic;
system('.\\models\\NOModel_full.exe');
toc;

%%
ne = 3000;
ni = 1000;
spike_som = load("outputs//model_full//spike.txt");
index = spike_som(:,2)<ne;
spike_E = spike_som(index,:);
index = logical((spike_som(:,2)>=ne).*(spike_som(:,2)<=ne+ni));
spike_I = spike_som(index,:);

scatter(spike_E(:,1)*1000, spike_E(:,2)+1,10,'.','r');
hold on;
scatter(spike_I(:,1)*1000, spike_I(:,2)+1,10,'.','b');
hold on;
xlim([0,3000]);
title('Raster-Full-n=4000-t=3');
set(gcf,'Position',[10,10,2000,300]);
saveas(gcf,'Raster-Full-n=4000-t=3.png');
%%
spike_som = load("outputs//model_full//spike.txt");
spike = zeros(10000,ne+ni);
for i=1:size(spike_som,1)
    neuron_i = spike_som(i,2);
    spike(1,neuron_i+1) = spike(1,neuron_i+1) + 1;
    spike(spike(1,neuron_i+1)+1,neuron_i+1) = spike_som(i,1)*1000;
end

param.ne = ne;
param.ni = ni;
param.duration = 3000;
param.sdbin = 2.5;

param.spectrogram_timewindow = 200;
param.frequency_range = [5,100];

res.spike = spike;
sd = spikedensity(res, param);
spectrogram(sd.e, param);
title("Spec-Full-n=4000-t=3");
set(gcf,'Position',[10,10,2000,300]);
saveas(gcf,'Spec-Full-n=4000-t=3.png');