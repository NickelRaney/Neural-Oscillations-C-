function [] = rasterplot_SOM(res,param)
% A function to show scatterplots and save them in folder output

ne = param.ne;
ni = param.ni;
nsom = param.nsom;
duration = param.duration;
figure
for i=1:ne
times = res.spike(:,i);
num   = size(times, 1);
scatter(times, i*ones(num, 1),10,'.','r');
hold on
end

for i=(ne+1):(ne+ni)
times = res.spike(:,i);
%times = times(times > 1000);
%times = times(times < 3000);
num   = size(times, 1);
scatter(times, i*ones(num, 1),10,'.','b');
hold on
end

for i=(ne+ni+1):(ne+ni+nsom)
times = res.spike(:,i);
%times = times(times > 1000);
%times = times(times < 3000);
num   = size(times, 1);
scatter(times, i*ones(num, 1),10,'.','g');
hold on
end
xlim([0, 3000]);
%title(name);
yticks([100 200 300 300+ni]);
%xticks([1000 1500 2000 2500 3000]);
%set(gca, 'xtick',[]);
ylabel('Index');
set(gcf,'Position',[10,10,1500,300]);
set(gca,'fontsize',11);
end


