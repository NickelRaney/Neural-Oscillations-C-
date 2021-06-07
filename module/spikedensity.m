function [sd] = spikedensity(res,param)
ne = param.ne;
ni = param.ni;
spike = res.spike;
bin=param.sdbin;
duration_time = param.duration;
t=duration_time*1000/bin;
sd.e=zeros(1,t+1);
sd.i=zeros(1,t+1);
for i=1:ne
    for j=2:spike(1,i)+1
        index = ceil(spike(j,i)*1000/bin);
        if index <= t+1
        sd.e(index) = sd.e(index) +1;
        end
    end
end

for i=ne+1:ne+ni
    for j=2:spike(1,i)+1
        index = ceil(spike(j,i)*1000/bin);
        if index <=t+1
        sd.i(index) = sd.i(index) +1;
        end
    end
end

sd.e=sd.e/((ne)*bin/1000);
sd.i=sd.i/((ni)*bin/1000);