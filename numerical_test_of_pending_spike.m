%% numerical experiment of pending spikes
%%
spike=exprnd(1,1,3000);
var_spike=zeros(1,201);
for t=0.05:0.05:10
    var_spike(round(t/0.05))=var(spike<t);
end

% var_spike(41:end)=var_spike(41:end)+0.5*var_spike(1:161);
plot(var_spike);
hold on
%%
plot(var_spike2);
hold on
%%
init=20*randn(1,3000)+30;
var_init=var(init)
var_spike=zeros(1,200);
spike=exprnd(1,400,3000);
spike=[spike;exprnd(1,200,3000)+2];
for t=0.05:0.05:10
    var_spike(round(t/0.05))=var(init+sum(spike<t));
end
plot(var_spike-var_init);
%% numerical test of pending I spike in another coordinate.

init_e=randn(1,10000)*sqrt(30)+50;
init_i=pic(init_e);
%%

spike=[exprnd(4.5,20,10000);exprnd(4.5,20,10000)+2];
var_spike_i=zeros(1,200);
mean_spike_i=zeros(1,200);
var_spike_e=zeros(1,200);
mean_spike_e=zeros(1,200);
d=init_e;
for t=0.05:0.05:10
    
    index = sum(spike > t-0.05 & spike<= t,1);
    d=update(d,index);
    var_spike_e(round(t/0.05))=var(d);
    mean_spike_e(round(t/0.05))=mean(d);
    var_spike_i(round(t/0.05))=var(init_i-sum(spike<t,1));
    mean_spike_i(round(t/0.05))=mean(init_i-sum(spike<t,1));
end

%%
init_i=randn(1,10000)*sqrt(30)+10;
init_e=pec(init_i);
histogram(init_e);
%%
mean_c=zeros(1,200);
var_c=zeros(1,200);
for i=1:200
    [mean_c(i),var_c(i)]=pecd(mean_spike_i(i),var_spike_i(i));
end

%%
[m,v]=picd(40.32,61.27)
[m2,v2]=picd(40.32,61.271)
v2-v
%%
function d=update(d,index)
for i=1:length(d)
    for j=1:index(i)
        d(i)=d(i)-(d(i)+66)/166;
    end
end
end


%%
function [m1,v1]=picd(m,v)
if v==0
    v1=0;
    m1=pic(m);
else
    x=m-3*sqrt(v):6*sqrt(v)/50:m+3*sqrt(v);
    y=1/sqrt(2*pi*v)*exp(-(x-m).^2/2/v)*6*sqrt(v)/50;
    x=pic(x');
    m1=y*x;
    v1=y*x.^2-m1^2;
end
end

function [m1,v1]=pecd(m,v)
if v==0
    v1=0;
    m1=pec(m);
else
    x=m-3*sqrt(v):6*sqrt(v)/20:m+3*sqrt(v);
    y=1/sqrt(2*pi*v)*exp(-(x-m).^2/2/v)*6*sqrt(v)/20;
    x=pec(x');
    m1=y*x;
    v1=y*x.^2-m1^2;
end
end

function y=pic(x)
y=166*log((66+x)/166)+100;
end

function x=pec(y)
x=exp((y-100)/166)*166-66;
end

