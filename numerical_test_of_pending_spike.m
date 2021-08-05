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

init_e=randn(1,100000)*sqrt(0.1)+50;
init_i=pic(init_e);
%%

spike_e=[exprnd(1.3,10,100000)];
spike_i=[exprnd(4.5,50,100000)+1];
var_spike_e=zeros(1,400);
mean_spike_e=zeros(1,400);
d=init_e;
for t=0.05:0.05:20
    index = sum(spike_i > t-0.05 & spike_i<= t,1);
    d=update(d,index);
    index2 = sum(spike_e > t-0.05 & spike_e<= t,1);
    d=d+index2;
    var_spike_e(round(t/0.05))=var(d);
    mean_spike_e(round(t/0.05))=mean(d);
end

%%
m=mean(init_e);
v=var(init_e);

dt=0.01;
mean_c=zeros(1,20000);
var_c=zeros(1,20000);
dv_c2=zeros(1,20000);
lmv_c2=zeros(1,20000);
count=1;
lm=log(m+66);
lmv=log(v+(m+66)^2);
m=m+66;
for t=0.01:0.01:1

    
    m=m+10*(exp(-t/1.3))/1.3*dt;
    v=v+10*(2/1.3*exp(-2*t/1.3)-exp(-t/1.3)/1.3)*dt;
    
    lmv=log(v+m^2);
    lm=log(m);
    
    lmv_c2(count)=lmv;
    dv_c2(count)=dv_i;
    mean_c(count)=m-66;
    var_c(count)=v;
    count=count+1;
end
 
 for t=1.01:0.01:20
    b=t-1;
    lm=lm+50*(-1/4.5*exp(-b/4.5)/166)/(-(1-exp(-b/4.5))/166+1)*dt;
    lmv=lmv+50*(-1/4.5*exp(-b/4.5)+1/4.5*exp(-b/4.5)*(1-1/166)^2)/(exp(-b/4.5)+(1-exp(-b/4.5))*(1-1/166)^2)*dt;
    
    dv_i=exp(lmv)-exp(2*lm)-v;
    
    m=exp(lm)+10*(exp(-t/1.3))/1.3*dt;
    v=v+10*(2/1.3*exp(-2*t/1.3)-exp(-t/1.3)/1.3)*dt+dv_i;
    
    lmv=log(v+m^2);
    lm=log(m);
    
    lmv_c2(count)=lmv;
    dv_c2(count)=dv_i;
    mean_c(count)=m-66;
    var_c(count)=v;
    count=count+1;
 end
 
%  for t=4.51:0.01:20
%      
%     b=t-1;
%     lm=lm+50*(-1/4.5*exp(-b/4.5)/166)/(-(1-exp(-b/4.5))/166+1)*dt;
%     lmv=lmv+50*(-1/4.5*exp(-b/4.5)+1/4.5*exp(-b/4.5)*(1-1/166)^2)/(exp(-b/4.5)+(1-exp(-b/4.5))*(1-1/166)^2)*dt;
%     a=t-4.5;
%     lm=lm+5*(-1/4.5*exp(-a/4.5)/166)/(-(1-exp(-a/4.5))/166+1)*dt;
%     lmv=lmv+5*(-1/4.5*exp(-a/4.5)+1/4.5*exp(-a/4.5)*(1-1/166)^2)/(exp(-a/4.5)+(1-exp(-a/4.5))*(1-1/166)^2)*dt;
%     
%     dv_i=exp(lmv)-exp(2*lm)-v;
%     m=exp(lm)+10*(exp(-t/1.3))/1.3*dt;
%     v=v+10*(2/1.3*exp(-2*t/1.3)-exp(-t/1.3)/1.3)*dt+dv_i;
%     
%     lmv=log(v+m^2);
%     lm=log(m);
%     
%     mean_c(count)=m-66;
%     var_c(count)=v;
%     count=count+1;
%  end
 %%
subplot(1,2,1)
plot(mean_c(1:5:2000))
hold on
plot(mean_spike_e)
subplot(1,2,2)
plot(var_c(1:5:2000))
hold on
plot(var_spike_e)
 %%
m=mean(init_e);
v=var(init_e);
hi=1;
he=1;
dt=0.01;
mean_c=zeros(1,20000);
var_c=zeros(1,20000);
dv_c1=zeros(1,20000);
lmv_c1=zeros(1,20000);
count=1;
lm=log(m+66);
lmv=log(v+(m+66)^2);
m=m+66;
for t=0:0.01:4.5
    m=m-m/166*hi/4.5*dt;
    lmv=lmv+(-1/4.5*exp(-t/4.5)+1/4.5*exp(-t/4.5)*(1-1/166)^2)/(exp(-t/4.5)+(1-exp(-t/4.5))*(1-1/166)^2)*dt;
    dv_i=exp(lmv)-m^2-v;
   
    m=m+he/1.3*dt;
    v=v+(2/1.3*exp(-2*t/1.3)-exp(-t/1.3)/1.3)*dt+dv_i;
    
    he=he-he/1.3*dt;
    hi=hi-hi/4.5*dt;
    
    lmv=log(v+m^2);
    
    lmv_c1(count)=lmv;
    dv_c1(count)=dv_i;
    mean_c(count)=m-66;
    var_c(count)=v;
    count=count+1;
end

 for t=4.501:0.01:20
     
    lm=lm+(-1/4.5*exp(-t/4.5)/166)/(-(1-exp(-t/4.5))/166+1)*dt;
    lmv=lmv+(-1/4.5*exp(-t/4.5)+1/4.5*exp(-t/4.5)*(1-1/166)^2)/(exp(-t/4.5)+(1-exp(-t/4.5))*(1-1/166)^2)*dt;
    a=t-4.5;
    lm=lm+(-1/4.5*exp(-a/4.5)/166)/(-(1-exp(-a/4.5))/166+1)*dt;
    lmv=lmv+(-1/4.5*exp(-a/4.5)+1/4.5*exp(-a/4.5)*(1-1/166)^2)/(exp(-a/4.5)+(1-exp(-a/4.5))*(1-1/166)^2)*dt;
    
    dv_i=exp(lmv)-exp(2*lm)-v;
    m=exp(lm)+(exp(-t/1.3))/1.3*dt;
    v=v+(2/1.3*exp(-2*t/1.3)-exp(-t/1.3)/1.3)*dt+dv_i;
    
    lmv=log(v+m^2);
    lm=log(m);
    
    mean_c(count)=m-66;
    var_c(count)=v;
    count=count+1;
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
[v]=picd(50,30)
[v2]=picd(50,30+0.01)
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

function dv=ct(dv_ic,m,v,dm)
    
    p=(picd(m+0.0001,v)-picd(m,v))/0.0001;
    
    q=(picd(m,v+0.0001)-picd(m,v))/0.0001;
    dv=(dv_ic-p*dm)/q;
end
function [v1]=picd(m,v)
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

