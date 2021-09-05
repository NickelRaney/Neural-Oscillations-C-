function [res]=ode_full(param,esi_e,esi_i)
%This is the ode version of full model. 
ne = param.ne;
ni = param.ni;
p_ee = param.p_ee;
p_ie = param.p_ie;
p_ei = param.p_ei;
p_ii = param.p_ii;
see = param.s_ee;
sie = param.s_ie;
sei = param.s_ei;
sii = param.s_ii;
tau_ee = param.tau_ee;
tau_ie = param.tau_ie;
tau_ei = param.tau_ei;
tau_ii = param.tau_ii;
dt = param.gridsize;

ref_e = param.tau_re/dt-1;
ref_i = param.tau_ri/dt-1;

M = param.M;
Mr = param.Mr;
duration = param.duration*1000;

lambda_e = param.lambda_e/1000;
lambda_i = param.lambda_i/1000;


esi_e=exprnd(1/lambda_e,lambda_e*(duration*1.05),ne); %external spike interval
esi_i=exprnd(1/lambda_i,lambda_i*(duration*1.05),ni); %external spike interval
while min(sum(esi_e))<duration || min(sum(esi_i))<duration
    esi_e=exprnd(1/lambda_e,lambda_e*(duration*1.05),ne); %external spike interval
    esi_i=exprnd(1/lambda_i,lambda_i*(duration*1.05),ni); %external spike interval
end



ex_e=zeros(duration/dt,ne);
ex_i=zeros(duration/dt,ni);
%convert esi to the external effect at each time point.
qe=exp(-dt/tau_ee);
qi=exp(-dt/tau_ie);
for i=1:ne
    t=esi_e(1,i);
    count=1;
    while t<duration
        ind=ceil(t/dt);
        ex_e(ind,i)=ex_e(ind,i)+exp((t-ind*dt)/tau_ee)/tau_ee;
        count=count+1;
        t=t+esi_e(count,i);
    end
end
for i=1:ni
    t=esi_i(1,i);
    count=1;
    while t<duration
        ind=ceil(t/dt);
        ex_i(ind,i)=ex_i(ind,i)+exp((t-ind*dt)/tau_ie)/tau_ie;
        count=count+1;
        t=t+esi_i(count,i);
    end
end
for i=2:size(ex_e,1)
    ex_e(i,:)=ex_e(i,:)+ex_e(i-1,:)*qe;
    ex_i(i,:)=ex_i(i,:)+ex_i(i-1,:)*qi;
end



he=zeros(2,300);
hi=zeros(2,100);


ve=zeros(1,ne);
vi=zeros(1,ni);
refc_e=zeros(1,ne);%ref clock
refc_i=zeros(1,ni);

sl=duration;
spike_e=zeros(sl,ne);
spike_i=zeros(sl,ni);

nearray=[0:ne-1];
niarray=[0:ni-1];

res.HE = zeros(ceil(duration/dt)+1,ne+ni);
res.HI = zeros(ceil(duration/dt)+1,ne+ni);
res.VE = zeros(ceil(duration/dt)+1,ne);
res.VI = zeros(ceil(duration/dt)+1,ni);

res.spikecount_e = zeros(ceil(duration/dt)+1,1);
res.spikecount_i = zeros(ceil(duration/dt)+1,1);

for step=2:duration/dt
    rind_e= refc_e==0; %ref index
    rind_i= refc_i==0;
    refc_e(~rind_e)=refc_e(~rind_e)-1;
    refc_i(~rind_i)=refc_i(~rind_i)-1;
    
    ve(rind_e)=ve(rind_e)+(7+see*he(1,rind_e)/tau_ee-sei*he(2,rind_e)/tau_ei.*(ve(rind_e)+Mr)/(M+Mr))*dt;
    vi(rind_i)=vi(rind_i)+(7+sie*hi(1,rind_i)/tau_ie-sii*hi(2,rind_i)/tau_ii.*(vi(rind_i)+Mr)/(M+Mr))*dt;
    
    he=he.*[exp(-dt/tau_ee);exp(-dt/tau_ei)];
    hi=hi.*[exp(-dt/tau_ie);exp(-dt/tau_ii)];
    
    sind_e = ve>100; %spike index
    sind_i = vi>100; %spike index
    spikecount_e=sum(sind_e);
    spikecount_i=sum(sind_i);
    res.spikecount_e(step-1)=spikecount_e;
    res.spikecount_i(step-1)=spikecount_i;
    
    if spikecount_e>0
        if ref_e<0
            ve(sind_e)=0;
        else
            ve(sind_e)=0;
            refc_e(sind_e)=refc_e(sind_e)+ref_e;
        end
        he(1,:)=he(1,:)+binornd(spikecount_e,p_ee,1,ne);
        hi(1,:)=hi(1,:)+binornd(spikecount_e,p_ie,1,ni);
        spike_e(1,sind_e)=spike_e(1,sind_e)+1; %write down the spike time
        spikeind=nearray(sind_e)*sl+spike_e(1,sind_e)+1;
        spike_e(round(spikeind))=step*dt;
    end
    if spikecount_i>0
         if ref_i<0
            vi(sind_i)=vi(sind_i)-100;
         else
            vi(sind_i)=0;
            refc_i(sind_i)=refc_i(sind_i)+ref_i;
         end
        he(2,:)=he(2,:)+binornd(spikecount_i,p_ei,1,ne);
        hi(2,:)=hi(2,:)+binornd(spikecount_i,p_ii,1,ni);
        spike_i(1,sind_i)=spike_i(1,sind_i)+1; %write down the spike time
        spikeind=niarray(sind_i)*sl+spike_i(1,sind_i)+1;
        spike_i(round(spikeind))=step*dt;
    end
    res.VE(step,:)=ve(:);
    res.VI(step,:)=vi(:);
    res.HE(step,:)=[he(1,:),hi(1,:)];
    res.HI(step,:)=[he(2,:),hi(2,:)];
    res.spikecount_e(step) = spikecount_e;
    res.spikecount_i(step) = spikecount_i;
end
res.spike=[spike_e,spike_i];
res.spike(2:end,:)=res.spike(2:end,:)/1000;
end
