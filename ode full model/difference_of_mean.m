
spikecount_e=res_lif.spikecount_e;
spikecount_i=res_lif.spikecount_i;
mfe_init=[];
i=1;
while i<length(spikecount_e)
    if spikecount_e(i)~=0
        mfe_init=[mfe_init,i];
        i=i+50;
    else
        i=i+1;
    end
end
%%

dm=zeros(1,length(mfe_init));
for i=1:length(mfe_init)
    ind=mfe_init(i)-3;
    dm(i)=mean(a.VE(ind,:))-mean(a.VI(ind,:));
end

%%

dm=zeros(21,400);
res=cell(21,400);
for s=1:21
    s
    d = -10 +s;
    init.ve=sqrt(75)*randn(1,300)+d;
    init.vi=sqrt(80)*randn(1,100);
    m=max(max(init.ve),max(init.vi))-100;
    init.ve=init.ve-m-1;
    init.vi=init.vi-m-1;
    init.he=zeros(2,300);
    init.hi=zeros(2,100);
    param2.duration=0.05;
    for trial=1:400
        res_lif=ode_full(param2,init);
        mfe_init=[];
        i=1;
        while res_lif.spikecount_e(i)==0
            i=i+1;
        end
        i=i+80;
        while res_lif.spikecount_e(i)==0
            i=i+1;
        end
        j=1;
        while res_lif.spikecount_i(j)==0
            j=j+1;
        end
        j=j+80;
        while res_lif.spikecount_i(j)==0
            j=j+1;
        end
        ind=min(i,j);
        dm(s,trial)=mean(res_lif.VE(ind,:))-mean(res_lif.VI(ind,:));
        res{s,trial}=res_lif;
    end
end

%%
<<<<<<< Updated upstream
m=mean(dm,2);
v=sqrt(var(dm')/400);
errorbar([-10:10],m,v);
xlabel('m_{VE}-m_{VI} before MFE');
ylabel('m_{VE}-m_{VI} after MFE');
=======

>>>>>>> Stashed changes
