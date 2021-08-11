
spikecount_e=a.spikecount_e;
spikecount_i=a.spikecount_i;
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
