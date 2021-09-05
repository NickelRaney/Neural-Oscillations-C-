addpath('module')

%%
param.ne = 300;
param.ni = 100;
param.M        = 100;
param.Mr       = 66;
param.p_ee     = 0.15;
param.p_ie     = 0.5;
param.p_ei     = 0.5;
param.p_ii     = 0.4;
param.s_ee     = 5;
param.s_ie     = 2;
param.s_ei     = 4.91;
param.s_ii     = 4.91;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.5;
param.s_ii     = 4.91*0.4;

param.lambda_e = 7000;
param.lambda_i = 7000;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.tau_re=0;
param.tau_ri=0;
param.duration = 1;
param.delta_time = 1;
param.delta_time2 = 0.1;
param.factor_Leak=0;
param.LeakE = 20;
param.LeakI = 16.7;
param.factor_Leak = inf;
param.gridsize=0.1;
%%
ind=zeros(2,121);
ind(1,:)=kron([0:0.5:5],ones(1,11));
ind(2,:)=kron(ones(1,11),[0:0.5:5]);
for k=12:length(ind)
    i=ind(1,k);
    j=ind(2,k);
    param2=param;  
    param2.tau_re=i;
    param2.tau_ri=j;
    name=['tau_re=',num2str(i),'-tau_ri=',num2str(j),'-P=1'];
    disp(name);
    [res,dm,res_lif]=dm_plot(param2,name);
end

    function [res,dm,res_lif]=dm_plot(param,name)
    ve = zeros(1,300);
    vi = zeros(1,100);
    he = zeros(2,300);
    hi = zeros(2,100);
    res_lif=model_LIF3(param,ve,vi,he,hi);
    
    dm=zeros(45,24);
    res=cell(45,24);
    index=zeros(45,24);
    
    xrange=[-20:24];
    for d=1:numel(xrange)
        disp(['d= ', num2str(d)]);
        param2=param;
        param2.duration=0.08;
        pardm=zeros(1,24);
        parres=cell(1,24);
        parind=zeros(1,24);
        parfor trial=1:24
            disp(['Trial: ', num2str(trial)]);
            ve=sqrt(75)*randn(1,300)+xrange(d);
            vi=sqrt(80)*randn(1,100);
            m=max(max(ve),max(vi))-100;
            ve=ve-m-1;
            vi=vi-m-1;
            he=zeros(2,300);
            hi=zeros(2,100);
            
            res_lif_small=model_LIF3(param2,ve,vi,he,hi);
            i=1;
            while res_lif_small.spikecount_e(i)==0 && i<length(res_lif_small.spikecount_e)
                i=i+1;
            end
            i=i+170;
            while res_lif_small.spikecount_e(i)==0 && i<length(res_lif_small.spikecount_e)
                i=i+1;
            end
            j=1;
            while res_lif_small.spikecount_i(j)==0 && j<length(res_lif_small.spikecount_i)
                j=j+1;
            end
            j=j+170;
            while res_lif_small.spikecount_i(j)==0 && j<length(res_lif_small.spikecount_i)
                j=j+1;
            end
            ind=min(i,j)
            
            pardm(trial)=mean(res_lif_small.VE(ind,:))-mean(res_lif_small.VI(ind,:));
            parres{trial}=res_lif_small;
            parind(trial)=ind;
        end
        dm(d,:)=pardm;
        res(d,:)=parres;
        index(d,:)=parind;
    end
    
    subplot(2,5,[1,2,6,7]);
    xc=xrange'*ones(1,24);
    scatter(xc(:),dm(:),'.')
    grid on
    hold on
    plot([-30:30],[-30:30]);
    xlim([-30 30]);
    ylim([-30 30]);
    xticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
    yticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
    xlabel('difference of mean in the last MFE');
    ylabel('difference of mean in the next MFE');
    title(['LIF ', name]);
    subplot(2,5,[3,4,5]);
    rasterplot(res_lif,param);
    title('rasterplot')
    subplot(2,5,[8,9,10]);
    plot(mean(res_lif.VE(2:end,:),2)-mean(res_lif.VI(2:end,:),2))
    title('difference of mean')
    set(gcf,'Position',[10,10,1200,400]);
    saveas(gcf ,['figure/LIF/Delta_mean/tau_ref',name,'.fig']);
    close(gcf);
    end