function res = model_ode_complete(param,init)
%this function delete the effect of pending spike on variance
%and add effect of connecting probability on variance.
%this is the ode model of gamma
duration= param.duration*1000;
delta_time=param.delta_time;
dt=param.dt;
ne=param.ne;
ni=param.ni;
s_ee=param.s_ee;
s_ie=param.s_ie;
s_ei=param.s_ei;
s_ii=param.s_ii;
p_ee=param.p_ee;
p_ie=param.p_ie;
p_ei=param.p_ei;
p_ii=param.p_ii;

lambda_e=param.lambda_e;
lambda_i=param.lambda_i;

Mr=param.Mr;
M=param.M;

%peak_e&i will record the state of each peak using 3 variables, it can
%record 10 peaks maximally. The first row is the mean of each peak, the
%second is the variance while the third is the share of each peak in total
%population.

if isempty(init)
    peak_e=zeros(3,10);
    peak_i=zeros(3,10);
    peak_e(3,1)=1;
    peak_i(3,1)=1;
    h=zeros(1,4);
    npe=1;
    npi=1;
    index_e=1;
    index_i=1;
else
    peak_e=init.peak_e;
    peak_i=init.peak_i;
    h=init.h;
    npe=init.npe;
    npi=init.npi;
    index_e=init.index_e;
    index_i=init.index_i;
end



%npe&i record the number of peaks.



tau=zeros(1,4);
tau(1)=param.tau_ee;
tau(2)=param.tau_ie;
tau(3)=param.tau_ei;
tau(4)=param.tau_ii;
%in order: ee,ie,ei,ii


t_temp=0;

res.peak_e=zeros(duration/delta_time,30);
res.peak_i=zeros(duration/delta_time,30);
res.npi = zeros(duration/delta_time,1);
res.npe = zeros(duration/delta_time,1);
res.h=zeros(duration/delta_time,4);
res.index_e=zeros(duration/delta_time,1);
res.index_i=zeros(duration/delta_time,1);





record=zeros(duration/dt,8);




%kernel of mean and variance
t=tau(1)*4:-dt:dt;
lk_ee=round(tau(1)*4/dt);
mk_ee=exp(-t/tau(1))/tau(1)*s_ee;
vk_ee=(2/tau(1)*exp(-2*t/tau(1))-exp(-t/tau(1))/tau(1))*s_ee^2;


t=tau(2)*4:-dt:dt;
lk_ie=round(tau(2)*4/dt);
mk_ie=exp(-t/tau(2))/tau(2)*s_ie;
vk_ie=(2/tau(2)*exp(-2*t/tau(2))-exp(-t/tau(2))/tau(2))*s_ie^2;


t=tau(3)*4:-dt:dt;
lk_ei=round(tau(3)*4/dt);
mk_ei=(-1/tau(3)*exp(-t/tau(3))*s_ei/(M+Mr))./(-(1-exp(-t/tau(3)))*s_ei/(M+Mr)+1);
vk_ei=(-exp(-t/tau(3))/tau(3)+exp(-t/tau(3))/tau(3)*(1-s_ei/(M+Mr))^2)./(exp(-t/tau(3))+(1-exp(-t/tau(3)))*(1-s_ei/(M+Mr))^2);


t=tau(4)*4:-dt:dt;
lk_ii=round(tau(4)*4/dt);
mk_ii=(-1/tau(4)*exp(-t/tau(4))*s_ii/(M+Mr))./(-(1-exp(-t/tau(4)))*s_ii/(M+Mr)+1);
vk_ii=(-exp(-t/tau(4))/tau(4)+exp(-t/tau(4))/tau(4)*(1-s_ii/(M+Mr))^2)./(exp(-t/tau(4))+(1-exp(-t/tau(4)))*(1-s_ii/(M+Mr))^2);

t=0;
lmax=round(max(tau)*4/dt);

spikecount_e=zeros(lmax+duration/dt,1);
spikecount_i=zeros(lmax+duration/dt,1);

%pending spike variance


ntime=lmax;

while t < duration
    t = t+dt;
    t_temp = t_temp+dt;
    update();
    
    if t_temp > 0.99*delta_time
        index1=round(t/delta_time);
        res.peak_e(index1,:)=peak_e(:)';
        res.peak_i(index1,:)=peak_i(:)';
        res.npi(index1) = npi;
        res.npe(index1) = npe;
        res.index_e(index1)=index_e;
        res.index_i(index1)=index_i;

        t_temp=0;
    end
end
res.spikecount_e=spikecount_e;
res.spikecount_i=spikecount_i;
res.record=record;




    function update()
        %峰的运动一共分为三种情况。对于每个峰，我们认为它只有3sigma宽。3sigma以外的分布不考虑。
        %case1的含义是只有一个峰，且右边界（mean+3*sigma）<M
        %当峰逐渐右移，右边界超过M，来到case2。该状态代表一部分神经元放电，如果下一个dt里有神经元放电（dh>0）,则保持case2。
        %case2会一直持续到HI抑制强于向右推的力量，峰会向左移动，一旦dh<0，生成一个新的小峰，并更新到case3，保持case3一直到dh重新大于0。
        %当主峰全部过去（peak第一列对应的峰），所有小峰按照相同的方差均值合并成一个主峰。
        %注意！所有小峰都是完整的高斯，但是主峰不是！主峰在case2和3里会是一个右边截断一部分的高斯！通过peak（3,1）即主峰剩余的百分比可以计算
        %主峰此刻截断到哪里了。详见phinv函数（phi函数的反函数，phi函数是gauss的CDF）。
        ntime=ntime+1;
        
        dpsm_ee=mk_ee*spikecount_e(ntime-lk_ee:ntime-1);
        dpsv_ee=vk_ee*spikecount_e(ntime-lk_ee:ntime-1);
        dpsm_ei=mk_ei*spikecount_i(ntime-lk_ei:ntime-1);
        dpsv_ei=vk_ei*spikecount_i(ntime-lk_ei:ntime-1);
        
        lm=log(peak_e(1,1:npe)+Mr);
        lmv=log(peak_e(2,1:npe)+(peak_e(1,1:npe)+Mr).^2);

        lm=lm+dpsm_ei*dt;
        lmv=lmv+dpsv_ei*dt;
        
        dv_i=exp(lmv)-exp(2*lm)-peak_e(2,1:npe);
        
        record(ntime-lmax,1:4)=[dpsm_ee,dpsv_ee,dpsm_ei,dpsv_ei];
        
        peak_e(1,1:npe)=exp(lm)+(lambda_e+dpsm_ee)*dt-Mr;
        peak_e(2,1:npe)=peak_e(2,1:npe)+(lambda_e+dpsv_ee)*dt+dv_i;
        
        switch index_e
            case 1
                spikecount_e(ntime)=0;
                if peak_e(1,1)+3*sqrt(peak_e(2,1))>M
                    dh = peak_e(3,1) - phi(M-peak_e(1,1),peak_e(2,1));
                    index_e = 2;
                    npe = npe + 1;
                    [m_new, v_new] = newpeak(peak_e(2,1), phi(M-peak_e(1,1),peak_e(2,1)), peak_e(3,1));
                    peak_e(1,npe) = m_new;
                    peak_e(2,npe) = v_new;
                    peak_e(3,npe) = dh;
                    peak_e(3,1) = peak_e(3,1) - dh;
                    spikecount_e(ntime)=dh*ne;
                end
                
            case 2
                dh=peak_e(3,1)-phi(M-peak_e(1,1),peak_e(2,1));
                if dh>0
                    [m_new, v_new] = newpeak(peak_e(2,1), phi(M-peak_e(1,1),peak_e(2,1)), peak_e(3,1));
                    % Merge the small area and the left peak
                    m_old = peak_e(1,npe);
                    peak_e(1,npe) = (m_new * dh + peak_e(1,npe) *peak_e(3,npe))/(dh + peak_e(3,npe));
                    peak_e(2,npe) = ((m_new^2 + v_new)*dh + (m_old^2 + peak_e(2, npe))*peak_e(3, npe))/(dh+peak_e(3,npe)) - peak_e(1,npe)^2;
                    peak_e(3,npe)=peak_e(3,npe)+dh;
                    peak_e(3,1)=peak_e(3,1)-dh;
                    
                    if peak_e(3,1)<0.0001
                        %                         [m_new,v_new]=newpeak(peak_e(2,1),peak_e(3,1),peak_e(3,npe+1)+peak_e(3,1));
                        %                         peak_e(1,1)=m_new;
                        %                         peak_e(2,1)=v_new;
                        %
                        %                         peak_e(3,1)=peak_e(3,npe+1)+peak_e(3,1);
                        %                         peak_e(3,npe+1)=0;
                        m_new = sum(peak_e(3,1:npe).*peak_e(1,1:npe));
                        peak_e(2,1) = sum(peak_e(3,1:npe).*(peak_e(2,1:npe)+peak_e(1,1:npe).^2))-m_new^2;
                        %                         (2,1) = peak_e(2,1)/3;
                        peak_e(1,1) = m_new;
                        peak_e(3,1) = 1;
                        peak_e(:,2:10) = 0;
                        npe = 1;
                        index_e = 1;
                    end
                else
                    %                     npe=npe+1;
                    %                     [m_new,v_new]=newpeak(peak_e(2,1),peak_e(3,1),peak_e(3,npe)+peak_e(3,1));
                    %                     peak_e(1,npe)=m_new;
                    %                     peak_e(2,npe)=v_new;
                    index_e = 3;
                end
                
                spikecount_e(ntime)=max(dh,0)*ne;
            case 3
                dh=peak_e(3,1)-phi(M-peak_e(1,1),peak_e(2,1));
                if dh>0
                    npe = npe + 1;
                    [m_new, v_new] = newpeak(peak_e(2,1), phi(M-peak_e(1,1),peak_e(2,1)), peak_e(3,1));
                    peak_e(1,npe) = m_new;
                    peak_e(2,npe) = v_new;
                    peak_e(3,npe) = dh;
                    peak_e(3,1) = peak_e(3,1) - dh;
                    index_e = 2;
                end
                spikecount_e(ntime)=max(dh,0)*ne;
        end
        if peak_e(3,2)*peak_e(2,2)~=0
            if (peak_e(1,1)-sqrt(peak_e(2,1))<peak_e(1,2)+sqrt(peak_e(2,2))) || (peak_e(1,2)+3*sqrt(peak_e(2,2))>M)
                [~,v]=newpeak(peak_e(2,1),0.0001,peak_e(3,1));
                peak_e(2,1)=v;
                m_new=sum(peak_e(3,1:2).*peak_e(1,1:2))/(peak_e(3,1)+peak_e(3,2));
                v_new=sum(peak_e(3,1:2).*(peak_e(2,1:2)+peak_e(1,1:2).^2))/(peak_e(3,1)+peak_e(3,2))-m_new^2;
                peak_e(1,1)=m_new;
                peak_e(2,1)=v_new;
                peak_e(3,1)=peak_e(3,1)+peak_e(3,2);
                peak_e(:,2:9)=peak_e(:,3:10);
                npe=npe-1;
                if npe==1
                    index_e = 1;
                end
            end
        end
        
        dpsm_ie=mk_ie*spikecount_e(ntime-lk_ie:ntime-1);
        dpsv_ie=vk_ie*spikecount_e(ntime-lk_ie:ntime-1);
        dpsm_ii=mk_ii*spikecount_i(ntime-lk_ii:ntime-1);
        dpsv_ii=vk_ii*spikecount_i(ntime-lk_ii:ntime-1);
        
        lm=log(peak_i(1,1:npi)+Mr);
        lmv=log(peak_i(2,1:npi)+(peak_i(1,1:npi)+Mr).^2);
        lm=lm+dpsm_ii*dt;
        lmv=lmv+dpsv_ii*dt;
        
        dv_i=exp(lmv)-exp(2*lm)-peak_i(2,1:npi);
        record(ntime-lmax,5:8)=[dpsm_ie,dpsv_ie,dpsm_ii,dpsv_ii];
        
        peak_i(1,1:npi)=exp(lm)+(lambda_i+dpsm_ie)*dt-Mr;
        peak_i(2,1:npi)=peak_i(2,1:npi)+(lambda_i+dpsv_ie)*dt+dv_i;
        
        switch index_i
            case 1
                spikecount_i(ntime)=0;
                if peak_i(1,1)+3*sqrt(peak_i(2,1))>M
                    dh=peak_i(3,1)-phi(M-peak_i(1,1),peak_i(2,1));
                    
                    index_i = 2;
                    h(3)=h(3)+dh*ni*p_ei;
                    h(4)=h(4)+dh*ni*p_ii;
                    npi = npi+1;
                    [m_new,v_new]=newpeak(peak_i(2,1),phi(M-peak_i(1,1),peak_i(2,1)),peak_i(3,1));
                    peak_i(1, npi) = m_new;
                    peak_i(2, npi) = v_new;
                    peak_i(3, npi) = dh;
                    peak_i(3,1) = peak_i(3,1) - dh;
                    spikecount_i(ntime)=dh*ni;
                end
                
            case 2
                dh=peak_i(3,1)-phi(M-peak_i(1,1),peak_i(2,1));
                if dh>0
                    [m_new, v_new] = newpeak(peak_i(2,1), phi(M-peak_i(1,1),peak_i(2,1)), peak_i(3,1));
                    % Merge the small area and the left peak
                    m_old = peak_i(1,npi);
                    peak_i(1,npi) = (m_new * dh + peak_i(1,npi) *peak_i(3,npi))/(dh + peak_i(3,npi));
                    peak_i(2,npi) = ((m_new^2 + v_new)*dh + (m_old^2 + peak_i(2, npi))*peak_i(3, npi))/(dh+peak_i(3,npi)) - peak_i(1,npi)^2;
                    peak_i(3,npi)=peak_i(3,npi)+dh;
                    peak_i(3,1)=peak_i(3,1)-dh;

                    if peak_i(3,1)<0.0001
                        %                         [m_new,v_new]=newpeak(peak_i(2,1),peak_i(3,1),peak_i(3,npi+1)+peak_i(3,1));
                        %                         peak_i(1,1)=m_new;
                        %                         peak_i(2,1)=v_new;
                        %
                        %                         peak_i(3,1)=peak_i(3,npi+1)+peak_i(3,1);
                        %                         peak_i(3,npi+1)=0;
                        m_new = sum(peak_i(3,1:npi).*peak_i(1,1:npi));
                        peak_i(2,1) = sum(peak_i(3,1:npi).*(peak_i(2,1:npi)+peak_i(1,1:npi).^2))-m_new^2;
                        
                        %                         peak_i(2,1)=peak_i(2,1)/3;
                        
                        peak_i(1,1) = m_new;
                        peak_i(3,1) = 1;
                        peak_i(:,2:10) = 0;
                        npi = 1;
                        index_i = 1;
                    end
                else
                    %                     npi=npi+1;
                    %                     [m_new,v_new]=newpeak(peak_i(2,1),peak_i(3,1),peak_i(3,npi)+peak_i(3,1));
                    %                     peak_i(1,npi)=m_new;
                    %                     peak_i(2,npi)=v_new;
                    index_i = 3;
                end

                spikecount_i(ntime)=max(dh,0)*ni;
            case 3
                dh=peak_i(3,1)-phi(M-peak_i(1,1),peak_i(2,1));
                if dh>0
                    npi = npi+1;
                    [m_new,v_new]=newpeak(peak_i(2,1),phi(M-peak_i(1,1),peak_i(2,1)),peak_i(3,1));
                    peak_i(1,npi) = m_new;
                    peak_i(2,npi) = v_new;
                    peak_i(3,npi) = dh;
                    peak_i(3,1) = peak_i(3,1) - dh;

                    index_i = 2;
                end
                spikecount_i(ntime)=max(dh,0)*ni;
        end
        if peak_i(3,2)*peak_i(2,2)~=0
            if (peak_i(1,1)-sqrt(peak_i(2,1))<peak_i(1,2)+sqrt(peak_i(2,2))) || (peak_i(1,2)+3*sqrt(peak_i(2,2))>M)
                [~,v]=newpeak(peak_i(2,1),0.0001,peak_i(3,1));
                peak_i(2,1)=v;
                m_new=sum(peak_i(3,1:2).*peak_i(1,1:2))/(peak_i(3,1)+peak_i(3,2));
                v_new=sum(peak_i(3,1:2).*(peak_i(2,1:2)+peak_i(1,1:2).^2))/(peak_i(3,1)+peak_i(3,2))-m_new^2;
                peak_i(1,1)=m_new;
                peak_i(2,1)=v_new;
                peak_i(3,1)=peak_i(3,1)+peak_i(3,2);
                peak_i(:,2:9)=peak_i(:,3:10);
                npi=npi-1;
                if npi==1
                    index_i = 1;
                end
            end
        end
    end


    function h = phi(x,v)
        t
        % Return the CDF of the Gaussian N(0,sqrt(v)).
        h=0.5*(1+erf(x/sqrt(v)/sqrt(2)));
    end

    function x = phinv(v,h)
        % The inverse of CDF of the Gaussian N(0,sqrt(v)).
        % note: assume original distribution has zero mean so that there is
        % no m variable.
        h=min(h,0.9999);
        h=max(h,0.0001);
        x = erfinv((2*h-1))*sqrt(2)*sqrt(v);
    end

    function [m,v]=newpeak(v,a,b)
        %this function calculates the mean and variance of new peak. a<b are the ratios of the region.
        a=phinv(v,a);
        b=phinv(v,b);
        x=[a+(b-a)/200:(b-a)/100:b-(b-a)/200];
        % Calculate the new mean
        m=(1/sqrt(2*pi*v)*exp(-x.^2/2/v).*((b-a)/100))*(x-a)';
        % Calculate the new variance
        v=(1/sqrt(2*pi*v)*exp(-x.^2/2/v).*((b-a)/100))*((x-a).^2)'-m^2;
    end
end