function res = model_ode(param)
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

lambda_e=param.lambda_e;
lambda_i=param.lambda_i;

Mr=param.Mr;
M=param.M;

%peak_e&i will record the state of each peak using 3 variables, it can
%record 10 peaks maximally. The first row is the mean of each peak, the
%second is the variance while the third is the share of each peak in total
%population.

peak_e=zeros(3,10);
peak_i=zeros(3,10);

%npe&i record the number of peaks.

npe=1;
npi=1;
peak_e(3,1)=1;
peak_i(3,1)=1;
h=zeros(1,4);
tau=zeros(1,4);
tau(1)=param.tau_ee;
tau(2)=param.tau_ie;
tau(3)=param.tau_ei;
tau(4)=param.tau_ii;
%in order: ee,ie,ei,ii

t=0;
t_temp=0;

res.peak_e=zeros(duration/delta_time,30);
res.peak_i=zeros(duration/delta_time,30);
res.npi = zeros(duration/delta_time,1);
res.npe = zeros(duration/delta_time,1);
res.h=zeros(duration/delta_time,4);
res.index_e=zeros(duration/delta_time,1);
res.index_i=zeros(duration/delta_time,1);

index_e=1;
index_i=1;
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
        res.h(index1,:)=h(:)';
        t_temp=0;
    end
end

    function update()
        %峰的运动一共分为三种情况。对于每个峰，我们认为它只有3sigma宽。3sigma以外的分布不考虑。
        %case1的含义是只有一个峰，且右边界（mean+3*sigma）<M
        %当峰逐渐右移，右边界超过M，来到case2。该状态代表一部分神经元放电，如果下一个dt里有神经元放电（dh>0）,则保持case2。
        %case2会一直持续到HI抑制强于向右推的力量，峰会向左移动，一旦dh<0，生成一个新的小峰，并更新到case3，保持case3一直到dh重新大于0。
        %当主峰全部过去（peak第一列对应的峰），所有小峰按照相同的方差均值合并成一个主峰。
        %注意！所有小峰都是完整的高斯，但是主峰不是！主峰在case2和3里会是一个右边截断一部分的高斯！通过peak（3,1）即主峰剩余的百分比可以计算
        %主峰此刻截断到哪里了。详见phinv函数（phi函数的反函数，phi函数是gauss的CDF）。
        h_temp=h;
        switch index_e
            case 1
                M_e=(Mr+peak_e(1,1:npe))./(M+Mr);
                V_e(1,:)=-2*peak_e(2,1:npe)./(M+Mr);
                V_e(2,:)=((peak_e(1,1:npe)+Mr)/(Mr+M)).^2+peak_e(2,1:npe)/(Mr+M)^2;
                peak_e(2,1:npe)=peak_e(2,1:npe)+(lambda_e+1/tau(1)*h_temp(1)*s_ee^2+1/tau(3)*V_e(1,:)*h_temp(3)*s_ei+1/tau(3)*V_e(2,:)*h_temp(3)*s_ei^2)*dt;
                peak_e(1,1:npe)=peak_e(1,1:npe)+(lambda_e+1/tau(1)*h_temp(1)*s_ee-1/tau(3)*M_e*h_temp(3)*s_ei)*dt;
                if peak_e(1,1)+3*sqrt(peak_e(2,1))>M
                    index_e = 2;
                end
            case 2
                M_e=(Mr+peak_e(1,1:npe))./(M+Mr);
                V_e(1,:)=-2*peak_e(2,1:npe)./(M+Mr);
                V_e(2,:)=((peak_e(1,1:npe)+Mr)/(Mr+M)).^2+peak_e(2,1:npe)/(Mr+M)^2;
                peak_e(2,1:npe)=peak_e(2,1:npe)+(lambda_e+1/tau(1)*h_temp(1)*s_ee^2+1/tau(3)*V_e(1,:)*h_temp(3)*s_ei+1/tau(3)*V_e(2,:)*h_temp(3)*s_ei^2)*dt;
                peak_e(1,1:npe)=peak_e(1,1:npe)+(lambda_e+1/tau(1)*h_temp(1)*s_ee-1/tau(3)*M_e*h_temp(3)*s_ei)*dt;
                dh=peak_e(3,1)-phi(M-peak_e(1,1),peak_e(2,1));
                if dh>0
                    peak_e(3,npe+1)=peak_e(3,npe+1)+dh;
                    peak_e(3,1)=peak_e(3,1)-dh;
                    h(1)=h(1)+dh*ne;
                    h(2)=h(2)+dh*ne;
                    if peak_e(3,1)<0.0013
                        [m_new,v_new]=newpeak(peak_e(2,1),peak_e(3,1),peak_e(3,npe+1)+peak_e(3,1));
                        peak_e(1,1)=m_new;
                        peak_e(2,1)=v_new;
                        
                        peak_e(3,1)=peak_e(3,npe+1)+peak_e(3,1);
                        peak_e(3,npe+1)=0;
                        m_new = sum(peak_e(3,1:npe).*peak_e(1,1:npe));
                        peak_e(2,1) = sum(peak_e(3,1:npe).*(peak_e(2,1:npe)+peak_e(1,1:npe).^2))-m_new^2;
                        
                        peak_e(2,1) = peak_e(2,1)/3;
                        
                        peak_e(1,1) = m_new;
                        peak_e(3,1) = 1;
                        peak_e(:,2:10) = 0;
                        npe = 1;
                        index_e = 1;
                    end
                else
                    npe=npe+1;
                    [m_new,v_new]=newpeak(peak_e(2,1),peak_e(3,1),peak_e(3,npe)+peak_e(3,1));
                    peak_e(1,npe)=m_new;
                    peak_e(2,npe)=v_new;
                    index_e = 3;
                end
            case 3
                M_e=(Mr+peak_e(1,1:npe))./(M+Mr);
                V_e(1,:)=-2*peak_e(2,1:npe)./(M+Mr);
                V_e(2,:)=((peak_e(1,1:npe)+Mr)/(Mr+M)).^2+peak_e(2,1:npe)/(Mr+M)^2;
                peak_e(2,1:npe)=peak_e(2,1:npe)+(lambda_e+1/tau(1)*h_temp(1)*s_ee^2+1/tau(3)*V_e(1,:)*h_temp(3)*s_ei+1/tau(3)*V_e(2,:)*h_temp(3)*s_ei^2)*dt;
                peak_e(1,1:npe)=peak_e(1,1:npe)+(lambda_e+1/tau(1)*h_temp(1)*s_ee-1/tau(3)*M_e*h_temp(3)*s_ei)*dt;
                dh=peak_e(3,1)-phi(M-peak_e(1,1),peak_e(2,1));
                if dh>0
                    peak_e(3,npe+1)=peak_e(3,npe+1)+dh;
                    peak_e(3,1)=peak_e(3,1)-dh;
                    h(1)=h(1)+dh*ne;
                    h(2)=h(2)+dh*ne;
                    index_e = 2;
                end
        end
        if peak_e(3,2)*peak_e(2,2)~=0
            if (peak_e(1,1)-sqrt(peak_e(2,1))<peak_e(1,2)+sqrt(peak_e(2,2))) || (peak_e(1,2)+3*sqrt(peak_e(2,2))>M)
                [m,v]=newpeak(peak_e(2,1),0.0001,peak_e(3,1));
                peak_e(2,1)=v;
                m_new=sum(peak_e(3,1:2).*peak_e(1,1:2))/(peak_e(3,1)+peak_e(3,2));
                v_new=sum(peak_e(3,1:2).*(peak_e(2,1:2)+peak_e(1,1:2).^2))/(peak_e(3,1)+peak_e(3,2))-m_new^2;
                peak_e(1,1)=m_new;
                peak_e(2,1)=v_new/3;
                peak_e(3,1)=peak_e(3,1)+peak_e(3,2);
                peak_e(:,2:9)=peak_e(:,3:10);
                npe=npe-1;
                if npe==1
                    index_e = 1;
                end
            end
        end
        switch index_i
            case 1
                M_i=(Mr+peak_i(1,1:npi))./(M+Mr);
                V_i(1,:)=-2*peak_i(2,1:npi)./(M+Mr);
                V_i(2,:)=((peak_i(1,1:npi)+Mr)/(Mr+M)).^2+peak_i(2,1:npi)/(Mr+M)^2;
                peak_i(2,1:npi)=peak_i(2,1:npi)+(lambda_i+1/tau(2)*h_temp(2)*s_ie^2+1/tau(4)*V_i(1,:)*h_temp(4)*s_ii+1/tau(4)*V_i(2,:)*h_temp(4)*s_ii^2)*dt;
                peak_i(1,1:npi)=peak_i(1,1:npi)+(lambda_i+1/tau(2)*h_temp(2)*s_ie-1/tau(4)*M_i*h_temp(4)*s_ii)*dt;
                if peak_i(1,1)+3*sqrt(peak_i(2,1))>M
                    index_i = 2;
                end
            case 2
                M_i=(Mr+peak_i(1,1:npi))./(M+Mr);
                V_i(1,:)=-2*peak_i(2,1:npi)./(M+Mr);
                V_i(2,:)=((peak_i(1,1:npi)+Mr)/(Mr+M)).^2+peak_i(2,1:npi)/(Mr+M)^2;
                peak_i(2,1:npi)=peak_i(2,1:npi)+(lambda_i+1/tau(2)*h_temp(2)*s_ie^2+1/tau(4)*V_i(1,:)*h_temp(4)*s_ii+1/tau(4)*V_i(2,:)*h_temp(4)*s_ii^2)*dt;
                peak_i(1,1:npi)=peak_i(1,1:npi)+(lambda_i+1/tau(2)*h_temp(2)*s_ie-1/tau(4)*M_i*h_temp(4)*s_ii)*dt;
                dh=peak_i(3,1)-phi(M-peak_i(1,1),peak_i(2,1));
                if dh>0
                    peak_i(3,npi+1)=peak_i(3,npi+1)+dh;
                    peak_i(3,1)=peak_i(3,1)-dh;
                    h(3)=h(3)+dh*ni;
                    h(4)=h(4)+dh*ni;
                    if peak_i(3,1)<0.0013
                        [m_new,v_new]=newpeak(peak_i(2,1),peak_i(3,1),peak_i(3,npi+1)+peak_i(3,1));
                        peak_i(1,1)=m_new;
                        peak_i(2,1)=v_new;
                        
                        peak_i(3,1)=peak_i(3,npi+1)+peak_i(3,1);
                        peak_i(3,npi+1)=0;
                        m_new = sum(peak_i(3,1:npi).*peak_i(1,1:npi));
                        peak_i(2,1) = sum(peak_i(3,1:npi).*(peak_i(2,1:npi)+peak_i(1,1:npi).^2))-m_new^2;
                        
                        peak_i(2,1)=peak_i(2,1)/3;
                        
                        peak_i(1,1) = m_new;
                        peak_i(3,1) = 1;
                        peak_i(:,2:10) = 0;
                        npi = 1;
                        index_i = 1;
                    end
                else
                    npi=npi+1;
                    [m_new,v_new]=newpeak(peak_i(2,1),peak_i(3,1),peak_i(3,npi)+peak_i(3,1));
                    peak_i(1,npi)=m_new;
                    peak_i(2,npi)=v_new;
                    index_i = 3;
                end
            case 3
                M_i=(Mr+peak_i(1,1:npi))./(M+Mr);
                V_i(1,:)=-2*peak_i(2,1:npi)./(M+Mr);
                V_i(2,:)=((peak_i(1,1:npi)+Mr)/(Mr+M)).^2+peak_i(2,1:npi)/(Mr+M)^2;
                peak_i(2,1:npi)=peak_i(2,1:npi)+(lambda_i+1/tau(2)*h_temp(2)*s_ie^2+1/tau(4)*V_i(1,:)*h_temp(4)*s_ii+1/tau(4)*V_i(2,:)*h_temp(4)*s_ii^2)*dt;
                peak_i(1,1:npi)=peak_i(1,1:npi)+(lambda_i+1/tau(2)*h_temp(2)*s_ie-1/tau(4)*M_i*h_temp(4)*s_ii)*dt;
                dh=peak_i(3,1)-phi(M-peak_i(1,1),peak_i(2,1));
                if dh>0
                    peak_i(3,npi+1)=peak_i(3,npi+1)+dh;
                    peak_i(3,1)=peak_i(3,1)-dh;
                    h(3)=h(3)+dh*ni;
                    h(4)=h(4)+dh*ni;
                    index_i = 2;
                end
        end
        if peak_i(3,2)*peak_i(2,2)~=0
            if (peak_i(1,1)-sqrt(peak_i(2,1))<peak_i(1,2)+sqrt(peak_i(2,2))) || (peak_i(1,2)+3*sqrt(peak_i(2,2))>M)
                [m,v]=newpeak(peak_i(2,1),0.0001,peak_i(3,1));
                peak_i(2,1)=v;
                m_new=sum(peak_i(3,1:2).*peak_i(1,1:2))/(peak_i(3,1)+peak_i(3,2));
                v_new=sum(peak_i(3,1:2).*(peak_i(2,1:2)+peak_i(1,1:2).^2))/(peak_i(3,1)+peak_i(3,2))-m_new^2;
                peak_i(1,1)=m_new;
                peak_i(2,1)=v_new/3;
                peak_i(3,1)=peak_i(3,1)+peak_i(3,2);
                peak_i(:,2:9)=peak_i(:,3:10);
                npi=npi-1;
                if npi==1
                    index_i = 1;
                end
            end
        end
        h = h-(h_temp./tau)*dt;
        
    end


    function h = phi(m,v)
        h=0.5*(1+erf(m/sqrt(v)/sqrt(2)));
    end

    function x = phinv(v,h)
        % note: assume original distribution has zero mean so that there is
        % no m variable.
        h=min(h,0.9999);
        h=max(h,0.0001);
        x = erfinv((2*h-1))*sqrt(2)*sqrt(v);
    end

    function [m,v]=newpeak(v,a,b)
        %this function calculates the mean and variance of new peak. a<b.
        a=phinv(v,a);
        b=phinv(v,b);
        x=[a+(b-a)/200:(b-a)/100:b-(b-a)/200];
        m=(1/sqrt(2*pi*v)*exp(-x.^2/2/v).*((b-a)/100))*(x-a)';
        v=(1/sqrt(2*pi*v)*exp(-x.^2/2/v).*((b-a)/100))*((x-a).^2)'-m^2;
    end
end