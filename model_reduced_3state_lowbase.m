function [res] = model_reduced_3state_lowbase(s,param,P)
c                                 = zeros(4,param.ne + param.ni);
c(1,1:param.ne)                   = param.lambda_e;
c(1,param.ne+1:param.ne+param.ni) = param.lambda_i;
c(2,1:param.ne)                   = param.tau_ee;
c(2,param.ne+1:param.ne+param.ni) = param.tau_ie;
c(3,:)                            = param.tau_i;



if s == false
    s=zeros(3,param.ne+param.ni);
end
%s is state matrix. The three rows indicate the triplet(V,H_e,H_i)for
%each neuron. -1 is low_base, 0 is base, 1 is gate.

m=zeros(3,param.ne+param.ni);
m(1,:)=c(1,:);
res.spike=zeros(2000,param.ne+param.ni);
duration_time = param.duration;
t=0;
i=1;
res.V_e = [];
res.V_i = [];
res.H_ie = [];
res.H_ii = [];
res.H_ee = [];
res.H_ei = [];
res.N_LBE = [];
res.N_LBI = [];
res.N_BE = [];
res.N_BI = [];
res.t = [];
time_check = 1;
time_delta = 0;
res.t_H_e = [];
res.t_H_i = [];
res.ratio = [];

while t<= duration_time
    N_LBE                             = sum(s(1,1:param.ne)<0);
    N_LBI                             = sum(s(1,param.ne+1:param.ne+param.ni)<0);
    N_BE                              = sum(s(1,1:param.ne)==0);
    N_BI                              = sum(s(1,param.ne+1:param.ne+param.ni)==0);

    is_spike = 0;
    %disp(["iteration: ", step]);
    m(2:3,:)=c(2:3,:)./s(2:3,:);
    a=-log(1-rand(3,param.ne+param.ni)).*m;
    %exponential distribution random variable
    min_a=min(min(a));
    [x,y]=find(a==min_a);
    %the position of the minimum decides the next operation.
    t = t+a(x,y);
    res.delta_t(i) = a(x,y);
    i = i+1;
    if x == 1 %external input operates
        if s(1,y) == -1 && y<=param.ne
            if rand(1)< P.P_LBE_Ex(N_LBE+1,N_BE+1) 
                s(1,y)=0;
            end
        elseif s(1,y) == 0 && y<= param.ne
            if rand(1)< P.P_BE_Ex(N_LBE+1,N_BE+1) 
                s(1,y)=1;
            end
        elseif s(1,y)==1 && y<=param.ne
            if rand(1)<P.P_GE_Ex(N_LBE+1,N_BE+1)
                is_spike = 1;
                s(1,y) = -1;
            end
        elseif s(1,y) == -1 && y>param.ne
            if rand(1)< P.P_LBI_Ex(N_LBI+1,N_BI+1) 
                s(1,y)=0;
            end
        elseif s(1,y)==0 && y>param.ne
            if rand(1)<P.P_BI_Ex(N_LBI+1,N_BI+1) 
                s(1,y)=1;
            end
        elseif s(1,y)==1 && y>param.ne
            if rand(1)<P.P_GI_Ex(N_LBI+1,N_BI+1) 
                is_spike = 1;
                s(1,y) = -1;
            end
        end
    elseif x == 2 %H_e operates
        s(x,y) = s(x,y)-1;
        if s(1,y) == -1 && y<=param.ne
            if rand(1)< P.P_LBE_E(N_LBE+1,N_BE+1) 
                s(1,y) = 0;
            end
        elseif s(1,y) == 0 && y<=param.ne
            if rand(1)<P.P_BE_E(N_LBE+1,N_BE+1) 
                s(1,y) = 1;
            end
        elseif s(1,y)==1 && y<= param.ne
            if rand(1)<P.P_GE_E(N_LBE+1,N_BE+1) 
                s(1,y)= -1;
                is_spike = 1;
            end
        elseif s(1,y) == -1 && y>param.ne
            if rand(1)< P.P_LBI_E(N_LBI+1,N_BI+1) 
                s(1,y) = 0;
            end
        elseif s(1,y)==0 && y> param.ne
            if rand(1)<P.P_BI_E(N_LBI+1,N_BI+1) 
                s(1,y) = 1;
            end
        elseif s(1,y)==1 && y> param.ne
            if rand(1)< P.P_GI_E(N_LBI+1,N_BI+1) 
                s(1,y)= -1;
                is_spike = 1;
            end
        end
    elseif x==3 %H_i operates
        s(x,y)=s(x,y)-1;
        if s(1,y) == 1 && y <= param.ne
            if rand(1)< P.P_GE_I(N_LBE+1,N_BE+1)
                s(1,y)=0;
            end
        elseif s(1,y) == 0 && y <= param.ne
            if rand(1)< P.P_BE_I(N_LBE+1,N_BE+1)
                s(1,y)=-1;
            end  
        elseif s(1,y)==1 && y > param.ne
            if rand(1)< P.P_GI_I(N_LBI+1,N_BI+1)
                s(1,y)=0;
            end
        elseif s(1,y) == 0 && y > param.ne
            if rand(1)< P.P_BI_I(N_LBI+1,N_BI+1)
                s(1,y)=-1;
            end 
        end
    end
    
    if is_spike == 1
        if y<=param.ne %an excitatory spike
            %decide the postsynaptic neurons
            p=rand(1,param.ne+param.ni);
            p(1:param.ne)=-p(1:param.ne) + param.p_ee; p(param.ne+1:param.ne+param.ni)=-p(param.ne+1:param.ne+param.ni)+param.p_ie;
            p=sign(abs(p)+p);
            s(2,:)=s(2,:)+p;
            res.t_H_e = [res.t_H_e [s(2,y), s(3,y)]'];
        else %an inhibitory spike
            %decide the postsynaptic neurons
            p=rand(1,param.ne+param.ni);
            p(1:param.ne)=-p(1:param.ne)+param.p_ei; p(param.ne+1:param.ne+param.ni)=-p(param.ne+1:param.ne+param.ni)+param.p_ii;
            p=sign(abs(p)+p);
            s(3,:)=s(3,:)+p;
            res.t_H_i = [res.t_H_i [s(2,y), s(3,y)]'];
        end
        res.spike(1,y)=res.spike(1,y)+1; %write down the spike time
        res.spike(res.spike(1,y)+1,y)=t/1000;
    end
    time_delta = time_delta + a(x,y); 

        
    if time_delta >= time_check
        res.t    = [res.t t];
        res.V_e  = [res.V_e, s(1,1:param.ne)];
        res.V_i  = [res.V_i, s(1,param.ne+1: param.ne+ param.ni)];
        res.H_ee = [res.H_ee s(2,1:param.ne)];
        res.H_ei = [res.H_ei s(3,1:param.ne)];
        res.H_ie = [res.H_ie s(2,param.ne+1:param.ne+param.ni)];
        res.H_ii = [res.H_ii s(3,param.ne+1:param.ne+param.ni)];
        res.N_LBE = [res.N_LBE N_LBE];
        res.N_BE = [res.N_BE N_BE];
        res.N_LBI = [res.N_LBI N_LBI];
        res.N_BI = [res.N_BI N_BI];
        time_delta = 0;
    end
end

end

