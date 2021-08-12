function [res] = model_L(param)
%this is the new version of model_full imitating cpp code
ne = param.ne;
ni = param.ni;
p_ee = param.p_ee;
p_ie = param.p_ie;
p_ei = param.p_ei;
p_ii = param.p_ii;
s_ee = param.s_ee;
s_ie = param.s_ie;
s_ei = param.s_ei;
s_ii = param.s_ii;
tau_ee = 1000/param.tau_ee;
tau_ie = 1000/param.tau_ie;
tau_ei = 1000/param.tau_ei;
tau_ii = 1000/param.tau_ii;
tau_ri = 1000/param.tau_ri;
tau_re = 1000/param.tau_re;
gLeak = -1 + (1 / 2.71828)^(1 / param.factor_Leak);
M = param.M;
Mr = param.Mr;
duration_time=param.duration;
c = zeros(1, 10);
c(1) = ne * param.lambda_e;
c(2) = ni * param.lambda_i;
if param.factor_Leak==inf
    c(3)=0;
    c(4)=0;
else
c(3) = ne * 1000 * param.factor_Leak/ param.LeakE;
c(4) = ni * 1000 * param.factor_Leak/ param.LeakI;
end
delta_time = param.delta_time/1000;
%c is the clock variable, with elements representing 1.E drive 2.I drive
%3.E leak 4.I leak 5.HEE 6.HEI 7.HIE 8.HII 9.E ref 10.I ref

res.spike=zeros(2000,ne+ni);
VE = zeros(1,ne);
VI = zeros(1,ni);
HEE = zeros(1,100000);
HIE = zeros(1,100000);
HEI = zeros(1,100000);
HII = zeros(1,100000);

Eref = zeros(1, ne+1);
Iref = zeros(1, ni+1);

awakeE = ones(1, ne);
awakeI = ones(1, ni);
spikeE = 0;
spikeI = 0;
t = 0;
t_temp  = t;
res.spike = zeros(ne+ni, 100000);
res.HE = zeros(ceil(duration_time/delta_time)+1,ne+ni);
res.HI = zeros(ceil(duration_time/delta_time)+1,ne+ni);
res.VE = zeros(ceil(duration_time/delta_time)+1,ne);
res.VI = zeros(ceil(duration_time/delta_time)+1,ni);
res.time =[];
index1 = 1;
while t < duration_time
    delta = log(1-rand(1)) / sum(c);
    t = t-delta;
    t_temp = t_temp - delta;
    if t_temp > delta_time
        HEE_temp = HEE(2: HEE(1)+1);
        HIE_temp = HIE(2: HIE(1)+1);
        HEI_temp = HEI(2: HEI(1)+1);
        HII_temp = HII(2: HII(1)+1);
        res.HE(index1,1:(ne)) = histcounts(HEE_temp,1:(ne+1));
        res.HE(index1,(ne+1):(ne+ni)) = histcounts(HIE_temp,1:(ni+1));
        res.HI(index1,1:ne) = histcounts(HEI_temp,1:(ne+1));
        res.HI(index1,(ne+1):(ne+ni)) = histcounts(HII_temp,1:(ni+1));
        res.VE(index1,:) = VE - (Mr+1)*(1-awakeE);
        res.VI(index1,:) = VI - (Mr+1)*(1-awakeI);
        res.time = [res.time, t*1000];
        index1  = index1 +1;
        t_temp = 0;
    end
    
    index = find_index(c);      %type of events
    
    switch index
        case 1
            whichhit = ceil(rand(1)*ne);
            if awakeE(whichhit) == 1
                VE(whichhit) = VE(whichhit) + 1;
                if VE(whichhit) >= M
                    spikeE = 1;
                end
            end
        case 2
            whichhit = ceil(rand(1)*ni);
            if awakeI(whichhit) == 1
                VI(whichhit) = VI(whichhit) + 1;
                if VI(whichhit) >= M
                    spikeI = 1;
                end
            end
        case 3
            whichhit = ceil(rand(1)*ne);
            if awakeE(whichhit) == 1
                Leak = gLeak * VE(whichhit);
                VE(whichhit) = VE(whichhit) + real2int(Leak);
            end
        case 4
            whichhit = ceil(rand(1)*ni);
            if awakeI(whichhit) == 1
                Leak = gLeak * VI(whichhit);
                VI(whichhit) = VI(whichhit) + real2int(Leak);
            end
        case 5
            index = ceil(rand(1)*HEE(1))+1;
            whichhit = HEE(index);
            HEE(index) = HEE(HEE(1)+1);
            HEE(HEE(1)+1) = 0;
            HEE(1)=HEE(1)-1;
            if awakeE(whichhit) == 1
                VE(whichhit) = VE(whichhit) + real2int(s_ee);
                if VE(whichhit) >= M
                    spikeE = 1;
                end
            end
            c(5) = c(5)-tau_ee;
        case 6
            index = ceil(rand(1)*HIE(1))+1;
            whichhit = HIE(index);
            HIE(index) = HIE(HIE(1)+1);
            HIE(HIE(1)+1) = 0;
            HIE(1)=HIE(1)-1;
            if awakeI(whichhit) == 1
                VI(whichhit) = VI(whichhit) + real2int(s_ie);
                if VI(whichhit) >= M
                    spikeI = 1;
                end
            end
            c(6) = c(6)-tau_ie;
        case 7
            index = ceil(rand(1)*HEI(1))+1;
            whichhit = HEI(index);
            HEI(index) = HEI(HEI(1)+1);
            HEI(HEI(1)+1) = 0;
            HEI(1)=HEI(1)-1;
            if awakeE(whichhit) == 1
                VE(whichhit) = VE(whichhit) - real2int((VE(whichhit)+Mr)*s_ei/(M + Mr));
                if VE(whichhit) < -Mr
                    VE(whichhit) = Mr;
                end
            end
            c(7) = c(7)-tau_ei;
        case 8
            index = ceil(rand(1)*HII(1))+1;
            whichhit = HII(index);
            HII(index) = HII(HII(1)+1);
            HII(HII(1)+1) = 0;
            HII(1)=HII(1)-1;
            if awakeI(whichhit) == 1
                VI(whichhit) = VI(whichhit) - real2int((VI(whichhit)+Mr)*s_ii/(M + Mr));
                if VI(whichhit) < -Mr
                    VI(whichhit) = Mr;
                end
            end
            c(8) = c(8)-tau_ii;
        case 9
            index = ceil(rand(1)*Eref(1))+1;
            whichhit = Eref(index);
            Eref(index) = Eref(Eref(1)+1);
            Eref(Eref(1)+1) = 0;
            Eref(1)=Eref(1)-1;
            awakeE(whichhit) = 1;
            c(9) = c(9) - tau_re;
        case 10
            index = ceil(rand(1)*Iref(1))+1;
            whichhit = Iref(index);
            Iref(index) = Iref(Iref(1)+1);
            Iref(Iref(1)+1) = 0;
            Iref(1)=Iref(1)-1;
            awakeI(whichhit) = 1;
            c(10) = c(10) - tau_ri;
    end
    if spikeE == 1
        spikeE = 0;
        VE(whichhit) = 0;
        awakeE(whichhit) = 0;
        Eref(Eref(1)+2)=whichhit;
        Eref(1) = Eref(1)+1;
        c(9) = c(9) + tau_re;
        p = -rand(1, ne) + p_ee;
        HEE_new = find(p>0);
        HEE(HEE(1)+2:HEE(1)+1+length(HEE_new)) = HEE_new;
        HEE(1)=HEE(1)+length(HEE_new);
        p = -rand(1, ni) + p_ie;
        HIE_new = find(p>0);
        HIE(HIE(1)+2:HIE(1)+1+length(HIE_new)) = HIE_new;
        HIE(1)=HIE(1)+length(HIE_new);
        c(5) = tau_ee * HEE(1);
        c(6) = tau_ie* HIE(1);
        res.spike(1,whichhit)=res.spike(1,whichhit)+1; %write down the spike time
        res.spike(res.spike(1,whichhit)+1,whichhit)=t;
    end
    if spikeI == 1
        spikeI = 0;
        VI(whichhit) = 0;
        awakeI(whichhit) = 0;
        Iref(Iref(1)+2)=whichhit;
        Iref(1) = Iref(1)+1;
        c(10) = c(10) + tau_ri;
        p = -rand(1, ne) + p_ei;
        HEI_new = find(p>0);
        HEI(HEI(1)+2:HEI(1)+1+length(HEI_new)) = HEI_new;
        HEI(1)=HEI(1)+length(HEI_new);
        p = -rand(1, ni) + p_ii;
        HII_new = find(p>0);
        HII(HII(1)+2:HII(1)+1+length(HII_new)) = HII_new;
        HII(1)=HII(1)+length(HII_new);
        c(7) = tau_ei * HEI(1);
        c(8) = tau_ii * HII(1);
        res.spike(1,whichhit+ne)=res.spike(1,whichhit+ne)+1; %write down the spike time
        res.spike(res.spike(1,whichhit+ne)+1,whichhit+ne)=t;
    end
end
res.HE = res.HE(1:index1-1,:);
res.HI = res.HI(1:index1-1,:);
res.VE = res.VE(1:index1-1,:);
res.VI = res.VI(1:index1-1,:);
end


function index = find_index(c)
l=length(c);
s = sum(c);
index = 1;
temp = rand(1);
frac = c(1) / s;
while frac < temp && index <= l
    index = index + 1;
    frac = frac + c(index) / s;
end
end


function a = real2int(b)
f = floor(b);
if rand(1) < b - f
    a = f + 1;
else
    a = f;
end
end