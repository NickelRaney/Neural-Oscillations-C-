%% Parameters
ne = 300;
ni = 100;
ve = zeros(1, ne);
vi = zeros(1, ni);
he = ones(2, ne);
he(1,:) = he(1,:) * 400;
he(2,:) = he(2,:) * 300;
hi = ones(2, ni);
hi(1,:) = hi(1,:) * 300;
hi(2,:) = hi(2,:) *200;
tau_ee = 1.4;
tau_ie = 1.2;
tau_ei = 4.5;
tau_ii = 4.5;
lambda_e = 7;
lambda_i = 7;
Mr = 66;
M = 100;
s_ee1 = 5*0.15;
s_ie1 = 2*0.5;
s_ei1 = 4.91*0.42;
s_ii1 = 4.91*0.4;
dt = 0.001;

%% LIF
duration = 100;
var_lif = zeros(ceil(duration/dt)+1, 2);
mean_lif = zeros(ceil(duration/dt)+1, 2);
he1 = he;
hi1 = hi;
for step=2:duration/dt+1
    ve=ve+(s_ee1*he1(1,:)/tau_ee-s_ei1*he1(2,:)/tau_ei.*(ve+Mr)/(M+Mr))*dt;
    vi=vi+(s_ie1*hi1(1,:)/tau_ie-s_ii1*hi1(2,:)/tau_ii.*(vi+Mr)/(M+Mr))*dt;
    he1=he1.*[exp(-dt/tau_ee);exp(-dt/tau_ei)];
    hi1=hi1.*[exp(-dt/tau_ie);exp(-dt/tau_ii)];
    var_lif(step,1) = var(ve);
    var_lif(step,2) = var(vi);
    mean_lif(step,1) = mean(ve);
    mean_lif(step,2) = mean(vi);
end

%% MIF
VE = zeros(1,ne);
VI = zeros(1,ni);
duration =1000;
HEE = zeros(1,500000);
HIE = zeros(1,500000);
HEI = zeros(1,500000);
HII = zeros(1,500000);
alpha = 1;
s_ee = s_ee1 / alpha;
s_ei = s_ei1 / alpha;
s_ie = s_ie1 / alpha;
s_ii = s_ii1 / alpha;
var_mif1 = zeros(ceil(duration/dt)+1,2);
mean_mif1 = zeros(ceil(duration/dt)+1,2);

for i =1:ne
    index_ee = HEE(1);
    index_ei = HEI(1);
    HEE(index_ee+2: index_ee+he(1,i)*alpha+1) = i;
    HEI(index_ei+2: index_ei+he(2,i)*alpha+1) = i;
    HEE(1) = HEE(1)+he(1,i)*alpha;
    HEI(1) = HEI(1)+he(2,i)*alpha;
end

for i =1:ni
    index_ie = HIE(1);
    index_ii = HII(1);
    HIE(index_ie+2: index_ie+hi(1,i)*alpha+1) = i;
    HII(index_ii+2: index_ii+hi(2,i)*alpha+1) = i;
    HIE(1) = HIE(1) + hi(1,i)*alpha;
    HII(1) = HII(1) + hi(2,i)*alpha;
end

t=0;
t_ee = 1/tau_ee;
t_ei = 1/tau_ei;
t_ie = 1/tau_ie;
t_ii = 1/tau_ii;
c(1) = ne*lambda_e*1;
c(2) = ne*lambda_i*1;
c(3) = HEE(1)*t_ee;
c(4) = HEI(1)*t_ei;
c(5) = HIE(1)*t_ie;
c(6) = HIE(1)*t_ii;
t_count = 0;
step = 2;

while t< duration
    delta = log(1-rand(1))/sum(c);
    t = t-delta;
    t_count = t_count -delta;
    if t_count >dt
        var_mif1(step,1) = var(VE);
        var_mif1(step,2) = var(VI);
        mean_mif1(step,1) = mean(VE);
        mean_mif1(step,2) = mean(VI);
        step = step + 1;
        t_count = 0;
    end 
    index = find_index(c);
    switch index
        case 1
            whichhit = ceil(rand(1)*ne);
             VE(whichhit) = VE(whichhit) ;
        case 2
            whichhit = ceil(rand(1)*ni);
            VI(whichhit) = VI(whichhit) ;
        case 3
            index = ceil(rand(1)*HEE(1))+1;
            whichhit = HEE(index);
            HEE(index) = HEE(HEE(1)+1);
            HEE(HEE(1)+1) = 0;
            HEE(1)=HEE(1)-1;
            VE(whichhit) = VE(whichhit) + s_ee;
            c(3) = c(3)-t_ee;
       case 5
            index = ceil(rand(1)*HIE(1))+1;
            whichhit = HIE(index);
            HIE(index) = HIE(HIE(1)+1);
            HIE(HIE(1)+1) = 0;
            HIE(1)=HIE(1)-1;
            VI(whichhit) = VI(whichhit) + s_ie;
            c(4) = c(4)-tau_ie;
        case 4
            index = ceil(rand(1)*HEI(1))+1;
            whichhit = HEI(index);
            HEI(index) = HEI(HEI(1)+1);
            HEI(HEI(1)+1) = 0;
            HEI(1)=HEI(1)-1;
            VE(whichhit) = VE(whichhit) - (VE(whichhit)+Mr)*s_ei/(M + Mr);
            if VE(whichhit) < -Mr
                    VE(whichhit) = -Mr;
            end
            c(5) = c(5)-tau_ei;
        case 6
            index = ceil(rand(1)*HII(1))+1;
            whichhit = HII(index);
            HII(index) = HII(HII(1)+1);
            HII(HII(1)+1) = 0;
            HII(1)=HII(1)-1;
            VI(whichhit) = VI(whichhit) - (VI(whichhit)+Mr)*s_ii/(M + Mr);
            if VI(whichhit) < -Mr
                    VI(whichhit) = -Mr;
            end
            c(6) = c(6)-tau_ii;
    end
end

%% MIF alpha =100
duaration = 1000;
VE = zeros(1,ne);
VI = zeros(1,ni);
duration =100;
HEE = zeros(1,500000);
HIE = zeros(1,500000);
HEI = zeros(1,500000);
HII = zeros(1,500000);
alpha = 200;
s_ee = s_ee1 / alpha;
s_ei = s_ei1 / alpha;
s_ie = s_ie1 / alpha;
s_ii = s_ii1 / alpha;
var_mif2 = zeros(ceil(duration/dt)+1,2);
mean_mif2 = zeros(ceil(duration/dt)+1,2);
t_mif2 = zeros(ceil(duration/dt)+1,1);
for i =1:ne
    index_ee = HEE(1);
    index_ei = HEI(1);
    HEE(index_ee+2: index_ee+he(1,i)*alpha+1) = i;
    HEI(index_ei+2: index_ei+he(2,i)*alpha+1) = i;
    HEE(1) = HEE(1)+he(1,i)*alpha;
    HEI(1) = HEI(1)+he(2,i)*alpha;
end

for i =1:ni
    index_ie = HIE(1);
    index_ii = HII(1);
    HIE(index_ie+2: index_ie+hi(1,i)*alpha+1) = i;
    HII(index_ii+2: index_ii+hi(2,i)*alpha+1) = i;
    HIE(1) = HIE(1) + hi(1,i)*alpha;
    HII(1) = HII(1) + hi(2,i)*alpha;
end

t=0;
t_ee = 1/tau_ee;
t_ei = 1/tau_ei;
t_ie = 1/tau_ie;
t_ii = 1/tau_ii;
c(1) = ne*lambda_e*1;
c(2) = ne*lambda_i*1;
c(3) = HEE(1)*t_ee;
c(4) = HEI(1)*t_ei;
c(5) = HIE(1)*t_ie;
c(6) = HII(1)*t_ii;
t_count = 0;
step = 2;


while t< duration
    delta = log(1-rand(1))/sum(c);
    t = t-delta;
    t_count = t_count -delta;
    if t_count >dt
        var_mif2(step,1) = var(VE);
        var_mif2(step,2) = var(VI);
        mean_mif2(step,1) = mean(VE);
        mean_mif2(step,2) = mean(VI);
        t_mif2(step) = t;
        step = step + 1;
        t_count = 0;
    end 
    index = find_index(c);
    switch index
        case 1
            whichhit = ceil(rand(1)*ne);
             VE(whichhit) = VE(whichhit) ;
        case 2
            whichhit = ceil(rand(1)*ni);
            VI(whichhit) = VI(whichhit)  ;
        case 3
            index = ceil(rand(1)*HEE(1))+1;
            whichhit = HEE(index);
            HEE(index) = HEE(HEE(1)+1);
            HEE(HEE(1)+1) = 0;
            HEE(1)=HEE(1)-1;
            VE(whichhit) = VE(whichhit) + s_ee;
            c(3) = c(3)-t_ee;
       case 5
            index = ceil(rand(1)*HIE(1))+1;
            whichhit = HIE(index);
            HIE(index) = HIE(HIE(1)+1);
            HIE(HIE(1)+1) = 0;
            HIE(1)=HIE(1)-1;
            VI(whichhit) = VI(whichhit) + s_ie;
            c(5) = c(5)-t_ie;
        case 4
            index = ceil(rand(1)*HEI(1))+1;
            whichhit = HEI(index);
            HEI(index) = HEI(HEI(1)+1);
            HEI(HEI(1)+1) = 0;
            HEI(1)=HEI(1)-1;
            VE(whichhit) = VE(whichhit) - (VE(whichhit)+Mr)*s_ei/(M + Mr);
            if VE(whichhit) < -Mr
                    VE(whichhit) = -Mr;
            end
            c(4) = c(4)-t_ei;
        case 6
            index = ceil(rand(1)*HII(1))+1;
            whichhit = HII(index);
            HII(index) = HII(HII(1)+1);
            HII(HII(1)+1) = 0;
            HII(1)=HII(1)-1;
            VI(whichhit) = VI(whichhit) - (VI(whichhit)+Mr)*s_ii/(M + Mr);
            if VI(whichhit) < -Mr
                    VI(whichhit) = -Mr;
            end
            c(6) = c(6)-t_ii;
    end
end
mean_mif2 = mean_mif2(1:step-1,:);
var_mif2 = var_mif2(1:step-1,:);
t_mif2 = t_mif2(1:step-1,:);


%%
t= 0:dt:duration;
subplot(2,2,1);
plot(t, mean_lif(:,1));
hold on;
plot(t_mif2, mean_mif2(:,1));
legend('LIF', 'MIF alpha=200');
xlabel('Time (ms)');
ylabel('mean(VE)');
title('VE');
subplot(2,2,2);
plot(t, mean_lif(:,2));
hold on;
plot(t_mif2, mean_mif2(:,2));
legend('LIF', 'MIF alpha=200');
xlabel('Time (ms)');
ylabel('mean(VI)');
title('VI');

t= 0:dt:duration;

subplot(2,2,3);
plot(t, var_lif(:,1));
hold on;
plot(t_mif2, var_mif2(:,1));
legend('LIF','MIF alpha=200');
ylabel('var(VE)');
xlabel('Time (ms)');

subplot(2,2,4);
plot(t, var_lif(:,2));
hold on;
plot(t_mif2, var_mif2(:,2));
legend('LIF','MIF alpha=200');
ylabel('var(VI)');
xlabel('Time (ms)');
sgtitle('MIF vs LIF');