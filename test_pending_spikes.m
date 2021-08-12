%% Pending E spike
dt = 0.01;
duration = 30;
init_var = 10;
init_mean = 10;
num_neuron = 3000;
HE = 3;
SE = 5;
tauE = 1.4;
d = rand(1, num_neuron)*sqrt(init_var) + init_mean;
var0 = var(d);
m20 = mean(d.^2);
m0 = mean(d);

% Simulation
pending_E = exprnd(tauE, HE, num_neuron);
times = 0:dt:duration;
num_t = size(times,2);
d_varS = zeros(num_t, 1);
d_varS(1) = 0;
d_m2S = zeros(num_t,1);
d_m2S(1)=0;
m_S = zeros(num_t,1);
for i=2:num_t
    bool_index = pending_E > times(i-1) & pending_E <= times(i);
    sum_E = sum(bool_index,1);
    d = d + sum_E*SE;
    d_varS(i) = var(d) - var0;
    d_m2S(i) = mean(d.^2);
    m_S(i) = mean(d);
end

% Estimation
var_temp = var0;
d_varE = zeros(num_t,1);
d_varE(1) = 0;
d_m2E = zeros(num_t,1);
d_m2E(1) = 0;
m2_temp = m20;
m_temp =m0;
m_E = zeros(num_t,1);
for i =2:num_t
    var_temp = var_temp + HE*dt*(-1/tauE*exp(-times(i)/tauE) + 2/tauE*exp(-2*times(i)/tauE))*SE^2;
    d_varE(i) = var_temp - var0;
    gradientmE =HE*1/tauE*exp(-(times(i))/tauE)*SE;
    m_temp = m_temp + dt*gradientmE;
    gradientm2E =  HE*(-1/tauE*exp(-times(i)/tauE) + 2/tauE*exp(-2*times(i)/tauE))*SE^2+...
                   2*m_temp*gradientmE;
    m2_temp = m2_temp + dt*gradientm2E;
    d_m2E(i) =m2_temp;
    m_E(i) = m_temp;

end

subplot(1,2,1);
plot(times, m_S);
hold on;
plot(times, m_E);
xlabel('ms');
ylabel('Mean');
title('Pending E&I effects');

subplot(1,2,2);
plot(times, d_m2S);
hold on;
plot(times, d_m2E);
xlabel('ms');
ylabel('\Delta Var(V)');
title('Pending E effects');
legend('Simulation', 'Estimation');

%% Pending I spikes
dt = 0.0005;
duration = 150;
init_var = 10;
init_mean = 40;
num_neuron = 3000;
HI = 6;
SI = 4.91;
tauI = 4.5;
d = randn(1, num_neuron)*sqrt(init_var) + init_mean;
var0 = var(d);
m0 = mean(d);

% Simulation
pending_I = exprnd(tauI, HI, num_neuron);
pending_I(3:4,:) = pending_I(3:4,:) + 30;
pending_I(5:6,:) = pending_I(5:6,:) + 60;
times = 0:dt:duration;
num_t = size(times,2);
d_varS = zeros(num_t, 1);
d_varS(1) = 0;
m_S = zeros(num_t,1);
m_S(1) = m0;
for i=2:num_t
    bool_index = pending_I > times(i-1) & pending_I <= times(i);
    sum_I = sum(bool_index,1);
    d = (d+66).*(1-SI/166).^(sum_I)-66;
    d_varS(i) = var(d) - var0;
    m_S(i) = mean(d);
end

% Estimation
var_temp = var0;
m_temp = m0+66;
m2_temp = m_temp^2 + var_temp;
d_varE = zeros(num_t,1);
d_varE(1) = 0;
m_E = zeros(num_t,1);
m_E(1) = m0;
for i =2:num_t
    i
    t = times(i);
    %m_temp = m_temp + HI*(m_temp +66)*(-1/tauI*exp(-times(i)/tauI))*SI/166*dt;
%     a = 1+HI*dt*((1/tauI*exp(-times(i)/tauI))*(1-SI/166)^2 -1/tauI*exp(-times(i)/tauI));
%     b = HI*dt*(2/tauI*exp(-2*times(i)/tauI) - 1/tauI*exp(-t/tauI))*(SI/166)^2*(10+66)^2;
%     var_temp = (var_temp - b/(1-a))^HI + b/(1-a);
    eff = SI/166;
     gradient = -2*1/tauI*exp(-t/tauI)*eff/(1-(1-exp(-t/tauI))*eff);
     gradient = gradient -(t>30)* 2* 1/tauI*exp(-(t-30)/tauI)*eff/(1-(1-exp(-(t-30)/tauI))*eff);
     gradient = gradient -(t>60)* 2* 1/tauI*exp(-(t-60)/tauI)*eff/(1-(1-exp(-(t-60)/tauI))*eff);
     m_temp = exp(log(m_temp) + dt*gradient);
     m_E(i) = m_temp - 66;
     gradient1 = 2*(1/tauI*exp(-t/tauI)*(1-eff)^2 - 1/tauI*exp(-t/tauI))/((1-exp(-t/tauI))*(1-eff)^2 + exp(-t/tauI));
     gradient1 = gradient1 + (t>30)* 2*(1/tauI*exp(-(t-30)/tauI)*(1-eff)^2 - 1/tauI*exp(-(t-30)/tauI))/((1-exp(-(t-30)/tauI))*(1-eff)^2 + exp(-(t-30)/tauI));
     gradient1 = gradient1 + (t>60)* 2*(1/tauI*exp(-(t-60)/tauI)*(1-eff)^2 - 1/tauI*exp(-(t-60)/tauI))/((1-exp(-(t-60)/tauI))*(1-eff)^2 + exp(-(t-60)/tauI));
     m2_temp = exp(log(m2_temp) + dt*gradient1);
    d_varE(i) = m2_temp -m_temp^2 - var0;
end

%%
subplot(1,2,1);
plot(times, m_S);
hold on;
plot(times, m_E);
xlabel('ms');
ylabel('Mean');
title('Pending I effects');

subplot(1,2,2);
plot(times, d_varS);
hold on;
plot(times, d_varE);
xlabel('ms');
ylabel('\Delta Var(V)');
title('Pending I effects');
legend('Simulation', 'Estimation');

%% Pending E & I spikes
dt = 0.001;
duration = 150;
init_var = 10;
init_mean = 10;
num_neuron = 3000;
HE = [20 20 30 40];
tE = [0 10 25 40];
HI = [10 50 60];
tI = [10 20 30];
SI = 4.91;
SE = 5;
tauI = 4.5;
tauE = 1.4;
d = randn(1, num_neuron)*sqrt(init_var) + init_mean;
var0 = var(d);
m0 = mean(d);

% Simulation
pending_E = exprnd(tauE, sum(HE), num_neuron);
pending_I = exprnd(tauI, sum(HI), num_neuron);
indexE = 1;
indexI = 1;
for i =1:size(HE,2)
    pending_E(indexE:indexE+HE(i)-1,:) = pending_E(indexE:indexE+HE(i)-1,:) + tE(i);
    indexE = indexE + HE(i);
end
for i =1:size(HI,2)
    pending_I(indexI:indexI+HI(i)-1,:) = pending_I(indexI:indexI+HI(i)-1,:) + tI(i);
    indexI = indexI + HI(i);
end

for j =1:size(HI)
end
% pending_I(3:4,:) = pending_I(3:4,:) + 30;
% pending_I(5:6,:) = pending_I(5:6,:) + 60;
times = 0:dt:duration;
num_t = size(times,2);
d_varS = zeros(num_t, 1);
d_varS(1) = 0;
m_S = zeros(num_t,1);
m_S(1) = m0;

%%
for i=2:num_t
    bool_index = pending_I > times(i-1) & pending_I <= times(i);
    sum_I = sum(bool_index,1);
    d = (d+66).*(1-SI/166).^(sum_I)-66;
    bool_index = pending_E > times(i-1) & pending_E <= times(i);
    sum_E = sum(bool_index,1);
    d = d + SE * sum_E;
    d_varS(i) = var(d) - var0;
    m_S(i) = mean(d);
end

%% Estimation
var_temp = var0;
m_temp = m0+66;
m2_temp = m_temp^2 + var_temp;
d_varE = zeros(num_t,1);
d_varE(1) = 0;
m_E = zeros(num_t,1);
m_E(1) = m0;
for i =2:num_t
    i
    t = times(i);
    %m_temp = m_temp + HI*(m_temp +66)*(-1/tauI*exp(-times(i)/tauI))*SI/166*dt;
%     a = 1+HI*dt*((1/tauI*exp(-times(i)/tauI))*(1-SI/166)^2 -1/tauI*exp(-times(i)/tauI));
%     b = HI*dt*(2/tauI*exp(-2*times(i)/tauI) - 1/tauI*exp(-t/tauI))*(SI/166)^2*(10+66)^2;
%     var_temp = (var_temp - b/(1-a))^HI + b/(1-a);
     eff = SI/166;
     gradientmI = 0;
     gradientmE = 0;
     gradientm2I = 0;
     gradientm2E = 0;
     for j =1:size(HE,2)
              gradientmE = gradientmE+ (t-tE(j)>0)*HE(j)*1/tauE*exp(-(t-tE(j))/tauE)*SE;
              gradientm2E = gradientm2E + (t-tE(j)>0)*HE(j)*(-1/tauE*exp(-(t-tE(j))/tauE) + 2/tauE*exp(-2*(t-tE(j))/tauE))*SE^2;
     end
     
     for j =1:size(HI,2)
              gradientmI = gradientmI+ (t-tI(j)>0)*m_temp*HI(j)*(-1/tauI*exp(-(t-tI(j))/tauI)*eff)/(1-(1-exp(-(t-tI(j))/tauI))*eff);
             gradientm2I =gradientm2I+ (t-tI(j)>0)*m2_temp*HI(j)*(1/tauI*exp(-(t-tI(j))/tauI)*(1-eff)^2 - 1/tauI*exp(-(t-tI(j))/tauI))/((1-exp(-(t-tI(j))/tauI))*(1-eff)^2 + exp(-(t-tI(j))/tauI));
     end
    gradientm2E =  gradientm2E + 2*m_temp*gradientmE;
     m_temp = m_temp + dt*(gradientmI + gradientmE);
     m2_temp = m2_temp + dt*(gradientm2I + gradientm2E);
     d_varE(i) = m2_temp -m_temp^2 - var0;
     m_E(i) = m_temp - 66;
end

%%
subplot(1,2,1);
plot(times, m_S);
hold on;
plot(times, m_E);
xlabel('ms');
ylabel('Mean');
title('Pending E&I effects');

subplot(1,2,2);
plot(times, d_varS);
hold on;
plot(times, d_varE);
xlabel('ms');
ylabel('\Delta Var(V)');
title('Pending E&I effects');
legend('Simulation', 'Estimation');

