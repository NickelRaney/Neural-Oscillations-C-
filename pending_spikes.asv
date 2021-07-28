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

% Simulation
pending_E = exprnd(tauE, HE, num_neuron);
times = 0:dt:duration;
num_t = size(times,2);
d_varS = zeros(num_t, 1);
d_varS(1) = 0;
for i=2:num_t
    bool_index = pending_E > times(i-1) & pending_E <= times(i);
    sum_E = sum(bool_index,1);
    d  = d + sum_E*SE;
    d_varS(i) = var(d) - var0;
end

% Estimation
var_temp = var0;
d_varE = zeros(num_t,1);
d_varE(1) = 0;
for i =2:num_t
    var_temp = var_temp + HE*dt*(-1/tauE*exp(-times(i)/tauE) + 2/tauE*exp(-2*times(i)/tauE))*SE^2;
    d_varE(i) = var_temp - var0;
end


plot(times, d_varS);
hold on;
plot(times, d_varE);
xlabel('ms');
ylabel('\Delta Var(V)');
title('Pending E effects');
legend('Simulation', 'Estimation');

%% Pending I spikes
dt = 0.001;
duration = 150;
init_var = 40;
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

%%

