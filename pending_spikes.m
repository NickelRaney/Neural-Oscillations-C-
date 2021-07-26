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
duration = 100;
init_var = 10;
init_mean = 10;
num_neuron = 3000;
HI = 1;
SI = 4.91;
tauI = 4.5;
d = randn(1, num_neuron)*sqrt(init_var) + init_mean;
var0 = var(d);
m0 = mean(d);

% Simulation
pending_I = exprnd(tauI, HI, num_neuron);
times = 0:dt:duration;
num_t = size(times,2);
d_varS = zeros(num_t, 1);
d_varS(1) = 0;
for i=2:num_t
    bool_index = pending_I > times(i-1) & pending_I <= times(i);
    sum_I = sum(bool_index,1);
    d = (d+66).*(1-SI/166).^(sum_I)-66;
    d_varS(i) = var(d) - var0;
end

% Estimation
var_temp = var0;
m_temp = m0;
d_varE = zeros(num_t,1);
d_varE(1) = 0;
for i =2:num_t
    t = times(i);
    m_temp = m_temp + HI*(m_temp +66)*(-1/tauI*exp(-times(i)/tauI))*SI/166*dt;
    a = 1+HI*dt*((1/tauI*exp(-times(i)/tauI))*(1-SI/166)^2 -1/tauI*exp(-times(i)/tauI));
    b = HI*dt*(2/tauI*exp(-2*times(i)/tauI) - 1/tauI*exp(-t/tauI))*(SI/166)^2*(10+66)^2;
    var_temp = (var_temp - b/(1-a))^HI + b/(1-a);
%     var_temp = var_temp + HI*dt*(var_temp*((1/tauI*exp(-times(i)/tauI))*(1-SI/166)^2 -1/tauI*exp(-times(i)/tauI))+...
%         (2/tauI*exp(-2*times(i)/tauI) - 1/tauI*exp(-t/tauI))*(SI/166)^2*(10+66)^2);
    d_varE(i) = var_temp - var0;
end


plot(times, d_varS);
hold on;
plot(times, d_varE);
xlabel('ms');
ylabel('\Delta Var(V)');
title('Pending I effects');
legend('Simulation', 'Estimation');