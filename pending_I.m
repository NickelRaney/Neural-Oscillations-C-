%Effect of I, changing Var
num_n = 10000;
V = zeros(num_n,5);
V(:,1) = randn(num_n,1)*sqrt(10)+20;
V(:,2) = randn(num_n,1)*sqrt(20)+20;
V(:,3) = randn(num_n,1)*sqrt(30)+20;
V(:,4) = randn(num_n,1)*sqrt(40)+20;
V(:,5) = randn(num_n,1)*sqrt(50)+20;
init_var = var(V,1);
pending_I = [exprnd(4.5,num_n,5),exprnd(4.5,num_n,5),exprnd(4.5,num_n,5)+4.5];
pending_E = [exprnd(1,num_n,5)+1];
% pending_I = [exprnd(4.5,num_n,5)];

t_max = max(max(pending_I));
t = 0:0.01:20;
d_var = zeros(size(t,2),5);
d_var(1,:) = 0;
for i = 2:size(t,2)
    for j =1:5
        for k=1:3
            index = pending_I(:,j+(k-1)*5) > t(i-1) & pending_I(:,j+(k-1)*5)<= t(i);
            V(index, j) = V(index, j) - 2.*(V(index,j)+66)/166;
        end
        for k=1
            index = pending_E(:,j+(k-1)*5) > t(i-1) & pending_E(:,j+(k-1)*5)<= t(i);
            V(index, j) = V(index, j) +1;
        end
        d_var(i,j) = var(V(:,j))- init_var(j);
    end
end
plot(t, d_var(:,1));
hold on;
plot(t, d_var(:,2));
hold on;
plot(t, d_var(:,3));
hold on;
plot(t, d_var(:,4));
hold on;
plot(t, d_var(:,5));
xlabel('ms');
ylabel('\Delta V');
legend('m=20 v=10','m=20 v=20','m=20 v=30','m=20 v=40','m=20 v=50');
title('Effect of pending I spike on variance');
saveas(gcf,'figure/ode/EffectOfI1.png');
saveas(gcf,'figure/ode/EffectofI1.fig');


%% ODE simulation
x=V(:,1);
v=zeros(5,1);
dv=zeros(5,1);
m=zeros(5,1);
dm=zeros(5,1);
v(1)=var(x);
m(1)=mean(x);
for t=0.01:0.01:1
    v(1:3)=v(1:3)+dv(1:3)*0.01;
    m(1:3)=m(1:3)+dm(1:3)*0.01;
    c1=f(t,0,4.5);
    dv(3)=c1(1)*dv(2)+c1(2)*v(2)+2*c1(3)*(m(2)+66)*dm(2)...
        +c1(4)*(m(2)+66).^2;
    dm(3)=c1(5)*dm(2)+c1(6)*m(2)-exp(t/4.5)/4.5*2/166*66;
    
end



%% Effect of I, changing mean
num_n = 20000;
V = zeros(num_n,5);
V(:,1) = randn(num_n,1)*sqrt(100)+10;
V(:,2) = randn(num_n,1)*sqrt(100)+20;
V(:,3) = randn(num_n,1)*sqrt(100)+30;
V(:,4) = randn(num_n,1)*sqrt(100)+40;
V(:,5) = randn(num_n,1)*sqrt(100)+50;
init_var = var(V,1);
pending_I = exprnd(4.5,num_n,5);
t_max = max(max(pending_I));
t = 0:0.01:t_max;
d_var = zeros(size(t,2),5);
d_var(1,:) = 0;
for i = 2:size(t,2)
    for j =1:5
    index = pending_I(:,j) > t(i-1) & pending_I(:,j)<= t(i);
    V(index, j) = V(index, j) - (0.4*4.91).*(V(index,j)+66)/166;
    d_var(i,j) = var(V(:,j))- init_var(j);
    end
end
figure;
plot(t, d_var(:,1));
hold on;
plot(t, d_var(:,2));
hold on;
plot(t, d_var(:,3));
hold on;
plot(t, d_var(:,4));
hold on;
plot(t, d_var(:,5));
xlabel('ms');
ylabel('\Delta V');
legend('m=10 v=100','m=20 v=100','m=30 v=100','m=40 v=100','m=50 v=100');
title('Effect of pending I spike on variance');
saveas(gcf,'figure/ode/EffectOfI2.png');
saveas(gcf,'figure/ode/EffectofI2.fig');

%%
init_var = 10;
m = 20;
num_n = 3000;
HI = 5;
tauI = 4.5;
dt = 0.01;
SI = 4.91;

d = randn(1,num_n)*sqrt(init_var) + m;
init_V = var(d);

% Real simulation
pending_I = exprnd(tauI, HI, num_n);
t = max(max(pending_I));
time_points = 0:dt:t;
size_t =size(time_points,2);
delta_V_real = zeros(size_t,1);
delta_V_real(1) = 0;
for i =2:size_t
    bool_temp = pending_I > time_points(i-1) & pending_I <=time_points(i);
    power_index = sum(bool_temp,1);
    d = (d+66).*(1-SI/166).^(power_index)-66;
    V = var(d);
    delta_V_real(i) = V-init_V; 
end

% Estimation
V=init_V;
delta_V_esti = zeros(size_t,1);
for i =2:size_t
    V = V*(1-2*HI*SI*dt/(166*tauI));
    HI = HI - HI*dt/tauI;
    delta_V_esti(i) = V - init_V;
end

plot(time_points, delta_V_real);
hold on;
plot(time_points, delta_V_esti);
xlabel('ms');
ylabel('\Delta V');
title('Simulation versus Estimation');
legend('Simulation','Estimation');


%%

function c=f(t,tn,tau)
    
    a=2/166;
    b=1-exp(-t/tau);
    f=b*(1-a)^2+1-b;
    g=b(1-b)*a^2;
    h=1-b*a;
    df=(1-a)^2/tau*exp(-(t-tn)/tau)-exp(-(t-tn)/tau)/tau;
    dg=a^2*(2/tau*exp(-2*(t-tn)/tau)-1/tau*exp(-(t-tn)/tau));
    dh=a*exp(-(t-tn)/tau)/tau;
    c=[f,df,g,dg,h,dh];
end
