%Effect of I, changing Var
num_n = 10000;
V = zeros(num_n,5);
V(:,1) = randn(num_n,1)*sqrt(10)+20;
V(:,2) = randn(num_n,1)*sqrt(20)+20;
V(:,3) = randn(num_n,1)*sqrt(30)+20;
V(:,4) = randn(num_n,1)*sqrt(40)+20;
V(:,5) = randn(num_n,1)*sqrt(50)+20;
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
