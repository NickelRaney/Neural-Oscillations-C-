<<<<<<< Updated upstream
save_bool = true;
param1 = param;
param1.factor_Leak = inf;
beta = 30;
alpha = 30;
param1.M = 100*beta;
param1.Mr = 66*beta;
extra_name = 'p=1-M=3000-ns=30-ref=0';
param1.p_ee     = 1;
param1.p_ie     = 1;
param1.p_ei     = 1;
param1.p_ii     = 1;
% param1.tau_ee   = 1.4;
% param1.tau_ie   = 1.2;
% param1.tau_ei   = 4.5;
% param1.tau_ii   = 4.5;
param1.s_ee     = 5*0.15/alpha*beta;
param1.s_ie     = 2*0.5/alpha *beta;
param1.s_ei     = 4.91*0.5/alpha*beta;
param1.s_ii     = 4.91*0.4/alpha*beta;
param1.ns_ee    = alpha;
param1.ns_ie    = alpha;
param1.ns_ei    = alpha;
param1.ns_ii    = alpha;
param1.tau_ri   = 0;
param1.tau_re   = 0;
param1.lambda_e = 7000*beta;
param1.lambda_i = 7000*beta;
bar = 500;
tic;
res = model_L(param1);
toc;

%%
subplot(4,1,1);
    xlim([2000,3000]);
for i =2:4
    subplot(4,1,i);
    xlim([2,3]);
end
=======
for i=1:4
    subplot(4,1,i);
    xlim([1500,2000]);
end
>>>>>>> Stashed changes
