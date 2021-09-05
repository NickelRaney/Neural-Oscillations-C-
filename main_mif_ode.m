%% Setting paths
addpath('module');

%% ode2

param.ne       = 300;
param.ni       = 100;
param.s_ee     = 5*0.15;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.435;
param.s_ii     = 4.91*0.4;

param.M        = 100;
param.Mr       = 66;
param.lambda_e = 7;
param.lambda_i = 7;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.duration = 1;
param.delta_time = 0.5;
param.dt = 0.01;

%%
tic;
res=model_ode2(param);
toc;

%% ode3
param.ne       = 300;
param.ni       = 100;
param.s_ee     = 5*0.157;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.4;
param.s_ii     = 4.91*0.4;
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;

param.M        = 100;
param.Mr       = 66;
param.lambda_e = 7;
param.lambda_i = 7;
param.tau_ee   = 1.4;
param.tau_ie   = 1.2;
param.tau_ei    = 4.5;
param.tau_ii    = 4.5;
param.duration = 0.0188;
param.delta_time = 0.1;
param.dt = 0.01;
%%
init=[];
tic;    
res=model_ode_complete(param,init);
toc;

%%
s_path='./';
ode_video(param,res,s_path)

%%
t_track=136;
init.peak_e=reshape(res.peak_e(t_track,:),3,10);
init.peak_i=reshape(res.peak_i(t_track,:),3,10);
init.npe=res.npe(t_track);
init.npi=res.npi(t_track);
init.h=res.h(t_track,:);
init.index_e=res.index_e(t_track);
init.index_i=res.index_i(t_track);
init.peak_e(2,1)=257;
init.peak_i(2,2)=191;
%%
tic;
res2=model_ode3(param,init);
toc;
%%
s_path='./';
ode_video(param,res2,s_path)
%%
param.p_ee     = 1;
param.p_ie     = 1;
param.p_ei     = 1;
param.p_ii     = 1;
param.lambda_e = 7000;
param.lambda_i = 7000;
param.tau_re=0;
param.tau_ri=0;
param.s_ee     = 5*0.157;
param.s_ie     = 2*0.5;
param.s_ei     = 4.91*0.4;
param.s_ii     = 4.91*0.4;
param.ns_ee=1;
param.ns_ie=1;
param.ns_ei=1;
param.ns_ii=1;
param.factor_Leak=0;
param.LeakE = 20;
param.LeakI = 16.7;
param.factor_Leak = inf;
param.delta_time = 0.1;

tic;
res_mif=model_L(param);
toc;

%%

s_path='./';
compare_video(param,res,res_mif,s_path)

%%

range=[1:300];
subplot(2,2,1);
plot(sum(res_mif.HE(range,1:300),2))
hold on;
plot(res.h(range,1)*300*100);
xlabel('HEE');
legend('MIF','ode');

subplot(2,2,2);
plot(sum(res_mif.HE(range,301:400),2))
hold on;
plot(res.h(range,2)*100*100);
xlabel('HIE');

subplot(2,2,3);
plot(sum(res_mif.HI(range,1:300),2))
hold on;
plot(res.h(range,3)*300*100);   
xlabel('HEI');

subplot(2,2,4);
plot(sum(res_mif.HI(range,301:400),2))
hold on;
plot(res.h(range,4)*100*100);
xlabel('HII');

%%

mVE_ode = zeros(300,1);
mVI_ode = zeros(300,1);
npe = res.npe;
npi = res.npi;
for i =1:300
    peak_e = res.peak_e(i, :);
    peak_i = res.peak_i(i, :);
    peak_e = reshape(peak_e, 3, 10);
    peak_i = reshape(peak_i, 3, 10);
    me = peak_e(1,1);
    ve = peak_e(2,1);
    re = peak_e(3,1);
    mi = peak_i(1,1);
    vi = peak_i(2,1);
    ri = peak_i(3,1);
    xe = erfinv((2*re-1))*sqrt(2)*sqrt(ve);
    xi = erfinv((2*ri-1))*sqrt(2)*sqrt(vi);
    mVE_ode(i) = -sqrt(ve/2/pi)*exp(-xe^2/2/ve) + me;
    mVI_ode(i) = -sqrt(vi/2/pi)*exp(-xi^2/2/vi) + mi;
    if npe(i)>=2
        mVE_ode(i) = mVE_ode(i)*peak_e(3,1) + sum(peak_e(1,2:npe(i)).*peak_e(3,2:npe(i)));
    end
    if npi(i)>=2
        mVI_ode(i) = mVI_ode(i)*peak_i(3,1) + sum(peak_i(1,2:npi(i)).*peak_i(3,2:npi(i)));
    end
end

range = [1:300];
subplot(1,2,1);
plot(mean(res_mif.VE(range,:),2));
hold on;
plot(mVE_ode);
subplot(1,2,2);
plot(mean(res_mif.VI(range,:),2));
hold on;
plot(mVI_ode);