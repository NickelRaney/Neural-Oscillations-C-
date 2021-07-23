%% numerical experiment of pending spikes
%%
spike_t=exprnd(1,1,3000);
var_spike=zeros(1,201);
for t=0.05:0.05:10
    var_spike(round(t/0.05))=var(spike_t<t);
end

var_spike(41:end)=var_spike(41:end)+0.5*var_spike(1:161);
plot(var_spike*400);
hold on
%%
plot(var_spike2);
hold on
%%
init=20*randn(1,3000)+30;
var_init=var(init)
var_spike=zeros(1,200);
spike_t=exprnd(1,400,3000);
spike_t=[spike_t;exprnd(1,200,3000)+2];
for t=0.05:0.05:10
    var_spike(round(t/0.05))=var(init+sum(spike_t<t));
end
plot(var_spike-var_init);
