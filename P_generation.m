% Generate P from C++ code
bar.e            = 40;
bar.i            = 40;
% param.ne         = 300;
% param.ni         = 100;
% param.s_ee       = 5;
% param.s_ie       = 2;
% param.s_ei       = 5;
% param.s_ii       = 5;
param.ne         = 75;
param.ni         = 25;
param.s_ee       = 20;
param.s_ie       = 8;
param.s_ei       = 20;
param.s_ii       = 20;
param.M          = 100;
param.Mr         = 66;


mp=load('membrane_potential_sample.txt');
res.V_e=mp(:,1:param.ne)';
[nr,nc]=size(res.V_e);
res.V_e=reshape(res.V_e,1,nr*nc);

res.V_i=mp(:,param.ne+1:param.ne+param.ni)';
[nr,nc]=size(res.V_i);
res.V_i=reshape(res.V_i,1,nr*nc);

rare_res.V_e=[];
rare_res.V_i=[];

P = P_generation_statistics(res, rare_res, param, bar);
%%
fop = fopen('Probability.txt','wt');
for i =1:param.ne+1
fprintf(fop, ' %s', num2str(P.P_BE_Ex(i)));
end
fprintf(fop, '\n' );
for i =1:param.ne+1
fprintf(fop, ' %s', num2str(P.P_GE_Ex(i)));
end
fprintf(fop, '\n' );
for i =1:param.ne+1
fprintf(fop, ' %s', num2str(P.P_BE_E(i)));
end
fprintf(fop, '\n' );
for i =1:param.ne+1
fprintf(fop, ' %s', num2str(P.P_GE_E(i)));
end
fprintf(fop, '\n' );
for i =1:param.ne+1
fprintf(fop,  ' %s', num2str(P.P_GE_I(i)));
end
for i =1:param.ni+1
fprintf(fop, ' %s', num2str(P.P_BI_Ex(i)));
end
fprintf(fop, '\n' );
for i =1:param.ni+1
fprintf(fop, ' %s', num2str(P.P_GI_Ex(i)));
end
fprintf(fop, '\n' );
for i =1:param.ni+1
fprintf(fop, ' %s', num2str(P.P_BI_E(i)));
end
fprintf(fop, '\n' );
for i =1:param.ni+1
fprintf(fop, ' %s', num2str(P.P_GI_E(i)));
end
fprintf(fop, '\n' );
for i =1:param.ni+1
fprintf(fop,  ' %s', num2str(P.P_GI_I(i)));
end
back = fclose( fop ) ;