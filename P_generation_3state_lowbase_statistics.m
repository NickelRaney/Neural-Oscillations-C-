%the input res.V_e is a n*ne matrix
function [P] = P_generation_3state_lowbase_statistics(res, param, bar)
bar_low_e   = bar.low_e;
bar_high_e  = bar.high_e;
bar_low_i   = bar.low_i;
bar_high_i  = bar.high_i;

V_e     = res.VE;
V_i     = res.VI;
ne      = param.ne;
ni      = param.ni;
Mr      = param.Mr;
M       = param.M;
s_ee    = ceil(param.s_ee);
s_ie    = ceil(param.s_ie);
s_ei    = param.s_ei;
s_ii    = param.s_ii;
s_ei_high    = ceil((66*s_ei+bar_high_e*166)/(166-s_ei) - bar_high_e);
s_ii_high    = ceil((66*s_ii+bar_high_i*166)/(166-s_ii) - bar_high_i);
s_ei_low     = ceil((66*s_ei+bar_low_e*166)/(166-s_ei) - bar_low_e);
s_ii_low     = ceil((66*s_ii+bar_low_i*166)/(166-s_ii) - bar_low_i);

P.P_LBE_Ex= zeros(ne+1, ne+1, ni+1);
P.P_BE_Ex = zeros(ne+1, ne+1, ni+1);
P.P_GE_Ex = zeros(ne+1, ne+1, ni+1);
P.P_LBE_E = zeros(ne+1, ne+1, ni+1);
P.P_BE_E  = zeros(ne+1, ne+1, ni+1);
P.P_GE_E  = zeros(ne+1, ne+1, ni+1);
P.P_GE_I  = zeros(ne+1, ne+1, ni+1);
P.P_BE_I  = zeros(ne+1, ne+1, ni+1);

P.P_LBI_Ex = zeros(ni+1, ni+1, ne+1);
P.P_BI_Ex = zeros(ni+1, ni+1, ne+1);
P.P_GI_Ex = zeros(ni+1, ni+1, ne+1);
P.P_LBI_E = zeros(ni+1, ni+1, ne+1);
P.P_BI_E  = zeros(ni+1, ni+1, ne+1);
P.P_GI_E  = zeros(ni+1, ni+1, ne+1);
P.P_GI_I  = zeros(ni+1, ni+1, ne+1);
P.P_BI_I  = zeros(ni+1, ni+1, ne+1);

N_LBE=sum(V_e<=bar_low_e,2);
N_LBI=sum(V_i<=bar_low_i,2);
max_N_LBE = max(N_LBE);
max_N_LBI = max(N_LBI);


PDF_e_temp=zeros(1,Mr+M+1);
PDF_i_temp=zeros(1,Mr+M+1);

for i=1:max_N_LBE+1
    V_e_temp = V_e(N_LBE==(i-1),:);
    V_i_temp = V_i(N_LBE==(i-1),:);
    N_BE = sum((V_e_temp>bar_low_e).*(V_e_temp<=bar_high_e),2);
    max_N_BE = max(N_BE);
    for j=1:max_N_BE+1
        V_e_temp2 = V_e_temp(N_BE==(j-1),:);
        V_i_temp2 = V_i_temp(N_BE==(j-1),:);
        N_GI = sum(V_i_temp2 >= bar_low_i, 2);
        max_N_GI = max(N_GI);
        for j1=1:max_N_GI+1
            V_e_temp3 = V_e_temp2(N_GI ==j1-1,:);
        for k = 1:Mr+M+1
            PDF_e_temp(k) = sum(sum(V_e_temp3==k-2-Mr));
        end
        PDF_e_temp = PDF_e_temp/sum(PDF_e_temp);
        
        P.P_LBE_Ex(i,j,j1) = PDF_e_temp(bar_low_e+Mr+2)/sum(PDF_e_temp(1: bar_low_e+Mr+2));
        P.P_BE_Ex(i,j,j1) = PDF_e_temp(bar_high_e+Mr+2)/sum(PDF_e_temp(bar_low_e+Mr+3: bar_high_e+Mr+2));
        P.P_GE_Ex(i,j,j1) = PDF_e_temp(M+Mr+1)/sum(PDF_e_temp(bar_high_e+Mr+3: M+Mr+1));
        
        P.P_LBE_E(i,j,j1) = sum(PDF_e_temp((bar_low_e+Mr+2-s_ee+1):bar_low_e+Mr+2))/sum(PDF_e_temp(1: bar_low_e+Mr+2));
        P.P_BE_E(i,j,j1) = sum(PDF_e_temp((bar_high_e+Mr+2-s_ee+1):bar_high_e+Mr+2))/sum(PDF_e_temp(bar_low_e+Mr+3: bar_high_e+Mr+2));
        P.P_GE_E(i,j,j1) = sum(PDF_e_temp((Mr+M+2-s_ee):M+Mr+1))/sum(PDF_e_temp(bar_high_e+Mr+3:M+Mr+1));

        P.P_GE_I(i,j,j1) = sum(PDF_e_temp(bar_high_e+Mr+3:bar_high_e+Mr+s_ei_high+2))/sum(PDF_e_temp(bar_high_e+Mr+3:M+Mr+1));
        P.P_BE_I(i,j,j1) = sum(PDF_e_temp(bar_low_e+Mr+3:bar_low_e+Mr+s_ei_low+2))/sum(PDF_e_temp(bar_low_e+Mr+3: bar_high_e+Mr+2));
        end
    end
end
% P.P_LBE_Ex=P_extend_2d(P.P_LBE_Ex);
% P.P_BE_Ex=P_extend_2d(P.P_BE_Ex);
% P.P_GE_Ex=P_extend_2d(P.P_GE_Ex);
% P.P_LBE_E=P_extend_2d(P.P_LBE_E);
% P.P_BE_E=P_extend_2d(P.P_BE_E);
% P.P_GE_E=P_extend_2d(P.P_GE_E);
% P.P_GE_I=P_extend_2d(P.P_GE_I);
% P.P_BE_I=P_extend_2d(P.P_BE_I);

for i=1:max_N_LBI+1
    V_i_temp = V_i(N_LBI==(i-1),:);
    V_e_temp = V_e(N_LBI==(i-1),:);
    N_BI = sum((V_i_temp>bar_low_i).*(V_i_temp<=bar_high_i),2);
    max_N_BI = max(N_BI);
    for j=1:max_N_BI+1
        V_i_temp2 = V_i_temp(N_BI==(j-1),:);
        V_e_temp2 = V_e_temp(N_BI==(j-1),:);
        N_GE = sum(V_e_temp2 >= bar_low_e, 2);
        max_N_GE = max(N_GE);
        for j1=1:max_N_GE+1
        V_i_temp3 = V_i_temp2(N_GE ==j1-1,:);

        for k = 1:Mr+M+1
            PDF_i_temp(k) = sum(sum(V_i_temp3==k-2-Mr));
        end
        PDF_i_temp = PDF_i_temp/sum(PDF_i_temp);
        
        P.P_LBI_Ex(i,j,j1) = PDF_i_temp(bar_low_i+Mr+2)/sum(PDF_i_temp(1: bar_low_i+Mr+2));
        P.P_BI_Ex(i,j,j1) = PDF_i_temp(bar_high_i+Mr+2)/sum(PDF_i_temp(bar_low_i+Mr+3: bar_high_i+Mr+2));
        P.P_GI_Ex(i,j,j1) = PDF_i_temp(M+Mr+1)/sum(PDF_i_temp(bar_high_i+Mr+3: M+Mr+1));
        
        P.P_LBI_E(i,j,j1) = sum(PDF_i_temp((bar_low_i+Mr+2-s_ie+1):bar_low_i+Mr+2))/sum(PDF_i_temp(1: bar_low_i+Mr+2));
        P.P_BI_E(i,j,j1) = sum(PDF_i_temp((bar_high_i+Mr+2-s_ie+1):bar_high_i+Mr+2))/sum(PDF_i_temp(bar_low_i+Mr+3: bar_high_i+Mr+2));
        P.P_GI_E(i,j,j1) = sum(PDF_i_temp((Mr+M+2-s_ie):M+Mr+1))/sum(PDF_i_temp(bar_high_i+Mr+3:M+Mr+1));
        
        P.P_GI_I(i,j,j1) = sum(PDF_i_temp(bar_high_i+Mr+3:bar_high_i+Mr+s_ii_high+2))/sum(PDF_i_temp(bar_high_i+Mr+3:M+Mr+1));
        P.P_BI_I(i,j,j1) = sum(PDF_i_temp(bar_low_i+Mr+3:bar_low_i+Mr+s_ii_low+2))/sum(PDF_i_temp(bar_low_i+Mr+3: bar_high_i+Mr+2));
        end
    end
end
% P.P_LBI_Ex=P_extend_2d(P.P_LBI_Ex);
% P.P_BI_Ex=P_extend_2d(P.P_BI_Ex);
% P.P_GI_Ex=P_extend_2d(P.P_GI_Ex);
% P.P_LBI_E=P_extend_2d(P.P_LBI_E);
% P.P_BI_E=P_extend_2d(P.P_BI_E);
% P.P_GI_E=P_extend_2d(P.P_GI_E);
% P.P_GI_I=P_extend_2d(P.P_GI_I);
% P.P_BI_I=P_extend_2d(P.P_BI_I);

P.P_LBI_Ex=P_extend_3d(P.P_LBI_Ex);
P.P_BI_Ex=P_extend_3d(P.P_BI_Ex);
P.P_GI_Ex=P_extend_3d(P.P_GI_Ex);
P.P_LBI_E=P_extend_3d(P.P_LBI_E);
P.P_BI_E=P_extend_3d(P.P_BI_E);
P.P_GI_E=P_extend_3d(P.P_GI_E);
P.P_GI_I=P_extend_3d(P.P_GI_I);
P.P_BI_I=P_extend_3d(P.P_BI_I);

P.P_LBE_Ex=P_extend_3d(P.P_LBE_Ex);
P.P_BE_Ex=P_extend_3d(P.P_BE_Ex);
P.P_GE_Ex=P_extend_3d(P.P_GE_Ex);
P.P_LBE_E=P_extend_3d(P.P_LBE_E);
P.P_BE_E=P_extend_3d(P.P_BE_E);
P.P_GE_E=P_extend_3d(P.P_GE_E);
P.P_GE_I=P_extend_3d(P.P_GE_I);
P.P_BE_I=P_extend_3d(P.P_BE_I);
end



