%the input res.V_e is a n*ne matrix
function [P] = P_generation_3state_lowbase_statistics(res, param, bar)
bar_low_e   = bar.low_e;
bar_high_e  = bar.high_e;

bar_low_i   = bar.low_i;
bar_high_i  = bar.high_i;
V_e     = res.V_e;
V_i     = res.V_i;
ne      = param.ne;
ni      = param.ni;
Mr      = param.Mr;
M       = param.M;
s_ee    = param.s_ee;
s_ie    = param.s_ie;
s_ei    = param.s_ei;
s_ii    = param.s_ii;
s_ei    = ceil((66*s_ei+bar_e*166)/(166-s_ei) - bar_e);
s_ii    = ceil((66*s_ii+bar_i*166)/(166-s_ii) - bar_i);

P.P_LBE_Ex = zeros(ne+1, ne+1);
P.P_BE_Ex = zeros(ne+1, ne+1);
P.P_GE_Ex = zeros(ne+1, ne+1);
P.P_LBE_E  = zeros(ne+1, ne+1);
P.P_BE_E  = zeros(ne+1, ne+1);
P.P_GE_E  = zeros(ne+1, ne+1);
P.P_GE_I  = zeros(ne+1, ne+1);
P.P_BE_I  = zeros(ne+1, ne+1);

P.P_LBI_Ex = zeros(ni+1, ni+1);
P.P_BI_Ex = zeros(ni+1, ni+1);
P.P_GI_Ex = zeros(ni+1, ni+1);
P.P_LBI_E = zeros(ni+1, ni+1);
P.P_BI_E  = zeros(ni+1, ni+1);
P.P_GI_E  = zeros(ni+1, ni+1);
P.P_GI_I  = zeros(ni+1, ni+1);
P.P_BI_I  = zeros(ni+1, ni+1);

N_LBE=sum(ve<=bar_low_e,2);
N_LBI=sum(vi<=bar_low_i,2);
max_N_LBE = max(N_LBE);
max_N_LBI = max(N_LBI);


PDF_e_temp=zeros(1,Mr+M+1);
PDF_i_temp=zeros(1,Mr+M+1);

for i=1:max_N_LBE+1
    V_e_temp = V_e(N_LBE==(i-1),:);
    N_BE = sum((V_e_temp>bar_low_e&V_e_temp<=bar_high_e),2);
    max_N_BE = max(N_BE);
    for j=1:max_N_BE+1
        V_e_temp2 = V_e_temp(N_BE==(j-1),:);
        for k = 1:Mr+M+1
            PDF_e_temp(k) = sum(sum(V_e_temp2==k-2-Mr));
        end
        PDF_e_temp = PDF_e_temp/sum(PDF_e_temp);
        
        P.P_LBE_Ex(i,j) = PDF_e_temp(bar_low_e+Mr+2)/sum(PDF_e_temp(1: bar_low_e+Mr+2));
        P.P_BE_Ex(i,j) = PDF_e_temp(bar_high_e+Mr+2)/sum(PDF_e_temp(bar_low_e+Mr+3: bar_high_e+Mr+2));
        P.P_GE_Ex(i,j) = PDF_e_temp(M+Mr+1)/sum(PDF_e_temp(bar_high_e+Mr+3: M+Mr+1));
        
        P.P_LBE_E(i,j) = sum(PDF_e_temp((bar_low_e+Mr+2-s_ee+1):bar_low_e+Mr+2))/sum(PDF_e_temp(1: bar_low_e+Mr+2));
        P.P_BE_E(i.j) = sum(PDF_e_temp((bar_high_e+Mr+2-s_ee+1):bar_high_e+Mr+2))/sum(PDF_e_temp(bar_low_e+Mr+3: bar_high_e+Mr+2));
        P.P_GE_E(i,j) = sum(PDF_e_temp((Mr+M+2-s_ee):M+Mr+1))/sum(PDF_e_temp(bar_high_e+Mr+3:M+Mr+1));
        
        P.P_GE_I(i,j) = sum(PDF_e_temp(bar_high_e+Mr+3:bar_high_e+Mr+s_ei+2))/sum(PDF_e_temp(bar_high_e+Mr+3:M+Mr+1));
        P.P_BE_I(i,j) = sum(PDF_e_temp(bar_low_e+Mr+3:bar_low_e+Mr+s_ei+2))/sum(PDF_e_temp(bar_low_e+Mr+3: bar_high_e+Mr+2));
        
    end
    if isempty(find(P.P_BE_Ex(i,:)~=0&isnan(P.P_BE_Ex(i,:))==0))||isempty(find(P.P_GE_Ex(i,:)~=0&isnan(P.P_GE_Ex(i,:))==0))...
            ||isempty(find(P.P_BE_E(i,:)~=0&isnan(P.P_BE_E(i,:))==0))||isempty(find(P.P_GE_E(i,:)~=0&isnan(P.P_GE_E(i,:))==0))||...
            isempty(find(P.P_GE_I(i,:)~=0&isnan(P.P_GE_I(i,:))==0))||isempty(find(P.P_LBE_Ex(i,:)~=0&isnan(P.P_LBE_Ex(i,:))==0))...
            ||isempty(find(P.P_LBI_Ex(i,:)~=0&isnan(P.P_LBI_Ex(i,:))==0))||isempty(find(P.P_BE_I(i,:)~=0&isnan(P.P_BE_I(i,:))==0))
        P.P_LBE_Ex(i,:)=0;
        P.P_BE_Ex(i,:)=0;
        P.P_GE_Ex(i,:)=0;
        P.P_LBE_E(i,:)=0;
        P.P_BE_E(i,:)=0;
        P.P_GE_E(i,:)=0;
        P.P_GE_I(i,:)=0;
        P.P_BE_I(i,:)=0;
        break;
    end
    
    position=max(find(P.P_LBE_Ex(i,:)~=0));
    P.P_LBE_Ex(i,position+1:ne+1)=P.P_LBE_Ex(i,position);
    position=min(find(P.P_LBE_Ex(i,:)~=0&isnan(P.P_LBE_Ex(i,:))==0));
    P.P_LBE_Ex(i,1:position-1)=P.P_LBE_Ex(i,position);
    position=find(isnan(P.P_LBE_Ex(i,:))==1|P.P_LBE_Ex(i,:)==0);
    for j=1:length(position)
        P.P_LBE_Ex(i,position(j))=P.P_LBE_Ex(i,position(j)-1);
    end
    
    position=max(find(P.P_BE_Ex(i,:)~=0));
    P.P_BE_Ex(i,position+1:ne+1)=P.P_BE_Ex(i,position);
    position=min(find(P.P_BE_Ex(i,:)~=0&isnan(P.P_BE_Ex(i,:))==0));
    P.P_BE_Ex(i,1:position-1)=P.P_BE_Ex(i,position);
    position=find(isnan(P.P_BE_Ex(i,:))==1|P.P_BE_Ex(i,:)==0);
    for j=1:length(position)
        P.P_BE_Ex(i,position(j))=P.P_BE_Ex(i,position(j)-1);
    end
    
    
    position=max(find(P.P_GE_Ex(i,:)~=0));
    P.P_GE_Ex(i,position+1:ne)=P.P_GE_Ex(i,position);
    position=min(find(P.P_GE_Ex(i,:)~=0&isnan(P.P_GE_Ex(i,:))==0));
    P.P_GE_Ex(i,1:position-1)=P.P_GE_Ex(i,position);
    position=find(isnan(P.P_GE_Ex(i,:))==1|P.P_GE_Ex(i,:)==0);
    for j=1:length(position)
        P.P_GE_Ex(i,position(j))=P.P_GE_Ex(i,position(j)-1);
    end
    
    position=max(find(P.P_LBE_E(i,:)~=0));
    P.P_LBE_E(i,position+1:ne+1)=P.P_LBE_E(i,position);
    position=min(find(P.P_LBE_E(i,:)~=0&isnan(P.P_LBE_E(i,:))==0));
    P.P_LBE_E(i,1:position-1)=P.P_LBE_E(i,position);
    position=find(isnan(P.P_LBE_E(i,:))==1|P.P_LBE_E(i,:)==0);
    for j=1:length(position)
        P.P_LBE_E(i,position(j))=P.P_LBE_E(i,position(j)-1);
    end
    
    position=max(find(P.P_BE_E(i,:)~=0));
    P.P_BE_E(i,position+1:ne)=P.P_BE_E(i,position);
    position=min(find(P.P_BE_E(i,:)~=0&isnan(P.P_BE_E(i,:))==0));
    P.P_BE_E(i,1:position-1)=P.P_BE_E(i,position);
    position=find(isnan(P.P_BE_E(i,:))==1|P.P_BE_E(i,:)==0);
    for j=1:length(position)
        P.P_BE_E(i,position(j))=P.P_BE_E(i,position(j)-1);
    end
    
    position=max(find(P.P_GE_E(i,:)~=0));
    P.P_GE_E(i,position+1:ne)=P.P_GE_E(i,position);
    position=min(find(P.P_GE_E(i,:)~=0&isnan(P.P_GE_E(i,:))==0));
    P.P_GE_E(i,1:position-1)=P.P_GE_E(i,position);
    position=find(isnan(P.P_GE_E(i,:))==1|P.P_GE_E(i,:)==0);
    for j=1:length(position)
        P.P_GE_E(i,position(j))=P.P_GE_E(i,position(j)-1);
    end
    
    position=max(find(P.P_GE_I(i,:)~=0));
    P.P_GE_I(i,position+1:ne)=P.P_GE_I(i,position);
    position=min(find(P.P_GE_I(i,:)~=0&isnan(P.P_GE_I(i,:))==0));
    P.P_GE_I(i,1:position-1)=P.P_GE_I(i,position);
    position=find(isnan(P.P_GE_I(i,:))==1|P.P_GE_I(i,:)==0);
    for j=1:length(position)
        P.P_GE_I(i,position(j))=P.P_GE_I(i,position(j)-1);
    end
    
    position=max(find(P.P_BE_I(i,:)~=0));
    P.P_BE_I(i,position+1:ne+1)=P.P_BE_I(i,position);
    position=min(find(P.P_BE_I(i,:)~=0&isnan(P.P_BE_I(i,:))==0));
    P.P_BE_I(i,1:position-1)=P.P_BE_I(i,position);
    position=find(isnan(P.P_BE_I(i,:))==1|P.P_BE_I(i,:)==0);
    for j=1:length(position)
        P.P_BE_I(i,position(j))=P.P_BE_I(i,position(j)-1);
    end
end

for i=1:max_N_LBI+1
    V_i_temp = V_i(N_LBI==(i-1),:);
    N_BI = sum((V_i_temp>bar_low_i&V_i_temp<=bar_high_i),2);
    max_N_BI = max(N_BI);
    for j=1:max_N_BI+1
        V_i_temp2 = V_i_temp(N_BI==(j-1),:);
        for k = 1:Mr+M+1
            PDF_i_temp(k) = sum(sum(V_i_temp2==k-2-Mr));
        end
        PDF_i_temp = PDF_i_temp/sum(PDF_i_temp);
        
        P.P_LBI_Ex(i,j) = PDF_i_temp(bar_low_i+Mr+2)/sum(PDF_i_temp(1: bar_low_i+Mr+2));
        P.P_BI_Ex(i,j) = PDF_i_temp(bar_high_i+Mr+2)/sum(PDF_i_temp(bar_low_i+Mr+3: bar_high_i+Mr+2));
        P.P_GI_Ex(i,j) = PDF_i_temp(M+Mr+1)/sum(PDF_i_temp(bar_high_i+Mr+3: M+Mr+1));
        
        P.P_LBI_E(i,j) = sum(PDF_i_temp((bar_low_i+Mr+2-s_ie+1):bar_low_i+Mr+2))/sum(PDF_i_temp(1: bar_low_i+Mr+2));
        P.P_BI_E(i.j) = sum(PDF_i_temp((bar_high_i+Mr+2-s_ie+1):bar_high_i+Mr+2))/sum(PDF_i_temp(bar_low_i+Mr+3: bar_high_i+Mr+2));
        P.P_GI_E(i,j) = sum(PDF_i_temp((Mr+M+2-s_ie):M+Mr+1))/sum(PDF_i_temp(bar_high_i+Mr+3:M+Mr+1));
        
        P.P_GI_I(i,j) = sum(PDF_i_temp(bar_high_i+Mr+3:bar_high_i+Mr+s_ii+2))/sum(PDF_i_temp(bar_high_i+Mr+3:M+Mr+1));
        P.P_BI_I(i,j) = sum(PDF_i_temp(bar_low_i+Mr+3:bar_low_i+Mr+s_ii+2))/sum(PDF_i_temp(bar_low_i+Mr+3: bar_high_i+Mr+2));
        
    end
    if isempty(find(P.P_BI_Ex(i,:)~=0&isnan(P.P_BI_Ex(i,:))==0))||isempty(find(P.P_GI_Ex(i,:)~=0&isnan(P.P_GI_Ex(i,:))==0))...
            ||isempty(find(P.P_BI_E(i,:)~=0&isnan(P.P_BI_E(i,:))==0))||isempty(find(P.P_GI_E(i,:)~=0&isnan(P.P_GI_E(i,:))==0))||...
            isempty(find(P.P_GI_I(i,:)~=0&isnan(P.P_GI_I(i,:))==0))||isempty(find(P.P_LBI_Ex(i,:)~=0&isnan(P.P_LBI_Ex(i,:))==0))...
            ||isempty(find(P.P_LBI_Ex(i,:)~=0&isnan(P.P_LBI_Ex(i,:))==0))||isempty(find(P.P_BI_I(i,:)~=0&isnan(P.P_BI_I(i,:))==0))
        P.P_LBI_Ex(i,:)=0;
        P.P_BI_Ex(i,:)=0;
        P.P_GI_Ex(i,:)=0;
        P.P_LBI_E(i,:)=0;
        P.P_BI_E(i,:)=0;
        P.P_GI_E(i,:)=0;
        P.P_GI_I(i,:)=0;
        P.P_BI_I(i,:)=0;
        break;
    end
    
    position=max(find(P.P_LBI_Ex(i,:)~=0));
    P.P_LBI_Ex(i,position+1:ne+1)=P.P_LBI_Ex(i,position);
    position=min(find(P.P_LBI_Ex(i,:)~=0&isnan(P.P_LBI_Ex(i,:))==0));
    P.P_LBI_Ex(i,1:position-1)=P.P_LBI_Ex(i,position);
    position=find(isnan(P.P_LBI_Ex(i,:))==1|P.P_LBI_Ex(i,:)==0);
    for j=1:length(position)
        P.P_LBI_Ex(i,position(j))=P.P_LBI_Ex(i,position(j)-1);
    end
    
    position=max(find(P.P_BI_Ex(i,:)~=0));
    P.P_BI_Ex(i,position+1:ne+1)=P.P_BI_Ex(i,position);
    position=min(find(P.P_BI_Ex(i,:)~=0&isnan(P.P_BI_Ex(i,:))==0));
    P.P_BI_Ex(i,1:position-1)=P.P_BI_Ex(i,position);
    position=find(isnan(P.P_BI_Ex(i,:))==1|P.P_BI_Ex(i,:)==0);
    for j=1:length(position)
        P.P_BI_Ex(i,position(j))=P.P_BI_Ex(i,position(j)-1);
    end
    
    
    position=max(find(P.P_GI_Ex(i,:)~=0));
    P.P_GI_Ex(i,position+1:ne)=P.P_GI_Ex(i,position);
    position=min(find(P.P_GI_Ex(i,:)~=0&isnan(P.P_GI_Ex(i,:))==0));
    P.P_GI_Ex(i,1:position-1)=P.P_GI_Ex(i,position);
    position=find(isnan(P.P_GI_Ex(i,:))==1|P.P_GI_Ex(i,:)==0);
    for j=1:length(position)
        P.P_GI_Ex(i,position(j))=P.P_GI_Ex(i,position(j)-1);
    end
    
    position=max(find(P.P_LBI_E(i,:)~=0));
    P.P_LBI_E(i,position+1:ne+1)=P.P_LBI_E(i,position);
    position=min(find(P.P_LBI_E(i,:)~=0&isnan(P.P_LBI_E(i,:))==0));
    P.P_LBI_E(i,1:position-1)=P.P_LBI_E(i,position);
    position=find(isnan(P.P_LBI_E(i,:))==1|P.P_LBI_E(i,:)==0);
    for j=1:length(position)
        P.P_LBI_E(i,position(j))=P.P_LBI_E(i,position(j)-1);
    end
    
    position=max(find(P.P_BI_E(i,:)~=0));
    P.P_BI_E(i,position+1:ne)=P.P_BI_E(i,position);
    position=min(find(P.P_BI_E(i,:)~=0&isnan(P.P_BI_E(i,:))==0));
    P.P_BI_E(i,1:position-1)=P.P_BI_E(i,position);
    position=find(isnan(P.P_BI_E(i,:))==1|P.P_BI_E(i,:)==0);
    for j=1:length(position)
        P.P_BI_E(i,position(j))=P.P_BI_E(i,position(j)-1);
    end
    
    position=max(find(P.P_GI_E(i,:)~=0));
    P.P_GI_E(i,position+1:ne)=P.P_GI_E(i,position);
    position=min(find(P.P_GI_E(i,:)~=0&isnan(P.P_GI_E(i,:))==0));
    P.P_GI_E(i,1:position-1)=P.P_GI_E(i,position);
    position=find(isnan(P.P_GI_E(i,:))==1|P.P_GI_E(i,:)==0);
    for j=1:length(position)
        P.P_GI_E(i,position(j))=P.P_GI_E(i,position(j)-1);
    end
    
    position=max(find(P.P_GI_I(i,:)~=0));
    P.P_GI_I(i,position+1:ne)=P.P_GI_I(i,position);
    position=min(find(P.P_GI_I(i,:)~=0&isnan(P.P_GI_I(i,:))==0));
    P.P_GI_I(i,1:position-1)=P.P_GI_I(i,position);
    position=find(isnan(P.P_GI_I(i,:))==1|P.P_GI_I(i,:)==0);
    for j=1:length(position)
        P.P_GI_I(i,position(j))=P.P_GI_I(i,position(j)-1);
    end
    
    position=max(find(P.P_BI_I(i,:)~=0));
    P.P_BI_I(i,position+1:ne+1)=P.P_BI_I(i,position);
    position=min(find(P.P_BI_I(i,:)~=0&isnan(P.P_BI_I(i,:))==0));
    P.P_BI_I(i,1:position-1)=P.P_BI_I(i,position);
    position=find(isnan(P.P_BI_I(i,:))==1|P.P_BI_I(i,:)==0);
    for j=1:length(position)
        P.P_BI_I(i,position(j))=P.P_BI_I(i,position(j)-1);
    end
end



