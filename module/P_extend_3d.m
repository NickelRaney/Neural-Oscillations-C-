function [P] = P_extend_3d(P)
[n1, n2, n3] = size(P);
P_ind = zeros(n1+2, n2+2, n3+2);
P_ind_temp=(P~=0).*(1-isnan(P));
P(isnan(P)) = 0;
ke = zeros(3,3,3);
ke(2,2,1) = 1;
ke(1,2,2) = 1;
ke(2,1,2) = 1;
ke(2,3,2) = 1;
ke(3,2,2) = 1;
ke(2,2,3) = 1;
P_ind(2:end-1, 2:end-1, 2:end-1)=(1-P_ind_temp).*convn(P_ind_temp,ke,'same');
pos_temp=find(P_ind~=0);
P_ind(2:end-1, 2:end-1, 2:end-1)=P_ind_temp;
P_temp=zeros(n1+2, n2+2, n3+2);
P_temp(2:end-1, 2:end-1, 2:end-1)=P;
P=P_temp;
pos=zeros(1,n1*n2*n3+1);
pos2 = zeros(1, n1*n2*n3+1);
pos(2:length(pos_temp)+1)=pos_temp;
pos(1)=length(pos_temp);
P_ind(1,:,:)= -1;
P_ind(end,:,:)= -1;
P_ind(:,1,:)= -1;
P_ind(:,end,:)= -1;
P_ind(:,:,1)= -1;
P_ind(:,:,end)= -1;

while 1
    for i =2:pos(1)+1
        P(pos(i)) = P_expand(pos(i));
        P_ind(pos(i)) = 1;
    end
    pos=pos2;
    pos2(:)=0;
    if pos(1)==0
        break;
    end
end
P=P(2:end-1,2:end-1,2:end-1);

function [s]=P_expand(position)
    mean1 = P(position+1)*(P_ind(position+1)>0) + P(position-1)*(P_ind(position-1)>0);
    mean1 = mean1/((P_ind(position+1)>0)*(P_ind(position-1)>0) + 1);
    mean2 = P(position+n1+2)*(P_ind(position+n1+2)>0) + P(position-n1-2)*(P_ind(position-n1-2)>0);
    mean2 = mean2/((P_ind(position+n1+2)>0)*(P_ind(position-n1-2)>0) + 1);
    mean3 = P(position+(n1+2)*(n2+2))*(P_ind(position+(n1+2)*(n2+2))>0) + P(position-(n1+2)*(n2+2))*(P_ind(position-(n1+2)*(n2+2))>0);
    mean3 = mean3/((P_ind(position+(n1+2)*(n2+2))>0)*(P_ind(position-(n1+2)*(n2+2))>0) + 1);
    a = [mean1 mean2 mean3];
    weight = [0.4 0.4 0.4];
    s = (weight*a')/(weight*(a>0)');
    instack(position+1);
    instack(position-1);
    instack(position+n1+2);
    instack(position-n1-2);
    instack(position+(n1+2)*(n2+2));
    instack(position-(n1+2)*(n2+2));
end

function instack(position)
    if P_ind(position)==0
        P_ind(position)=-1;
        pos2(1)=pos2(1)+1;
        pos2(pos2(1)+1)=position;
    end
end
end

