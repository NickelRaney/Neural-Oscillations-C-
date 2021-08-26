function P = P_extend_2d(P)
[r,c]=size(P);
P_ind=zeros(r+2,c+2);
P_ind_temp=(P~=0).*(1-isnan(P));

% pos=find(P_ind_temp==0);
% P(pos)=0.5;


P_ind(2:end-1,2:end-1)=(1-P_ind_temp).*conv2(P_ind_temp,[0,1,0;1,0,1;0,1,0],'same');
P_temp=zeros(r+2,c+2);
P_temp(2:end-1,2:end-1)=P;
P=P_temp;
pos=zeros(1,r*c+1);
pos2=zeros(1,r*c+1);
pos_temp=find(P_ind~=0);
pos(2:length(pos_temp)+1)=pos_temp;
pos(1)=length(pos_temp);
P_ind(2:end-1,2:end-1)=P_ind_temp;
P_ind(1,:)=16;
P_ind(end,:)=16;
P_ind(:,1)=16;
P_ind(:,end)=16;
r=r+2;

while 1
    index=zeros(1,r*c+1);
    index(2:pos(1)+1)= P_ind(pos(2:pos(1)+1)-1) + 2 * P_ind(pos(2:pos(1)+1)+1)...
        + 4*P_ind(pos(2:pos(1)+1)-r) + 8 * P_ind(pos(2:pos(1)+1)+r);
    index=mod(index,16);
    P_ind(pos(2:pos(1)+1))=1;
    for i=2:pos(1)+1
        %up=1 down=2 left=4 right=8
        switch index(i)
            case 1
                s=P(pos(i)-1);
                instack(pos(i)+1);
                instack(pos(i)-r);
                instack(pos(i)+r);
            case 2
                s=P(pos(i)+1);
                instack(pos(i)-1);
                instack(pos(i)-r);
                instack(pos(i)+r);
            case 3
                s=0.5*P(pos(i)-1)+0.5*P(pos(i)+1);
                instack(pos(i)-r);
                instack(pos(i)+r);
            case 4
                s=P(pos(i)-r);
                instack(pos(i)-1);
                instack(pos(i)+1);
                instack(pos(i)+r);
            case 5
                s=0.6*P(pos(i)-r)+0.4*P(pos(i)-1);
                instack(pos(i)+1);
                instack(pos(i)+r);
            case 6
                s=0.6*P(pos(i)-r)+0.4*P(pos(i)+1);
                instack(pos(i)-1);
                instack(pos(i)+r);
            case 7
                s=0.6*P(pos(i)-r)+0.2*P(pos(i)-1)+0.2*P(pos(i)+1);
                instack(pos(i)+r);
            case 8
                s=P(pos(i)+r);
                instack(pos(i)-1);
                instack(pos(i)-r);
                instack(pos(i)+1);
            case 9
                s=0.4*P(pos(i)-1)+0.6*P(pos(i)+r);
                instack(pos(i)-r);
                instack(pos(i)+1);
            case 10
                s=0.4*P(pos(i)+1)+0.6*P(pos(i)+r);
                instack(pos(i)-1);
                instack(pos(i)-r);
            case 11
                s=0.2*P(pos(i)-1)+0.2*P(pos(i)+1)+0.6*P(pos(i)+r);
                instack(pos(i)-r);
            case 12
                s=0.5*P(pos(i)+r)+0.5*P(pos(i)-r);
                instack(pos(i)-1);
                instack(pos(i)+1);
            case 13
                s=0.4*P(pos(i)-1)+0.3*P(pos(i)+r)+0.3*P(pos(i)-r);
                instack(pos(i)+1);
            case 14
                s=0.4*P(pos(i)+1)+0.3*P(pos(i)+r)+0.3*P(pos(i)-r);
                instack(pos(i)-1);
            case 15
                s=0.2*P(pos(i)-1)+0.2*P(pos(i)+1)+0.3*P(pos(i)+r)+0.3*P(pos(i)-r);
        end
        P(pos(i))=s;
    end
    pos=pos2;
    pos2(:)=0;
    if pos(1)==0
        break;
    end
end
P=P(2:end-1,2:end-1);


function instack(position)
    if P_ind(position)==0
        P_ind(position)=16;
        pos2(1)=pos2(1)+1;
        pos2(pos2(1)+1)=position;
    end
end

end
