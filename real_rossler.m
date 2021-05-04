%real Rossler simulation
%%
a=[0.2,0.2,5.7];
tspan=[0:0.1:1200];
x0=[0,0,0];%假设t=0时s1=0,s2=0,如果不是自己改
[t,x]=ode45(@(t,x) rossler(a,x),tspan,x0,[]);

plot3(x(:,1),x(:,2),x(:,3))
xlabel('x');
ylabel('y');
zlabel('z');

%%
pm=[];
for i=1:12000
    if x(i,1)<0 && x(i+1,1)>0 && x(i,2)<0;
        coor=zeros(1,2);
        coor(1)=0;
        coor(2)=(-x(i,1)/(x(i+1,1)-x(i,1)))*x(i+1,2)+(x(i+1,1)/(x(i+1,1)-x(i,1)))*x(i,2);
        pm=[pm;coor];
    end
end

%%
scatter(pm(:,1),pm(:,2));
for j = 1:10:size(pm,1)-1
    PlotLineArrow(gca, [pm(j,1), pm(j+1,1)], [pm(j,2), pm(j+1,2)], 'b', 'r');
end
%%
scatter(-pm(1:203,2),-pm(2:204,2));
%%
function dx=rossler(a,x)
    c=100;
    dx=zeros(3,1);
    dx(1)=-x(2)-x(3)+c*(0.5-rand(1));
    dx(2)=x(1)+a(1)*x(2)+c*(0.5-rand(1));
    dx(3)=a(2)+x(3)*(x(1)-a(3))+c*(0.5-rand(1));
end
