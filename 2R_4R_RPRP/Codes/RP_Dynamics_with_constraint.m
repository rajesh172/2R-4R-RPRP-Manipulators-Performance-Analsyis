%% Initialization
clear 
close all
clc
%%
%syms theta1 theta2 theta3 l1 l2 l3   
T=0:1/1E3:5;
x0=[0;0;1.2;0];

[t,x]=ode23s(@ode_solver_script_RP,T,x0);
clc

%%
xo=x(:,3).*cos(x(:,1));
yo=x(:,3).*sin(x(:,1));

dx=cos(x(:,1)).*x(:,4) - sin(x(:,1)).*x(:,3).*x(:,2);
dy=sin(x(:,1)).*x(:,4) + cos(x(:,1)).*x(:,3).*x(:,2);


Lmax=2;
Lmin=1;
m=1;
g=9.81;

for i=1:1:length(xo)
    G=1E10;
    ll=sqrt(xo(i)^2 + yo(i)^2);
    Kmax=G*0; Kmin=G*0;
    if(ll>Lmax)
    Kmax=G;
    end
    if(ll<Lmin)
    Kmin=G;
    end
    dL=x(i,4);
    L=x(i,3);
    Th=x(i,1);
    dTh=x(i,2);

TE(i,:)=(Kmax*Lmax^2)/2 + (Kmin*Lmin^2)/2 + (m*dL^2)/2 + (Kmax*L^2)/2 + (Kmin*L^2)/2 + (m*L^2*dTh^2)/2 - Kmax*Lmax*L - Kmin*Lmin*L + g*m*sin(Th)*L;
end

%%
figure
plot(T,TE)
ylim([TE(2)-5 TE(2)+5])
xlabel('Time','interpreter','latex','fontsize',18)
ylabel('Total Energy','interpreter','latex','fontsize',18)
grid minor
set(gca,'Fontsize',18);
%%

animate_rp(x,t)

%%

function Out= ode_solver_script_RP(t,x)

%% Input parameters
m=1; g=9.81; Lmx=2; Lmn=1; br=0; bl=0;

tau1=0;
tau2=0;

G=1E10;
Kmx=G*0; Kmn=G*0;

%% Equation of motion 
M11=m*x(3)^2;
M12=0;
M21=0;
M22=m;


if(x(3)>Lmx)
    Kmx=G;
end
if(x(3)<Lmn)
    Kmn=G;
end

H1 =2*m*x(3)*x(2)*x(4) - br*x(2);
H2 = -m*x(2)*x(2)*x(3) + Kmx*(x(3)-Lmx) + Kmn*(x(3)-Lmn)- bl*x(2);

G1=m*g*x(3)*cos(x(1));
G2=m*g*sin(x(1));

T=[tau1; tau2];

M=[M11 M12;M21 M22];

HG = [H1 + G1; H2 + G2];

%% Equation in terms of acceleration

ddth = (inv(M)) * (T  - HG )    ;

OP=zeros(4,1);

%% Output
OP(1)=x(2); %Intergretion of velocity will give the position for theta 1
OP(2)=ddth(1);%Intergretion of acceleration will give the velocity for theta 1
OP(3)=x(4); %Intergretion of velocity will give the position  for theta 2
OP(4)=ddth(2); %Intergretion of acceleration will give the velocity for theta 2


Out=OP; % Output

end


function []=animate_rp(x,t)


figure('units','normalized','outerposition',[0 0 1 1])

for i=1:floor(length(x)/1000):length(x)
    
theta1=x(i,1);
dtheta1=x(i,2); 
L=x(i,3); 
dL=x(i,4);
l1=L; %Input the l length
% Homogeneus transformation matrix
H01 = [cos(theta1) -sin(theta1) 0 l1*cos(theta1);sin(theta1) cos(theta1) 0 l1*sin(theta1);0 0 1 0;0 0 0 1]; %Frame 0 to 1 tranformation
O=[0,0];                   %Joint 1 position
P1=[H01(1,4) H01(2,4)];    %Joint 2 position
plot(P1(1),P1(2),'ok','LineWidth',5)
hold on
plot(0,0,'ok','LineWidth',10)
xlim([-2.5 2.5])
ylim([-2.5 2.5])
axis square
grid minor
plot([0 P1(1)], [0 P1(2)],'r','LineWidth',5)

% title(strcat('Time = ',num2str(t(i,1))),'Interpreter','latex')
xlabel('X axis (m)','Interpreter','latex')
ylabel('Y axis (m)','Interpreter','latex')
set(gca,'FontSize',18)
pmx=2*cos(0:1/100:2*pi);
pmy=2*sin(0:1/100:2*pi);
plot(pmx,pmy,'--k')
pmx=1*cos(0:1/100:2*pi);
pmy=1*sin(0:1/100:2*pi);
plot(pmx,pmy,'--k')
hold off

drawnow
end
end



