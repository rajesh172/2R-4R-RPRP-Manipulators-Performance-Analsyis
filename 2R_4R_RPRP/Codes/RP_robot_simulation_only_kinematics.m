%%
clear;close all;clc
%%
dt=1/400;
p = 1/2;
t=p+dt:-dt:0;
theta_d = p*pi - t*pi;
r_d_i = 1:1/length(theta_d):2;
r_d_f = 2:-1/length(theta_d):1;
r_d = [r_d_i r_d_f];

x = [];
y = [];
k = 1;
m1 = 1;
m2 = 1;
l1 = 1;
index = 1:1:length(t)*length(t)-2;

for j=1:1:length(theta_d)
    for i=1:1:length(theta_d)
        if (rem(j,2)~= 0)

        x(k) = r_d(i).*cos(theta_d(j));
        y(k) = r_d(i).*sin(theta_d(j));
        k = k+1;
        else
        x(k) = r_d(i+length(theta_d)).*cos(theta_d(j));
        y(k) = r_d(i+length(theta_d)).*sin(theta_d(j));
        k = k+1;
        end

    end
end

dx=(diff(x)./dt);
dy=(diff(y)./dt);
v=[dx;dy];
KX=[1 0;0 0];
%%
ath=0;
K=[cosd(ath) 0;sind(ath) 0];
kxx=K(1,1);
kxy=K(1,2);
kyx=K(2,1);
kyy=K(2,2);
q=[0;0];

for k=1:1:length(x)-1
    theta=q(1,k);r=q(2,k);
    J = jaco_2(theta,r);
    
    
    
    [xo,yo] = fwd_kin2(q(:,k));
    Err=[x(k);y(k)]-[xo;yo];
    q(:,k+1)=q(:,k) + 1*pinv(J)*(v(:,k) + 1e2*Err)*dt ;
end
animate_2r(q,K,x,y)
r = q(2,:);
theta = q(1,:);
theta_dot = diff(q(1,:));
r_dot = diff(q(2,:));
theta_ddot = diff(theta_dot);
r_ddot = diff(r_dot);
theta_dot(end) = [];
theta(end) = [];
theta_end = [];
r_dot(end) = [];
r(end) = [];
r(end) = [];
% tau1 = m1*l1*l1*theta_ddot/3 + m2*l1*l1*theta_ddot/3 + m2*r.*r.*theta_ddot + 2*m2.*r.*r_dot.*theta_dot + 2*m2*l1.*r.*theta_ddot + 2*m2*l1*r_dot.*theta_dot;
% tau2 = m2*r_ddot - m2*r.*theta_dot.*theta_dot - m2*l1*theta_dot.*theta_dot;
% Energy = tau1 + tau2;
energy = m1*l1*l1*theta_dot.*theta_dot/6 + m2*l1*l1*theta_dot.*theta_dot/2 + m2*r.*r.*theta_dot.*theta_dot/2 + m2*r.*theta_dot.*theta_dot*l1 + m2*r_dot.*r_dot/2;
plot(index,energy)


%%

function J = jaco_2(theta,r)
l1=1;
J=[-(l1+r)*sin(theta),cos(theta);(l1+r)*cos(theta),sin(theta)];
end

function [x,y] = fwd_kin2(q)
l1=1;
x=(l1+q(2,:))*cos(q(1,:));
y=(l1+q(2,:))*sin(q(1,:));
end
    
function []=animate_2r(q,K,xt,yt)
xi=q;
figure('WindowState','maximized')
%pause(7)
c=1;
for i=1:5:length(xi)
    theta=xi(1,i);
    r=xi(2,i);
    l1=1; %Input the l length
    
    % Homogeneus transformation matrix
    H01 = [-sin(theta) 0 cos(theta) 0;cos(theta) 0 sin(theta) 0;0 0 1 0;0 0 0 1]; %Frame 0 to 1 tranformation
    H12 = [1 0 0 0;0 1 0 0;0 0 1 l1+r;0 0 0 1]; %Frame 1 to 2 tranformation
    H02=H01*H12;      %Frame 0 to 2 tranformation

%     P1=[l1*cos(theta) l1*sin(theta)];
%     P2=[(l1+r)*cos(theta) (l1+r)*sin(theta)];
    P1=[H01(1,4) H01(2,4)];
    P2=[H02(1,4) H02(2,4)];
    P2x(c,:)=P2(1,1);
    P2y(c,:)=P2(1,2);
    plot(xt,yt,'--r')
    hold on
    plot(P1(1),P1(2),'ok','LineWidth',1)

    plot(P2x,P2y,'k')
    plot(P2(1),P2(2),'ok','LineWidth',1)
    % J=[-l1*sin(theta1)-l2*sin(theta2+theta1)-l3*sin(theta3+theta2+theta1)-l4*sin(theta4+theta3+theta2+theta1),-l2*sin(theta2+theta1)-l3*sin(theta3+theta2+theta1)-l4*sin(theta4+theta3+theta2+theta1),-l3*sin(theta1+theta2+theta3)-l4*sin(theta4+theta3+theta2+theta1),-l4*sin(theta4+theta3+theta2+theta1);l1*cos(theta1)+l2*cos(theta2+theta1)+l3*cos(theta3+theta2+theta1)+l4*cos(theta4+theta3+theta2+theta1),l2*cos(theta2+theta1)+l3*cos(theta3+theta2+theta1)+l4*cos(theta4+theta3+theta2+theta1),l3*cos(theta1+theta2+theta3)+l4*cos(theta4+theta3+theta2+theta1),l4*cos(theta4+theta3+theta2+theta1)];
   
    plot(0,0,'ok','LineWidth',3)
%     plot(Xof,Yof,'r')
%     plot(Xofk,Yofk,'m')
    xlim([-4 4])
    ylim([-4 4])
    axis square; 
    grid minor
    plot([0 P1(1)], [0 P1(2)],'g','LineWidth',3)
    plot([P1(1) P2(1)], [P1(2) P2(2)],'r','LineWidth',3)
    xlabel('X axis (m)','Interpreter','latex')
    ylabel('Y axis (m)','Interpreter','latex')
    set(gca,'FontSize',18)
    drawnow
    hold off
    c=c+1;
end
end