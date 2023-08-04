clear all;
close all;
clc;
%%
m1=1; m2=1;l1=1; l2=1; g=9.81;
dt=1/1000;
t=0:dt:2+dt;
r = 2;
x_ = ones(size(t));
y_ = linspace(-sqrt(3),sqrt(3),length(t));

fx_ = ones(size(t));
fy_ = zeros(size(t));

dx_=diff(x_)./dt;
dy_=diff(y_)./dt;

xd = [t' x_(1,:)' y_(1,:)'];

fd = [t' fx_(1,:)' fy_(1,:)'];

dx = [t(1,1:end-1)' dx_' dy_'];

sim("Hybrid_control_final.slx");
theta1 = ans.theta.Data(:,1);
theta2 = ans.theta.Data(:,2);
fx = ans.force.Data(:,1);
fy = ans.force.Data(:,2);
q=[theta1';theta2'];
f=[fx';fy'];

animate_2r(q,x_,y_,f)

function []=animate_2r(q,xt,yt,f)
xi=q;
figure('WindowState','maximized')
%pause(7)
c=1;
for i=1:15:length(xi)
    theta1=xi(1,i);
    theta2=xi(2,i);
    f1 = f(1,i);
    f2 = f(2,i);
    
    l1=1; %Input the l length
    l2=1; %Input the l length
   
    % Homogeneus transformation matrix
    H01 = [cos(theta1) -sin(theta1) 0 l1*cos(theta1);sin(theta1) cos(theta1) 0 l1*sin(theta1);0 0 1 0;0 0 0 1]; %Frame 0 to 1 tranformation
    H12 = [cos(theta2) -sin(theta2) 0 l2*cos(theta2);sin(theta2) cos(theta2) 0 l2*sin(theta2);0 0 1 0;0 0 0 1]; %Frame 1 to 2 tranformation
    % H23 = [cos(theta3) -sin(theta3) 0 l3*cos(theta3);sin(theta3) cos(theta3) 0 l3*sin(theta3);0 0 1 0;0 0 0 1]; %Frame 1 to 2 tranformation
    H02=H01*H12;      %Frame 0 to 2 tranformation
    % H03=H01*H12*H23;

    P1=[H01(1,4) H01(2,4)];
    P2=[H02(1,4) H02(2,4)];
    % P3=[H03(1,4) H03(2,4)];
    P2x(c,:)=P2(1,1);
    P2y(c,:)=P2(1,2);
    plot(xt,yt,'--y')
    hold on
    plot(P1(1),P1(2),'ok','LineWidth',1)

    plot(P2x,P2y,'k')
    plot(P2(1),P2(2),'ok','LineWidth',5)
    % quiver(P2(1),P2(2),f1,f2,'linewidth',1,'color','r')
    quiver(P2(1),P2(2),1,0,'linewidth',1,'color','r')
    % plot(P3(1),P3(2),'ok','LineWidth',1)

    % J=[-l1*sin(theta1)-l2*sin(theta2+theta1)-l3*sin(theta3+theta2+theta1),-l2*sin(theta2+theta1)-l3*sin(theta3+theta2+theta1),-l3*sin(theta1+theta2+theta3);l1*cos(theta1)+l2*cos(theta2+theta1)+l3*cos(theta3+theta2+theta1),l2*cos(theta2+theta1)+l3*cos(theta3+theta2+theta1),l3*cos(theta1+theta2+theta3)];
   
    plot(0,0,'ok','LineWidth',3)
%     plot(Xof,Yof,'r')
%     plot(Xofk,Yofk,'m')
    xlim([-3.5 3.5])
    ylim([-3.5 3.5])
    axis square; 
    grid minor
    plot([0 P1(1)], [0 P1(2)],'g','LineWidth',2)
    plot([P1(1) P2(1)], [P1(2) P2(2)],'b','LineWidth',2)
   % plot([P2(1) P3(1)], [P2(2) P3(2)],'g','LineWidth',2)
    xlabel('X axis (m)','Interpreter','latex')
    ylabel('Y axis (m)','Interpreter','latex')
    set(gca,'FontSize',18)
    drawnow
    hold off
    c=c+1;
end
end



