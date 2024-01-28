clear all;
clc;

%% Trajectory
dt=1/1000;
t=0:dt:2+dt;
theta = t*pi + pi/4;
r = sqrt(2);
x = r.*cos(theta);
y = r.*sin(theta);

% dt=1/1000;
% t=0:dt:2-dt;
% r = 2;
% x = ones(size(t));
% y = linspace(-sqrt(3)/2,sqrt(3)/2,length(t));

dx=diff(x)./dt;
dy=diff(y)./dt;
v=[dx;dy];

q=[0;pi/2];
dq = [];

for k=1:1:length(x)-1
    th1=q(1,k);th2=q(2,k);
    J = jaco_2(th1,th2);
    [xo,yo] = fwd_kin2(q(:,k));
    Err=[x(k);y(k)]-[xo;yo];
    q(:,k+1)=q(:,k) + 1*pinv(J)*(v(:,k) + 80*Err)*dt ;
    dq(:,k) = pinv(J)*v(:,k);
end


qd = [t' q'];


% dq1 = diff(q(1,:))/dt;
% dq2 = diff(q(2,:))/dt;

ddq1 = diff(dq(1,:))/dt;
ddq2 = diff(dq(2,:))/dt;

% dqd = [t(1,1:end-1)' dq1' dq2'];

dqd = [t(1,1:end-1)' dq'];

ddqd = [t(1,1:end-2)' ddq1' ddq2'];

sim("Adaptive_Control_2R_2021a_final.slx");

theta1 = ans.theta1.Data(1,:);
theta2 = ans.theta2.Data(1,:);

q1=[theta1;theta2];

animate_2r(q1,x,y)

function []=animate_2r(q1,xt,yt)
xi=q1;
figure('WindowState','maximized')
%pause(7)
c=1;
for i=1:15:length(xi)
    theta1=xi(1,i);
    theta2=xi(2,i);
    
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
    %quiver(P2(1),P2(2),1,0,'linewidth',1,'color','r')
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

function J = jaco_2(th1,th2)
l1=1;l2=1;
J=[-l1*sin(th1),-l2*sin(th2);l1*cos(th1),l2*cos(th2)];
end

function [x,y] = fwd_kin2(q)
l1=1;l2=1;
x=l1*cos(q(1,:)) + l2*cos(q(2,:));
y=l1*sin(q(1,:)) + l2*sin(q(2,:));
end
