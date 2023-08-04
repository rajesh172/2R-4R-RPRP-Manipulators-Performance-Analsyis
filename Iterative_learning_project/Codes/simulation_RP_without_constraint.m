clear 
close all
clc
%%
%syms theta1 theta2 theta3 l1 l2 l3   
T = 0:1/10:10;
[t,x]=ode45('ode_solver_RP',T,[0,1,0,2]);

m1=0.5; m2=0.4;  l1=2; g=9.81;

th1=x(:,1); th1_d=x(:,2);r1=x(:,3); r1_d=x(:,4) ; %Joint position and velocities


%%
c=1;
for i=1:1:length(x)

theta1=x(i,1);
dtheta1=x(i,2);
r_1 = x(i,3);
dr_1 = x(i,4);
% theta2=x(i,5); 
% dtheta2=x(i,6);
% r_2 = x(i,7);
% dr_2 = x(i,8);

% Homogeneus transformation matrix
H01 = [-sin(theta1) 0 cos(theta1) 0;cos(theta1) 0 sin(theta1) 0;0 1 0 0;0 0 0 1]; %Frame 0 to 1 tranformation
H12 = [-1 0 0 0; 0 0 1 0; 0 1 0 (l1+r_1) ; 0 0 0 1]; %Frame 1 to 2 tranformation
%H23 = [-cos(theta2-theta1) 0 -sin(theta2-theta1) 0; -sin(theta2-theta1) 0 cos(theta2-theta1) 0; 0 1 0 0; 0 0 0 1]; %Frame 2 to 3 transformation
%H34 = [1 0 0 0; 0 1 0 0; 0 0 1 l2+r_2; 0 0 0 1]; % Frame 3 to 4 transformation
H02=H01*H12;      %Frame 0 to 2 tranformation     
%H03 = H02*H23;    %Frame 0 to 3 tranformation 
%H04 = H03*H34;    %Frame 0 to 4 tranformation 

O=[0,0];                   %Joint 1 position
P1=[H01(1,4) H01(2,4)];    %Joint 2 position
P2=[H02(1,4) H02(2,4)];    %Joint 3 position
%P3=[H03(1,4) H03(2,4)];    %Joint 4 position
%P4=[H04(1,4) H04(2,4)];    % End effector Position

P2x(c,:) = P2(1,1);
P2y(c,:) = P2(1,2);

% subplot(211)
plot(P1(1)+l1*cos(theta1),P1(2)+l1*sin(theta1),'oc','LineWidth',3)
hold on
plot(P2(1),P2(2),'ok','LineWidth',3)
%plot(P3(1)+l2*cos(theta2-theta1),P3(2)+l2*sin(theta2-theta1),'or','LineWidth',3)
%plot(P4(1),P4(2),'ob','LineWidth',3)
plot(P2x,P2y,'b')

plot(0,0,'ok','LineWidth',2)
xlim([-5 5])
ylim([-5 5])

grid minor
plot([0 P1(1)+l1*cos(theta1)], [0 P1(2)+l1*sin(theta1)],'b','LineWidth',2)
plot([P1(1)+l1*cos(theta1) P2(1)], [P1(2)+l1*sin(theta1) P2(2)],'m','LineWidth',2)
%plot([P2(1) P3(1)+l2*cos(theta2-theta1)], [P2(2) P3(2)+l2*sin(theta2-theta1)],'g','LineWidth',2)
%plot([P3(1)+l2*cos(theta2-theta1) P4(1)], [P3(2)+l2*sin(theta2-theta1) P4(2)],'r','LineWidth',2)
xlabel('X axis (m)','Interpreter','latex')
ylabel('Y axis (m)','Interpreter','latex')
set(gca,'FontSize',18)
drawnow
hold off
c=c+1;
end

