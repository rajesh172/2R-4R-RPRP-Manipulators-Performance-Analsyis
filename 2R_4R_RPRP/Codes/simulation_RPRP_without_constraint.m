clear 
close all
clc
%%
%syms theta1 theta2 theta3 l1 l2 l3   
T = 0:1/10:100;
[t,x]=ode45('ode_solver_without_constraint',T,[0,1,0,0,0,0,0,0]);

m1=0.5; m2=0.4; m3=1; m4=0.6;  l1=1; l2=1; g=-9.81;

th1=x(:,1); th1_d=x(:,2);r1=x(:,3); r1_d=x(:,4);th2=x(:,5); th2_d=x(:,6);r2=x(:,7); r2_d=x(:,8) ; %Joint position and velocities


%%
for i=1:1:length(th1)

% Equation for kinetic energy
L1 = 1/2*(m1*l1^2)*(th1_d(i))^2;
L2 = 1/2*m2*((l1+r1(i))^2)*(th1_d(i)^2) + 1/2*m2*(r1_d(i)^2);
L3 = 1/2*m3*(r1_d(i)^2) + m3*l2*r1_d(i)*th2_d(i)*sin(th1(i)-th2(i)) + 1/2*m3*(l2)^2*(th2_d(i)^2) + 1/2*m3*((l1+r1(i))^2)*(th1_d(i)^2) + m3*l2*(l1+r1(i))*th1_d(i)*th2_d(i)*cos(th1(i)-th2(i));
L4 = 1/2*m4*(r1_d(i)^2+r2_d(i)^2) + m4*(l2+r2(i))*r1_d(i)*th2_d(i)*sin(th1(i)-th2(i)) + 1/2*m4*(l2)^2*(th2_d(i)^2) + 1/2*m4*((l1+r1(i))^2)*(th1_d(i)^2) + m4*(l2+r2(i))*(l1+r1(i))*th1_d(i)*th2_d(i)*cos(th1(i)-th2(i)) + m4*r2(i)*(th2_d(i)^2)*(l2+ r2(i)/2) + m4*r1_d(i)*r2_d(i)*cos(th1(i)-th2(i)) - m4*(l1+r1(i))*r2_d(i)*th1_d(i)*sin(th1(i)-th2(i));

KE(i,1) = L1 + L2 + L3 + L4;

% Equation for potential energy
PE(i,1) = 0; % Assuming P.E to be zero
end

%Total energy
TE=KE+PE;

figure('units','normalized','outerposition',[0 0 1 1])
c=1;
for i=1:1:length(x)
    
theta1=x(i,1);
dtheta1=x(i,2);
r_1 = x(i,3);
dr_1 = x(i,4);
theta2=x(i,5); 
dtheta2=x(i,6);
r_2 = x(i,7);
dr_2 = x(i,8);

% Homogeneus transformation matrix
H01 = [-sin(theta1) 0 cos(theta1) 0;cos(theta1) 0 sin(theta1) 0;0 1 0 0;0 0 0 1]; %Frame 0 to 1 tranformation
H12 = [-1 0 0 0; 0 0 1 0; 0 1 0 (l1+r_1) ; 0 0 0 1]; %Frame 1 to 2 tranformation
H23 = [-cos(theta2-theta1) 0 -sin(theta2-theta1) 0; -sin(theta2-theta1) 0 cos(theta2-theta1) 0; 0 1 0 0; 0 0 0 1]; %Frame 2 to 3 transformation
H34 = [1 0 0 0; 0 1 0 0; 0 0 1 l2+r_2; 0 0 0 1]; % Frame 3 to 4 transformation
H02=H01*H12;      %Frame 0 to 2 tranformation     
H03 = H02*H23;    %Frame 0 to 3 tranformation 
H04 = H03*H34;    %Frame 0 to 4 tranformation 

O=[0,0];                   %Joint 1 position
P1=[H01(1,4) H01(2,4)];    %Joint 2 position
P2=[H02(1,4) H02(2,4)];    %Joint 3 position
P3=[H03(1,4) H03(2,4)];    %Joint 4 position
P4=[H04(1,4) H04(2,4)];    % End effector Position

P4x(c,:) = P4(1,1);
P4y(c,:) = P4(1,2);

% subplot(211)
plot(P1(1)+l1*cos(theta1),P1(2)+l1*sin(theta1),'oc','LineWidth',3)
hold on
plot(P2(1),P2(2),'ok','LineWidth',3)
plot(P3(1)+l2*cos(theta2-theta1),P3(2)+l2*sin(theta2-theta1),'or','LineWidth',3)
plot(P4x,P4y,'b')
plot(P4(1),P4(2),'ob','LineWidth',3)
plot(0,0,'ok','LineWidth',2)
xlim([-25 25])
ylim([0 50])

grid minor
plot([0 P1(1)+l1*cos(theta1)], [0 P1(2)+l1*sin(theta1)],'b','LineWidth',2)
plot([P1(1)+l1*cos(theta1) P2(1)], [P1(2)+l1*sin(theta1) P2(2)],'m','LineWidth',2)
plot([P2(1) P3(1)+l2*cos(theta2-theta1)], [P2(2) P3(2)+l2*sin(theta2-theta1)],'g','LineWidth',2)
plot([P3(1)+l2*cos(theta2-theta1) P4(1)], [P3(2)+l2*sin(theta2-theta1) P4(2)],'r','LineWidth',2)
xlabel('X axis (m)','Interpreter','latex')
ylabel('Y axis (m)','Interpreter','latex')
set(gca,'FontSize',18)
drawnow
hold off
% title(strcat('Time = ',num2str(t(i,1))),'Interpreter','latex')

% pause(0.000000000000000000000000000000000000000001);

% subplot(212)
% plot(t(i,1),TE(i,1),'ok','LineWidth',1)
% hold on
% title("Total Energy",'Interpreter','latex')
% xlabel('Time (s)','Interpreter','latex')
% ylabel('Energy  (J)','Interpreter','latex')
% ylim([-40 40])
% set(gca,'FontSize',18)
% grid minor


% subplot(222)
% plot(t(i,1),KE(i,1),'or','LineWidth',1)
% hold on
% title("Kinetic Energy",'Interpreter','latex')
% xlabel('Time (s)','Interpreter','latex')
% ylabel('Energy (J) ','Interpreter','latex')
% ylim([-40 40])
% set(gca,'FontSize',18)
% grid minor
% 
% subplot(224)
% plot(t(i,1),PE(i,1),'ob','LineWidth',1)
% hold on
% title("Potential Energy",'Interpreter','latex')
% xlabel('Time (s)','Interpreter','latex')
% ylabel('Energy  (J) ','Interpreter','latex')
% ylim([-40 40])
% set(gca,'FontSize',18)
% grid minor
% set(gca)
% F(i) = getframe(gcf) ;

c=c+1;
end

% % create the video writer with 30 fps
%   writerObj = VideoWriter('animation_Q1_a.avi');
%   writerObj.FrameRate = 30;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);