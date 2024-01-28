clear all
close all
clc
%% ODE solver
 T=0:1/1000:2;
[t,x]=ode45('ode_solver_flip_velocity',T,[0,0.5,0,0.8,0,1,0,2]); 

m1=0.5; m2=0.4; m3=1; m4=0.6;  l1=1; l2=1; g=9.81;l1_max=1;l1_min=0;l2_max=1;l2_min=0;

th1=x(:,1); th1_d=x(:,2);r1=x(:,3); r1_d=x(:,4);th2=x(:,5); th2_d=x(:,6);r2=x(:,7); r2_d=x(:,8) ; %Joint position and velocities

%% Energy Calculation
%G = 200;
%G1 = 200;
for i=1:1:length(th1)
% if r1(i)>l1_max;
%     k1_max=G;
% end
% 
% 
% if r1(i)<=l1_max;
%     k1_max=0;
% end
% 
% 
% if r1(i)>=l1_min;
%     k1_min=0;
% end
% 
% 
% if r1(i)<l1_min;
%     k1_min=G;
% end
% 
% 
% if r2(i)>l2_max;
%     k2_max=G1;
% end
% if r2(i)<=l2_max;
%     k2_max=0;
% end
% if r2(i)>=l2_min;
%     k2_min=0;
% end
% if r2(i)<l2_min;
%     k2_min=G1;
% end

% Equation for kinetic energy
L1 = 1/2*(m1*l1^2)*(th1_d(i))^2;
L2 = 1/2*m2*((l1+r1(i))^2)*(th1_d(i)^2) + 1/2*m2*(r1_d(i)^2);
L3 = 1/2*m3*(r1_d(i)^2) + m3*l2*r1_d(i)*th2_d(i)*sin(th1(i)-th2(i)) + 1/2*m3*(l2)^2*(th2_d(i)^2) + 1/2*m3*((l1+r1(i))^2)*(th1_d(i)^2) + m3*l2*(l1+r1(i))*th1_d(i)*th2_d(i)*cos(th1(i)-th2(i));
L4 = 1/2*m4*(r1_d(i)^2+r2_d(i)^2) + m4*(l2+r2(i))*r1_d(i)*th2_d(i)*sin(th1(i)-th2(i)) + 1/2*m4*(l2)^2*(th2_d(i)^2) + 1/2*m4*((l1+r1(i))^2)*(th1_d(i)^2) + m4*(l2+r2(i))*(l1+r1(i))*th1_d(i)*th2_d(i)*cos(th1(i)-th2(i)) + m4*r2(i)*(th2_d(i)^2)*(l2+ r2(i)/2) + m4*r1_d(i)*r2_d(i)*cos(th1(i)-th2(i)) - m4*(l1+r1(i))*r2_d(i)*th1_d(i)*sin(th1(i)-th2(i));
% L1, L2, L3 and L4 are the individual K.E of the 4 masses w.r.t origin.
% L1 = 1/2*(m1*l1^2)*(th1_d(i))^2;
% L2 = 1/2*m2*((l1+r1(i))^2)*(th1_d(i)^2) + 1/2*m2*(r1_d(i)^2);
% L3 = m3/2*((sin(th1(i))*r1_d(i) + cos(th1(i))*th1_d(i)*(l1 + r1(i)) + l2*cos(th2(i))*th2_d(i))^2 + (l2*sin(th2(i))*th2_d(i) - cos(th1(i))*r1_d(i) + sin(th1(i))*th1_d(i)*(l1 + r1(i)))^2);
% L4 = m4/2*((sin(th1(i))*r1_d(i) + sin(th2(i))*r2_d(i) + cos(th1(i))*th1_d(i)*(l1 + r1(i)) + cos(th2(i))*th2_d(i)*(l2 + r2(i)))^2 + (cos(th1(i))*r1_d(i) + cos(th2(i))*r2_d(i) - sin(th1(i))*th1_d(i)*(l1 + r1(i)) - sin(th2(i))*th2_d(i)*(l2 + r2(i)))^2);
KE(i,1) = L1 + L2 + L3 + L4;
%KE(i,1) = 0;

% Equation for potential energy

PE_g(i,1) = g*l1*sin(th1(i))*(m1+m2+m3+m4) +g*r1(i)*sin(th1(i))*(m2+m3+m4) + g*l2*sin(th2(i))*(m3+m4) + g*r2(i)*sin(th2(i))*(m4) ; % Assuming P.E to be zero
%PE_s(i,1) = 0.5*k1_max*((l1_max-r1(i))^2) + 0.5*k1_min*((l1_min-r1(i))^2) + 0.5*k2_max*((l2_max-r2(i))^2) + 0.5*k2_min*((l2_min-r2(i))^2);
PE(i,1) = PE_g(i,1); % + PE_s(i,1);


end

%Total energy
TE=KE+PE;

figure
plot(TE)


%% Display The Results

%Ploting energies
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(211);
plot(t,KE,'r','LineWidth',1.5);
title("Kinetic Energy",'Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Energy (J) ','Interpreter','latex');
ylim([0 150]);
set(gca,'FontSize',18);
grid minor;


subplot(212);
plot(t,PE,'b','LineWidth',1.5);
title("Potential Energy",'Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Energy  (J) ','Interpreter','latex');
ylim([-150 150]);
set(gca,'FontSize',18);
grid minor;
set(gca);
% saveas(gcf,'Q1_a_KE_PE1.png')

figure('units','normalized','outerposition',[0 0 1 1]);
plot(t,TE,'LineWidth',1.5);
title("Total Energy",'Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Energy  (J)','Interpreter','latex');
ylim([0 10]);
set(gca,'FontSize',18);
grid minor;
% saveas(gcf,'Q1_a_TE1.png')