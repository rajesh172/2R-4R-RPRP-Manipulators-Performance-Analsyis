%%
clear;close all;clc

%%
dt=1/1000;
p = 200;
t=0:dt:p+dt;
theta = p*pi - t*pi;
r = 4*theta/(p*pi);
x = r.*cos(theta);
y = r.*sin(theta);

dx=diff(x)./dt;
dy=diff(y)./dt;
v=[dx;dy];

%%
q=[0;1;0;1];
manipulability_index = [];
Global_Manipulability_index = 0;
condition_no_index = [];
Local_conditioning_index = [];
weighted_Frobenius_norm = [];
harmonic_mean = [];
stochastic_index = [];
task_space_control_index = [];
dynamic_manipulability= [];
minimum_singular_value_index = [];
global_Isotropy_Index = [];
structural_length_index = [];
velocity_index = [];
error_index = [];
Product_manipulability_condition_no = [];
Direct_selective_index = [];
kinematic_performance_index = [];

for k=1:1:length(x)-1
    
    th1 = q(1,k); r1 = q(2,k); th2 = q(3,k); r2 = q(4,k);
    
    J = jaco_RPRP(th1,r1,th2,r2);
    [xo,yo] = fwd_kin_RPRP(q(:,k));
    [xp,yp] = fwd_kin_RPRP(q(:,k) + 0.05*q(:,k));
    [xm,ym] = fwd_kin_RPRP(q(:,k) - 0.05*q(:,k));

    manipulability_index(1,k) = manipul_index(J);
    Global_Manipulability_index = Global_Manipulability_index + Glo_manipul_index(manipulability_index(1,k),r(1,k),dt,p);
    condition_no_index(1,k) = Cond_no_index(J);
    Local_conditioning_index(1,k) = 1/(condition_no_index(1,k));
    weighted_Frobenius_norm(1,k) = Frob_norm(J);
    harmonic_mean(1,k) = h_mean(J);
    stochastic_index(1,k) = stochastic(J);
    task_space_control_index(1,k) = task_control(J);
    dynamic_manipulability(1,k) = Dynamic_manipul(th1,r1,th1+th2,r2,J);
    minimum_singular_value_index(1,k) = MSV(J);
    global_Isotropy_Index(1,k) = Global_Iso(J);
    structural_length_index(1,k) = Struc_len_index(p,r1,r2);
    velocity_index(1,k) = vel_index(J);
    Product_manipulability_condition_no(1,k) = pro_manipul_cond(manipulability_index(1,k),condition_no_index(1,k));
    Direct_selective_index(1,k) = direct_selec_index(manipulability_index(1,k));
    kinematic_performance_index(1,k) = kin_perform_index(J);
    error_index(1,k) = err_index(xo,yo,xp,yp,xm,ym);

    Err=[x(k);y(k)]-[xo;yo];  
    
    if (q(2,k)<=1 && q(2,k)>=0 && q(4,k)<=1 && q(4,k)>=0)
        q(:,k+1)=q(:,k) + 1*pinv(J)*(v(:,k) + 1000*Err)*dt;

    elseif ((q(2,k)>1 || q(2,k)<0) && q(4,k)<=1 && q(4,k)>=0)
        q(:,k+1)=q(:,k) + [1;0;1;1].*(pinv(J)*(v(:,k) + 1000*Err)*dt);

    elseif ((q(4,k)>1 || q(4,k)<0) && q(2,k)<=1 && q(2,k)>=0)
        q(:,k+1)=q(:,k) + [1;1;1;0].*(pinv(J)*(v(:,k) + 1000*Err)*dt);

    else
        q(:,k+1)=q(:,k) + [1;0;1;0].*(pinv(J)*(v(:,k) + 1000*Err)*dt);
    end
end

animate_4r(q,x,y)
%%
x(end)=[];
y(end)=[];

ti = -4:0.003:4;
[x1,y1] = meshgrid(ti,ti);

%%
F = TriScatteredInterp(x',y',manipulability_index');
z1 = F(x1,y1);
figure('WindowState','maximized')
mesh(x1,y1,z1);
zlim([0 6.5]);
legend('Manipulability index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Manipulability_index_RPRP.png')

%%
F = TriScatteredInterp(x',y',condition_no_index');
z1 = F(x1,y1);
figure('WindowState','maximized')
mesh(x1,y1,z1);
legend('condition no index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Condition_Number_RPRP.png')

%%
F = TriScatteredInterp(x',y',Local_conditioning_index');
z1 = F(x1,y1);
figure('WindowState','maximized')
mesh(x1,y1,z1);
zlim([0 1]);
legend('Local conditioning index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Local_Conditioning_index_RPRP.png')

%%
F = TriScatteredInterp(x',y',weighted_Frobenius_norm');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 3]);
legend('Weighted Frobenis norm RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Weighted_Frobenius_Norm_RPRP.png')

%%
F = TriScatteredInterp(x',y',harmonic_mean');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 1.45]);
legend('Harmonic mean manipulability index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Harmonic_Mean_Manipulability_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',stochastic_index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 4.2]);
legend('Stochastic manipulability index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Stochastic_Manipulability_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',task_space_control_index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 12]);
legend('Task Space Control index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Task_Space_Control_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',dynamic_manipulability');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 25]);
legend('Dynamic manipulability index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Dynamic_Manipulability_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',minimum_singular_value_index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 1.8]);
legend('Minimum Singular Value Index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Minimum_Singular_Value_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',global_Isotropy_Index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 1]);
legend('Global Isotropy Index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Global_Isotropy_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',structural_length_index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0.1 0.2]);
legend('Structural Length_Index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Structural_Length_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',velocity_index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 11]);
legend('Velocity Index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Velocity_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',Product_manipulability_condition_no');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 31]);
legend('Product Manipulability Condition Number RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Product_Manipulability_Condition_Number_RPRP.png')

%%
F = TriScatteredInterp(x',y',Direct_selective_index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
legend('Direct Selective Index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Direct_Selective_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',kinematic_performance_index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
legend('Kinematic Performance Index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Kinematic_Performance_Index_RPRP.png')

%%
F = TriScatteredInterp(x',y',error_index');
z1 = F(x1,y1);
figure('WindowState','maximized');
mesh(x1,y1,z1);
zlim([0 12]);
legend('Error Index RPRP');
xlabel("x");ylabel("y"),zlabel('Index');
saveas(gcf,'Error_Index_RPRP.png')

%%
function J = jaco_RPRP(th1,r1,th2,r2)
l1=1;l2=1;

J = [-(l1+r1)*sin(th1) - (l2+r2)*sin(th1+th2), cos(th1), - (l2+r2)*sin(th1+th2), cos(th1+th2); (l1+r1)*cos(th1) + (l2+r2)*cos(th1+th2), sin(th1), (l2+r2)*cos(th1+th2), sin(th1+th2)];
end

%%
function [x,y] = fwd_kin_RPRP(q)
l1=1;l2=1;

x = (l1+q(2,:))*cos(q(1,:)) + (l2+q(4,:))*cos(q(1,:)+q(3,:));
y = (l1+q(2,:))*sin(q(1,:)) + (l2+q(4,:))*sin(q(1,:)+q(3,:));
end

%%
function I = manipul_index(J)
I = sqrt(abs(det(J*J')));

end

%%
function I = Glo_manipul_index(meu,r,dt,p)
A = meu*r*dt*pi; % where dt*pi is d(theta)
B = p*pi;
I = A/B;
end

%%
function I = Cond_no_index(J)
% N = [1/2,0;0,1/2];
% I = (1/2)*sqrt(trace(J*N*J')*trace(inv(J)*N*inv(J')));
I = norm(J)*norm(pinv(J));
end

%%
function I = Frob_norm(J)
n = 4;  % n denotes no. of degree of freedom of manipulator
I = sqrt(trace(J*J')/n);
end

%% 
function I = h_mean(J)
I = 1/sqrt(trace(inv(J*J')));
end

%% 
function I= stochastic(J)
n =4;  % n denotes no. of degree of freedom of manipulator
t= 2;  % t is no. of degree of freedom for task space
if det(J*J')~=0
    I = sqrt(n*t)*h_mean(J);
else
    I = 0;
end
end

%% 
function I = task_control(J)
I = sqrt(det(J*J')*Cond_no_index(J));
end

%%
function I= Dynamic_manipul(th1,r1,th2,r2,J)
m1=0.5; m2=0.5; m3=0.5; m4=0.5;  l1=1; l2=1;

%M matrix or Inertia Matrix
M11 = (m1+m2+m3+m4)*(l1^2) + (m2+m3+m4)*(r1)^2 +2*(m2+m3+m4)*l1*r1;
M12 = 0;
M13 = m4*cos(th1-th2)*r1*r2 + (m3+m4)*l2*(l1)*cos(th1-th2) + (m3+m4)*l2*r1*cos(th1-th2) + m4*l1*r2*cos(th1-th2);
M14 = -m4*sin(th1-th2)*(l1+r1);

M21 = 0;
M22 = (m2+m3+m4);
M23 = m3*l2*sin(th1-th2) + m4*(l2)*sin(th1-th2) + m4*r2*sin(th1-th2);
M24 = m4*cos(th1-th2);

M31 = m4*cos(th1-th2)*(l1+r1)*(l2+r2) + m3*l2*(l1+r1)*cos(th1-th2);
M32 = m3*l2*sin(th1-th2) + m4*(l2+r2)*sin(th1-th2);
M33 = m3*(l2^2) + m4*(l2^2) + 2*m4*r2*(l2+ r2/2);
M34 = 0;

M41 = -m4*sin(th1-th2)*(l1+r1);
M42 = m4*cos(th1-th2);
M43 = 0;
M44 = m4;
M = [M11 M12 M13 M14 ; M21 M22 M23 M24 ; M31 M32 M33 M34 ; M41 M42 M43 M44];
I = sqrt(det(J*(inv(M*M'))*J'));
end

%%
function I = MSV(J)
I = min(svd(J));
end

%%
function I = Global_Iso(J)
sigma_min = min(svd(J));
sigma_max = max(svd(J));
I = sigma_min/sigma_max;
end

%%
function I = Struc_len_index(p,r1,r2)
I = (2+r1+r2)/sqrt(p*pi);
end

%%
function I = vel_index(J)
x_dot  = J*[1;1;1;1];
I = sqrt(x_dot'*x_dot);
end

%%
function I = pro_manipul_cond(m,k)
I = m*k;
end

%%
function I = direct_selec_index(m) 
I = 1/m;
end

%%
function I = kin_perform_index(J)
I = (1/2)*sqrt(trace(J*J')*trace(inv(J*J')));

end

%%
function I = err_index(xo,yo,xp,yp,xm,ym)
err = [xo;yo;xo;yo] - [xp;yp;xm;ym];
I = sqrt(err'*err);
end

%%
function [] = animate_4r(q,xt,yt)
xi=q;
figure('WindowState','maximized')
%pause(7)
c=1;
for i=1:10:length(xi)
    theta1 = xi(1,i);
    lr1 = xi(2,i);
    theta2 = xi(3,i);
    lr2 = xi(4,i);

    %Input the l length
    l1 = 1;
    l2 = 1;

    % Homogeneus transformation matrix
    H01 = [cos(theta1) -sin(theta1) 0 l1*cos(theta1);sin(theta1) cos(theta1) 0 l1*sin(theta1);0 0 1 0;0 0 0 1]; %Frame 0 to 1 tranformation
    H12 = [1 0 0 lr1;0 1 0 0;0 0 1 0;0 0 0 1]; %Frame 1 to 2 tranformation
    H23 = [cos(theta2) -sin(theta2) 0 l2*cos(theta2);sin(theta2) cos(theta2) 0 l2*sin(theta2);0 0 1 0;0 0 0 1]; %Frame 2 to 3 tranformation
    H34 = [1 0 0 lr2;0 1 0 0;0 0 1 0;0 0 0 1]; %Frame 3 to 4 tranformation
    H02 = H01*H12;      %Frame 0 to 2 tranformation
    H03 = H02*H23;
    H04 = H03*H34;

    P1 = [H01(1,4) H01(2,4)];
    P2 = [H02(1,4) H02(2,4)];
    P3 = [H03(1,4) H03(2,4)];
    P4 = [H04(1,4) H04(2,4)];

    P4x(c,:) = P4(1,1);
    P4y(c,:) = P4(1,2);

    plot(xt,yt,'--r')
    hold on
    plot(P1(1),P1(2),'ok','LineWidth',1)

    plot(P4x,P4y,'k')
    plot(P2(1),P2(2),'ok','LineWidth',1)
    plot(P3(1),P3(2),'ok','LineWidth',1)
    plot(P4(1),P4(2),'ok','LineWidth',1)
    
    plot(0,0,'ok','LineWidth',3)
    xlim([-5 5])
    ylim([-5 5])
    axis square; 
    grid minor
    plot([0 P1(1)], [0 P1(2)],'g','LineWidth',2)
    plot([P1(1) P2(1)], [P1(2) P2(2)],'b','LineWidth',2)
    plot([P2(1) P3(1)], [P2(2) P3(2)],'y','LineWidth',2)
    plot([P3(1) P4(1)], [P3(2) P4(2)],'c','LineWidth',2)
    xlabel('X axis (m)','Interpreter','latex')
    ylabel('Y axis (m)','Interpreter','latex')
    set(gca,'FontSize',18)
    drawnow
    hold off
    c=c+1;
end
end
