clear all;
close all;
clc;
%%
%syms theta1 theta2 theta3 l1 l2 l3   
% T = 0:1/1000:5;
m1=1; m2=1;l1=1; l2=1; g=9.81;
dt=1/1000;
t=0:dt:2+dt;
theta = 2*pi - t*pi;
r = 2;
x_ = r.*cos(theta);
y_ = r.*sin(theta);

dx=diff(x_)./dt;
dy=diff(y_)./dt;
v=[dx;dy];

%%
q=[0.01;0];
for k=1:1:length(x_)-1
    th1=q(1,k);th2=q(2,k);
    J = jaco_2(th1,th2);
    [xo,yo] = fwd_kin2(q(1,k),q(2,k));
    Err=[x_(k);y_(k)]-[xo;yo];
    q(:,k+1)=q(:,k) + 1*pinv(J)*(v(:,k) + 1e3*Err)*dt ;    
end
q1d = [t' q(1,:)'];
q2d = [t' q(2,:)'];
dq1d = diff(q(1,:))./dt;
dq2d = diff(q(2,:))./dt;
q1d_d = [t(1,1:end-1)' dq1d'];
q2d_d = [t(1,1:end-1)' dq2d'];

function J = jaco_2(th1,th2)
l1=1;l2=1;
J=[-l1*sin(th1) - l2*sin(th1+th2),-l2*sin(th1+th2);l1*cos(th1)+l2*cos(th1+th2),l2*cos(th1+th2)];
end

function [x,y] = fwd_kin2(q1,q2)
l1=1;l2=1;
x=l1*cos(q1(1,:)) + l2*cos(q1(1,:)+ q2(1,:));
y=l1*sin(q1(1,:)) + l2*sin(q1(1,:)+ q2(1,:));
end