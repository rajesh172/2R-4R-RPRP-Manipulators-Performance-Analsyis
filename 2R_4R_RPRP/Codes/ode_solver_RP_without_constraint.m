function Out= ode_solver_RP_without_constraint(t,x)


%% Input parameters
m1=0.5; m2=0.4;  l1=2; g=9.81;



%% Equation of motion 

%M matrix or Inertia Matrix
M11 = (m1+m2)*(l1^2) + (m2)*(x(3))^2 +2*(m2)*l1*x(3);
M12 = 0;
%M13 = m4*cos(x(1)-x(5))*(x(3))*(x(7)) + (m3+m4)*l2*(l1)*cos(x(1)-x(5)) + (m3+m4)*l2*x(3)*cos(x(1)-x(5)) + m4*l1*x(7)*cos(x(1)-x(5));
%M14 = -m4*sin(x(1)-x(5))*(l1+x(3));

M21 = 0;
M22 = (m2);
%M23 = m3*l2*sin(x(1)-x(5)) + m4*(l2)*sin(x(1)-x(5)) + m4*x(7)*sin(x(1)-x(5));
%M24 = m4*cos(x(1)-x(5));

% M31 = m4*cos(x(1)-x(5))*(l1+x(3))*(l2+x(7)) + m3*l2*(l1+x(3))*cos(x(1)-x(5));
% M32 = m3*l2*sin(x(1)-x(5)) + m4*(l2+x(7))*sin(x(1)-x(5));
% M33 = m3*(l2^2) + m4*(l2^2) + 2*m4*x(7)*(l2+ x(7)/2);
% M34 = 0;
% 
% M41 = -m4*sin(x(1)-x(5))*(l1+x(3));
% M42 = m4*cos(x(1)-x(5));
% M43 = 0;
% M44 = m4;


%H and G matrix or Coriolisis 

H1 =  (l1 + x(3))*(2*m2*x(4)*x(2)) + l1*g*cos(x(1))*(m1+m2)+m2*g*x(3)*cos(x(1))  ;
 
H2 =  - m2*(x(3)+l1)*x(2)^2  + g*sin(x(1))*(m2);
 
 

M = [M11 M12 ; M21 M22];

HG = [H1 ; H2];


T=[0; 0;]; %Torques and Forces

%% Equation in terms of acceleration

th_dd = (inv(M)) * (T  - HG )    ;% Joint  accelaration

OP=zeros(4,1);

%% Output
OP(1)=x(2); % It is theta1_dot
OP(2)=th_dd(1);% It is theta1_double_dot
OP(3)=x(4); % It is r1_dot
OP(4)=th_dd(2); % It is r1_double_dot
% OP(5) = x(6);% It is theta2_dot
% OP(6)=th_dd(3);% It is theta2_double_dot
% OP(7) = x(8);% It is r2_dot
% OP(8)=th_dd(4);% It is r2_double_dot

Out=OP; % Output

end