clear all;
clc;
close all;
g = 9.81;
m = 2;
l = 1; % Length of string
b = 5; %Damping factor
lambda = 0.9; %% Learning factor
t_final = 8;
t_step = 0.1;
t_in = linspace(0,t_final,t_final/t_step +1)'; %% Time as a column vector
C = [1,0];
B = [0;t_step/(m*l*l)];

    


theta_d =  tanh(t_in); % Desired theta
Tau_in = [t_in 0*t_in];
Tau_ff = [t_in 0*t_in];
%simTau = Tau_in;
theta = zeros(size(t_in)); % Creating theta as a vector of 0
theta_dot = zeros(size(t_in)); % Creating theta_dot as a vector of 0

for k = 1:length(t_in)-1
    Tau_in(k+1,2) = Tau_in(k,2) + lambda*tanh(theta_d(k,1) - theta(k,1)); % Iterative learning equation
    A = [1,t_step;-(m*g*l*cos(theta_d(k,1))*t_step/(m*l*l)),((-b/(m*l*l))*t_step+1)];
    Tau_ff(k+1,2) =   inv(C*A*B)*(theta_d(k+1)-C*A*A*[theta(k,1);theta_dot(k,1)]);
   
    simTau = Tau_in;
    simTau_ff = Tau_ff;
    sim("pendulum_ILC_ff.slx");
    
   
    theta_dot = ans.simX.Data(:,1);
    theta = ans.simX.Data(:,2);
end


figure
hold on
plot(t_in,theta_d)
plot(t_in,theta)
legend('theta_d','theta')
hold off



