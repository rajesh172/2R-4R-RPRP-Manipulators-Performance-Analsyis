clear all;
clc;
close all;
g = 9.81;
m = 2;
l = 1; % Length of string
b = 5; %Damping factor
lambda = 0.8; %% Learning factor
t_final = 8;
t_step = 0.1;
t_in = linspace(0,t_final,t_final/t_step +1)'; %% Time as a column vector

theta_d =  tanh(t_in); % Desired theta
Tau_in = [t_in 0*t_in];
simTau = Tau_in;
theta = zeros(size(t_in)); % Creating theta as a vector of 0
for k = 1:length(t_in)-1
    Tau_in(k+1,2) = Tau_in(k,2) + lambda*(theta_d(k,1) - theta(k,1)); % Iterative learning equation
    
    sim("Pendulum_dynamics_2021a.slx");
    simTau = Tau_in;
   
    
    theta = ans.simX.Data(:,2);
end
%theta(end) = [];
figure
hold on
plot(t_in,theta_d)
plot(t_in,theta)
legend('theta_d','theta')
hold off

