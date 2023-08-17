clc;
clear;
close all;

% Thrust Co-efficient
roh = 1.225;W = 70000;Omega = 27;R = 8.18;
C_t = W/(roh*(pi*R^2)*(Omega*R)^2);

% Rec_Variables
Rec_mu = [];
Rec_theta_0 = [];
Rec_theta_1c = [];
Rec_theta_1s = [];
Rec_lambda = [];
Rec_phi_s = [];
Rec_alpha_s = [];

%% Updating Advance Ratio
for mu = 0:0.03:0.45

% Initial guess
theta_0= deg2rad(10);
theta_1c= deg2rad(1.5);
theta_1s = deg2rad(-8);
alpha_s=deg2rad(-5);
phi_s=deg2rad(-3);

lambda=sqrt(C_t/2);
if(mu==0) % Hover Condition
    lambda=sqrt(C_t/2);
end
if(mu~=0)
    lambda=mu*tan(alpha_s)+(C_t)*1.15/(2*sqrt(mu^2+lambda^2));
end
Res=Residue(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s);
error =norm(Res);
delta_t=zeros(6,1);
i=1;
while error>1e-2
%------------------- Updating Control Inputs -------------------------
    theta_0 = theta_0 + delta_t(1);
    theta_1c= theta_1c+ delta_t(2);
    theta_1s= theta_1s+ delta_t(3);
    lambda  = lambda  + delta_t(4);   
    phi_s   = phi_s   + delta_t(5); 
    alpha_s = alpha_s + delta_t(6);
%--------------------- Control Inputs --------------------------------
    Res=Residue(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s);
    Jacb=jacobian(theta_0, theta_1c, theta_1s,lambda,mu,phi_s,alpha_s);
    delta_t=-(Jacb\Res);
    i=i+1;
    error=norm(Res);
end

% Recording Variables
Rec_theta_0 = [Rec_theta_0, theta_0];
Rec_theta_1c = [Rec_theta_1c, theta_1c];
Rec_theta_1s = [Rec_theta_1s, theta_1s];
Rec_lambda = [Rec_lambda, lambda];
Rec_phi_s = [Rec_phi_s, phi_s];
Rec_alpha_s = [Rec_alpha_s, alpha_s];
Rec_mu = [Rec_mu , mu];
end

plot_AS5_q2;
