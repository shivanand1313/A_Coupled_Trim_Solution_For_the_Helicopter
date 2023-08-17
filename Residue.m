
function Res = Residue( theta_0, theta_1c, theta_1s, lamda,mu, phi_s, alpha_s)

W = 70000;
R = 8.18;
rho = 1.225;
Omega=27;
v_tip=Omega*R;
h = 1.83;
l_t= 9.75;
f = 1.85;
Xcg = -0.6;
Ycg = 0;
M_xF=0;
M_yF=0;
kh = 1.15;
theta_fp = 0;

%%
fnm = Helicopter_Moment(theta_0, theta_1c, theta_1s, lamda,mu, phi_s, alpha_s);
T=fnm(3);H=fnm(1);Y=fnm(2);
Q=fnm(6);MX=fnm(4);MY=fnm(5);

%%
D=0.5*rho*(mu*v_tip)^2*f;
Yf=0;
theta_fp=0;
Ct=T/(rho*(pi*R^2)*v_tip^2);

%% ---------------------------------------- Calculating Residuals ----------------------------------------
R1=W-T*cos(alpha_s)*cos(phi_s)+Y*sin(phi_s)-H*sin(alpha_s)*cos(phi_s)+Yf*sin(phi_s)+D*sin(theta_fp);
R2=D*cos(theta_fp)+H*cos(alpha_s)-T*sin(alpha_s);
R3=(Y+Yf)*cos(phi_s)+T*cos(alpha_s)*sin(phi_s)+H*sin(alpha_s)*sin(phi_s);
R4=MY+M_yF-W*(Xcg*cos(alpha_s)-h*sin(alpha_s))-D*(Xcg*sin(alpha_s)+h*cos(alpha_s));
R5=MX+M_xF+Yf*h+W*(h*sin(phi_s)-Ycg*cos(phi_s));

%% Inflow Ratio Calculation
if mu==0
    lambbda_0=sqrt(Ct/2);
lambda_c = lambbda_0;
else
% Fixed Point Iteration Algorithm
lambbda_0 = sqrt(Ct/2);  %initial inflow ratio value for Hover
for iteration = 1:1:3  % No. of iteration
     lam = (mu*tan(alpha_s)) + Ct/(2*sqrt(mu^2 + lambbda_0^2)) ;
     lambbda_0 = lam ;
     iteration = iteration + 1 ;
end
lambda_c = lam;
end
%%

R6=lamda-lambda_c;
% R7 = Q - Yf*l_t;
%% ----------------------------------------Calculating Residuals----------------------------------------
Res=[R1;R2;R3;R4;R5;R6];
