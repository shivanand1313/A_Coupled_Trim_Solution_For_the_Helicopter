function fnm = Helicopter_Moment(theta_0, theta_1c, theta_1s, lambda,mu, phi_s, alpha_s)
%% Given_Data;
roh = 1.225;    %Density of air, kg/m3
Nb = 4;         %Number of blades
R = 8.18;       %Blade radius, m
C = 0.46;       %Blade chord, m
Cd0 = 0.01;     %Profile Drag coefficient
Cl_alpha = 5.73;%Lift curve slope,
v_B = 1.04 ;    %Blade flap frequency, per rev
Gamma = 8.0;    %Lock number
W = 70000;      %Weight of the aircraft, N
theta_tw = 0;   %Blade twist rate
I_beta = (roh*Cl_alpha*C*(R^4))/Gamma;
A = pi*R^2;     %DisK Area m2


Omega = 27;%Rotor angular speed, rad/sec

%% Harmonic Balance
beta_0_HB = (Gamma/(v_B^2))*(theta_0*(1+mu^2)/8 + theta_tw*(1+5*mu^2/6)/10 + mu*theta_1s/6 - lambda/6);
P = [1 (1+(mu^2)/2)*(Gamma/(8*(v_B^2 - 1))); -(1-(mu^2)/2)*(Gamma/(8*(v_B^2 - 1))) 1];
Q = (Gamma/(v_B^2 - 1))*[ (theta_1c/8)*(1+(mu^2)/2) - (mu/6)*beta_0_HB;(theta_1s/8)*(1-(mu^2)/2) ...
    + (mu/3)*theta_0 - (mu/4)*lambda + (mu^2/4)*theta_1s + mu*theta_tw/4];
Ra = P\Q;
beta_1c_HB = Ra(1);
beta_1s_HB = Ra(2);

beta_0 = beta_0_HB; 
beta_1c = beta_1c_HB;
beta_1s = beta_1s_HB;

%% Gaussian;
weight= [0.171324492 0.360761573 0.467913935 0.467913935 0.360761573 0.171324492];
node=[-0.9324695142031520278123,-0.661209386466264513661,-0.2386191860831969086305,....
    0.238619186083196908631,0.661209386466264513661,0.9324695142031520278123];
n = 10; % span devided into 10 radial segments
dy = R/n; % width of each segment

%% Main Code 
% Use Dimensional Equation for calculation of U_t, U_p, U_r
% Check Lambda error Code. 

% Recording variable
    Rec_i_dS_z = [];
    Rec_i_dn_z = [];
    Rec_F_x = [];
    Rec_F_y = [];
    Rec_F_z = [];
    Rec_M_x = [];
    Rec_M_y = [];
    Rec_M_z = [];
count = 0;
for k = 0:1:Nb

% Recording variable    
    Rec_dM_num = [];
    Rec_dM_ana = [];
    Rec_T_dS_x = [];
    Rec_T_dS_r = [];
    Rec_T_dS_z = [];
    Rec_T_dm_x = [];
    Rec_T_dm_z = [];
    Rec_T_dm_y = [];
    Rec_N_dS_z = [];
    Rec_T_F_x = [];
    Rec_T_F_y = [];
    Rec_T_F_z = [];
    Rec_T_M_x = [];
    Rec_T_M_y = [];
    Rec_T_M_z = [];
    
%% Implementing logic
DD = 60;


for psi = ( pi*k /2):(pi/DD):((2*pi)+( pi*k /2))

        count = count + 1;


    M_n = 0;    % Initialising moments for blade for each PSI location
    D_dS_x_i = 0;
    D_dS_r_i = 0;
    D_dS_z_i = 0;
    T_dS_x = 0;
    T_dS_y = 0;
    T_dS_z = 0;
    DN_dS_z_i = 0;
    T_N_dS_z  = 0;
    T_dn_f = 0;
    T_dn_l = 0;
    T_dn_t = 0;

for i = 1:1:n % span divided into n number of radial segments
      rev = 6;
      [beta,beta_str,beta_2str] = Newmarks_beta(psi,rev,mu,alpha_s,phi_s,lambda,theta_0,theta_1c,theta_1s);
      beta = beta(count,1);
      beta_str = beta_str(count,1);
      beta_2str = beta_2str(count,1);
    % blade flap angle and its derivatives w.r.t psi
%     beta = beta_0 + beta_1c*cos(psi) + beta_1s*sin(psi);
%     beta_str = -beta_1c*sin(psi) + beta_1s*cos(psi);
%     beta_2str = -beta_1c*cos(psi) - beta_1s*sin(psi);
     
    % Gaussian Quadrature limits
    a= dy*(i-1); % integration lower limit
    b= dy*i;     % integration upper limit
    
    M_i = 0;    % Initialising for each segment of gaussian Quadrature
    dS_x_i = 0;    
    dS_r_i = 0;
    dS_z_i = 0;
    N_dS_z_i = 0;
    dn_f_i = 0;
    dn_l_i = 0;
    dn_t_i = 0;

    for j = 1:1:6 % 6-Point Gaussian Quadrature
        t = ((b-a)/2)*node(j) + ((b+a)/2); % change of variable
        
        %% Gaussian integration of 6 segments

    % -----------------------------Dimensional_calculation-----------------------------------------

    %  Vilocity for Dimensional calculation
    U_t = t*Omega + (mu*sin(psi)*Omega*R); 
    U_p = (lambda*cos(beta)*Omega*R) + (t*beta_str*Omega) + (mu*cos(psi)*sin(beta)*Omega*R);
    U_r = mu*Omega*R*cos(psi);
    
    % Angles required for Calculations
    theta_psi_t = theta_0 + theta_tw*(t/R) + theta_1c*cos(psi) + theta_1s*sin(psi);
    phi = atan(U_p/U_t);
    Tau = atan(U_r/U_t);

    % Calculation of Aerodynamic Forces
    d_Lift = 0.5*roh*(U_t^2)*C*Cl_alpha*(theta_psi_t-phi);
    d_Drag = 0.5*roh*(U_t^2)*C*Cd0;

     % Aerodynamic Forces on Blade elenent
    dF_z = (d_Lift*cos(phi) - d_Drag*sin(phi))*cos(beta);
    dF_x = d_Lift*sin(phi) + d_Drag*cos(phi)*cos(Tau);
    dF_r = -d_Lift*sin(beta) + d_Drag*sin(Tau);

    % Shear loads on Blade element
    dS_z = dF_z + ((3*roh*Cl_alpha*C*(R^2)*(Omega^2)*beta_2str*t))/(Gamma);
    dS_x = dF_x;
    dS_r = -(beta*dF_r) + ((3*roh*Cl_alpha*C*(R^2)*(Omega^2)*t)/(Gamma));

    % Non-Dimensional Shear loads on Blade element
    N_dS_z = dS_z;

    % moments on Blade element
    dn_f = (v_B^2 - 1)*I_beta*(Omega^2)*beta ;
    dn_l = t*dF_x ;
    dn_t = 0; % d_M_x
    %% 

        % Shear loads on Blade element
        dS_x_i = dS_x_i +  weight(j)*0.5*(b-a)*dS_x;    % numerically integrating the Shear loads from each segment
        dS_r_i = dS_r_i +  weight(j)*0.5*(b-a)*dS_r;
        dS_z_i = dS_z_i +  weight(j)*0.5*(b-a)*dS_z;

        % Non-dimensional rotating frame vertical shear force
        N_dS_z_i = N_dS_z_i +  weight(j)*0.5*(b-a)*N_dS_z;

        % Moments on Blade element
        dn_f_i = dn_f_i +  weight(j)*0.5*(b-a)*dn_f;    % numerically integrating the moments from each segment over t (gaussian function)
        dn_l_i = dn_l_i +  weight(j)*0.5*(b-a)*dn_l;
        dn_t_i = dn_t_i +  weight(j)*0.5*(b-a)*dn_t;

    end
    % numerically integrating to whole blade
    T_dS_x = T_dS_x + dS_x_i;
    T_dS_y = T_dS_y + dS_r_i;
    T_dS_z = T_dS_z + dS_z_i;
    DN_dS_z_i = N_dS_z_i/(roh*((R*Omega)^2)*C*R); % ---------------- > DOUBT
    T_N_dS_z = T_N_dS_z + DN_dS_z_i;

    T_dn_f = T_dn_f + dn_f_i;
    T_dn_l = T_dn_l + dn_l_i;
    T_dn_t = T_dn_t + dn_t_i;
     
end 

    % Recording the Root Shear loads on Blade
    Rec_T_dS_x = [Rec_T_dS_x ; T_dS_x];
    Rec_T_dS_r = [Rec_T_dS_r ; T_dS_y];
    Rec_T_dS_z = [Rec_T_dS_z ; T_dS_z];

    % Recording the moments on Blade
    Rec_T_dm_x = [Rec_T_dm_x ; T_dn_f];     T_dm_x = T_dn_f ;
    Rec_T_dm_y = [Rec_T_dm_y ; T_dn_t];     T_dm_y = T_dn_t ;
    Rec_T_dm_z = [Rec_T_dm_z ; T_dn_l];     T_dm_z = (-T_dn_l) ;
    
    % The variation of non-dimensional rotating frame vertical shear force
    Rec_N_dS_z = [Rec_N_dS_z ; T_N_dS_z];

    %% Transforming to fixed frame for each blade
    F_x = T_dS_y*cos(psi) + T_dS_x*sin(psi);
    F_y = T_dS_y*sin(psi) - T_dS_x*cos(psi);
    F_z = T_dS_z;

    M_x = T_dm_x*sin(psi) + T_dm_y*cos(psi);
    M_y = -T_dm_x*cos(psi) + T_dm_y*sin(psi);
    M_z = T_dm_z;
    
    %%
    % Recording the Fixed Frame on Blade for each psi location
    Rec_T_F_x = [Rec_T_F_x ; F_x];
    Rec_T_F_y = [Rec_T_F_y ; F_y];
    Rec_T_F_z = [Rec_T_F_z ; F_z];
    
    % Recording the Fixed Frame moments on Blade for each psi location
    Rec_T_M_x = [Rec_T_M_x ; M_x];
    Rec_T_M_y = [Rec_T_M_y ; M_y];
    Rec_T_M_z = [Rec_T_M_z ; M_z];


end

% Recording the Non-dimensional rotating frame vertical shear force
RR_N_dS = Rec_N_dS_z';
Rec_i_dS_z = [Rec_i_dS_z; RR_N_dS];

% Recording the Forces by each Blade
RR_F_x = Rec_T_F_x';
Rec_F_x = [Rec_F_x; RR_F_x];
RR_F_y = Rec_T_F_y';
Rec_F_y = [Rec_F_y; RR_F_y];
RR_F_z = Rec_T_F_z';
Rec_F_z = [Rec_F_z; RR_F_z];

% Recording the Moments by each Blade
RR_M_x = Rec_T_M_x';
Rec_M_x = [Rec_M_x; RR_M_x];
RR_M_y = Rec_T_M_y';
Rec_M_y = [Rec_M_y; RR_M_y];
RR_M_z = Rec_T_M_z';
Rec_M_z = [Rec_M_z; RR_M_z];


end

% Adding up Forces at each psi location for all 4 Blades
Sum_F_x = Rec_F_x(1,:) + Rec_F_x(2,:) + Rec_F_x(3,:) + Rec_F_x(4,:);
Sum_F_y = Rec_F_y(1,:) + Rec_F_y(2,:) + Rec_F_y(3,:) + Rec_F_y(4,:);
Sum_F_z = Rec_F_z(1,:) + Rec_F_z(2,:) + Rec_F_z(3,:) + Rec_F_z(4,:);

% Adding up Moments at each psi location for all 4 Blades
Sum_M_x = Rec_M_x(1,:) + Rec_M_x(2,:) + Rec_M_x(3,:) + Rec_M_x(4,:);
Sum_M_y = Rec_M_y(1,:) + Rec_M_y(2,:) + Rec_M_y(3,:) + Rec_M_y(4,:);
Sum_M_z = Rec_M_z(1,:) + Rec_M_z(2,:) + Rec_M_z(3,:) + Rec_M_z(4,:);

% Total Force for all 4 Blades
a = length(Rec_T_dS_z);
One = ones(a,1);
% Rotor Thrust T
T_Thrust = (1/a).*((Sum_F_z*One)) ;
Thrust = T_Thrust.*One;
% Rotor Drag
T_Rotor_Drag = (1/a).*((Sum_F_x*One));
Rotor_Drag = T_Rotor_Drag.*One;
% Rotor Side Force
T_Rotor_SF = (1/a).*((Sum_F_y*One));
Rotor_SF = T_Rotor_SF.*One;

% Total Moments for all 4
% Rotor Torque Mz
T_Rotor_Torque = (-1/a).*((Sum_M_z*One)) ;
Rotor_Torque = T_Rotor_Torque.*One;
% Rotor Roll Moment Mx
T_Rotor_Roll = (1/a).*((Sum_M_x*One)) ;
Rotor_Roll = T_Rotor_Roll.*One;
% Rotor Pitch Moment My
T_Rotor_Pitch = (1/a).*((Sum_M_y*One)) ;
Rotor_Pitch = T_Rotor_Pitch.*One;

fnm  = [ T_Rotor_Drag, T_Rotor_SF, T_Thrust, T_Rotor_Roll, T_Rotor_Pitch, T_Rotor_Torque ] ;



end