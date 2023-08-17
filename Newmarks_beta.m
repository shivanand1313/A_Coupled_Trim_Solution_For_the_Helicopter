function [beta,beta_str,beta_2str] = Newmarks_beta(~,rev,mu,alpha_s,phi_s,lambda,theta_0,theta_1c,theta_1s)

v_B = 1.04 ;     %Blade flap frequency, per rev
Gamma = 8.0;    %Lock number
theta_tw = 0;   %Blade twist rate

%% Calculated data
d_psi = deg2rad(3); % Delta psi steps

% Newmark's Parameters
beta_N = (1/4);
gamma_N = (1/2);
%zer=rev*2*pi/d_psi + 1;
%% variable arrays declarations
M_beta = zeros(721,1); % initial conditions for Newmark's Integration
beta_0 = zeros(721,1);
beta_00 = zeros(721,1);
beta_000 = zeros(721,1);

%% Numerical Integration
x=0;
    for psi=0:d_psi:(rev*2*pi) % defined over whole azimuth by number of revolutions
         x=x+1;

        % implementing Newmark's Algorithm
        % M_beta value
        M_beta(x) = (1/8+ (mu/3)*sin(psi) + ...
            ((mu^2)/4)*(sin(psi)^2))*(theta_0+...
            theta_1c*cos(psi)+theta_1s*sin(psi))+...
            theta_tw*((1/10)+(mu^2)/6*(sin(psi)^2)+mu/4*sin(psi))...
            - (lambda*cos(beta_0(x)))*((1/6)+mu/4*sin(psi)) - ...
            beta_00(x)*((1/8)+(mu/6)*sin(psi))...
            - mu*sin(beta_0(x))*cos(psi)*((1/6)+(mu/4)*sin(psi));

        % Calculating initial value of Beta_2str
        if x==1
            beta_000(x) = Gamma*M_beta(x) - (v_B^2)*beta_0(x);
        end

        % After M_beta-(n) and  beta_2str-(n) is known, we find beta-(n+1), beta_str-(n+1), beta_2str-(n+1)
        beta_2s =(Gamma*M_beta(x)-(v_B^2)*...
            (d_psi*beta_00(x)+ beta_0(x)+(d_psi^2)*(1-2*beta_N)*beta_000(x)/2))/...
            (1+(v_B^2)*(d_psi^2)*beta_N);
        beta_000(x+1) = beta_2s;

        beta_s = beta_00(x)+ d_psi*((1-gamma_N)...
            *beta_000(x)+gamma_N*beta_000(x+1));
        beta_00(x+1) = beta_s;

        beta_o = beta_0(x)+ d_psi*beta_00(x)...
            + ((d_psi^2)/2)*((1-2*beta_N)*beta_000(x)...
            + 2*beta_N*beta_000(x+1));
        beta_0(x+1) = beta_o;

    end
    
xx = 480; yy = 600 ;
beta_0deg = beta_0(xx:yy);
beta_90deg = beta_0((xx+30):(yy+30));
beta_180deg = beta_0((xx+60):(yy+60));
beta_270deg = beta_0((xx+90):(yy+90));
beta_360deg = beta_0((xx+120):(yy+120));
beta = [beta_0deg; beta_90deg; beta_180deg; beta_270deg; beta_360deg] ;     

beta_str_0deg = beta_00(xx:yy);
beta_str_90deg = beta_00((xx+30):(yy+30));
beta_str_180deg = beta_00((xx+60):(yy+60));
beta_str_270deg = beta_00((xx+90):(yy+90));
beta_str_360deg = beta_00((xx+120):(yy+120));
beta_str = [beta_str_0deg; beta_str_90deg; beta_str_180deg; beta_str_270deg; beta_str_360deg] ;

beta_2str_0deg = beta_000(xx:yy);
beta_2str_90deg = beta_000((xx+30):(yy+30));
beta_2str_180deg = beta_000((xx+60):(yy+60));
beta_2str_270deg = beta_000((xx+90):(yy+90));
beta_2str_360deg = beta_000((xx+120):(yy+120));
beta_2str = [beta_2str_0deg; beta_2str_90deg; beta_2str_180deg; beta_2str_270deg; beta_2str_360deg] ;

end