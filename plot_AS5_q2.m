%% Plot
close all;

figure(1)
plot(Rec_mu,rad2deg(Rec_theta_0),'r-.')
hold on
plot(Rec_mu,rad2deg(Rec_theta_1c),'g-.')
plot(Rec_mu,rad2deg(Rec_theta_1s),'b-.')
xlabel('\mu')
ylabel('Control Input angles [deg]')
title(' Variation of Control angles vs \mu');
legend('\theta_0','\theta_1_c','\theta_1_s');
ylim([-30 40])

figure(2)
plot(Rec_mu,rad2deg(Rec_alpha_s),'r-.')
hold on
plot(Rec_mu,rad2deg(Rec_phi_s),'b-.')
xlabel('\mu')
ylabel('Vehicle shaft angles [deg]')
title('Variation of Vehicle shaft angles vs \mu');
legend('\alpha_s','\phi_s');


%%