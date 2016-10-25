function [orbit] = bin_CUEBIT(orbit)

% pdf for E ==========================================

E_min    = min(orbit.E(1,:));
E_max    = max(orbit.E(1,:));

% -------------------------------------------
NE       = 50;
dE       = (E_max - E_min)/(NE-1);

orbit.E_bins1  = E_min: dE :E_max;

orbit.PE1      = hist(orbit.E(1,:),orbit.E_bins1);
orbit.PE1      = orbit.PE1/ (sum(orbit.PE1) * dE);

% -------------------------------------------
NE       = 100;
dE       = (E_max - E_min)/(NE-1);

orbit.E_bins2  = E_min: dE :E_max;

orbit.PE2      = hist(orbit.E(1,:),orbit.E_bins2);
orbit.PE2      = orbit.PE2/ (sum(orbit.PE2) * dE);

figure;
hold on;
semilogy(orbit.E_bins1, orbit.PE1,'k');
semilogy(orbit.E_bins2, orbit.PE2,'r');
title(['Rcoil = ',num2str(orbit.Rcoil)]);

set(gca, 'YScale','log');



% pdf for theta ======================================

eps   = 0.05;
i_index = find(orbit.E(1,:)< (1-eps) * orbit.E(1,1))
orbit.E(:, i_index) = NaN;

orbit.Rmag = 9.8246E-01;
orbit.Zmag = -2.8423E-04;

orbit.E(4,:) = atan2(orbit.E(3,:) - orbit.Zmag,  orbit.E(2,:) - orbit.Rmag);

theta_min    = min(orbit.E(4,:));
theta_max    = max(orbit.E(4,:));


% -------------------------------------------
Ntheta       = 50;
dtheta       = (theta_max - theta_min)/(Ntheta-1);

orbit.theta_bins1  = theta_min: dtheta :theta_max;
orbit.Ptheta1      = hist(orbit.E(4,:),orbit.theta_bins1);
orbit.Ptheta1      = orbit.Ptheta1/ (sum(orbit.Ptheta1) * dtheta);

% -------------------------------------------
Ntheta       = 100;
dtheta       = (theta_max - theta_min)/(Ntheta-1);

orbit.theta_bins2  = theta_min: dtheta :theta_max;
orbit.Ptheta2      = hist(orbit.E(4,:),orbit.theta_bins2);
orbit.Ptheta2      = orbit.Ptheta2/ (sum(orbit.Ptheta2) * dtheta);

% -------------------------------------------
% losses versus time
Nt = 100;
orbit.E(:, Nt+1:length(orbit.E))   = NaN;
Ntheta             = 50;
theta_min    = min(orbit.E(4,:));
theta_max    = max(orbit.E(4,:));
dtheta       = (theta_max - theta_min)/(Ntheta-1);

orbit.theta_bins3  = theta_min: dtheta :theta_max;
orbit.Ptheta3      = hist(orbit.E(4,:),orbit.theta_bins1);
orbit.Ptheta3      = orbit.Ptheta3/ (sum(orbit.Ptheta3) * dtheta);




figure;
hold on;
plot(orbit.theta_bins1, orbit.Ptheta1,'k');
plot(orbit.theta_bins2, orbit.Ptheta2,'r');
plot(orbit.theta_bins3, orbit.Ptheta3,'b');
title(['Rcoil = ',num2str(orbit.Rcoil)]);

return


xlabel(['t (s)']);
ylabel(['E (keV)']);
title(['Rcoil = ',num2str(orbit.Rcoil)]);

figure;
hold on;

plot(orbit.E(2,:), orbit.E(3,:),'k.');
xlabel('R(m)');
ylabel('Z(m)');
daspect([1 1 1]);

eps     = 0.005;
i_uconf = find((orbit.E(1,1)< (1-eps) * orbit.E(1,:))|(orbit.E(1,1) > (1+eps) * orbit.E(1,:)));
plot(orbit.E(2,i_uconf), orbit.E(3,i_uconf),'r.','LineWidth',2);
title(['Rcoil = ',num2str(orbit.Rcoil)]);

return

figure;
hold on;

plot(orbit.E(2,:), orbit.E(1,:),'kx');
xlabel('R(m)');
ylabel('E(keV)');

figure;
hold on;

plot3(orbit.E(1,:), orbit.E(2,:),orbit.E(3,:))

end

