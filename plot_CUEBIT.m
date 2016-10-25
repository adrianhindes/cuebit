function plot_CUEBIT( orbit)

figure;
hold on;
plot(orbit.t, orbit.E(1,:),'k');
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

