% plot_solovev

NR = 128;
rgrid = 0:2/(NR-1):2;
zgrid = -2:2/(NR-1):2;
[rmesh, zmesh] = meshgrid(rgrid, zgrid);

psi0= 0.9;
gamma = 0.8;
rmag = 0.964;
rbean= 0.93;
psi = psi0 * ( gamma/8 * ((rmesh.^2 - rmag^2).^2  - rbean^4) + (1-gamma)/2 * rmesh.^2 .* zmesh.^2 );

Npsi = 20;
psi_min   = min(min(psi));
psi_surfs = [psi_min:abs(psi_min/(Npsi-1)):0];

figure;
hold on;
contour(rmesh, zmesh, psi, psi_surfs);
daspect([1 1 1])
% plot(gs2.rlcf, gs2.zlcf ,'r','LineWidth',2);

load orbit.dat
plot(orbit(:,1),orbit(:,2)),xlabel('R(m)'), ylabel('Z(m)')
plot(orbit(1:1000,1),orbit(1:1000,2)),xlabel('R(m)'), ylabel('Z(m)')

% ,xlabel('R(m)'), ylabel('Z(m)')