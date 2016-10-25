% plot_equilibrium

NR = 128;
rgrid = 0:2/(NR-1):2;
zgrid = -2:2/(NR-1):2;
[rmesh, zmesh] = meshgrid(rgrid, zgrid);

psi0= 0.9;
gamma = 0.73;
rmag = 1.047;
rbean= 1.017;
psi = psi0 * ( gamma/8 * ((rmesh.^2 - rmag^2).^2  - rbean^4) + (1-gamma)/2 * rmesh.^2 .* zmesh.^2 );

Npsi = 20;
psi_min   = min(min(psi));
psi_surfs = [psi_min:abs(psi_min/(Npsi-1)):0];

figure;
hold on;
contour(rmesh, zmesh, psi, psi_surfs);
daspect([1 1 1])

% ,xlabel('R(m)'), ylabel('Z(m)')