function plot_orbit( gs2)

figure;
hold on;

[C,h]=contour(gs2.rmesh,gs2.zmesh,gs2.psi,20,'k');
plot(gs2.rlcf, gs2.zlcf ,'r','LineWidth',2);

load orbit.dat
plot(orbit(:,1),orbit(:,2));
xlabel('R(m)');
ylabel('Z(m)');
daspect([1 1 1])

end

