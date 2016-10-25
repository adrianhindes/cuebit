
function plot_GS2(gs2)

% MJH 20/09/02
% plots gs2 information 


figure;
subplot('Position',[0.1  0.3  0.3  0.65])
hold on;

[C,h]=contour(gs2.rmesh,gs2.zmesh,gs2.psi,20,'k');


plot(gs2.rlcf, gs2.zlcf ,'b','LineWidth',2);

title(gs2.header2);
daspect([1 1 1]);
xlabel('R [m]');
ylabel('Z [m]');
hold off;

subplot('Position',[0.1  0.05 0.3  0.15])
hold on;
% izero=fix(gs2.NZBOX/2)+1;
% plot(gs2.rmesh(:,izero),gs2.PSI(:,izero));
% axis([0 gs2.RBOXLEN+gs2.RBOXLFT min(gs2.PSI(:,izero)) 1.1*max(gs2.PSI(:,izero))]);
% plot(min(gs2.bound(:,1)),gs2.PSILCF,'ro');
% plot(max(gs2.bound(:,1)),gs2.PSILCF,'ro');
% plot(gs2.RAXIS,gs2.PSIAXIS,'ro');
ylabel('\psi');
xlabel('R [m]');
hold off;

subplot('Position',[0.55  0.825  0.4  0.15])
plot(gs2.psi1D, gs2.q,'b');
ylabel('q');


subplot('Position',[0.55  0.625  0.4  0.15])
plot(gs2.psi1D, gs2.p);
ylabel('P [Pa]');

subplot('Position',[0.55  0.425  0.4  0.15])
plot(gs2.psi1D, gs2.f);
ylabel('f [T m]');

subplot('Position',[0.55  0.25  0.4  0.15])
plot(gs2.psi1D, gs2.pp);
ylabel('p'' [Pa/Wb]');


return;

