



%sigma(Wn);
pw = (5.2*s^3 + 1723*s^2 + 3954*s + 436.5)/...
    (s^3 + 206.4*s^2 + 60.17*s + 1.206);
lm=(0.05457*s^3 + 0.1872*s^2 + 0.4375*s + 1.401)/...
    (s^3 + 13.95*s^2 + 14.38*s + 112.7);

vm = simplify(pw*inv(1-lm)); 

% figure(1)
% sigma(1/(sqrt(rho))*HFB)
% hold on
% sigma(vm,'r',omega)
% hold off
% title('\sigma(1/sqrt(\rho)*H\PhiB) > pw/(1-lm)')
% legend('1/sqrt(\rho)*H\PhiB','pw/(1-lm)')
% grid




% sigma (pw)
% hold on
% sigma (Wp)
% hold on
% sigma (Wu)
% hold on
% sigma (lm)