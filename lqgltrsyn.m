function [K_ltr,svdKltr,W1] = lqgltrsyn(usys)

%% Sintesi controllore LQG/LTR

%=====================================================%
%Parto dal sistema nominale
An = usys.NominalValue.A;
Bn = usys.NominalValue.B;
Cn = usys.NominalValue.C;
Dn = usys.NominalValue.D;
Bn = Bn(:,1);


 s = tf('s');
 
 Gn = simplify(Cn/(s*eye(4)-An)*Bn);
 Gn = Gn(:,1);
 %tzero(Gn)
 omega=logspace(-4, 6, 100);
%=====================================================%


%=====================================================% 
%Definizione dei bound

pw = (5.2*s^3 + 1723*s^2 + 3954*s + 436.5)/...
    (s^3 + 206.4*s^2 + 60.17*s + 1.206);

lm=(0.05457*s^3 + 0.1872*s^2 + 0.4375*s + 1.401)/...
    (s^3 + 13.95*s^2 + 14.38*s + 112.7);

vm = simplify(pw*inv(1-lm)); 
%=====================================================% 

%=====================================================%
%Tuning di H e rho
H = 10^2*eye(4);
rho = 0.1;

Phi = inv(s*eye(4)-An); %Phi
HPhiB = H*Phi*Bn;

figure(5)
sigma(1/(sqrt(rho))*HPhiB)
hold on
sigma(vm,'r',omega)
hold off
title('\sigma(1/sqrt(\rho)*H\PhiB) > pw/(1-lm)')
legend('1/sqrt(\rho)*H\PhiB','pw/(1-lm)')
grid

%=====================================================%
%% LQR
Q = H'*H;
R = rho;
omega=logspace(-4, 6, 10000);


Kc = lqr(usys.NominalValue(:,1), Q, R);
%=====================================================% 
%Verifica dei vincoli
figure(6)
sigma(Kc*Phi*Bn)
hold on 
title('Kc*Phi*Bn~1/(sqrt(\rho))*H\PhiBn')
sigma(1/(sqrt(rho))*HPhiB, 'r', omega)
legend('Kc*Phi*Bn','1/(sqrt(\rho))*H\PhiBn')
hold off

%=====================================================% 

%=====================================================% 
%Verifiche sulla robustezza

% margini di stabilità: sigma(I+inv(T_LQ))>0.5 (-6db),
% è come avere abs(lm)<0.5.
%Sto ipotizzando che le incertezze siano in ingresso al modello
T_LQ = Kc*Phi*Bn; %matrice di trasferimento desiderata: Tc (apro nel compensatore)
figure(7)
sigma(1 + inv(T_LQ))
grid 
title('\sigma(I+inv(T_{lq}))>0.5')
legend('\sigma(I+inv(T_{lq})')
hold on
yline(-6, 'r--', omega)
hold off


figure(8)
sigma(T_LQ * inv(1+(T_LQ)))
hold on
sigma(1/lm,'r', omega)
title('Tlq*inv(1+Tlq)<1/lm')
legend('T_{lq}*inv(1+T_{lq})','1/lm')
hold off
%=====================================================%

B = usys.B;
Gamma = B * B.';
mu = 0.1 * eye(4);

%=====================================================%
%% LTR

q = [1e5 1e6 1e7 1e8 1e9 1e10];
[K_ltr,svdKltr,W1] = ltrsyn(usys.NominalValue(:,1), Kc, Gamma, mu, q);
%=====================================================%