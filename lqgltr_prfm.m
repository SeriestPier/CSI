%=====================================================%
%Analisi MU
%=====================================================%

%Dichiarazioni preliminari
G_uncd = G_unc(:,1);
T = K_ltr*G_unc.NominalValue(:,1);

ltr_sys_cl = lft(G_uncd,K_ltr,1,4);

looptransfer=loopsens(G_uncd,K_ltr);

omega=logspace(-4, 6, 100); 

ltr_clp=looptransfer.Ti; % '.Ti' funzione di sensibilità complementare

ltr_clp_g=ufrd(ltr_sys_cl, omega); %calcola la risposta in frequenza di un modello incerto.


%=====================================================%

%% Analisi stabilità nominale
figure(11)
sv = sigma(ltr_clp.NominalValue,omega);
sys = frd(sv(1,:),omega);
semilogx(sys,'r-')
grid
xlabel('Frequency (rad/sec)')
title('Nominal performance')
%=====================================================%

%% Analisi stabilità robusta
opt = robopt('Display','on');
[stabmarg,destabu,report,info] = robuststab(ltr_clp_g,opt);
report;

figure(12)
semilogx(info.MussvBnds(1,1),'g-',info.MussvBnds(1,2),'b--')
grid

title('Robust stability')
xlabel('Frequency (rad/sec)')
ylabel('mu')
legend('\mu-upper bound','\mu-lower bound')
%=====================================================%

%% Analisi prestazioni robuste
opt = robopt('Display','on');
[perfmarg,perfmargunc,report,info] = robustperf(ltr_clp,opt);
report;
figure(13)
semilogx(info.MussvBnds(1,1),'r-',info.MussvBnds(1,2),'b--')
grid
xlabel('Frequency (rad/sec)')
ylabel('mu')
title('Robust performance')
legend('\mu-upper bound','\mu-lower bound')

%Risp gradino
[WCG,WCU]=wcgain(ltr_clp)
risp_grad(looptransfer.PSi,14)
%risp_grad(ltr_sys_cl,14)
%=====================================================%