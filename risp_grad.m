%% Risposta al gradino
function risp_grad(sys,i)
figure(i)
step(gridureal(sys(:,1), 50), 'm--') %sostituisce 40 campioni equispaziato con i parametri incerti di sys
hold on
grid

%Plotto solo dal primo ingresso
step(sys.NominalValue(:, 1), 'c')
legend('50 campioni incerti', 'Risposta Nominale')
hold off
