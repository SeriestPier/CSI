%% Genero 20 punti per limitare superiormente le incertezze moltiplicative relative creando una risposta in frequenza W_delta 

% Fits a logarithmic frequency response with a stable,
% minimum phase transfer function
%
[freq,resp_db] = ginput(20);       % pick 30 points
for i = 1:20                       % Converts the logarithmic response
    resp(i) = 10^(resp_db(i)/20);  % to a magnitude response
end                                %
sys = frd(resp,freq);              % creates frd object
ord = 3;                           % third order fit
W = fitmagfrd(sys,ord);            % fits the frequency response
Wm = tf(W)                         % converting into transfer function form