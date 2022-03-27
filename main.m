%=====================================================%
% Pierpaolo Serio - Gennaio 2022
% Progetto 6 - Controllo robot su due ruote
%=====================================================%

clear all
clc

%% Inizializzazione del modello

mod_robot
olp_robot_2dof

%% CONTROLLORE MU

% Sintesi controllore MU
dms_robot_2dof

% Prestazioni controllore MU
dmu_robot_2dof

%% LQG/LTR Input

% Sintesi controllore LQG/LTR
[K_ltr,svdKltr,W1] = lqgltrsyn(G_unc);

% Prestazioni controllore LQG/LTR
lqgltr_prfm

%% LQG/LTR Output

[K_ltr,svdKltr,W1] = lqgltrsynout(G_unc);

% Prestazioni controllore LQG/LTR
lqgltr_prfm
