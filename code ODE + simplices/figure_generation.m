% figure generatio

% roughly 3 hours - with validation online
% needs file stored_FHN_start.mat
FHN_simplex_Hopf

% runtime 10 hours - with validation a posteriori
% no prerequisite file
Fish_simplex_Hopf  % Fish ODE

% 10 hours, needs stored_lorenz84.mat, derivatives_Lorenz84.m and
% fn_Lorenz84.m
lorenz84_simplex_Hopf 


% the DDE example: SI model - takes from 3 to 10 days depending on
% computer, very memory intensive
SI_simplex_Hopf

% Bautin bifurcation
Bautin_simplex_Hopf