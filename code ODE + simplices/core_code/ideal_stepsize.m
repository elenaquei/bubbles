function delta = ideal_stepsize(Imin, Imax, Ys, Z0s, Z1s, Z2s)
% function tau = ideal_stepsize(Imin, Imax, Ys, Z0s, Z1s, Z2s)
%
% computation of the ideal stepsize given the bounds and minimal and
% maximal validation radius
%
% INPUT
% Imin, Imax            radii
% Ys, Z0s, Z1s, Z2s     all the bonds, each a matrix with 4 columns
% OUTPUT
% tau                   ideal stepsize 

r_center = (Imin+Imax)/2;

Y0 = (Ys(:,4));
Ydelta = (Ys(:,2));

Z0 = (Z1s(:,4)+Z0s(:,4));
Zdelta = (Z1s(:,2)+Z0s(:,2));

Z2 = (Z2s(:,2));

tau = sqrt( - (Y0 + (Z0 -1)* r_center + Z2 * r_center^2)./(Ydelta +(Zdelta)*r_center));

delta = min((0.9 + tau)/2); % rescaling to "calm down" the stepsize

delta = min(delta, 1.2);
