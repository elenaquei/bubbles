function [DxxFnV] = second_der_Hamiltonian(x,v,q)
% function [DxxFnV] = second_der_Hamiltonian(x,v,q)
%
% according to the system in 
% Switching to nonhyperbolic cycles from codim 2 bifurcations of equilibria in ODEs
% Yu.A. Kuznetsov, H.G.E. Meijer, W. Govaerts, B. Sautois

a = 2;
b = 3;
c = 1;
d = -1;

X = x(1); Y = x(2); Z = x(3); U = x(4);

v1=v(1); v2=v(2); v3=v(3); v4 = v(4);

DxxFnV = [0 0 0 0 
    0 0 0 0 
    0 0 0 0 
    4*c*v1+6*d*X*v1   0  0  0]*q;

