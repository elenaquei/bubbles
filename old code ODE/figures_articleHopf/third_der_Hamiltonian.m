function DDDF = third_der_Hamiltonian(~, v1,v2,v3)
% function DDDF = third_der_Hamiltonian(x, v1,v2,v3)

% DDDF = D^3_x F(x) v1 v2 v3 
% directional derivative! 

a = 2;
b = 3;
c = 1;
d = -1;


DDDF=[0
    0
    0
    6*d*v1(1)*v2(1)]*v3(1);