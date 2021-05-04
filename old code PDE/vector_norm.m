function Norm = vector_norm(w)
% function Norm = vector_norm(w)
%
% INPUT
% w         2D vector
% GLOBAL 
% nu        (nu1,n2)>(1,1)
% OUTPUT
% Norm      norm of the vector w

global nu

if length(nu)==1
    nu_loc = [nu,nu];
else
    nu_loc = nu;
end

span_k1 = (size(w,1)-1)/2;
span_k2 = (size(w,2)-1)/2;

NU_K1 = nu_loc(1).^(abs(-span_k1:span_k1));
NU_K2 = nu_loc(2).^(abs(-span_k2:span_k2));

Norm = sum(sum( NU_K1 *abs(w) *NU_K2'));