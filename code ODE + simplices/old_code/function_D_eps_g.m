function g = function_D_eps_g(x,epsilon,gsplit,K)
if nargin<3
    gsplit = 1E-4;
    K = 30;
end
ix=intval(x);    iepsilon=intval(epsilon);
if ~contains_zero(iepsilon) || sup(abs(iepsilon))>gsplit
    g = -x*exp(-epsilon*x)/epsilon - (exp(-epsilon*x)-1)/epsilon^2;
else
    err = sup(K/factorial(K+1)*abs(ix)^(K+1)*abs(iepsilon)^(K-1));
    g =  ix^2/2 + ((ix.^(3:K)).*(2:K-1))*((-iepsilon).^(1:K-2).') + midrad(0,err);
end
   if ~isintval(x) && ~isintval(epsilon)
       g = mid(g);
   end
end