function g = function_g(x,epsilon,gsplit,K)
if nargin<3
    gsplit = 1E-4;
    K = 20;
end
ix=intval(x);    iepsilon=intval(epsilon);
if ~contains_zero(iepsilon) && sup(abs(iepsilon))>gsplit
    % Crude computation: if epsilon doesn't contain zero and we are far
    % from the removable singularity at epsilon=0, we use the direct
    % function evaluation and trust INTLAB to get a tight enclosure.
    g = 1/iepsilon*(exp(-iepsilon*ix)-1);
elseif sup(abs(ix*iepsilon))<1
    % This uses a geometric series bound to handle the case when |epsilon|
    % is small. Necessary, since epsilon is actually the amplitude
    % parameter and over a continuation it will cross zero continuously.
    % Some minor compute savings since if epsilon is very small, K can be
    % chosen to be equal to 4.
    base = inf(1/abs(ix*iepsilon));
    Kspecial = max(4,ceil(sup(log(1E14*abs(ix))/log(base))));
    g = -ix - (ix.^(2:Kspecial))*((-iepsilon).^(1:Kspecial-1)).' + midrad(0,1E-15);
else
    % In general, we can fall back to a Taylor series error bound by
    % treating epsilon as formally fixed, then getting a Taylor remainder
    % in terms of derivatives in x. Conditioning K to keep the error below
    % 1E-15 is tedious, so we don't.
    err = sup(1/factorial(K+1)*abs(ix)^(K+1)*abs(iepsilon)^K*exp(abs(iepsilon*ix)));
    g = -ix - (ix.^(2:K))*((-iepsilon).^(1:K-1)).' + midrad(0,err);
end
   if ~isintval(x) && ~isintval(epsilon)
       g = mid(g);
   end
end