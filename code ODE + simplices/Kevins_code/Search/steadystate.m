function u = steadystate(p,R0,ic)
if nargin<3
    ic = rand;
end
f = @(x)-x + R0*exp(-p*x)*x*(1-x);
Df = @(x) -1 + R0*(-p)*exp(-p*x)*x*(1-x) + R0*exp(-p*x)*(1-2*x);
u = ic;
while abs(f(u))>1E-15 & abs(f(u))<1E4
    u = u - f(u)/Df(u);
end
end

