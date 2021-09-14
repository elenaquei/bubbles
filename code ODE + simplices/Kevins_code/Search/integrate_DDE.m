function sol = integrate_DDE(p,R0,tau,ic)
if nargin < 4
    ic = rand;
end
sol = dde23(@(t,y,Z)DDEfun(t,y,Z,p,R0),tau,ic,[0,1E5]);
end

function xdot = DDEfun(t,y,Z,p,R0)
xdot = -y + R0*exp(-p*Z)*y*(1-y);
end