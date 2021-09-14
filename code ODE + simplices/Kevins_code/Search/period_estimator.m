function T = period_estimator(sol,k,ic,r)
sol_pe = ddesd(@(t,y,Z)DDEfun(t,y,Z,sol,k),@(t,y)DELAYS(t,y),[ic;deval(sol,sol.x(end)-r)],[sol.x(end)-r,sol.x(end)]);
T = sol_pe.y(1,end);
end

function xdot = DDEfun(t,y,Z,sol,k)
xdot = k*[(Z(2,1)-y(2))*(deval(sol,t)-y(2)); 
        -y(2)+deval(sol,t)];
end

function Z = DELAYS(t,y)
Z = t - y(1);
end