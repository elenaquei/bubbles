function [a0,a1,a2,phi,epsilon] = get_orbit(u,p,sol,T,N)
x = linspace(sol.x(end)-T,sol.x(end),10*N+1);
y = deval(sol,x);
phi = 2*pi/T;
s = linspace(-T/2,T/2,10*N+1);
a = get_fourier(s,y,phi,N,T/(10*N)) - u*[zeros(N,1);1;zeros(N,1)];
K = diag(-N:N);
epsilon = sqrt(real((K^2*a)'*a));
a0 = a/epsilon;
y1 = function_g((y-u)/epsilon*p,epsilon);
y2 = exp(-p*(y-u));
a1 = get_fourier(s,y1,phi,N,T/(10*N));
a2 = get_fourier(s,y2,phi,N,T/(10*N));
end