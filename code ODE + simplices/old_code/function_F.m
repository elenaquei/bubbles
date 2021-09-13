function F = function_F(X,Xs,tau,L,K)
% Some comments added for Elena - this code was not very commented.
N = ((size(X)*[1;0] - 7)/3 - 1)/2;
% Steady state
f = -X(1) + X(3)*exp(-X(2)*X(1))*X(1)*(1-X(1));
% Phase
Phase = -1i*Xs(8:8+2*N)'*(K*X(8:8+2*N));
% Amplitude
Amplitude = Xs(8:8+2*N)'*(K^2*X(8:8+2*N))-1;
% Boundary condition
BC = [sum(X(9+2*N:9+4*N)) - function_g(X(2)*sum(X(8:8+2*N)),X(5)) ; 
     sum(X(10+4*N:end)) - exp(-X(5)*X(2)*sum(X(8:8+2*N)))];
% Periodic orbit
u=X(1); p=X(2); R0=X(3); phi=X(4); epsilon=X(5); eta1=X(6); eta2=X(7);
a0=X(8:8+2*N);  a1=X(9+2*N:9+4*N);  a2=X(10+4*N:end);   % EPSILON = AMPLITUDE. PHI = OMEGA0. Fourier are a0,a1,a2 instead of z0,z1,z2.
dtau = diag(exp(-1i*tau*phi*(-N:N)));
vfield_0 = -a0 + R0*exp(-p*u)*u*(1-u)*(dtau*a1) + R0*exp(-p*u)*( ...
    Convo(dtau*a2,-epsilon*Convo(a0,a0)+(1-2*u)*a0) );
fPO = -1i*phi*[K*a0;K*a1;K*a2] + ...
    [vfield_0; -p*Convo(a2,vfield_0); -p*epsilon*Convo(a2,vfield_0)] + ...
    [zeros(2*N+1,1);zeros(N,1);eta1;zeros(N,1);zeros(N,1);eta2;zeros(N,1)];
% Combine everything
F = [f;Phase;Amplitude;BC;fPO];
end