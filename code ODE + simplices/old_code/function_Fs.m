function F = function_Fs(X,Xs,PHIs,tau,L,K)
N = ((size(X)*[1;0] - 7)/3 - 1)/2;
% Steady state
f = -X(1) + X(3)*exp(-X(2)*X(1))*X(1)*(1-X(1));
% Simplex
Simplex = PHIs.'*(X-Xs);
% Phase
Phase = (L.phase(X(8:8+2*N)))'*Xs(8:8+2*N);
% Amplitude
Amplitude = (L.amplitude(X(8:8+2*N)))'*Xs(8:8+2*N)-1;
% Boundary condition
BC = [sum(X(9+2*N:9+4*N)) - function_g(sum(X(8:8+2*N)),X(5)) ; 
     sum(X(10+4*N:end)) - exp(-X(5)*X(2)*sum(X(8:8+2*N)))];
% Periodic orbit
u=X(1); p=X(2); R0=X(3); phi=X(4); epsilon=X(5); eta1=X(6); eta2=X(7);
a0=X(8:8+2*N);  a1=X(9+2*N:9+4*N);  a2=X(10+4*N:end);
dtau = diag(exp(-1i*tau*phi*(-N:N)));
a2_square = Convo(a2,a2);
vfield_0 = -a0 + R0*exp(-p*u)*u*(1-u)*(dtau*a0) + R0*exp(-p*u)*( ...
    -epsilon*Convo(Convo(a0,a0),dtau*a2_square) + ...
    +(1-2*u)*Convo(a0,dtau*a2_square) );
fPO = -1i*phi*[K*a0;K*a1;K*a2] + ...
    reshape(vfield_0*[1,-p,-epsilon*p],3*(2*N+1),1) + ...
    [zeros(2*N+1,1);zeros(N,1);eta1;zeros(N,1);zeros(N,1);eta2;zeros(N,1)];
% Combine everything
F = [Simplex;f;Phase;Amplitude;BC;fPO];
end