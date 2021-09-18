function F = function_F(X,Xs,tau,L,K)
% Some comments added for Elena - this code was not very commented.
n_nodes = X_sol.n_nodes;
K= -n_nodes:n_nodes;

% n_nodes = ((size(X_sol)*[1;0] - 7)/3 - 1)/2;

if isa(X,'Xi_vector')
    % define scalar unknowns
    u=X_sol.scalar(1);
    p=X_sol.scalar(2);
    R0=X_sol.scalar(3);
    phi=X_sol.scalar(4);
    amplitude=X_sol.scalar(5);
    eta1=X_sol.scalar(6);
    eta2=X_sol.scalar(7);
    
    % Fourier unknowns
    z0=X_sol.vector(1,:);
    z1=X_sol.vector(2,:);
    z2=X_sol.vector(3,:);
else
    u=X(1); p=X(2); R0=X(3); phi=X(4); amplitude=X(5); eta1=X(6); eta2=X(7);
    z0=X(8:8+2*N);  z1=X(9+2*N:9+4*N);  z2=X(10+4*N:end);
end
% PHI = OMEGA0

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
dtau = diag(exp(-1i*tau*phi*(-N:N)));
vfield_0 = -z0 + R0*exp(-p*u)*u*(1-u)*(dtau*z1) + R0*exp(-p*u)*( ...
    Convo(dtau*z2,-amplitude*Convo(z0,z0)+(1-2*u)*z0) );
fPO = -1i*phi*[K*z0;K*z1;K*z2] + ...
    [vfield_0; -p*Convo(z2,vfield_0); -p*amplitude*Convo(z2,vfield_0)] + ...
    [zeros(2*N+1,1);zeros(N,1);eta1;zeros(N,1);zeros(N,1);eta2;zeros(N,1)];

% Combine everything
F = [f;Phase;Amplitude;BC;fPO];
end