function DF = DF_xivector(X_sol,X_approx,delay)

n_nodes = X_sol.nodes;
K= -n_nodes:n_nodes;

% n_nodes = ((size(X_sol)*[1;0] - 7)/3 - 1)/2;

% define scalar unknowns
u=X_sol.scalar(1); 
p=X_sol.scalar(2); 
R0=X_sol.scalar(3); 
phi=X_sol.scalar(4); 
epsilon=X_sol.scalar(5); 
eta1=X_sol.scalar(6); 
eta2=X_sol.scalar(7);

% Fourier unknowns
z0=X_sol.vector(1,:).';  
z1=X_sol.vector(2,:).';  
z2=X_sol.vector(3,:).';

K_mat = diag(-n_nodes:n_nodes);

% Steady state
DF_equilib = [-1 + R0*exp(-p*u)*(1-2*u) - p*R0*exp(-p*u)*u*(1-u),...
    -u*R0*exp(-p*u)*u*(1-u), exp(-p*u)*u*(1-u), zeros(1,length(X_sol)-3)];
% Phase
DF_phase = [zeros(1,7),-1i*X_approx(8:8+2*n_nodes).'.*K,zeros(1,2*(2*n_nodes+1))]; 
% Amplitude
DF_amplitude = [zeros(1,7),X_approx(8:8+2*n_nodes).'.*K.^2,zeros(1,2*(2*n_nodes+1))];
% Boundary conditions
DF_BC1 = [0, sum(z0)*exp(-epsilon*p*sum(z0)), 0, 0, ...
    -function_D_eps_g(p*sum(z0),epsilon), 0, 0, ...
    ones(1,2*n_nodes+1)*p*exp(-epsilon*p*sum(z0)), ...
    ones(1,2*n_nodes+1), zeros(1,2*n_nodes+1)] ;  
DF_BC2 = [0, epsilon*sum(z0)*exp(-epsilon*p*sum(z0)), 0, 0, ...
    p*sum(z0)*exp(-epsilon*p*sum(z0)), 0, 0, ...
    ones(1,2*n_nodes+1)*epsilon*p*exp(-epsilon*p*sum(z0)), zeros(1,2*n_nodes+1), ...
    ones(1,2*n_nodes+1)];
% Vector field : Du
dtau = diag(exp(-1i*delay*phi*K));
DfPO_u0 = -p*R0*exp(-p*u)*(dtau*z1)*u*(1-u)...  
    + R0*exp(-p*u)*(dtau*z1)*(1-2*u)...
    - p*R0*exp(-p*u)*Convo((dtau*z2),(-epsilon*Convo(z0,z0)+(1-2*u)*z0))...
    + R0*exp(-p*u)*Convo((dtau*z2),-2*z0);
DfPO_u = [DfPO_u0 ; -p*Convo(z2, DfPO_u0) ; -epsilon*p*Convo(z2, DfPO_u0)];
% Vector field: Dp
f0 = -z0 + R0*exp(-p*u)*(dtau*z1*u*(1-u)...
    + Convo(dtau*z2,-epsilon*Convo(z0,z0)+(1-2*u)*z0));
f0_p = -u*(f0+z0);
DfPO_p = [zeros(2*n_nodes+1,1);-Convo(z2,f0);-epsilon*Convo(z2,f0)]...
    + [f0_p; -p*Convo(z2,f0_p); -epsilon*p*Convo(z2,f0_p)];
% Vector field: DR0
f0_R0 = exp(-p*u)*(dtau*z1)*u*(1-u)...
    + exp(-p*u)*Convo(dtau*z2,-epsilon*Convo(z0,z0) + (1-2*u)*z0);
DfPO_R0 = [f0_R0;-p*Convo(z2,f0_R0);-epsilon*p*Convo(z2,f0_R0)];
% Vector field: Dphi
ddtau = -1i*delay*K_mat*dtau;
f0_phi = R0*exp(-p*u)*ddtau*z1*u*(1-u)...
    + R0*exp(-p*u)*Convo(ddtau*z2,-epsilon*Convo(z0,z0)+(1-2*u)*z0);
DfPO_phi = -1i*[K.'.*z0;K.'.*z1;K.'.*z2]...
    + [f0_phi; -p*Convo(z2,f0_phi);-epsilon*p*Convo(z2,f0_phi)];
% Vector field: Depsilon
f0_eps = R0*exp(-p*u)*Convo(dtau*z2,-Convo(z0,z0));
DfPO_eps = [f0_eps; -p*Convo(z2,f0_eps); -epsilon*p*Convo(z2,f0_eps)]...
    + [zeros(2*(2*n_nodes+1),1); -p*Convo(z2,f0)];
% Vector field: Deta1 and Deta2
DfPO_eta1 = [zeros(2*n_nodes+1,1);zeros(n_nodes,1);1;zeros(n_nodes,1);zeros(2*n_nodes+1,1)];
DfPO_eta2 = [zeros(2*n_nodes+1,1);zeros(2*n_nodes+1,1);zeros(n_nodes,1);1;zeros(n_nodes,1)];
% Vector field: Da0
e0 = [zeros(n_nodes,1);1;zeros(n_nodes,1)];
df0_a0_star = -e0 + R0*exp(-p*u)*Convo(dtau*z2,-2*epsilon*z0 + (1-2*u)*e0);
DfPO_a0 = [-1i*phi*K + Convo_mat(df0_a0_star);...
    -p*Convo_mat(Convo(z2,df0_a0_star));...
    -epsilon*p*Convo_mat(Convo(z2,df0_a0_star))];
% Vector field: Da1
df0_a1_star = R0*exp(-p*u)*dtau*u*(1-u);
DfPO_a1 = [df0_a1_star;...
    -1i*phi*K - p*Convo_mat(z2)*df0_a1_star;...
    -epsilon*p*Convo_mat(z2)*df0_a1_star];
% Vector field: Da2
df0_a2_star = R0*exp(-p*u)*(-epsilon*Convo(z0,z0)+(1-2*u)*z0);
DfPO_a2 = [Convo_mat(df0_a2_star)*dtau;...
    -p*Convo_mat(f0) - p*Convo_mat(Convo(z2,df0_a2_star))*dtau;...
    -1i*phi*K - epsilon*p*Convo_mat(f0)...
    - epsilon*p*Convo_mat(Convo(z2,df0_a2_star))*dtau];
% Gather all blocks
DF_vfield = [DfPO_u, DfPO_p, DfPO_R0,DfPO_phi,DfPO_eps,DfPO_eta1,...
    DfPO_eta2,DfPO_a0,DfPO_a1,DfPO_a2];
DF = [DF_equilib; DF_phase; DF_amplitude; DF_BC1; DF_BC2; DF_vfield];
end