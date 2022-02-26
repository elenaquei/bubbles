function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
    DalphaalphaxFn,DxxxFnV] = derivatives_Lorenz84(x,alpha,v,parameters)
%function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
%    DalphaalphaxFn,DxxxFnV] = derivatives_Lorenz84(x,alpha,v,parameters)
%
% 
% according to the system in 
% Switching to nonhyperbolic cycles from codim 2 bifurcations of equilibria in ODEs
% Yu.A. Kuznetsov, H.G.E. Meijer, W. Govaerts, B. Sautois

T = alpha;
F = 2;
alpha = 0.25;
beta =1 ;
G = 0.25;
delta = 1.04;
gamma = 0.987;

X = x(1); Y = x(2); Z = x(3); U = x(4);

DxFn =[ - alpha     -2*Y        -2*Z        -2*gamma*U
    Y-beta*Z        X-1         -beta*X     0
    beta*Y+Z        beta*X      X-1         0
    gamma*U         0           0           -delta+gamma*X];

Dalphafn = [0
    0
    0
    1];

if nargout<3
    return
end

v1=v(1); v2=v(2); v3=v(3); v4 = v(4);

DxxFnV = [0      -2*v2      -2*v3       -2*gamma*v4
    v2-beta*v3   v1         -beta*v1    0
    beta*v2+v3   beta*v1     v1         0
    gamma*v4     0           0          gamma*v1];

DalphaxFn = zeros(4);

DalphaxxFn=zeros(4);
DalphaalphaFn=zeros(4,1);
DalphaalphaxFn=zeros(4);

DxxxFnV=zeros(4);