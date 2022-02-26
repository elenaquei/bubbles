function [fn,DxFn,Dalphafn] = fn_Lorenz84(x,alpha)
% function [fn,DxFn,Dalphafn] = fn_Lorenz84(x,alpha)
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

fn = [ - Y^2 - Z^2 - alpha*X + alpha*F - gamma*U^2
    X*Y - beta*X*Z - Y + G
    beta*X*Y + X*Z - Z
    - delta*U + gamma*U*X + T];

if nargout<2 
    return
end

% actually never needed
DxFn =[ - alpha     -2*Y        -2*Z        -2*gamma*U
    Y-beta*Z        X-1         -beta*X     0
    beta*Y+Z        beta*X      X-1         0
    gamma*U         0           0           -delta+gamma*X];

Dalphafn = [0
    0
    0
    0];