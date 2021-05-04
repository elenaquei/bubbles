function Df = DKS(arg1, arg2, arg3)
% function y = DKS(nu, vector1D, G)
% OR function y = DKS(vector1D, G)
%       with nu_Hopf global variable
%
% INPUTS
% (nu        real - not necessary if exists global variable nu_Hopf)
% vector1D  Fourier series
% G         vector of length vector1D + 1
% OUTPUT
% Df        operator in Fourier series, derivative of the 
%           Kuramoto-Sivashinki stationary problem

% tested
global nu_Hopf

if nargin == 3
    nu= arg1;
    vector1D = arg2;
    G = arg3;
elseif nargin==2
    nu = nu_Hopf;
    vector1D = arg1;
    G = arg2;
end

node = (length(vector1D)-1)/2;

K = diag(-node:node);

x_first = vector1D(node+1:end);
x_second = vector1D(node+1:-1:1);

x_first = [x_first, 0*x_first(2:end)];
x_second = [x_second, 0*x_second(2:end)];

Df_y = - nu * K^4 + K^2 + 1i * 2 * K * toeplitz(x_first, x_second);

Df = Df_y;
Df(1+node,:)= G(:);
end