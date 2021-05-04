function y = KS(arg1, arg2, arg3, arg4)
% function y = KS(nu, vector1D, G, const_G)
% OR function y = KS(vector1D, G, const_G)
%       with nu_Hopf global variable
%
% INPUTS
% (nu        real - not necessary if exists global variable nu_Hopf)
% vector1D  Fourier series
% G         vector of length vector1D + 1
% const_G   real
% OUTPUT
% y         Fourier series, result of the Kuramoto-Sivashinki stationary
%           problem

% tested comparing to KS_loc and KS_all_inputs in test_space_KS
global nu_Hopf

if nargin == 4
    nu= arg1;
    vector1D = arg2;
    G = arg3;
    const_G = arg4;
elseif nargin==3
    nu = nu_Hopf;
    vector1D = arg1;
    G = arg2;
    const_G = arg3;
end

node = (length(vector1D) - 1)/2;

K = diag(-node:node);

x_conv_x = conv(vector1D,vector1D,'same');
y = - nu * K^4* vector1D(:) + K^2*vector1D(:) + 1i*K*x_conv_x(:);

if length(G) == 2*node+1
    y(node+1) =sum(G(:).*vector1D(:))+const_G;
elseif length(G) == 2*node+2
    y(node+1) =G(1)*nu + sum(G(2:end).*vector1D(:))+const_G;
else
    error('Length of G incompatible with other inputs')
end
