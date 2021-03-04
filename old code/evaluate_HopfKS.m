function Hopf_PDE = evaluate_HopfKS(x_PDE, G, const_G, Es, const_Es, flagFullSize)
% function Hopf_PDE = evaluate_HopfKS(x_PDE, G, const_G,Es,const_Es, flagFullSize)
%
% INPUT
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
% Es        PDE_vector or vector compatible with the previous PDE_vector,
%           continuation equation
% const_Es  contant in the Es equation
% G         matrix of 2 lines, each line compatible with the
%           PDE_vector
% const_G   2 elements vector, constants in the G equation
% flag      (DEFAULT 0) full size of convolutions or not
% OUTPUT
% Hopf_PDE  PDE_vector solution of the Hopf system in the Kuramoto Shivashinki

% WORKS (also linear parts)
if x_PDE.size_real ~= 3
    error('Input does not respect dimensional contraints - parameters')
end
if x_PDE.size_vector ~= 1
    error('Input does not respect dimensional contraints - vector')
end
if nargin<6 || isempty(flagFullSize)
    flagFullSize=0;
end

Hopf_PDE = x_PDE;

Hopf_PDE.parameters(1) =  compute_linear_func(PDE_vec2vec(x_PDE), Es.', const_Es);
Hopf_PDE.parameters(2:3) = compute_linear_func(PDE_vec2vec(x_PDE), G(1:2,:), const_G(1:2));

central_node_vec1 =  compute_linear_func(PDE_vec2vec(x_PDE), G(3,:), const_G(3));
central_node_vec2 = compute_linear_func(PDE_vec2vec(x_PDE), G(4,:), const_G(4));

if flagFullSize
    x_PDE = reshape(x_PDE, 2*x_PDE.node_space, 2* x_PDE.node_time);
    Hopf_PDE = reshape(Hopf_PDE, 2*Hopf_PDE.node_space, 2* Hopf_PDE.node_time);
end

y = x_PDE.vector1D;
z = squeeze(x_PDE.vector2D);
nu = x_PDE.parameters(1);
lambda = x_PDE.parameters(2);
a = x_PDE.parameters(3);

%nonlinear_y = -  derivative_space(Convo(y,y));
node = (length(y) - 1)/2;

K = diag(-node:node);

y_conv_y = conv(y,y,'same');
Hopf_PDE.vector1D  = (- nu * K^4* y(:) + K^2*y(:) + 1i*K*y_conv_y(:)).';

Hopf_PDE.vector1D(x_PDE.node_space+1) = central_node_vec1;

nonlinear_z =  lambda * derivative_space( 2 * Convo(z,y) + a * Convo(z,z));

vector2D_temp = - derivative_time(z) - (lambda*nu)*derivative_space(z,4) - lambda*derivative_space(z,2) + nonlinear_z;
%vector2D = zeros([1,size(vector2D_temp)]);
%vector2D = shiftdim(vector2D_temp,-1);

vector2D_temp(x_PDE.node_time+1,x_PDE.node_space+1) = central_node_vec2;

Hopf_PDE.vector2D = shiftdim(vector2D_temp,-1);

end

function x_prime = derivative_space(x,n)
if nargin ==1
    n = 1;
end
if min(size(x))==1 && size(x,1)~=1
    x = x.';
end
nodes = (size(x,2) - 1)/2;
K = -nodes:nodes;
K = repmat(K,size(x,1),1);
%K = diag(K(:));
x_prime = (1i*K).^n .* x;
end

function x_prime = derivative_time(x,n)
if nargin ==1
    n = 1;
end
if min(size(x))==1
    x_prime = x;
    return
end
nodes = (size(x,1) - 1)/2;
K = -nodes:nodes;
K = repmat(K(:),1,size(x,2));
%K = diag(K(:));
x_prime = (1i*K).^n .* x;
end