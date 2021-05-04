function [Yvector,fullY]=Y_bound_new(A,xBar,alpha)
% function Yvector=Y_bound(A,xBar,alpha)
%
% INPUT
% A      complex matrix, approximation of the inverse of the function of
%        the problem;
% xBar   Xi_vector, numerical solution of the ODE system;
% alpha  full_problem
%
% OUTPUT
% Yvector   a standard vector of length = size_vector+size_scalar, the
%           values of this vector come from the Y bound in the associated pdf (we
%           refer to rigorous numerics for analytic solution od D.E..... for the
%           notation, mainly)

global nu
global talkative

if ~isa(alpha,'full_problem')
    error('Wrong input')
end

%Yvector=zeros(xBar.size_vector+xBar.size_scalar,1);
F_bigXi=apply(alpha,xBar,1);
%F_bigXi=F_function(xBar,alpha,coefs_linear,0);  % 0 => maximum size
F_bigvec=Xi_vec2vec(F_bigXi);% maximum size         % Xi_vec2vec

A_big = extend_approximate_inverse(A,xBar.size_scalar, xBar.size_vector,xBar.nodes,xBar.nodes*alpha.vector_field.deg_vector);

fullY=vec2Xi_vec(A_big*F_bigvec,xBar.size_scalar,xBar.size_vector,xBar.nodes*alpha.vector_field.deg_vector);

Yvector = cnorm_Xi_vector(fullY,nu);

% stop everything if the result is too big
if any(Yvector>1)
    fprintf('Y computed, %d\n',max(Yvector));
%    error('Y is bigger than 1, no interval found')
elseif talkative>2
    fprintf('Y computed, %d\n',max(Yvector));
    fprintf('\n');
elseif talkative>1
    fprintf('Y computed, %d\n',max(Yvector));
    fprintf('\n');
end
return
end

