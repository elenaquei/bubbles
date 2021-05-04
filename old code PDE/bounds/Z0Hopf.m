function Z0 = Z0Hopf(A_dagger,X,G,Es)
% function Z0 = Z0Hopf(X,G,Es)
%
% INPUT
% X         PDE_vector
% G, Es     vectors
% OUTPUT
% Z0        vector

not_big = 0;
A = matrix_AHopf(X, G, Es, not_big);
%A_dagger =  Hopf_derivative(X,G,Es);
Z0 = norm_matrix(eye(size(A,1)) - A*A_dagger, X);