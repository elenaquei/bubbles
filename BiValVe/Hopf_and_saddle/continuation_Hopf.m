function [s, x_n,sol_2] = continuation_Hopf ( lambda0, x0, f, n_nodes, n_iter,...
    h, s, eigenvec, eigenval, sign_FLC, index_saddle, bool_saddle)
% function [s, x_n, x_0] = continuation_Hopf ( lambda0, x0, f, n_nodes, n_iter, h, s, eigenvec, eigenval, sign_FLC, index_saddle,bool_saddle)
%
% INPUTS
% lambda0       paramenter of the Hopf bifurcation
% x0            scalar values of the Hopf bifurcation
% f             vector field defining the problem
% n_nodes       integer, number of nodes to use for Fourier
% n_iter        number of validated continuation iterations to do (DEFAULT:
%               100)
% h             step size in the continuation code
% s             string, path where to save results
% eigenvec      approximation of the eigenvector at the Hopf bifurcation
% eigenval      approximation of the eigenvalue at the Hopf bifurcation

% index_saddle  variable of the Hopf bifurcation
%
% OUTPUTS
% s              string, path of saved solutions
% x_n            Xi_vector, last validated solution
global use_intlab

temp_intval = use_intlab;
use_intlab = 0;
if nargin < 11 || isempty(index_saddle)
    index_saddle = 1:(3+length(x0));
end
if nargin<12
    bool_saddle = 1; % also validate the existence of a saddle node
end

[sol_2, big_Hopf, x_dot_0] = Hopf_system_setup (lambda0, x0, f, n_nodes,...
    eigenvec, eigenval, sign_FLC, h);

% calling continuation with the Hopf boolean (guarantees that we keep using
% the same scalar conditions over the continuation)
bool_Hopf = 1;
bool_fancy_scalar = 0; % already included in the Hopf boolean

use_intlab = temp_intval;

[s, x_n] = continuation ( sol_2, big_Hopf, n_iter, h, x_dot_0,s, 10^-6, ...
    bool_Hopf, bool_fancy_scalar, bool_saddle, index_saddle);

end