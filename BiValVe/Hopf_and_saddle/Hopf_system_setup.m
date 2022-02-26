function [sol_2, big_Hopf, x_dot_0] = Hopf_system_setup (lambda0, x0, f, n_nodes,...
    eigenvec, eigenval, sign_FLC, step_size)
% function [sol_2, big_Hopf, x_dot_0] = Hopf_system_setup (lambda0, x0, f, n_nodes,...
%     eigenvec, eigenval, sign_FLC, step_size)
% INPUTS
% lambda0       paramenter of the Hopf bifurcation
% x0            scalar values of the Hopf bifurcation
% f             vector field defining the problem
% n_nodes       integer, number of nodes to use for Fourier
% n_iter        number of validated continuation iterations to do (DEFAULT:
%               100)
% h             step size in the continuation code
% eigenvec      approximation of the eigenvector at the Hopf bifurcation
% eigenval      approximation of the eigenvalue at the Hopf bifurcation

% index_saddle  variable of the Hopf bifurcation
%
% OUTPUTS
% s              string, path of saved solutions
% x_n            Xi_vector, last validated solution

if length(x0)<2
    error('Need at least 2 dimensions to have a Hopf bifurcation')
end

% transformation of the vector field into the bigger system needed for Hopf
big_Hopf = Taylor_series_Hopf(f,n_nodes);

% analytic contruction of the solution at the Hopf bifurcation
T_star = imag(eigenval);
a_star = (abs(step_size));
if nargin>9 && ~isempty(sign_FLC)
    lambda0 = lambda0+ sign_FLC*a_star^2;
end

% guess of an approximation of the solution (assuming the eigenvalues and
% eigenvectors at the Hopf bifurcation to give a good enough idea
y = zeros(length(x0),2*n_nodes+1);

y (:,n_nodes)= conj(eigenvec);
y (:,n_nodes+2)= eigenvec;

% rescaling: we want sum(sum(y.* K^2.*y)) = 1
K_mat = repmat(-n_nodes:n_nodes,size(y,1),1);
rescaling = sum(sum(K_mat.^2.*conj(y).*y));
y = 1/sqrt(abs(rescaling))*y;

% constructing the full solution of the blowed up system
sol = Xi_vector([1/T_star, lambda0, a_star, horiz(x0)],y);

% adding 2 scalar condition
% first scalar condition: integral condition
big_Hopf.scalar_equations =fancy_scalar_condition(sol,big_Hopf.scalar_equations,1);
big_Hopf.scalar_equations = F_update_Hopf(big_Hopf.scalar_equations,sol);

% what I care for is for the amplitude to increase, since at the moment
% it's negative
DF0 = derivative_to_matrix(derivative(big_Hopf,sol,0));
x_dot_0 = kernel(DF0);

[~,index] = max(abs(real(x_dot_0(1:sol.size_scalar))));
angle = atan( imag(x_dot_0(index))/real(x_dot_0(index)));
x_dot_0 = exp( - 1i * angle) * x_dot_0; % bringing x_dot_0 to be symmetric (by multiplication with the appropriate complex rotation)
x_dot_0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_0,sol)));

big_H_square = big_Hopf;
big_H_square.scalar_equations = big_H_square.scalar_equations.change_lin_coef_vector([x_dot_0.',0]);
sol_2 = Newton_2(sol,big_H_square);

if max(norm(sol-sol_2))>step_size
    % update scalar equations
    big_Hopf.scalar_equations =fancy_scalar_condition(sol_2,big_Hopf.scalar_equations,1);
    
    % second sclar condition: some constraint on the norm of z (the rescaled
    % periodic orbit)
    p1 = Xi_vec2vec(sol_2).';
    p1(1:sol_2.size_scalar) = 0;
    p1 = conj(p1);
    const_p1 = -sum(p1.*conj(p1));
    big_Hopf.scalar_equations = change_lin_coef_vector(big_Hopf.scalar_equations,[p1.';const_p1],2);
    
end

if x_dot_0(3)*sol_2.scalar(3)>0 % I want the ampitude to go through 0
    x_dot_0 = -x_dot_0;
end

end