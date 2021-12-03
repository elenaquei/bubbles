
function [x, par_alpha, eigenvec, eigenval] = Hopf_numerical(x, par_alpha, parameter, eigenvec, eigenval, rhs, phi)
% function [x, par_alpha, eigenvec, eigenval] = Hopf_numerical(x, par_alpha, parameter)
%
% INPUT
% x     vector of 3 elements with an approximation of the position of a
% Hopf bifurcation
% par_alpha     double storing an approximation of hte parameter of the
% Hopf bifurcation
% parameter     double, fixed parameter in the search of a Hopf bifurcation
%
% OUTPUT
% x         vector of 3 elements, numerical coordinates of the Hopf bifurcation
% alpha     double, numerical parameter of the Hopf bifurcation

fn = @(x,alpha) rhs(x, alpha, parameter);

% (1i*beta, phi1+1i*phi2) approximate eigenpair
Df=finite_diff_fn(fn,x,par_alpha);
[V,D] = eigs(Df);
if nargin<5
    eigenvalues = diag(D);
    [~,beta_index] = min(abs(real(eigenvalues)));
    beta_eig = imag(eigenvalues(beta_index));
    phi1 = real(V(:,beta_index));
    phi2 = imag(V(:,beta_index));
else
    beta_eig = imag(eigenval);
    phi1 = real(eigenvec);
    phi2 = imag(eigenvec);
end
X_loc = [par_alpha; beta_eig; x; phi1; phi2];

X_loc =newton_Hopf(fn,X_loc, phi);

par_alpha = X_loc(1);
eigenval = 1i*X_loc(2);
x = X_loc(2+(1:length(x)));
eigenvec = X_loc(2+length(x)+(1:length(x))) + 1i*X_loc(2+2*length(x)+(1:length(x)));
end