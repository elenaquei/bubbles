function [x1,x_dot_1,F_new, iter_Newton] =update(x0,h,x_dot_0,...
    F_not_square,maxiter,min_res, bool_Hopf, bool_fancy_scalar)
%function [x1,x_dot1,F_new, iter_Newton] =update(xBar0,delta,x_dot0,...
%   alpha,coefs_linear,maxiter,min_res, bool_Hopf, bool_fancy_scalar)
%
% INPUT
% xBar0             Xi_vector, approximation of a solution
% delta             step size
% x_dot0            vector, compatible with xBar0, direction of the curve
% alpha             full_problem, not square
% maxiter           maximum number of iterations for Newton (DEFAULT 10)
% min_res           minimum residual for Newton (DEFALUT 10*norm(apply(alpha,xBar0)) )
% bool_Hopf         system is Hopf (DEFAULT 0)
% bool_fancy_scalar integral condition (DEFAULT 0)
% OUTPUT
% xBar1         next approximation
% x_dot1        next tangent to the curve
% F_new         full_problem, square, at xBar1
% iter_Newton   number of iterations in Newton

if  nargin<5 || isempty(maxiter)
    maxiter=10;
end
if  nargin<6 || isempty(min_res)
    min_res= max(10^-14,10 * max(norm(apply(F_not_square,x0))));
end
if nargin<7 || isempty(bool_Hopf)
    bool_Hopf =0;
end
if nargin<8 || isempty(bool_fancy_scalar)
    bool_fancy_scalar  =0;
end

if bool_Hopf
    if x0.size_scalar~=x0.size_vector+3
        error('If you call for a Hopf validation, the number of scalar variables should be the number of vector variables + 3')
    end
end

% Predictor step: compute new solution with arch-length parametrisation
x_tilde_1 = x0 + h * vec2Xi_vec(x_dot_0,x0);

% update the scalar equation that needs updating
F_new_tilde = F_update(F_not_square, x_tilde_1, x_dot_0, bool_Hopf, bool_fancy_scalar);


% Corrector step: get a better approximation
[x1,iter_Newton] = Newton_2(x_tilde_1,F_new_tilde,maxiter,min_res);
x1 = x1.symmetrise; % let it be real

% update the scalar equation that needs updating
F_new = F_not_square;
if bool_Hopf
    F_new.scalar_equations = F_update_Hopf(F_new.scalar_equations,x1);
end
if bool_fancy_scalar
    F_new.scalar_equations = fancy_scalar_condition(x1, F_new.scalar_equations);
end

% % compute kernel and set it real
% x_dot_1 = kernel(derivative_to_matrix(derivative(F_new, x1, 0)));
% [~,index] = max(abs(real(x_dot_1(1:x0.size_scalar))));
% angle = atan( imag(x_dot_1(index))/real(x_dot_1(index)));
% x_dot_1 = exp( - 1i * angle) * x_dot_1;
% x_dot_1 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_1,x0)));
% if dot(x_dot_1,x_dot_0)<0
%     x_dot_1=-x_dot_1;
% end
%const = -Xi_vec2vec(x1).'*x_dot_1;
% addition of continuation equation
% F_new = F_not_square;

%F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[x_dot_1;const],F_new.scalar_equations.number_equations_lin+1);


[F_new,x_dot_1] = continuation_eq(F_new,x1,x_dot_0,h);

% some checks on the step just made
% if max(abs(DF1*x_dot_1)) > 10^-6
%     warning('inversion not cool')
% end

if norm(norm(x0-x1)) < h/100
    warning('Distance between consecutive points quite small, %e\n',norm(norm(x0-x1)))
end

return

end


function F_new = F_update(F_old, x_new, x_dot, bool_Hopf, bool_fancy)
% function F_new = F_update(F_old, x_new, x_dot, bool_Hopf, bool_fancy)
%
% update the scalar equations in F_old
% in particular, it updates the continuation equation, and, when requested
% the second Hopf equation and the first fancy integral phase condition
%
% INPUTS
% F_old         full_problem, initial problem
% x_new         Xi_vector, approximate solution
% x_dot         vector, approximate tangent
% bool_Hopf     bool for Hopf system
% bool_fancy    bool for fancy integral scalar equation
% OUTPUTS
% F_new         full_problem with fitting scalar coefficients

F_new = F_old;
% update the scalar equation that needs updating
if bool_Hopf
    F_new.scalar_equations = F_update_Hopf(F_old.scalar_equations,x_new);
end
if bool_fancy
    F_new.scalar_equations = fancy_scalar_condition(x_new, F_new.scalar_equations);
end
const = -Xi_vec2vec(x_new).'*conj(x_dot);
F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[conj(x_dot);const],x_new.size_scalar);
end


% function F_new = F_update_Hopf(F_old,x_new)
% % function F_new = F_update_Hopf(F_old,x_new)
% F_new = F_old;
% 
% F_new.scalar_equations =fancy_scalar_condition(x_new,F_old.scalar_equations,1);
%     
%     
% % % p1 = Xi_vec2vec(x_new).';
% % % p1(1:x_new.size_scalar) = 0;
% % 
% % p1 = Xi_vec2vec(x_new).';
% % p1(1:x_new.size_scalar) = 0;
% % % K = x_new.nodes;
% % % K_vec = -K:K;
% % % K_mat = [0*(1:x_new.size_scalar),repmat(K_vec, 1,x_new.size_vector)];
% % % p1 = K_mat.^2.*conj(p1); % NEW CONDITION
% % % const_p1 = -1;%sum(p1.*conj(p1));% NEW CONDITION
% % p1 = conj(p1);
% % const_p1 = -sum(p1.*conj(p1));
% % F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[p1.';const_p1],2);
% 
% F_new = F_update_Hopf(F_new,x_new);
% 
% % 
% % F_new = F_old;
% % p1 = Xi_vec2vec(x_new)';
% % p1(1:5) = 0;
% % const_p1 = -1;
% % F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[p1.';const_p1],2);
% end