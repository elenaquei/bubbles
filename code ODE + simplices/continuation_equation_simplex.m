function F_new = continuation_equation_simplex(F_new, x, new_eq)

old_coefs = extract_all_lin_coef(F_new.scalar_equations);

DF_x = derivative_to_matrix(derivative(F_new,x,0));
[Q,~] = qr(DF_x.');
U_old = Q(:,end-1:end);     % Orthonormalize???
U = good_rotation(U_old,x,old_coefs);

if nargin<3 || isempty(new_eq)
    new_eq = F_new.scalar_equations.number_equations_lin + [1:2];
end

lin_coef_vec1 = [ U(:,1).', - U(:,1).' * Xi_vec2vec(x)];

F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec1, new_eq(1)); 

lin_coef_vec2 = [ U(:,2).', - U(:,2).' * Xi_vec2vec(x)];
F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec2, new_eq(2)); 

if norm(lin_coef_vec1(1:end-1) - old_coefs(1,:))>10^-3 || ...
        norm(lin_coef_vec2(1:end-1) - old_coefs(2,:))>10^-3
    warning('Old scalar equation and new ones are far apart')
end

end


function rotated_coefs = good_rotation(coefs, x0, old_coefs)

for i = 1:2
    x_vec = coefs(:,i);
    [~,index] = max(abs(real(x_vec(1:x0.size_scalar))));
    angle = atan( imag(x_vec(index))/real(x_vec(index)));
    x_vec = exp( - 1i * angle) * x_vec;
    % bringing x_Xi_vec to be symmetric (by multiplication with the appropriate complex rotation)
    rotated_x_vec = Xi_vec2vec(symmetrise(vec2Xi_vec(x_vec,x0)));
    coefs(:,i) = rotated_x_vec;
end

apply_theta_rotation = @(theta) [sin(theta) cos(theta); cos(theta) sin(theta)] * coefs.';
theta_best = 0;
norm_best = norm(old_coefs - apply_theta_rotation(theta_best));
for theta = 0:0.2:2*pi
    norm_rotated = norm(old_coefs - apply_theta_rotation(theta));
    if norm_rotated < norm_best
        theta_best = theta;
        norm_best = norm_rotated;
    end
end

rotated_coefs = apply_theta_rotation(theta_best).';
end