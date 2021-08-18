function F_new = continuation_equation_simplex(F_new, x, F_old)

DF_x = derivative_to_matrix(derivative(F_new,x,0));
[Q,~] = qr(DF_x.');
U_old = Q(:,end-1:end);     % Orthonormalize???
U = good_rotation(U_old,x);

if nargin>2
    old_coefs1 = extract_lin_coef(F_old.scalar_equations,3);
    old_coefs2 = extract_lin_coef(F_old.scalar_equations,4);
    old_coefs = [old_coefs1,old_coefs2].';
    
    rotated_coefs = rotate_base_tangentspace(U,old_coefs);
    U = rotated_coefs;
end

new_eq = F_new.scalar_equations.number_equations_lin + [1:2];

lin_coef_vec1 = [ U(:,1).', - U(:,1).' * Xi_vec2vec(x)];

F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec1, new_eq(1));

lin_coef_vec2 = [ U(:,2).', - U(:,2).' * Xi_vec2vec(x)];
F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec2, new_eq(2));


end


function rotated_coefs = good_rotation(coefs, x0)
rotated_coefs = 0*coefs;
for i = 1:2
    x_vec = coefs(:,i);
    [~,index] = max(abs(real(x_vec(1:x0.size_scalar))));
    angle = atan( imag(x_vec(index))/real(x_vec(index)));
    x_vec = exp( - 1i * angle) * x_vec;
    % bringing x_Xi_vec to be symmetric (by multiplication with the appropriate complex rotation)
    rotated_x_vec = Xi_vec2vec(symmetrise(vec2Xi_vec(x_vec,x0)));
    rotated_coefs(:,i) = rotated_x_vec;
end
end

function rotated_coefs = rotate_base_tangentspace(coefs,old_coefs)

% 2D tangent space can be defined by a rotated base??
% FUCK, it can also be flipped!
flipped_coefs = [coefs(:,2),coefs(:,1)];

apply_theta_rotation = @(theta) [sin(theta) cos(theta); cos(theta) -sin(theta)] * coefs.';
apply_theta_rotation_flipped = @(theta) [sin(theta) cos(theta); cos(theta) -sin(theta)] * flipped_coefs.';
theta_best = 0;
norm_best = norm(old_coefs - apply_theta_rotation(theta_best));
for theta = 0:0.2:2*pi
    norm_rotated = norm(old_coefs - apply_theta_rotation(theta));
    if norm_rotated < norm_best
        theta_best = theta;
        norm_best = norm_rotated;
        flipped = 0;
    end
    
    norm_rotated_flipped = norm(old_coefs - apply_theta_rotation_flipped(theta));
    if norm_rotated_flipped < norm_best
        theta_best = theta;
        norm_best = norm_rotated_flipped;
        flipped = 1;
    end
end

if flipped
    rotated_coefs = apply_theta_rotation_flipped(theta_best).';
else
    rotated_coefs = apply_theta_rotation(theta_best).';
end


end