function F_new = continuation_equation_simplex(F_new, x, F_old)

DF_x = derivative_to_matrix(derivative(F_new,x,0));
[Q,~] = qr(DF_x.');
U_old = Q(:,end-1:end);     % Orthonormalize???
U = good_rotation(U_old,x);

new_eq = F_new.scalar_equations.number_equations_lin + [1:2];

if nargin>2
    old_coefs1 = extract_lin_coef(F_old.scalar_equations,new_eq(1));
    old_coefs2 = extract_lin_coef(F_old.scalar_equations,new_eq(2));
    old_coefs = [old_coefs1,old_coefs2].';
    
    rotated_coefs = rotate_base_tangentspace(U,old_coefs);
    U = rotated_coefs;
end

lin_coef_vec1 = [ U(:,1).', - U(:,1).' * Xi_vec2vec(x)];

F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec1, new_eq(1));

lin_coef_vec2 = [ U(:,2).', - U(:,2).' * Xi_vec2vec(x)];

F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec2, new_eq(2));

end


function rotated_coefs = good_rotation_unused(coefs, x0)
rotated_coefs = 0*coefs;
for i = 1:size(coefs,2)
    x_vec = coefs(:,i);
    [~,index] = max(abs(real(x_vec(1:x0.size_scalar))));
    angle = atan( imag(x_vec(index))/real(x_vec(index)));
    x_vec = exp( - 1i * angle) * x_vec;
    % bringing x_Xi_vec to be symmetric (by multiplication with the appropriate complex rotation)
    rotated_x_vec = Xi_vec2vec(symmetrise(vec2Xi_vec(x_vec,x0)));
    if norm(rotated_x_vec - x_vec)>10^-5
        warning('The rotation might not have been sufficient')
    end
    rotated_coefs(:,i) = rotated_x_vec;
end

end


function rotated_coefs = good_rotation(coefs, x0)

first_rotation = ones(1,8);
x = coefs(1:x0.size_scalar,1);
y = coefs(1:x0.size_scalar,2);
conj_prod = @(x,y) dot(x,y);

A_first = [imag(x), real(x), imag(y), real(y)];
f_first = @(c) [A_first*c;
    dot(c(1)+1i*c(2),c(1)+1i*c(2)) + dot(c(3)+1i*c(4),c(3)+1i*c(4)) - 1];
df_first = @(c) [A_first
    2*c(1:4).'];

first_rescaling = newton(first_rotation(1:4).', f_first, df_first);

coefs_rescaled = 0*coefs;
coefs_rescaled(:,1) = (first_rescaling(1)+1i*first_rescaling(2)) * coefs(:,1) + ...
     (first_rescaling(3)+1i*first_rescaling(4)) * coefs(:,2);
if abs(first_rescaling(1)+1i*first_rescaling(2)) > (first_rescaling(3)+1i*first_rescaling(4))
    coefs_rescaled(:,2) = coefs(:,2);
else
    coefs_rescaled(:,2) = coefs(:,1);
end

% Gram-Schmidt reduction on coefs_rescaled
proj = @(u,v) dot(u,v)/dot(u,u) * u;

coefs_rescaled(:,1) = coefs_rescaled(:,1) / norm(coefs_rescaled(:,1));

coefs_rescaled(:,2) = coefs_rescaled(:,2) - proj( coefs_rescaled(:,1), coefs_rescaled(:,2));
coefs_rescaled(:,2) = coefs_rescaled(:,2) / norm(coefs_rescaled(:,2));

rotated_coefs = coefs_rescaled;
for i = 2:size(coefs,2)
    x_vec = coefs_rescaled(:,i);
    [~,index] = max(abs(real(x_vec(1:x0.size_scalar))));
    angle = atan( imag(x_vec(index))/real(x_vec(index)));
    x_vec = exp( - 1i * angle) * x_vec;
    % bringing x_Xi_vec to be symmetric (by multiplication with the appropriate complex rotation)
    rotated_x_vec = Xi_vec2vec(symmetrise(vec2Xi_vec(x_vec,x0)));
    if norm(rotated_x_vec - x_vec)>10^-5
        warning('The rotation might not have been sufficient')
    end
    rotated_coefs(:,i) = rotated_x_vec;
end


% 
% 
% 
% 
% A = [imag(x), real(x), imag(y), real(y), zeros(x0.size_scalar,4);
%     zeros(x0.size_scalar,4),imag(x), real(x), imag(y), real(y)];
% 
% 
% f = @(c) [A*c;
%     dot(c(1)+1i*c(2),c(1)+1i*c(2)) + dot(c(3)+1i*c(4),c(3)+1i*c(4)) - 1
%     dot(c(5)+1i*c(6),c(5)+1i*c(6)) + dot(c(7)+1i*c(8),c(7)+1i*c(8)) - 1
%     ((c(1)+1i*c(2))*x + (c(3)+1i*c(4))*y)'* ( (c(5)+1i*c(6))*x + (c(7)+1i*c(8))*y)];
% 
% der_perpend = @(c) [ conj(x).'*((c(5)+1i*c(6))*x + (c(7)+1i*c(8))*y) 
%     -conj(x).'*((c(5)+1i*c(6))*x + (c(7)+1i*c(8))*y) 
%     conj(y).'* ((c(5)+1i*c(6))*x + (c(7)+1i*c(8))*y)
%     -conj(y).'*((c(5)+1i*c(6))*x + (c(7)+1i*c(8))*y)
%     (((c(1)-1i*c(2))*conj(x) + (c(3)-1i*c(4))*conj(y)).' * x)
%     (((c(1)-1i*c(2))*conj(x) + (c(3)-1i*c(4))*conj(y)).' * x)
%     (((c(1)-1i*c(2))*conj(x) + (c(3)-1i*c(4))*conj(y)).' * y)
%     (((c(1)-1i*c(2))*conj(x) + (c(3)-1i*c(4))*conj(y)).' * y)].';
% 
% df = @(c) [A
%     2*c(1:4).', 0*c(5:8).'
%     0*c(1:4).', 2*c(5:8).'
%     der_perpend(c)];
% 
% rescaling_coefs = newton(first_rotation.', f, df);
% 
% coefs_rescaled = 0*coefs;
% coefs_rescaled(:,1) = (rescaling_coefs(1)+1i*rescaling_coefs(2)) * coefs(:,1) + ...
%      (rescaling_coefs(3)+1i*rescaling_coefs(4)) * coefs(:,2);
% coefs_rescaled(:,2) = (rescaling_coefs(5)+1i*rescaling_coefs(6)) * coefs(:,1) + ...
%      (rescaling_coefs(7)+1i*rescaling_coefs(8)) * coefs(:,2);

end

function x = newton(x0, f, df)
x = x0;
for i = 1:100
    x = x - df(x)\f(x);
    if norm(f(x))<10^-10
        break
    end
end
end

function rotated_coefs = rotate_base_tangentspace(coefs,old_coefs)

% 2D tangent space can be defined by a rotated base
% be careful, it can also be flipped!
flipped_coefs = [coefs(:,2),coefs(:,1)];

apply_theta_rotation = @(theta) [sin(theta) cos(theta); cos(theta) -sin(theta)] * coefs.';
apply_theta_rotation_flipped = @(theta) [sin(theta) cos(theta); cos(theta) -sin(theta)] * flipped_coefs.';
theta_best = 0;
norm_best = norm(old_coefs - apply_theta_rotation(theta_best));
flipped = 0;
for theta = 0:0.2:2*pi
    
    norm_rotated = norm(old_coefs - apply_theta_rotation(theta));
    if norm_rotated < norm_best
        flipped = 0;
        theta_best = theta;
        norm_best = norm_rotated;
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
