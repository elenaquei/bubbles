function F_new = continuation_equation_simplex(F_new, x)

DF_x = derivative_to_matrix(derivative(F_new,x,0));
[Q,~] = qr(DF_x.');
U = Q(:,end-1:end);     % Orthonormalize???
U(:,1) = good_rotation(U(:,1),x);
U(:,2) = good_rotation(U(:,2),x);

new_eq = F_new.scalar_equations.number_equations_lin+[1:2];

lin_coef_vec1 = [ U(:,1).', - U(:,1).' * Xi_vec2vec(x)];

F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec1, new_eq(1)); 

lin_coef_vec2 = [ U(:,2).', - U(:,2).' * Xi_vec2vec(x)];
F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec2, new_eq(2)); 
end


function rotated_x_vec = good_rotation(x_vec, x0)

% [~,index] = max(abs(real(x_vec(1:x0.size_scalar))));
% angle = atan( imag(x_vec(index))/real(x_vec(index)));
% x_vec = exp( - 1i * angle) * x_vec;
% bringing x_Xi_vec to be symmetric (by multiplication with the appropriate complex rotation)
rotated_x_vec = Xi_vec2vec(symmetrise(vec2Xi_vec(x_vec,x0)));
end