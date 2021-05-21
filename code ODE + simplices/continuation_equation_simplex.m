function F_new = continuation_equation_simplex(F_new, x)

DF_x = derivative_to_matrix(derivative(F_new,x,0));
[Q,~] = qr(DF_x.');
U = Q(:,end-1:end);     % Orthonormalize???

new_eq = F_new.scalar_equations.number_equations_lin+[1:2];

lin_coef_vec1 = [ U(:,1).', - U(:,1).' * Xi_vec2vec(x)];

F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec1, new_eq(1)); 

lin_coef_vec2 = [ U(:,2).', - U(:,2).' * Xi_vec2vec(x)];
F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec2, new_eq(2)); 
end
