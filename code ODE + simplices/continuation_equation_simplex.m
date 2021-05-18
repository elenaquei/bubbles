function F_new = continuation_equation_simplex(F_new,x1)

DF_x = derivative_to_matrix(derivative(F_new,x1,0));
U = null(DF_x);

lin_coef_vec1 = [ U(:,1).', - dot(U(:,1).', Xi_vec2vec(x1))];
new_eq = F_new.scalar_equations.number_equations_lin+1;
F_new = change_lin_coef_vector(F_new.scalar_equations, lin_coef_vec1, new_eq); 


lin_coef_vec2 = [ U(:,2).', - dot(U(:,2).', Xi_vec2vec(x1))];
new_eq = F_new.scalar_equations.number_equations_lin+1;
F_new = change_lin_coef_vector(F_new.scalar_equations, lin_coef_vec2, new_eq); 
end
