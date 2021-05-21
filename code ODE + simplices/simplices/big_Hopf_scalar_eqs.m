function Hopf_scalar_eqs = big_Hopf_scalar_eqs(sol_Xi, numerical_Hopfs)
% function Hopf_scalar_eqs = big_Hopf_scalar_eqs(sol_Xi, numerical_Hopfs)
%
% numerical_Hopfs is a 2-cell, where each element of the cell is a data
% structure to store x, par_alpha, par_beta, eigenvec, eigenval
if sol_Xi.size_scalar ~= sol_Xi.size_vector + 4
    error('Did not extend for this case yet!')
end


% scalar equations will have the form
% <iKz, z> = 0 i.e. fancy scalar equation
% < K^2z, z > = 1 i.e. Hopf equation
% x = x_s i.e. Hopf center continuation equation
% <T,z> = c (2 equations), i.e. simplex continuation equation
% the simplex continuation equations ARE NOT INCLUDED in this code

Hopf_scalar_eqs = F_update_Hopf([],x_new);

vec_coefs = [0*Xi_vec2vec(sol_Xi), 0];
alpha = sol_Xi.scalar(2);
x1 = numerical_Hopfs{1}.x;
x2 = numerical_Hopfs{2}.x;
alpha1 = numerical_Hopfs{1}.alpha;
alpha2 = numerical_Hopfs{2}.alpha;
s = (alpha - alpha1) / (alpha2 - alpha1);
xs = x1 + s * ( x2 -x1);
for i = 1:sol_Xi.size_vector
    n_eq = 2+i;
    
    vec_coefs(4+i) = 1;
    vec_coefs(end) = -xs(i);
    
    Hopf_scalar_eqs = change_lin_coef_vector(Hopf_scalar_eqs, ...
    vec_coefs, n_eq); 
    
    vec_coefs = 0*vec_coefs;
end
