function F_new = F_update_Hopf(F_old,x_new)
% function F_new = F_update_Hopf(F_old,x_new) 

F_new = F_old;

F_new.scalar_equations =fancy_scalar_condition(x_new,F_old.scalar_equations,1);
    

p1 = Xi_vec2vec(x_new).';
p1(1:x_new.size_scalar) = 0;
p1 = conj(p1);
const_p1 = -sum(p1.*conj(p1));
%F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[p1.';const_p1],2);

p2 = Xi_vec2vec(x_new).';  
p2(1:x_new.size_scalar) = 0;
K = x_new.nodes;
K_vec = -K:K;
K_mat = [0*(1:x_new.size_scalar),repmat(K_vec, 1,x_new.size_vector)];
p2 = -1i*K_mat.*conj(p2); % NEW CONDITION
const_p2 = -sum(Xi_vec2vec(x_new).'.*p2) ;

% NEW CONDITION  
% const_p2 could also be -1, but unfit for Rossler, somehow

F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[p2.';const_p2],2);

end