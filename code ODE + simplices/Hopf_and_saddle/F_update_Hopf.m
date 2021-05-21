function s_new = F_update_Hopf(s_old,x_new)
% function F_new = F_update_Hopf(F_old,x_new) 

if ~isa(s_old, 'scalar_eq')
    all_eqs = s_old;
    s_old = s_old.scalar_equations;
    bool_scal = 0;
else
    bool_scal = 1;
end


if isempty(s_old)
    s_new =fancy_scalar_condition(x_new);
else
    s_new =fancy_scalar_condition(x_new,s_old,1);
end



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
p2 = K_mat.^2.*conj(p2); % NEW CONDITION
const_p2 = -sum(Xi_vec2vec(x_new).'.*p2) ; % -1

% NEW CONDITION  
% const_p2 could also be -1, but unfit for Rossler, somehow

s_new = change_lin_coef_vector(s_new,[p2.';const_p2],2);

if bool_scal==0
    temp = s_new;
    s_new = all_eqs;
    s_new.scalar_equations = temp;
end