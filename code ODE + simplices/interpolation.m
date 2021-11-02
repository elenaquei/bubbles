function alpha = interpolation(alpha0, alpha1, alpha2)

if isa(alpha0, 'float') || isa(alpha0, 'intval')
    alpha = interpolation_float(alpha0, alpha1, alpha2);
elseif isa(alpha0, 'Xi_vector')
    alpha = interpolation_Xivec(alpha0, alpha1, alpha2);
elseif isa(alpha0, 'scalar_eqs')
    alpha = interpolation_scal(alpha0, alpha1, alpha2);
elseif isa(alpha0, 'full_problem')
    alpha = interpolation_full_problem(alpha0, alpha1, alpha2);
else
    error('Not valid input')
end

end

function alpha = interpolation_full_problem(alpha0, alpha1, alpha2)
alpha = alpha0;
alpha.scalar_equations = interpolation_scal(alpha0.scalar_equations, alpha1.scalar_equations, alpha2.scalar_equations);
end

function alpha = interpolation_scal(alpha0, alpha1, alpha2)
alpha = alpha0;
for i =1:3
    alpha.linear_coef{i} = interpolation_float(alpha0.linear_coef{i},...
        alpha1.linear_coef{i}, alpha2.linear_coef{i});
end
end

function alpha = interpolation_Xivec(alpha0, alpha1, alpha2)
vec = interpolation_float( Xi_vec2vec(alpha0), Xi_vec2vec(alpha1), Xi_vec2vec(alpha2));
alpha= vec2Xi_vec(vec, alpha0);
end

function alpha = interpolation_float(a, b, c)
alpha = three_intval_real(real(a), real(b), real(c)) + 1i* three_intval_real(imag(a), imag(b), imag(c));
end

function val_intval = three_intval_real(a,b,c)
val_min = min(a, min(b,c));
val_max = max(a, max(b,c));
if isintval(val_min)
    val_min = inf(val_min);
end
if isintval(val_max)
    val_max = sup(val_max);
end
val_intval = infsup(val_min, val_max);
end