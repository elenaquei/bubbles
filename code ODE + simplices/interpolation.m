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
alpha.scalar_eqs = interpolation_scal(alpha0.scalar_eqs, alpha1.scalar_eqs, alpha2.scalar_eqs);
end

function alpha = interpolation_scal(alpha0, alpha1, alpha2)
alpha = alpha0;
for i =1:3
    alpha.lin_coef = interpolation_float(alpha0.lin_coef, alpha1.lin_coef, alpha2.lin_coef);
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