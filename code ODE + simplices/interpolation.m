function alpha = interpolation(alpha0, alpha1, varargin)

if nargin > 2
    alpha = interpolation( interpolation(alpha0, alpha1), varargin{:});
    return
end

if isa(alpha0, 'float') || isa(alpha0, 'intval')
    alpha = interpolation_float(alpha0, alpha1);
elseif isa(alpha0, 'Xi_vector')
    alpha = interpolation_Xivec(alpha0, alpha1);
elseif isa(alpha0, 'scalar_eqs')
    alpha = interpolation_scal(alpha0, alpha1);
elseif isa(alpha0, 'full_problem')
    alpha = interpolation_full_problem(alpha0, alpha1);
elseif isa(alpha0, 'cell')
    alpha = cell(size(alpha0));
    for i = 1:size(alpha0,1)
        for j = 1: size(alpha0,2)
            alpha{i,j} = interpolation(alpha0{i,j}, alpha1{i,j});
        end
    end
else
    error('Not valid input')
end

end

function alpha = interpolation_full_problem(alpha0, alpha1)
alpha = alpha0;
alpha.scalar_equations = interpolation_scal(alpha0.scalar_equations, alpha1.scalar_equations);
end

function alpha = interpolation_scal(alpha0, alpha1)
alpha = alpha0;
for i =1:3
    alpha.linear_coef{i} = interpolation_float(alpha0.linear_coef{i},...
        alpha1.linear_coef{i});
end
end

function alpha = interpolation_Xivec(alpha0, alpha1)
vec = interpolation_float( Xi_vec2vec(alpha0), Xi_vec2vec(alpha1));
alpha= vec2Xi_vec(vec, alpha0);
end

function alpha = interpolation_float(a, b)
if all(imag(a)==0) && all(imag(b)==0)
    alpha = three_intval_real(real(a), real(b));
else
    alpha = three_intval_real(real(a), real(b)) + 1i* three_intval_real(imag(a), imag(b));
end
end

function val_intval = three_intval_real(a,b)
val_min = min(a, b);
val_max = max(a, b);
if isintval(val_min)
    val_min = inf(val_min);
end
if isintval(val_max)
    val_max = sup(val_max);
end
val_intval = infsup(val_min, val_max);
end