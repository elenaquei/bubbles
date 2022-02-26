function x  = Newton(x0, f, df, norm_nu_loc, tries_max)
global display
if nargin<5
    tries_max = 100;
    if nargin<4
        norm_nu_loc = @(x) vector_norm(x);
    end
end
x = x0;
tries = 0;
while max(norm_nu_loc(f(x)))>10^-6 && tries<tries_max
    DF = df(x(:).');
    if display>1
        fprintf('Cond %e,    Norm %e\n',rcond(DF), norm_nu_loc(f(x)))
    end
    F = f(x);
    x = x(:) - DF\F(:);
    tries = tries +1;
end
if max(norm_nu_loc(f(x)))>10^-6
    error('Newton did not converge')
end
if any(isnan(x))
    error('Result contains NaN')
end
end