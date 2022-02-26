% IS_POLYNOMIAL_COEF
function bool = is_polynomial_coef(alpha)
% function bool = is_polynomial_coef(alpha)
%
% INPUT
% alpha         input to be tested
% OUTPUT
% bool          1 if alpha is of the class polynomila_coef
%               0 otherwise

bool =0;
try
    alpha.size_scalar;
    alpha.size_vector;
    alpha.deg_scalar;
    alpha.deg_vector;
    alpha.n_equations;
    alpha.n_terms;
    alpha.value; % cell{n_equations}(n_terms) -- coefficient of the term
    alpha.dot; % cell{n_equations}{n_terms}(variables,:) -- number of derivatives applied to the given variable
    alpha.power_vector; % cell{n_equations}{n_terms}(variables,:)
    alpha.delay; % cell{n_equations}{n_terms}(variables,:)
    alpha.power_scalar; % cell{n_equations}(scalar_variables,n_terms)
    alpha.monomial_per_term;
catch
    return
end
bool = 1;
end
% end IS_POLYNOMIAL_COEF