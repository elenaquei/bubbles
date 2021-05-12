function big_Hopf = Taylor_series_Hopf(small_Hopf,nodes)
% function big_Hopf = Taylor_series_Hopf(small_Hopf,nodes)
% 
% construction of the vector field of the periodic Hopf problem
%
% INPUT 
% small_Hopf   instance of polynomial_coefs, stands for f(x, lambda)
% nodes        positive integer, number of nodes of big_Hopf
% OUTPUT
% big_Hopf     instance of full_problem, stands for f(x, lambda) 
%                 y_dot - sum d^k f(x,l) y^k a^(|k|-1) T
% 
if nargin<2
    nodes = 2;
end


n_equations_lin = 1;
small_Hopf = set_zero_derivative(small_Hopf); % delete the derivatives, they are not necessaries

n_equations_vec = small_Hopf.n_equations;%+1;
n_scalar = small_Hopf.size_scalar + small_Hopf.size_vector + 2;
n_vector = small_Hopf.size_vector;

% linear scalar eqautions
lin_coefficients = cell(3,1);
lin_coefficients{1}= zeros(n_equations_lin,n_scalar);
lin_coefficients{2}=zeros(n_equations_lin,n_vector,2*nodes+1);
lin_coefficients{3}=zeros(n_equations_lin,1);

lin_coefficients{1}(1,4)= 0; % x1
lin_coefficients{2}(1,1,:) = 0; % y1(0)
% y1(0) = x1


value =cell(n_equations_vec,1);
power_scalar = cell(n_equations_vec,1);
power_vector =cell(n_equations_vec,1);
%monomial_per_term = cell(n_equations_vec,1);
%n_terms = ones(n_equations_vec,1);
n_terms = small_Hopf.n_terms;
% n_terms(1:end-1) = small_Hopf.n_terms -1;
% value{n_equations_vec} = 1;
% power_scalar{n_equations_vec} = [0,zeros(1,small_Hopf.size_scalar),1,zeros(1,n_vector)].';
% power_vector{n_equations_vec} = cell(1);
% power_vector{n_equations_vec}{1} = zeros(n_vector,1);

for i = 1:n_equations_vec %-1
    value{i} = zeros(n_terms(i),1);
    if any(small_Hopf.monomial_per_term{i} ~=1)
        error('Not a fitting Hopf problem')
    end
    %monomial_per_term{i} = ones(n_terms(i),1);
    j_tilde = 1;
    for j = 1:small_Hopf.n_terms(i)
        if sum(sum(abs(small_Hopf.dot{i}{j}))) ==1
            continue;
        elseif sum(sum(abs(small_Hopf.dot{i}{j}))) ==0
            % do
            value{i}(j_tilde) = small_Hopf.value{i}(j); 
            power_vector{i}{j_tilde} = zeros(n_vector,1);
            power_scalar{i}(:,j_tilde) = [0; small_Hopf.power_scalar{i}(:,j);0;small_Hopf.power_vector{i}{j}];
            j_tilde = j_tilde+1;
        else
            error('Not a Hopf problem')
        end
    end
end

scalar_Hopf = polynomial_coefs(n_scalar, n_vector, n_equations_vec, ...
                n_terms, value,power_scalar,power_vector,power_vector,power_vector);
eq_G = scalar_eq(n_equations_lin, n_equations_vec, n_scalar, n_vector, lin_coefficients, scalar_Hopf);


% define number_term, coefficients,power_Scal,power_Vec

number_term= zeros(n_vector,1);
power_Scal = cell(n_vector,1);
power_Vec= cell(n_vector,1);
coefficients = cell(n_vector,1);
derivative =  cell(n_vector,1);
for j = 1 : n_vector % equation we are considering
    counter = 0;
    Non0 = 0;
    for i = 1:small_Hopf.n_terms(j)
        if size(small_Hopf.power_vector{j}{i},2)~=1
            error('Not compatible with Hopf');
        end
        Non0 = Non0 + prod(small_Hopf.power_vector{j}{i}+1)-1;
    end
    Non0 = 1 + Non0;
    number_term(j) = Non0;
    power_Scal{j} = zeros(n_scalar, Non0);
    power_Vec{j} = cell(Non0,1);
    coefficients{j} = zeros(Non0, 1);
    derivative{j} = cell(Non0,1);
    
    for m=1:small_Hopf.n_terms(j)
        
        p = small_Hopf.power_vector{j}{m};
        
        for i = 1: prod(p+1)-1
            counter = counter+1;
            k = convert_from_N_to_vec(i,p); % works right
            
            c = -small_Hopf.value{j}(m);
            c_new= c * prod(factorial(p)./(factorial(p-k)))/prod(factorial((k)));% .*factorial(p-k))); % k over p
            %c*prod(prod(construct_product_matrix(n_vector,k,p)));
            %%% THERE IS A PROBLEM HERE 
            
            q_new = [1, small_Hopf.power_scalar{j}(m), sum(k)-1, horiz(p-k)]; % T, lambda,a, x
            % if all (k==0) && i == prod(p+1)
            %    q_new = [1, small_Hopf.power_scalar{j}(m), 0, horiz(p-k)];
            % end
            p_new = horiz(k).';
            
            if any(q_new<0)
                error('What the hell')
            elseif any(p_new<0)
                error('Here too')
            end
            
            power_Scal{j}(:,counter) = q_new;
            power_Vec{j}{counter}= p_new;
            derivative{j}{counter}= 0*p_new;
            coefficients{j}(counter) = c_new;
        end
        
    end
    if counter ~= Non0-1
        error('here')
    end
    % here to add y_dot_j
    counter = counter +1;
    power_Scal{j}(:,counter) = 0*q_new;
    power_Vec{j}{counter}= 0*p_new;
    derivative{j}{counter}= 0*p_new;
    derivative{j}{counter}(j)= 1;
    coefficients{j}(counter) = 1;
end


eq_F = polynomial_coefs(n_scalar, n_vector, n_vector, ...
                number_term, coefficients,power_Scal,power_Vec,derivative);

big_Hopf= full_problem(eq_G, eq_F);

end

function i = convert_from_N_to_vec ( i_iter , deg_f)
% function i = convert_from_N_to_vec ( i_iter , deg_f)

if i_iter > prod(deg_f+1)
    warning ('i_iter too big');
end
if i_iter <=0
    error('index must be strictly positive');
end

f = @(k,n,d)mod( floor( k / prod( d(1:n-1))),d(n));
d = deg_f +1;
i = 0*deg_f;
for n =1:length(deg_f)
    i(n) = f(i_iter,n,d);
end
end


% % help function 
% function K = construct_product_matrix(N, k,p)
% 
% K = ones(N, max(k));
% for i = 1:max(k)
%     index = find( (p-i+1).*(k-i)>0);
%     K(index,i) = p(index)-i+1;
% end
% end