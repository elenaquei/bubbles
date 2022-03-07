function [y, y_dot, F_cont] = reverse_hopf_gluing(s_second, s, f,rhs,bool_fancy_scalar,step_size_loc)
% function [y, y_dot, F_cont] = reverse_hopf_gluing(s_second, s, f,rhs, bool_fancy_scalar,step_size_loc)
%
% INPUT
% s             string, result of a Hopf continuation approach
% f             polynomial vector field, before Hopf transformation
% bool_fancy_scalar boolean, the phase condition is of the form <iKx,x>=0,
%               DEFAULT = 1
% OUTPUT
% y             Xi_vector, solution to f, such that the branch from x_n in  
%               s to y is validated under the appropriate change of
%               variables
% y_dot         vector of length compatible with y, tangent to the solution
%               branch
% F_new         full_problem, start of continuation
%
% ASSUMPTION: fixed point basically reached

if nargin <5 || isempty(bool_fancy_scalar)
    bool_fancy_scalar = 1;
end
global use_intlab
global talkative
temp_intval = use_intlab;
use_intlab = 0 ;

step_size = zeros(3); % will be overwritten by the load(s)

load(s_second)

if nargin<6 || isempty(step_size_loc)
    step_size = step_size(end)/100;
else
    step_size = step_size_loc;
end

x = x_n;

% deconstructing the blown-up solution into parts
n_vector = x.size_vector;
n_scalar = x.size_scalar;

if n_scalar~=2
    error('We are not ready for this')
end


% adding the standard continuation equation in F(x_0)
[F_0,x_dot_n] = continuation_eq(F,x_n,x_dot_n,step_size);

 
% the new solution and new equation
min_res = 10^-10;
bool_Hopf = 0;
[x_n1,x_dot_n1,F_1] =update_new(x_n,min(step_size,10^-3),x_dot_n,...
    F,[],min_res, bool_Hopf, bool_fancy_scalar);

[flag_a,~,Imax_start]=radii_polynomials(x_n1,F_1); 
[flag_b] = radii_polynomials_cont_new(x_n,x_n1,[],[],...
        F_0,F_1);
if flag_b*flag_a < 1
    error('Failed reverse Hopf gluing before starting') % works
end

fprintf('First step in the reverse Hopf gluing\n')

% New OPTION
new_cont_e_Xi_Vec = vec2Xi_vec(conj(x_dot_n),x_n);

new_cont_e_Xi_Vec.vector(:,:) = 0; 
%new_cont_e_Xi_Vec.scalar(n_scalar+1:end) = 0; 

new_cont_e = Xi_vec2vec(new_cont_e_Xi_Vec);

% general continuation equation in F
const = -Xi_vec2vec(x_n1).'*new_cont_e;


F_1.scalar_equations = change_lin_coef_vector(F_1.scalar_equations,...
    [new_cont_e;const],F_1.scalar_equations.number_equations_lin);
% the continuation equation is always the last one

%x_n1 = Newton_2(x_n1, F_1);

[flag]=radii_polynomials(x_n1,F_1); % works
if flag < 1
    error('Failed Hopf gluing before starting') % works 
end


% SECOND SECTION

% constructing the blown-up solution
center_approx = x_n1.vector(:,x_n1.nodes+1);

% better center
rhs_new = @(x)rhs(x(1),x(2),x(3),x(4),x_n1.scalar(2));
center = Newton_handle(rhs_new, center_approx);

periodic_orbit = x_n1.vector;
periodic_orbit(:,x_n1.nodes+1)=periodic_orbit(:,x_n1.nodes+1) - center;

% definition of amplitude
K_mat = repmat(-x_n1.nodes:x_n1.nodes,size(periodic_orbit,1),1);
amplitude_square = sum(sum(K_mat.^2.*conj(periodic_orbit).*periodic_orbit));
amplitude = sign(amplitude_square)*sqrt(abs(amplitude_square));
periodic_orbit = 1/amplitude*periodic_orbit;

new_scalar = [x_n1.scalar, amplitude, center']; % parameter, amplitude, fixed point
% blown-up solution
y = Xi_vector(new_scalar,periodic_orbit);

% passing to the standard continuation
% we have to add all the scalar equations.

%scal_eq = fancy_scalar_condition(y); 
% phase condition - corresponds to before the transform

big_Hopf = Taylor_series_Hopf(f,y.nodes);
big_Hopf = F_update_Hopf(big_Hopf,y);
% big_Hopf has 2 scalar equations and 4 polynomial equations


% calling continuation with the Hopf boolean (guarantees that we keep using
% the same scalar conditions over the continuation)
bool_Hopf = 1;
bool_saddle = 1; % also validate the existence of a saddle node
bool_fancy_scalar = 0; % already included in the Hopf boolean


F_cont =big_Hopf;

% adding parameter continuation to the F_0
F_cont0 = F_cont;

% retrieval of the coefficients of the linear part
coef_cont_eq = extract_lin_coef(F_1.scalar_equations, F_1.scalar_equations.number_equations_lin);
small_coef_cont_eq = zeros( y.size_scalar+y.size_vector*(2*y.nodes+1)+1,1);
small_coef_cont_eq(1:n_scalar) = coef_cont_eq(1:n_scalar);
small_coef_cont_eq(end) = const;
if any (coef_cont_eq(n_scalar+1:end))
    error('not compatible')
end

% add continuation equation in F
F_cont0.scalar_equations = change_lin_coef_vector(F_cont0.scalar_equations,...
    small_coef_cont_eq,F_cont0.scalar_equations.number_equations_lin+1);

 
% impose direction in the curve
%y_dot_0_approx = [x_dot_n1(1:n_scalar);x_dot_n1(x_n1.size_scalar+1:end)];
y_dot_0 = kernel(derivative_to_matrix(derivative(F_cont, y,0)));
[~,index] = max(abs(real(y_dot_0(1:y.size_scalar))));
angle = atan( imag(y_dot_0(index))/real(y_dot_0(index)));
y_dot_0 = exp( - 1i * angle) * y_dot_0;
y_dot_0 = Xi_vec2vec(symmetrise(vec2Xi_vec(y_dot_0,y)));
if y_dot_0(2)*x_dot_n1(2)<0
    y_dot_0=-y_dot_0;
end

% add continuation equation in F_cont1
[y_1,y_dot,F_cont1] =update_new(y,min(step_size,10^-3),y_dot_0,...
    F_cont,[],[], 0, bool_fancy_scalar);

use_intlab = temp_intval;
% some simple checks
[flag,Imin_end]=radii_polynomials(y,F_cont0); 
if flag < 1
   error('Failed Hopf gluing before starting, second step') 
end

[flag]=radii_polynomials(y_1,F_cont1); 
if flag < 1
   error('Failed Hopf gluing before starting, second step') 
end

if talkative
    disp('All checks for gluing completed: now validating the gluing')
end


% validating the segment
use_intlab = temp_intval ;
%flag = 1;%
[flag] = substepping(x_n,x_n1,F_0,F_1); 
if ~flag
    error('First step of reverse Hopf gluing failed')
else 
    fprintf('--------- First step of reverse Hopf gluing successful ----------\n')
end


fprintf('Second step in the Hopf gluing\n')
[flag] = substepping(y,y_1,F_cont0,F_cont1);
if ~flag
    error('Second step of Hopf gluing failed')
else 
    fprintf('--------- Second step of reverse Hopf gluing successful ----------\n')
end


norm_x = norm(y);
a1 = norm_x(3);
v_max = max(norm_x(end-x_n.size_vector:end));

% following equation 6.7 in Rigorous verification of Hopf bifurcations
% via desingularization and continuation
test_Imin = (1+a1+v_max+Imin_end)*Imin_end; 
if test_Imin > Imax_start
    error('Balls are not nested')
end


y = y_1;



% if Imin_cont > Imax || Imin > Imax_cont
%     warning('Are we sure it makes sense?')
% end

end

function f_period = add_period(f)
% function f_period = add_period(f)
%
% multiply all terms except the derivative by the inverse of the period
% (the first scalar unknown)
f_period = f;
f_period.size_scalar = f.size_scalar +1;
for i = 1:f.n_equations
    f_period.power_scalar {i} = zeros(size(f.power_scalar {i})+[1,0]);
    vector_index = 1:f.n_terms(i);
    for j =vector_index
        if any(f.dot{i}{j})
            index_dot =j;
        end
    end
    vector_index(index_dot) = [];
    f_period.power_scalar {i}(2:end,:) = f.power_scalar {i};
    f_period.power_scalar {i}(1,vector_index) = 1;
end

end



