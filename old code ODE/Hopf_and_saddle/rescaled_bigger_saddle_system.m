function [saddle_problem0,saddle_problem1,DxHyy] = rescaled_bigger_saddle_system...
    (normal_system0,normal_system1,x_tilde_0, x_tilde_1, step_size)
% function [saddle_problem0,saddle_problem1,DxHyy] = rescaled_bigger_saddle_system
%   (normal_system0,normal_system1,x_tilde_0, x_tilde_1, step_size)
%
% INPUT
% normal_system0        full_problem at x0, square problem, continuation equation
%                       included
% normal_system1        full_problem at x1, square problem, continuation equation
%                       included
% x_tilde_0             numerical solution at the beginning 
% x_tilde_1             numerical solution at the end (the one used for the
%                       continuation equation)
% step_size             rescaling of the arc length parameter
%
% OUTPUT
% saddle_problem0       full_problem, square problem, extended system to
%                       validated the saddle bifurcation, at the beginning
% saddle_problem1       full_problem, square problem, extended system to
%                       validated the saddle bifurcation, at the end
% DxHyy                 full_problem, 
%

% two options: either requesting two systems at s=0 and s=1, or requesting
% v_0,v_1,x_0,x1

% at the moment, let us not care of E_s but just of F and G, that request
% just one system.
% REMARK: the handling of F and G is done in exactly the same way! No need
% to duplicate the function ;)

if nargin<4
    error('Too few inputs')
%elseif nargin==4
%    rescaling=100;
end
normal_G0 = normal_system0.scalar_equations;
normal_G1 = normal_system1.scalar_equations;
normal_F = normal_system0.vector_field;
% check the vector field is the same
if normal_system0.vector_field ~= normal_system1.vector_field
    error('The vector field is not suppose to change from x0 to x1');
end


[calligraphic_G0,calligraphic_G1] = saddle_scalar_eq...
    (normal_G0,normal_G1,x_tilde_0, x_tilde_1, step_size);  % scalar polynomial coefficients
% G includes E_s therefore is doubled
[~,DxGyy] = saddle_expansion(normal_G0.polynomial_equations);

[calligraphic_F, DxFyy] = saddle_expansion(normal_F); % vector polynomial coefficients

saddle_problem0 = full_problem(calligraphic_G0,calligraphic_F);
saddle_problem1 = full_problem(calligraphic_G1,calligraphic_F);

zero_cell = calligraphic_G0.linear_coef;
zero_cell{1} = 0*zero_cell{1}(1:normal_G0.number_equations_lin,:);
zero_cell{2} = 0*zero_cell{2}(1:normal_G0.number_equations_lin,:,:);
zero_cell{3} = 0*zero_cell{3}(1:normal_G0.number_equations_lin);

DxxGyy = scalar_eq(normal_G0.number_equations_lin, normal_G0.number_equations_pol, calligraphic_G0.size_scalar, calligraphic_G0.size_vector, zero_cell, DxGyy);
DxHyy = full_problem(DxxGyy, DxFyy);
end

function [calligraphic_G0,calligraphic_G1] = saddle_scalar_eq(normal_G0,normal_G1,~, ~,step_size)
% function [calligraphic_G0,calligraphic_G1] = saddle_scalar_eq(normal_G0,normal_G1,~, ~, step_size) 
%
% constructing the polynomial coefficients structure necessary for the
% validation of the saddle node
%
% INPUT 
% normal_G0,normal_G1      polynomial_coefs, with no dots nor delays, at
%                          the beginning and at the end of the segment
% ~, ~     Xi_vector, numerical solutions at the beginning
%                   and at the end
% step_size                 rescaling factor
% OUTPUT
% calligrafic_G     polynomial_coefs, with no dots nor delays, fitting for
%                   the validation of the saddle node


if normal_G0.polynomial_equations~=normal_G1.polynomial_equations
 error('The polynomial part of the scalar equations must agree')
else
    normal_G = normal_G0.polynomial_equations;
end

M = normal_G.size_scalar;
N = normal_G.size_vector;
new_size_scalar = 3*M;
new_size_vector = 3*N; 

% --- E_s section
coef_0 = normal_G0.linear_coef;
coef_1 = normal_G1.linear_coef;
coef_Delta = coef_0;
%coef_Delta{1} = coef_1{1} - coef_0{1};
%coef_Delta{2} = coef_1{2} - coef_0{2};
%coef_Delta{3} = coef_1{3} - coef_0{3};
coef_Delta{1} = (coef_1{1} - coef_0{1})/step_size;%with the REFinement, it works
coef_Delta{2} = (coef_1{2} - coef_0{2})/step_size;
coef_Delta{3} = (coef_1{3} - coef_0{3})/step_size;


new_coef_1_and_2 =@(vec,vec_Delta) cat(1,cat(2,vec, 0*vec, 0*vec),...
    cat(2,vec_Delta, vec, 0*vec),...
    cat(2,0*vec, 2*vec_Delta, vec));
new_coef_3 = @(vec,vec_Delta)[vert(vec);vert(vec_Delta); 0*vec];

calligraphic_E0 = cell(3,1);
calligraphic_E0{1} = new_coef_1_and_2(coef_0{1},coef_Delta{1});
calligraphic_E0{2} = new_coef_1_and_2(coef_0{2},coef_Delta{2});
calligraphic_E0{3} = new_coef_3(coef_0{3},coef_Delta{3});
calligraphic_E1 = cell(3,1);
calligraphic_E1{1} = new_coef_1_and_2(coef_1{1},coef_Delta{1});
calligraphic_E1{2} = new_coef_1_and_2(coef_1{2},coef_Delta{2});
calligraphic_E1{3} = new_coef_3(coef_1{3},coef_Delta{3});
% --- end E_s section

polynomials = saddle_expansion(normal_G1.polynomial_equations);

calligraphic_G0 = scalar_eq(length(calligraphic_E0{3}), polynomials.n_equations, ...
    new_size_scalar, new_size_vector, calligraphic_E0, polynomials);
calligraphic_G1 = scalar_eq(length(calligraphic_E0{3}), polynomials.n_equations, ...
    new_size_scalar, new_size_vector, calligraphic_E1, polynomials);

end

function [new_pol,DxHyy] = saddle_expansion(old_pol)
% change the initial polynomial to fit the extended system to validate a
% saddle node
% INPUT 
% old_pol       instance of polynomial_coef
% OUTPUT
% new_pol       instance of polynomial_coef
% DxHyy         also
%
% new_pol = [Dx(olp_pol)y ; Dxx(old_pol) yy + Dx(old_pol)z; old_pol]
% corresponding to the extended system referred to in the pdf

%global rescaling_saddle 

n_equation = old_pol.n_equations;
value_new = cell(3*n_equation,1);
scalar_pow_new = cell(3*n_equation,1);
vector_pow_new = cell(3*n_equation,1);
dot_new = cell(3*n_equation,1);
n_terms_new = zeros(3*n_equation,1);


value_yy = cell(n_equation,1);
scalar_pow_yy = cell(n_equation,1);
vector_pow_yy = cell(n_equation,1);
dot_yy = cell(n_equation,1);
n_terms_yy = zeros(n_equation,1);


for i=1:n_equation
    local_value = old_pol.value{i};       % in functions referred as c
    local_scalar_power = old_pol.power_scalar{i}.';      % in functions referred as d1
    local_vector_power = matrification(old_pol.power_vector{i});   %% in functions referred as d2
    local_derivative =matrification(old_pol.dot{i});   % in functions referred as d3
    
    % add variables y and z
    [c_big, d1_big, d2_big, d3_big] = padding(local_value,local_scalar_power,local_vector_power,local_derivative);
    
    % compute fitting coefficients for various derivatives
    [c_y, d1_y, d2_y, d3_y] = Dx_y(c_big, d1_big, d2_big, d3_big);
    [c_yy, d1_yy, d2_yy, d3_yy] = Dx_y(c_y, d1_y, d2_y, d3_y);
    [c_z, d1_z, d2_z, d3_z] = Dx_z(c_big, d1_big, d2_big, d3_big);
    [c, d1, d2, d3] = sum_coef(c_yy, d1_yy, d2_yy,d3_yy,c_z, d1_z, d2_z,d3_z);
%     if length(rescaling_saddle)==2
%         r1 = rescaling_saddle(1);
%         r2 = rescaling_saddle(2);
%         [c, d1, d2, d3] = sum_coef((r2/r1)*c_yy, d1_yy, d2_yy,d3_yy,c_z, d1_z, d2_z,d3_z);
%     end
    
    % delete zero terms
    [c_yy, d1_yy, d2_yy, d3_yy] = contract (c_yy, d1_yy, d2_yy, d3_yy);
    [c_yy_and_z, d1_yy_and_z, d2_yy_and_z, d3_yy_and_z] = contract(c,d1,d2,d3);
    [c_y, d1_y, d2_y, d3_y] = contract(c_y, d1_y, d2_y, d3_y);
    
    % managing of storage
    d2_yy_and_z = inverse_matrification(d2_yy_and_z);
    d3_yy_and_z = inverse_matrification(d3_yy_and_z);
    
    d2_y = inverse_matrification(d2_y);
    d3_y = inverse_matrification(d3_y);
    
    d2_yy = inverse_matrification(d2_yy);
    d3_yy = inverse_matrification(d3_yy);
    
    % store in fitting format
    value_new{i+0*n_equation} = c_big;
    scalar_pow_new{i+0*n_equation} = d1_big.';
    vector_pow_new{i+0*n_equation} = inverse_matrification(d2_big);
    dot_new{i+0*n_equation} = inverse_matrification(d3_big);
    n_terms_new(i+0*n_equation) = length(local_value);
    
    value_new{i+n_equation} = c_y;
    scalar_pow_new{i+n_equation} = d1_y.';
    vector_pow_new{i+n_equation} = d2_y;
    dot_new{i+n_equation} = d3_y;
    n_terms_new(i+n_equation) = length(c_y);
    
    value_new{i+2*n_equation} = c_yy_and_z;
    scalar_pow_new{i+2*n_equation} = d1_yy_and_z.';
    vector_pow_new{i+2*n_equation} = d2_yy_and_z;
    dot_new{i+2*n_equation} = d3_yy_and_z;
    n_terms_new(i+2*n_equation) = length(c_yy_and_z);
    
    value_yy{i} = c_yy;
    scalar_pow_yy{i} = d1_yy.';
    vector_pow_yy{i} = d2_yy;
    dot_yy{i} = d3_yy;
    n_terms_yy(i) = length(c_yy);
end

if n_equation>0
    new_pol = polynomial_coefs(3*old_pol.size_scalar, 3*old_pol.size_vector, 3*n_equation, ...
                n_terms_new, value_new, scalar_pow_new, vector_pow_new, dot_new);
    DxHyy = polynomial_coefs(3*old_pol.size_scalar, 3*old_pol.size_vector, n_equation, ...
                n_terms_yy, value_yy, scalar_pow_yy, vector_pow_yy, dot_yy);
else
    % if there aren't, standard empty polynomial
    new_pol = polynomial_coefs(3*old_pol.size_scalar, 3*old_pol.size_vector, 0, ...
        [], [],[],[]);
    DxHyy = new_pol;
end

end

% function [c_new, d1_new, d2_new, d3_new] = saddle_managment(c,d1,d2,d3)
% % function [c_new, d1_new, d2_new, d3_new] = saddle_managment(c,d1,d2,d3)
% %
% % managing the changes applied to a single equation to change for the saddle
% % system
% if min(size(c))~=1 || size(d1,1) ~= length(c) || size(d2,1) ~= length(c) ...
%         || size(d3,1) ~=length(c) || size(d2,2)~=size(d3,2)
%     error('Inputs not accepted')
% end
% 
% [c_big, d1_big, d2_big, d3_big] = padding(c,d1,d2,d3);
% [c_y, d1_y, d2_y] = Dx_y(c_big,d1_big,d2_big,d3_big);
% [c_yy, d1_yy, d2_yy] = Dx_y(c_y,d1_y,d2_y,d3_y);
% [c_z, d1_z, d2_z] = Dx_z(c_big,d1_big,d2_big,d3_big);
% 
% [c, d1, d2, d3] = sum_coef(c_yy, d1_yy, d2_yy,c_z, d1_z, d2_z);
% 
% [c_new, d1_new, d2_new, d3_new] = contract(c,d1,d2,d3);
% end


function [c_new, d1_new, d2_new, d3_new] = sum_coef(c,d1,d2,d3,C,D1,D2,D3)
% function [c_new, d1_new, d2_new, d3_new] = sum_coef(c,d1,d2,d3,C,D1,D2,D3)
% 
% with c, d1,d2,d3 referring to a sum of monomials, and C, D1,D2,D3 as
% well, the output will refer to the sum of the two monomials
% INPUT
% c         real vector, length K
% d1        integer matrix, rows K, columns M
% d2, d3    integer matrix, rows K, columns N
% C         real vector, length K'
% D1        integer matrix, rows K', columns M
% D2, D3    integer matrix, rows K', columns N
% OUTPUT
% c_new     real vector, length K+K'
% d1_new    integer matrix, rows K+K', columns M
% d2, d3    integer matrix, rows K+K', columns N

c_new = [vert(c);vert(C)];
d1_new = [d1;D1];
d2_new = [d2;D2];
d3_new = [d3;D3];

end

function [c_big, d1_big, d2_big, d3_big] = padding(c,d1,d2,d3)
% function [c_big, d1_big, d2_big, d3_big] = padding(c,d1,d2,d3)
% consider the monomial sum S c(i) lambda^d1(i,:) x^d2(i,:) x_dot^d3(i,:) and compose the new
% powers to consider the monomial 
% S c(i) (lambda, alpha, beta)^d1(i,:) (x,y,z)^d2(i,:)
% (x_dot,y_dot,z_dot)^d3(i,:)
%
% INPUT
% c         real vector, length K
% d1        integer matrix, rows K, columns M
% d2, d3    integer matrix, rows K, columns N
% OUTPUT
% c_big             real vector, length K
% d1_big            integer matrix, rows K, columns 3M
% d2_big, d3_big    integer matrix, rows K, columns 3N

if any(size(d2)~=size(d3))
    error('Debugging problem: the length of the powers and of the derivatives should be the same')
    % i.e. d2 refers to the powers of x, while d3 to the powers of x_dot,
    % they should be of the same length
end
if size(d1,1) ~= length(c) || size(d2,1) ~= length(c)
    error('Dimensions are not matching')
end

c_big = c;
d1_big = [d1,0*d1,0*d1];
d2_big = [d2,0*d2,0*d2];
d3_big = [d3,0*d3,0*d3];
end

function [c_new, d1_new, d2_new, d3_new] = Dx_var_mon (c,d1,d2,d3,var)
% function [c_new, d1_new, d2_new, d3_new] = Dx_var_mon (c,d1,d2,d3,var)
% INPUTS
% data reffering to a monomial c (lambda, alpha, beta)^d1 (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3
% c        real scalar
% d1       integer vector, length 3M, scalar power
% d2       integer vector, length 3N, vector power
% d3       integer vector, length 3N, derivative power
% d3       integer matrix, columns 3N, rows K, derivative power
% var      integer, 2 or 3, for 2 Dx_y, for 3 Dx_z
% OUTPUT
% data referring to a sum of monomials D_x [  c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ] var
% c_new        real vector
% d1_new       integer matrix, columns 3M, rows M+2N, scalar power
% d2_new       integer matrix, columns 3N, rows M+2N, vector power
% d3_new       integer matrix, columns 3N, rows M+2N, derivative power

if length(c)~=1
    error('Call Dx_y instead, here we deal with each monomial at the time')
end
if mod(length(d1),3)~=0 ||  mod(length(d2),3)~=0 ||  mod(length(d3),3)~=0
   error('The length of the powers should be multiple of 3')
end
if length(d2)~=length(d3)
    error('Debugging problem: the length of the powers and of the derivatives should be the same')
    % i.e. d2 refers to the powers of x, while d3 to the powers of x_dot,
    % they should be of the same length
end
d1 = horiz(d1);
d2 = horiz(d2);
d3 = horiz(d3);

M = length(d1)/3;
N = length(d2)/3;
e_M = eye(M);
e_N = eye(N);
e_MN = eye(M,N);
e_NM = eye(N,M);

c_new = c*[vert(d1(1:M));
    vert(d2(1:N));
    vert(d3(1:N))];

first = @(d,N) d(1:N);
second = @(d,N) d(1+N:2*N);
third = @(d,N) d(2*N+1:3*N);

rep_mat = @(d,N) repmat(d,N,1);

big_mat = @(d,N_var) [ rep_mat( first(d,N_var),M),   rep_mat( second(d,N_var), M),   rep_mat( third(d,N_var), M);
           rep_mat( first(d,N_var),N),   rep_mat( second(d,N_var), N),   rep_mat( third(d,N_var), N);
           rep_mat( first(d,N_var),N),   rep_mat( second(d,N_var), N),   rep_mat( third(d,N_var), N)];

Der = @(row, column,N_var) [ -eye(M,N_var)*(row==1) , eye(M,N_var)*(row==1)*(column==2), eye(M,N_var)*(row==1)*(column==3);
    -eye(N,N_var)*(row==2), eye(N,N_var)*(row==2)*(column==2), eye(N,N_var)*(row==2)*(column==3)
    -eye(N,N_var)*(row==3), eye(N,N_var)*(row==3)*(column==2), eye(N,N_var)*(row==3)*(column==3)];
       
       
d1_new = big_mat(d1,M) + Der(1,var,M);
%[ rep_mat( first(d1,M),M) - e_M,   rep_mat( second(d1,M), M) + e_M,   rep_mat( third(d1,M), M);
%           rep_mat( first(d1,M),N),   rep_mat( second(d1,M), N),   rep_mat( third(d1,M), N);
%           rep_mat( first(d1,M),N),   rep_mat( second(d1,M), N),   rep_mat( third(d1,M), N)];
       
       
d2_new =  big_mat(d2,N) + Der(2,var,N);
%[ rep_mat( first(d2,N),M),   rep_mat( second(d2,N), M),   rep_mat( third(d2,N), M);
%           rep_mat( first(d2,N),N) - e_N,   rep_mat( second(d2,N), N) + e_N,   rep_mat( third(d2,N), N);
%           rep_mat( first(d2,N),N),   rep_mat( second(d2,N), N),   rep_mat( third(d2,N), N)];
       
d3_new =   big_mat(d3,N) + Der(3,var,N);
%[ rep_mat( first(d3,N),M),   rep_mat( second(d3,N), M),   rep_mat( third(d3,N), M);
%            rep_mat( first(d3,N),N),   rep_mat( second(d3,N), N),   rep_mat( third(d3,N), N);
%            rep_mat( first(d3,N),N) - e_N,   rep_mat( second(d3,N), N)+ e_N,   rep_mat( third(d3,N), N)];
end


function [c_new, d1_new, d2_new, d3_new] = Dx_var(c,d1,d2,d3,var)
% function [c_new, d1_new, d2_new, d3_new] = Dx_var(c,d1,d2,d3,var)
% INPUTS
% data referring to a sum of monomials sum [  c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ] 
% c        real vector, length K
% d1       integer matrix, columns 3M, rows K, scalar power
% d2       integer matrix, columns 3N, rows K, vector power
% d3       integer matrix, columns 3N, rows K, derivative power
% var      integer, 2 or 3, for 2 Dx_y, for 3 Dx_z
% OUTPUT
% data referring to a sum of monomials D_x [   sum (c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ) ] VAR
% c_new        real vector, length K'
% d1_new       integer matrix, columns 3M, rows K', scalar power
% d2_new       integer matrix, columns 3N, rows K', vector power
% d3_new       integer matrix, columns 3N, rows K', derivative power

K = length(c);
if size(d1,1)~= K || size(d2,1)~= K || size(d3,1)~= K
    error('The inputs are not consistent')
end
if mod(size(d1,2),3)~= 0 ||mod(size(d2,2),3)~= 0 || mod(size(d3,2),3)~= 0 
    error('The inputs are not consistent')
end
if size(d2,2)~= size(d3,2)
    error('Inputs are not consistent')
end
if var~=3 && var~=2
    error('Impossible')
end

M = size(d1,2)/3;
N = size(d2,2)/3;

c_new = zeros( K * (M +2*N), 1);
d1_new = zeros( K * (M +2*N), 3*M);
d2_new = zeros( K * (M +2*N), 3*N);
d3_new = zeros( K * (M +2*N), 3*N);

for i = 1:K
    indeces_out = (1:(M+2*N)) +(i-1)*(M+2*N);
    [c_new(indeces_out), d1_new(indeces_out,:), d2_new(indeces_out,:),...
        d3_new(indeces_out,:)] = Dx_var_mon ( c(i), d1(i,:), d2(i,:), d3(i,:),var);
end

end

function [c_new, d1_new, d2_new, d3_new] = Dx_y(c,d1,d2,d3)
% function [c_new, d1_new, d2_new, d3_new] = Dx_y(c,d1,d2,d3)
% INPUTS
% data referring to a sum of monomials sum [  c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ] 
% c        real vector, length K
% d1       integer matrix, columns 3M, rows K, scalar power
% d2       integer matrix, columns 3N, rows K, vector power
% d3       integer matrix, columns 3N, rows K, derivative power
% OUTPUT
% data referring to a sum of monomials D_x [   sum (c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ) ] y
% c_new        real vector, length K'
% d1_new       integer matrix, columns 3M, rows K', scalar power
% d2_new       integer matrix, columns 3N, rows K', vector power
% d3_new       integer matrix, columns 3N, rows K', derivative power
[c_new, d1_new, d2_new, d3_new] = Dx_var(c,d1,d2,d3,2);
end
function [c_new, d1_new, d2_new, d3_new] = Dx_z(c,d1,d2,d3)
% function [c_new, d1_new, d2_new] = Dx_z(c,d1,d2,d3)
% INPUTS
% data referring to a sum of monomials sum [  c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ] 
% c        real vector, length K
% d1       integer matrix, columns 3M, rows K, scalar power
% d2       integer matrix, columns 3N, rows K, vector power
% d3       integer matrix, columns 3N, rows K, derivative power
% OUTPUT
% data referring to a sum of monomials D_x [   sum (c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ) ] z
% c_new        real vector, length K'
% d1_new       integer matrix, columns 3M, rows K', scalar power
% d2_new       integer matrix, columns 3N, rows K', vector power
% d3_new       integer matrix, columns 3N, rows K', derivative power
[c_new, d1_new, d2_new, d3_new] = Dx_var(c,d1,d2,d3,3);
end


function [c_new, d1_new, d2_new, d3_new] = contract(c,d1,d2,d3)
% function [c_new, d1_new, d2_new, d3_new] = contract(c,d1,d2,d3)
% INPUTS
% data referring to a sum of monomials sum [  c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ], including zero terms
% c        real vector, length K
% d1       integer matrix, columns 3M, rows K, scalar power
% d2       integer matrix, columns 3N, rows K, vector power
% d3       integer matrix, columns 3N, rows K, derivative power
% OUTPUT
% data referring to a sum of monomials sum [  c (lambda, alpha, beta)^d1
% (x,y,z)^d2 (x_dot,y_dot,z_dot)^d3 ], EXCLUDING zero terms
% c_new        real vector, length K'
% d1_new       integer matrix, columns 3M, rows K', scalar power
% d2_new       integer matrix, columns 3N, rows K', vector power
% d3_new       integer matrix, columns 3N, rows K', derivative power


indeces_empty = (c ==0);
c_new = c;
c_new(indeces_empty) = [];
d1_new = d1;
d2_new = d2;
d3_new = d3;
d1_new(indeces_empty,:) = [];
d2_new(indeces_empty,:) = [];
d3_new(indeces_empty,:) = [];

A = [d1_new, d2_new, d3_new]; % all powers in one go, to test for summable elements
length_A = size(A,1);
for i=1:length_A-1
    index_i = [];
    for j = i+1:length_A
        if all(A(i,:)==A(j,:))
            index_i = [index_i, j];
        end
    end
    c_new( i ) = c_new(i) + sum(c_new(index_i));
    c_new(index_i ) = 0;
end
indeces_empty = (c_new ==0);
c_new(indeces_empty) = [];
d1_new(indeces_empty,:) = [];
d2_new(indeces_empty,:) = [];
d3_new(indeces_empty,:) = [];
end

function matrix_A = matrification(cell_A)
% INPUT
% cell_A     1D cell, such that each cell_A is a vector of the same length
% OUTPUT
% matrix_A   2D vector, matrix_A(i,j) = cell_A{i}(j)
if min(size(cell_A)) ~=1
    error('Requested a 1D cell')
end
I = max(size(cell_A));
length_A = length(cell_A{1});
matrix_A = zeros(I, length_A);

for i=1:I
    if length(cell_A{i})~=length_A
        error('The input cell needs to have elements of the same length')
    end
    if min(size(cell_A{i}))~=1
        error('elements of the cell needs to be vectors')
    end
    matrix_A(i,:) = horiz(cell_A{i}(:));
end

end

function cell_A = inverse_matrification(matrix_A)
% INPUT
% matrix_A   2D vector, matrix_A(i,j) = cell_A{i}(j)
% OUTPUT
% cell_A     1D cell, such that each cell_A is a vector of the same length
if length(size(matrix_A))>2
    error('The input needs to be a matrix')
end

cell_A = cell(size(matrix_A,1),1);
for i=1:size(matrix_A,1)
    cell_A{i} = matrix_A(i,:).';
end

end
