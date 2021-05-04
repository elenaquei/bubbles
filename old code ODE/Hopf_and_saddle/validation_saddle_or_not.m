function [saddle, index_saddle_confirmed,rad]=validation_saddle_or_not...
    (F_old,F_new, x0,x1,num_variable)
% function [saddle_confirmed,saddle_x0,saddle_x1,rad]=validation_not_saddle...
%     (F_old,F_new, x0,x1,validation_radius,num_variable)
% 
% INPUT
% F_old, F_new      full_problem, square
% x0, x1            Xi_vector, solutions to F_old and F_new respectively
%                   the segment is assumed to have been validated
% num_variable      vector of integers, optional, parameter(s) w.r.t. wich 
%                   we can have saddle
% OUTPUT
% index_saddle_confirmed  index of saddle parameter
% saddle_x0         refinement for the saddle
% saddle_x1         other end point
% rad               validation radus of the saddle 
global use_intlab 
global talkative

temp_use_intlab = use_intlab;
use_intlab = 0;

if ~exist('num_variable','var')|| isempty(num_variable) ...
        ||  all(num_variable== 0)
    num_variable=1:x0.size_scalar;
end
%saddle_x0=[];saddle_x1=[];rad=[];

F0 = F_old;
F1 = F_new;

x0_prime = der_s(x0,F0,F1,0);
x1_prime = der_s(x1,F0,F1,1);

% [bool,num_variable] = is_a_saddle(x0_prime, x1_prime, num_variable);
% if ~bool
%     index_saddle_confirmed = 0;
%     return
% end
% if talkative 
%     disp('Found a segment with a numerical saddle')
% end

% average the linear coefficients out, all except continuation equation
lin_coef_middle = cell(3,1);
lin_coef_middle{1} = (1/2)*(F1.scalar_equations.linear_coef{1} + F0.scalar_equations.linear_coef{1});
lin_coef_middle{2} = (1/2)*(F1.scalar_equations.linear_coef{2} + F0.scalar_equations.linear_coef{2});
lin_coef_middle{3} = (1/2)*(F1.scalar_equations.linear_coef{3} + F0.scalar_equations.linear_coef{3});
 
F0.scalar_equations.linear_coef{1}(1:end-1,:) = lin_coef_middle{1}(1:end-1,:);
F0.scalar_equations.linear_coef{2}(1:end-1,:,:) = lin_coef_middle{2}(1:end-1,:,:);
F0.scalar_equations.linear_coef{3}(1:end-1) = lin_coef_middle{3}(1:end-1);
F1.scalar_equations.linear_coef{1}(1:end-1,:) = lin_coef_middle{1}(1:end-1,:);
F1.scalar_equations.linear_coef{2}(1:end-1,:,:) = lin_coef_middle{2}(1:end-1,:,:);
F1.scalar_equations.linear_coef{3}(1:end-1) = lin_coef_middle{3}(1:end-1);

% no refinement
%n_subinterval = 10;
%for i = 1:n_subinterval
%    [x0,x1,x0_prime,x1_prime, F0, F1] = bisect_saddle(x0,x1,...
%        x0_prime,x1_prime, F0, F1,num_variable);
%end

h = max(norm(x0-x1));

x0_prime = h^-1 * x0_prime;
x1_prime = h^-1 * x1_prime;

x0_primeprime = h^-2 * der_ss(x0,x0_prime, F0, F1, 0);
x1_primeprime = h^-2 * der_ss(x1,x1_prime, F0, F1, 1);

big_x0 = merge(x0,x0_prime, x0_primeprime);
big_x1 = merge(x1,x1_prime, x1_primeprime);

[big_F0, big_F1] = rescaled_bigger_saddle_system(F0, F1, x0, x1, h*10);

big_x0_N = Newton_2(big_x0, big_F0, 30, 10^-10);
big_x1_N = Newton_2(big_x1, big_F1, 30, 10^-10);

if max(abs(big_x0_N.scalar(x0.size_scalar + num_variable)))<10^-10 || ...
        max(abs(big_x1_N.scalar(x0.size_scalar + num_variable)))<10^-10
    warning on
    warning('The validation will likely fail due to the first derivative being incredibly small')
    warning off
end

if length(num_variable)==1
    [big_x0_res,big_x1_res,big_F0_res,big_F1_res] = rescale_saddle...
        (big_x0_N,big_x1_N,big_F0,big_F1,num_variable);
    big_x0_res = Newton_2(big_x0_res, big_F0_res, 30, 10^-10);
    big_x1_res = Newton_2(big_x1_res, big_F1_res, 30, 10^-10);
else
    big_x0_res = big_x0_N;
    big_x1_res = big_x1_N;
    big_F0_res = big_F0;
    big_F1_res = big_F1;
end

use_intlab = temp_use_intlab;
[validated, saddle, index_saddle_confirmed,rad] = saddle_validation...
    (big_x0_res,big_x1_res,big_F0_res,big_F1_res,num_variable);

if ~ validated
    error('Could not validate the existence or absence of saddle node in this segment')
end

return

end
 
function [validated, saddle, index_saddle_confirmed,rad] = saddle_validation(saddle_x0, saddle_x1, saddle_problem0, saddle_problem1,num_variable)

size_scalar = saddle_x0.size_scalar/3;
y0_point = (saddle_x0.scalar(size_scalar+1:size_scalar*2));
y1_point = (saddle_x1.scalar(size_scalar+1:size_scalar*2));
 
z0_point = (saddle_x0.scalar(size_scalar*2+1:size_scalar*3));
z1_point = (saddle_x1.scalar(size_scalar*2+1:size_scalar*3));

[~,Imin_sad_y0]=radii_polynomials(saddle_x0,saddle_problem0);

[~,Imin_sad_y1]=radii_polynomials(saddle_x1,saddle_problem1);

if any(abs(y0_point(num_variable))<Imin_sad_y0) ||...
    any(abs(y1_point(num_variable))<Imin_sad_y1)
    error('Bounds for Y still too big! \n Impossible to validate the saddle because even end point derivative sign is not validated')
end


%[flag,I_min] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
%    saddle_problem0,saddle_problem1,[],Aold_saddle,Anew_saddle);
%if flag <=0
    [flag,I_min] = substepping(saddle_x0,saddle_x1,saddle_problem0,saddle_problem1);
    if flag <=0
        error('Could not validate saddle')
    end
%end





y0 = midrad(y0_point,Imin_sad_y0);
y1 = midrad(y1_point,Imin_sad_y1);
 
z0 = midrad(z0_point,I_min);
z1 = midrad(z1_point,I_min);

% option no saddle
index_saddle = find(sup(y0(num_variable).*y1(num_variable))<0,1);
if isempty(index_saddle)
    if all(inf(y0(num_variable).*y1(num_variable))>0)
        validated = 1;
        saddle = 0;
        index_saddle_confirmed = [];
        rad = [];
        return
    else
        validated = 0;
        saddle = 0;
        index_saddle_confirmed = [];
        rad = [];
        return
    end
else
    if all(inf(z0(num_variable(index_saddle)).*z1(num_variable(index_saddle)))>0)
        validated = 1;
        saddle = 1;
        index_saddle_confirmed = index_saddle;
        rad = I_min;
        return
    else
        validated = 0;
        saddle = 1;
        index_saddle_confirmed = index_saddle;
        rad = [];
        return
    end
end




end
 




function [ x_prime0 ] = der_s(x,F0,F1, begin)
% function [ x_prime0 ] = der_s(x,F0,F1, begin)
%
% numerical computation of the first derivative w.r.t. the archlength
%
% INPUTS
% F0, F1        full_problem, at the beginning and at the end of the
%               segment
% x0            Xi_vector, numerical solution
% begin         0 - refers to x0
%               1 - refers to x1
%
% OUTPUTS
% x_prime0      Xi_vector, approximate first derivatives w.r.t. the
%               archlength
if nargin==3
    begin =1;
end

if begin == 0   
    Dx0H = derivative_to_matrix(derivative(F0,x,0));
else
    Dx0H = derivative_to_matrix(derivative(F1,x,0));
end
    
x0_vec = vert(Xi_vec2vec(x));
v_0 = extract_all_lin_coef(F0.scalar_equations);
v_1 = extract_all_lin_coef(F1.scalar_equations);
v_Delta = v_1 - v_0;

c_0 = vert(F0.scalar_equations.linear_coef{3});
c_1 = vert(F1.scalar_equations.linear_coef{3});
c_Delta = c_1-c_0;

Ds0H = v_Delta*x0_vec + 0*v_0*x0_vec +c_Delta;
Ds0H = cat(1,Ds0H, zeros(length(v_0)-F0.scalar_equations.number_equations_lin,1));

x_prime0 = - Dx0H \ Ds0H;

x_prime0 = vec2Xi_vec(x_prime0,x);

x_prime0 = symmetrise(x_prime0);
end

function x_primeprime = der_ss(x,x_prime, F0, F1, begin)
% function x1_primeprime = der_ss(x1,x1_prime, F0, F1, begin)
%
% numerical computation of the second derivative w.r.t. the archlength
%
% INPUTS
% x0            Xi_vector, numerical solution
% x0_prime      Xi_vector, numerical approximation of first derivative
% F0, F1        full_problem, at the beginning and at the end of the
%               segment
% begin         0 - refers to x0
%               1 - refers to x1
%
% OUTPUTS
% x_primeprime  Xi_vector, approximate second derivatives w.r.t. the
%               archlength

if begin == 0
    DxH = derivative_to_matrix(derivative(F0,x,0));
else
    DxH = derivative_to_matrix(derivative(F1,x,0));
end

v_0 = extract_all_lin_coef(F0.scalar_equations);
v_1 = extract_all_lin_coef(F1.scalar_equations);
v_Delta = v_1 - v_0;

DsxH = v_Delta;
DsxH = cat(1,DsxH, zeros(length(v_0)-F0.scalar_equations.number_equations_lin,length(v_0)));

x_primeprime = -2* (DxH \ (DsxH* Xi_vec2vec(x_prime)));

x_primeprime = vec2Xi_vec(x_primeprime,x_prime);

x_primeprime = symmetrise(x_primeprime);
end


function [bool,num_var] = is_a_saddle(x0_prime, x1_prime, num_variable)
% function [bool,num_var] = is_a_saddle(x0_prime, x1_prime, num_variable)
% 
% given two derivatives, find if there is a saddle
%
% INPUT
% x0_prime, x1_prime        Xi_vectors, first derivative w.r.t. archlength
% num_variable              which scalar to take care off
% OUTPUT
% bool                      possibility of saddle 
%                           (first derivative changed sign)
% num_var                   index of the variable with the saddle node

lambda0_prime = x0_prime.scalar(num_variable);
lambda1_prime = x1_prime.scalar(num_variable);

if any(lambda0_prime.*lambda1_prime<-10^-25)
    bool = 1;
else
    bool = 0;
end

if bool && length(num_variable)>1
    num_var = find ((lambda0_prime.*lambda1_prime)<0);
else
    num_var = num_variable;
end
end

function big_x = merge(x,x_prime, x_primeprime)
% function big_x = merge(x,x_prime, x_primeprime)
% 
% INPUT
% three Xi_vector
% OUTPUT 
% one, merged Xi_vector
big_x = Xi_vector([x.scalar,x_prime.scalar,x_primeprime.scalar],...
    [x.vector;x_prime.vector;x_primeprime.vector]);

end


function [big_x0_res,big_x1_res,big_F0_res,big_F1_res] = rescale_saddle...
    (big_x0,big_x1,big_F0,big_F1,num_var)
% function [big_x0_res,big_x1_res,big_F0_res,big_F1_res] = rescale_saddle...
%     (big_x0,big_x1,big_F0,big_F1,radius,num_var)
%
% addign rescaling to the problem, in such a way that the first derivative
% and the second derivative are big enough for the validation

max_scale = max(norm(big_x0-big_x1));

n_scal = big_x0.size_scalar/3;
prime_index = n_scal + num_var;

lambda0_prime = big_x0.scalar(prime_index);
lambda1_prime = big_x1.scalar(prime_index);

rescale_derivative = min(max(max_scale./ abs(lambda0_prime-lambda1_prime)),10^7);
% pushing the limit of the rescaling

big_F0_res = rescale(big_F0, prime_index, rescale_derivative);
big_F1_res = rescale(big_F1, prime_index, rescale_derivative);

% rescale scalar variables in the same way
big_x0_res = big_x0;
big_x1_res = big_x1;
big_x0_res.scalar(prime_index) = rescale_derivative * big_x0.scalar(prime_index);
big_x1_res.scalar(prime_index) = rescale_derivative * big_x1.scalar(prime_index);

if abs(big_x0_res.scalar(prime_index)) < 10 ^-11 ||...
    abs(big_x1_res.scalar(prime_index)) < 10 ^-11
    warning('Probably validation will be impossible')
end

end

function [x0,x1,x0_prime,x1_prime, F0, F1] = bisect_saddle(x0,x1,...
        x0_prime,x1_prime, F0, F1,num_variable)
% function [x0,x1,x0_prime,x1_prime, F0, F1] = bisect_saddle(x0,x1,...
%         x0_prime,x1_prime, F0, F1,num_variable)
%
% bisectthe segment and chose on which side is the saddle node
x_middle = 0.5*(x1+x0);
F_middle = F0;
lin_coef_middle = cell(3,1);

lin_coef_middle{1} = (1/2)*(F1.scalar_equations.linear_coef{1} + F0.scalar_equations.linear_coef{1});
lin_coef_middle{2} = (1/2)*(F1.scalar_equations.linear_coef{2} + F0.scalar_equations.linear_coef{2});
lin_coef_middle{3} = (1/2)*(F1.scalar_equations.linear_coef{3} + F0.scalar_equations.linear_coef{3});
 
F_middle.scalar_equations.linear_coef{1} = lin_coef_middle{1};
F_middle.scalar_equations.linear_coef{2} = lin_coef_middle{2};
F_middle.scalar_equations.linear_coef{3} = lin_coef_middle{3};
 
x_middle = Newton_2(x_middle, F_middle);

x_middle_prime = der_s(x_middle, F0, F_middle, 1);

bool = is_a_saddle(x0_prime, x_middle_prime,num_variable);

if bool
    x1 = x_middle;
    F1 = F_middle;
else
    x0 = x_middle;
    F0 = F_middle;
end

end
