function [s,x1] = continuation_numerical( x0, F_not_square, n_iter, h, x_dot_0,s, min_res_N, bool_Hopf, bool_fancy_scalar)
% function [s,x1] = continuation_numerical( x0, F_not_square, n_iter, h, x_dot_0,s, min_res_N, bool_Hopf, bool_fancy_scalar)
% 
% numerical continuation, in the same (or very similar) setting to
% continuation (that is also validated on the go).
% INPUTS
% x0              Xi_vector, first approximation of the solution
% F_not_square    full_problem, underdermined, usually 1D of freedom,
% depending on the extra choices explained later
% n_iter          number of iterations of the continuation algorithm to be
% done (DEFAULT 10)
% h               initial stepsize (DEFAULT 10^-5) 
% x_dot_0         direction to follow (DEFAULT kernel(DF(x0))
% s               string, the results will be saved there (DEFAULT
% 'simulation_of_the_%s_at_%s',datestr(now,'ddmmmyy'),datestr(now,'HHMM'))
% min_res_N       requested residual in Newton (DEFAULT 10^-14)
% bool_Hopf       bool, is the system Hopf (DEFAULT 0) if the system is
% Hopf, the necessary restriction apply, in particular,
% F_not_square.scalar_eq will be overwritten and shouls have 3 linear
% scalar eqautions. The numebr of scalar unknowns should be 3+the number of
% vector unknown
% bool_fancy_scalar   bool, requesting the FIRST linear scalar equation to
% be \int_[0,2pi] x(t) dot x'(t) dt = 0 and updated at every iteration

global talkative
global use_intlab

use_intlab =0;
old_talkative = talkative;
talkative =0;

% test of inputs
if ~isa(x0,'Xi_vector')
    error('Invalid input solution')
end
if ~isa(F_not_square,'full_problem')
    error('Invalid input system')
end
if nargin<3 || isempty(n_iter)
    n_iter = 10;
elseif ~isa(n_iter,'double')
    warning('Invalid number of iterations, reset to the default value')
    n_iter = 10;
end
if nargin<4 || isempty(h)
    h = 10^-5;
elseif ~isa(h,'double')
    warning('Invalid number of iterations, reset to the default value')
    h = 10^-5;
end
if nargin<6 || isempty(s)
    s = sprintf('simulation_of_the_%s_at_%s',datestr(now,'ddmmmyy'),datestr(now,'HHMM'));
end
if nargin<7|| isempty(min_res_N)
    min_res_N = 10^-14;
end
if nargin<8 || isempty(bool_Hopf)
    bool_Hopf =0;
end
if nargin<9 || isempty(bool_fancy_scalar)
    bool_fancy_scalar = 0;
end

if bool_Hopf
    if x0.size_scalar~=x0.size_vector+3
        error('If you call for a Hopf validation, the number of scalar variables should be the number of vector variables + 3')
    end
end


% some useful variables
size_scalar = x0.size_scalar;
size_vector = x0.size_vector;
nodes = x0.nodes;

if F_not_square.scalar_equations.num_equations ~= size_scalar-1 &&( bool_Hopf && F_not_square.scalar_equations.num_equations ~= size_scalar-2)
    error('Invalid input problem: the number of scalar equations is inconsistent')
end

if F_not_square.vector_field.n_equations ~= size_vector
    error('Invalid input problem: the number of vector equations is inconsistent')
end

size_x = size_scalar+size_vector*(nodes*2+1);

if ~(nargin<5 || isempty(x_dot_0)) && (length(x_dot_0)~=size_x || min(size(x_dot_0)) ~= 1)
    warning('Invalid direction, default direction computed')
    x_dot_0 = [];
end


F_old = F_not_square;
% addition of some equations for Hopf
if bool_Hopf
    %DF0_0 = derivative_to_matrix(derivative(F_old,x0,0));
    %DF0_0(1,:)=0;
    p0 = zeros(size_x,1);%conj(kernel(DF0_0));
    p0(1:5) = 0;
    p0(6:(5+2*x0.nodes+1)) = 1;
    p0((5+2*x0.nodes+2):end) = 0;
    const_p0 = 0;
    F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[p0;const_p0],1);
    
    p1 = Xi_vec2vec(x0)';
    p1(1:5) = 0;
    const_p1 = -1;
    F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[p1.';const_p1],2);
    F_not_square = F_old;
end

% addition of some equations for fancy_scalar_equation
if bool_fancy_scalar
    F_old.scalar_equations = fancy_scalar_condition(x0,F_old.scalar_equations);
end


% last input test
if nargin<5 || isempty(x_dot_0) % definition of x_dot if not defined yet
    DF0 = derivative_to_matrix(derivative(F_old,x0,0));
    x_dot_0 = kernel(DF0);
    if min(size(x_dot_0)) ~= 1
        error('At the computed point, the derivative is singular')
    end
end
% finished test of inputs

% construction of old system
[~,index] = max(abs(real(x_dot_0(1:x0.size_scalar))));
angle = atan( imag(x_dot_0(index))/real(x_dot_0(index)));
x_dot_0 = exp( - 1i * angle) * x_dot_0; % bringing x_dot_0 to be symmetric (by multiplication with the appropriate complex rotation)
x_dot_0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_0,x0)));
const = -Xi_vec2vec(x0).'*x_dot_0;
F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[x_dot_0;const],F_old.scalar_equations.number_equations_lin+1); % addition of the continuation equation


% preallocate storage
norm_x = zeros(x0.size_scalar+x0.size_vector,n_iter);
% stored_x = zeros(x0.size_scalar+x0.size_vector*(2*x0.nodes+1),30);
% stored_xdot = stored_x;

fprintf('Start of numerical continuation\n')


for i =1 : n_iter
    % compute new step (and system to solve at the new step) with a
    % predictor-corrector algorithm
    if old_talkative>0
        if mod(i, 50)==0
            fprintf('Iteration %i\n',i);
        end
    end
    [x1,x_dot_1,~, iter_Newton] =update(x0,h,x_dot_0,...
        F_not_square,[],min_res_N, bool_Hopf, bool_fancy_scalar);
    
    if iter_Newton<2
        h = 2*h;
    elseif iter_Newton >3
        h = 0.5*h;
    end
    
    % iteration
    x_dot_0 = x_dot_1;
    x0 = x1;
    
    %storage
    norm_x(:,i)    = vert(norm(x0));
%     if i > 220 && i <250
%         stored_x(:,i - 220)    = Xi_vec2vec(x0);
%         stored_xdot(:,i - 220) = x_dot_0;
%     end
end
fprintf('End of numerical continuation\n')

% save to location the elements that are interesting
save(s, 'norm_x','x1','x_dot_1'); % 'stored_x', 'stored_xdot',
talkative = old_talkative;
return