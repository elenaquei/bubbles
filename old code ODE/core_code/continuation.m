function [s, x_n] = continuation ( x0, F_not_square, n_iter, h, x_dot_0,s, ...
    min_res_N, bool_Hopf,bool_fancy_scalar, bool_saddle,num_var_saddle)
% function [s, x_n] = continuation ( x0, F_not_square, n_iter, h, x_dot_0,s,...
%   min_res_N, bool_Hopf,bool_fancy_scalar, bool_saddle,num_var_saddle)
%
% INPUT
% x0             Xi_vector, solution at one point
% F_not_square   full_problem, underdetermined
% n_iter         integer, number of interations (DEFAULT: 10)
% h              real, initial step size (DFAULT: 10^-5)
% x_dot_0        complex vector of length M+N*(nodes*2+1), kernel of
%                DF(x0) (DEFAULT: kernel(DF(x0))
% min_res_N      double, Newton's residual ( DEFAULT: 10^-14 )
% bool_Hopf      bool, if 1 appended scalar equation during continuation
%                ||sol||=1 and y_1(0)=0 (DEFAULT: 0)
% bool_fancy_scalar   bool, requesting the FIRST linear scalar equation to
%                be \int_[0,2pi] x(t) dot x'(t) dt = 0 and updated at every
%                iteration (DEFAULT: 0)
% bool_saddle    bool, if 1 check for saddle nodes after each segment
%                (numerically) and if found one, validate it (DEFAULT: 0)
%
% OUTPUT
% s              string, path of saved solutions
% x_n            Xi_vector, last validated solution
%

global talkative
global use_intlab
global nu

temp_use_intval = use_intlab;
use_intlab =0;

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
    warning('Invalid atep size, reset to the default value')
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
if bool_Hopf
    if x0.size_scalar~=x0.size_vector+3
        error('If you call for a Hopf validation, the number of scalar variables should be the number of vector variables + 3')
    end
end
if nargin<9 || isempty(bool_fancy_scalar)
    bool_fancy_scalar = 0;
end
if nargin<10 || isempty(bool_saddle)
    bool_saddle = 0;
end

if bool_saddle+bool_Hopf>0 && nargin<11
    if bool_Hopf
        num_var_saddle = 2;
        index_amplitude = x0.size_scalar - x0.size_vector;
    else
        
        num_var_saddle = 1:x0.size_scalar;
        if length(num_var_saddle)>1
            warning('Impossible to rescale saddle node validation system due to searching for saddle nodes w.r.t. multiple variables')
        end
    end
end
if bool_Hopf
    
    index_amplitude = x0.size_scalar - x0.size_vector;
end


x_0 = x0; % used in the saving at the end

% data coming from the previous iteration
previous_iter.Y = [];
previous_iter.Z0 = [];
previous_iter.Z1 = [];
previous_iter.Z2 = [];

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


% addition of some equations for Hopf
F_old = F_not_square;

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

% addition of the 10th July 2018
if h<0
    x_dot_0=-x_dot_0;
    h=-h;
end
% finished test of inputs

% construction of old system
[~,index] = max(abs(real(x_dot_0(1:x0.size_scalar))));
angle = atan( imag(x_dot_0(index))/real(x_dot_0(index)));
x_dot_0 = exp( - 1i * angle) * x_dot_0; % bringing x_dot_0 to be symmetric (by multiplication with the appropriate complex rotation)
x_dot_0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_0,x0)));
%const = -Xi_vec2vec(x0).'*x_dot_0;
%F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[x_dot_0;const],F_old.scalar_equations.number_equations_lin+1); % addition of the continuation equation

[F_old,x_dot_0] = continuation_eq(F_old,x0,x_dot_0, 10^-14);


DF0 = derivative_to_matrix(derivative(F_old,x0,0));
Aold = inv(DF0);

% preallocate storage
step_size = zeros(n_iter,1);
number_modes = zeros(n_iter+1,1);
number_modes(1) = x0.nodes;
norm_x = zeros(x0.size_scalar+x0.size_vector,n_iter);
Interval = zeros(2,n_iter);
Z0_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Z1_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Z2_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Y_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);

saddle_x0_stored =cell(1);
saddle_x1_stored =cell(1);
saddle_y0 =[];
saddle_y1 =[];
saddle_z0 =[];
saddle_z1 =[];
saddle_index_var =[];

% validating old system % need a good starting point
x0 = Newton_2(x0,F_old,[],min_res_N);

% symmetrising initial condition and phase conditions
x0 = x0.symmetrise;
lin_coef_old = extract_all_lin_coef(F_old.scalar_equations);
for i = 1:F_old.scalar_equations.number_equations_lin
    new_lin_coef_i = Xi_vec2vec(symmetrise(vec2Xi_vec(lin_coef_old(i,:),x0)));
    new_lin_coef_i(end+1) = F_old.scalar_equations.linear_coef{3}(i);
    F_old.scalar_equations = change_lin_coef_vector(...
        F_old.scalar_equations,...
        new_lin_coef_i ,...
        i);
end
lin_coef_old = extract_all_lin_coef(F_not_square.scalar_equations);
for i = 1:F_not_square.scalar_equations.number_equations_lin
    new_lin_coef_i = Xi_vec2vec(symmetrise(vec2Xi_vec(lin_coef_old(i,:),x0)));
    new_lin_coef_i(end+1) = F_not_square.scalar_equations.linear_coef{3}(i);
    F_not_square.scalar_equations = change_lin_coef_vector(...
        F_not_square.scalar_equations,...
        new_lin_coef_i ,...
        i);
end

% - % - %
% update fancy scalar condition, x_dot_0 and continuation equation?

%try
if talkative >1
    fprintf('Validation of starting point\n');
end
use_intlab = temp_use_intval;
bool = radii_polynomials(x0,F_old,DF0,Aold);
use_intlab = 0;
%catch
if ~bool
    if talkative >1
        fprintf('Starting point require extra precision, running Newton\n');
    end
    x0 = Newton_2(x0,F_old,[],min_res_N);
    x0 = x0.symmetrise;
    
    % set up some requesed elements for the validation
    DF0 = derivative_to_matrix(derivative(F_old,x0,0));
    Aold = inv(DF0);
    %try
    % pointwise validation
    bool = radii_polynomials(x0,F_old,DF0,Aold);
    if ~bool
        error('Initial point not good enough, even after Newton. Consider changing the number of nodes and or nu');
    end
end

norm_x(:,1)    = vert(norm(x0));
norm_x(1:x0.size_scalar,1) = x0.scalar;

if talkative
    fprintf('Start of continuated validation\n')
end

for i =1 : n_iter
    % compute new step (and system to solve at the new step) with a
    % predictor-corrector algorithm
    [x1,x_dot_1,F_new] =update_new(x0,h,x_dot_0,...
        F_not_square,[],min_res_N, bool_Hopf,bool_fancy_scalar);
    
    % set up some requesed elements for the validation
    DF1 = derivative_to_matrix(derivative(F_new,x1,0));
    Anew = inv(DF1);
    
    use_intlab = temp_use_intval;
    
    % validate
    [flag,Imin,Imax,previous_iter,Yvector,Z0vector,Z1vector,Z2vector,new_step] = ...
        radii_polynomials_cont_new(x0,x1,DF0,DF1,...
        F_old,F_new,previous_iter,Aold,Anew);
    
    use_intlab = 0;
    
    if flag <1 % if the validation did not work, redo the previous steps
        tries =1;
        while flag<1 && tries<10
            previous_iter = [];
            h = h* new_step;
            tries = tries+1;
            if talkative >0
                fprintf('could not validate the %i-th interval, h decreased to %e\n',i, h)
            end
            % compute new step (and system to solve at the new step) with a
            % predictor-corrector algorithm
            use_intlab = 0;
            [x1,x_dot_1,F_new] =update_new(x0,h,x_dot_0,...
                F_not_square,[],min_res_N, bool_Hopf, bool_fancy_scalar);
            
            % set up some requesed elements for the validation
            DF1 = derivative_to_matrix(derivative(F_new,x1,0));
            Anew = inv(DF1);
            
            use_intlab = temp_use_intval;
            % validate
            [flag,Imin,Imax,~,Yvector,Z0vector,Z1vector,Z2vector,new_step] = radii_polynomials_cont_new(x0,x1,DF0,DF1,...
                F_old,F_new,[],Aold,Anew);
            
            use_intlab = 0;
        end
        if flag<1
            error('Could not validate the interval')
        end
    end
    % vadlidation succeeded or program crashed
    if bool_Hopf
        amplitude_0 = x0.scalar(index_amplitude);
        amplitude_1 = x1.scalar(index_amplitude);
    end
    if bool_saddle || ( bool_Hopf && amplitude_0*amplitude_1<0)
        
        
        %if numerical_check
        use_intlab = temp_use_intval;
        %[saddle_confirmed,saddle_x0,saddle_x1]=validation_saddle_new...
        %    (F_old,F_new, x0,x1,num_var_saddle);
        [saddle_confirmed, index_saddle_confirmed,rad]=validation_saddle_or_not...
            (F_old,F_new, x0,x1,num_var_saddle);
        use_intlab = 0;
        if saddle_confirmed                                                   %-%
            % save the initial solution, plus the derivatives of the
            % scalar parts, as well as the index of the variable that
            % has a saddle node
            if bool_Hopf
                disp('THE HOPF BIFURCATION HAS BEEN CONFIRMED!')
            else
                disp('THE SADDLE HAS BEEN CONFIRMED!')
            end
            %pause(3)
            
            % by the way, now we are confident of its uniqueness
            if saddle_confirmed>1
                warning('There are %i saddles in this segment', saddle_confirmed)
            end
            n_saddle = length(saddle_index_var)+1;
            saddle_iter(n_saddle) = i;
            saddle_x0_stored{n_saddle} = x0;
            saddle_x1_stored{n_saddle} = x1;
            saddle_index_var(n_saddle) = index_saddle_confirmed;
            %saddle_interval(n_saddle) = rad;
        end
    end
    
    
    
    % test on the smoothness condition
    flag_smoothness = smoothness_condition(x0, x_dot_0, x1,x_dot_1, Imin);
    if flag_smoothness == 0
        bool_refine_anyway = 1;
        [flag,Imin] = substepping(x0,x1,F_old,F_new,[], bool_refine_anyway);
        if ~flag
            error('Refining was problematic')
        end
        flag_smoothness = smoothness_condition(x0, x_dot_0, x1,x_dot_1, Imin);
        if ~flag_smoothness
            error('At iteration %i, the curve is not guaranteed to be smooth\n',i)
        end
    end
    
    % check if too many or too little nodes
    %
    % introduce 2 counters: too many, too little
    % if was already too many, new solution smaller
    % if too many, last mode(s) of new step imposed zero
    % if too little, add mode to new step asap
    % check for "ideal" number of nodes
    first_big_mode = find(log10(sum(abs(x1.vector),1))>-14.5,1);
    ideal_small_modes = max(5,0.3*x0.nodes);
    
    if first_big_mode <= 0.9*ideal_small_modes
        add_node = floor(-first_big_mode+ideal_small_modes);
        decrease_node = 0;
    elseif first_big_mode>1.1*ideal_small_modes
        decrease_node = ceil(-first_big_mode+ideal_small_modes);
        add_node = 0;
    else
        add_node =0;
        decrease_node = 0;
    end
    
    % additional check: if previous_iter.Z1-1 too small, add some nodes
    if ~isempty(previous_iter)
        %&& 1-previous_iter.Z1_norm-max(Z0vector)<0.2
        
        %previous_iter=struct('Y',new_iter_Y,'Z1',new_iter_Z1,'Z1_extrema',Z1s(:,end),...
        %'Z0_extrema',Z0s(:,end),'Y_extrema',Ys(:,end),'Z2',Z2);
        delta_small = 0.1;
        delta_big = 0.2;
        test_modes = @(delta) 1-previous_iter.Z0_extrema-previous_iter.Z1_extrema-...
            sqrt(delta^2 + 4*previous_iter.Y_extrema.*previous_iter.Z2);
        if any(test_modes(delta_small)<0)
            %if nu^x1.nodes>1.2
            %    warning('Consider decreasing nu')
            %end
            if add_node<=2
                decrease_node = 0;
                add_node = 2;
            end
        end
        if any(test_modes(delta_big)<0)
            decrease_node = 0;
        end
    end
    
    
    if add_node
        if talkative >0
            fprintf('At iteration %i, the number of nodes is increased by %i. New number of modes %i\n', i, add_node, x0.nodes+add_node);
        end
        x1 = reshape_Xi(x1, x1.nodes+add_node);
        x_dot_1 = Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot_1,x0), ...
            x0.nodes+add_node));
        F_new.scalar_equations = reshape(F_new.scalar_equations, x0.nodes+add_node);
        F_not_square.scalar_equations = reshape(F_not_square.scalar_equations, x0.nodes+add_node);
        DF1 = derivative_to_matrix(derivative(F_new,x1,0));
        Anew = inv(DF1);
        previous_iter =[];
        new_step = new_step*0.9;
    end
    if decrease_node && ~bool_Hopf
        if talkative >0
            fprintf('At iteration %i, the number of nodes is decreased by %i. New number of modes %i\n', i, abs(decrease_node), x0.nodes+decrease_node);
            fprintf('\n Leap iteration to decrease the number of nodes \n')
        end
        % necessary temporary elements for "one point validation"
        x1_big = x1;
        Anew_big = Anew;
        DF1_big = DF1;
        F_new_big = F_new;
        x1 = reshape_Xi(x1, x1.nodes+decrease_node);
        x_dot_1 = Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot_1,x0),...
            x0.nodes+decrease_node));
        F_new.scalar_equations = reshape(F_new.scalar_equations, x0.nodes+decrease_node);
        F_not_square.scalar_equations = reshape(F_not_square.scalar_equations, x0.nodes+decrease_node);
        % following node with less modes
        
        % one point validation
        x1_small = reshape_Xi(x1, x0.nodes);
        F_new_small = F_new;
        F_not_square_small = F_not_square ;
        F_new_small.scalar_equations = reshape(F_new.scalar_equations, x0.nodes);
        F_not_square_small.scalar_equations = reshape(F_not_square.scalar_equations, x0.nodes);
        DF1_small = derivative_to_matrix(derivative(F_new_small,x1_small,0));
        Anew_small = inv(DF1_small);
        
        use_intlab = temp_use_intval;
        % validate
        [flag] = radii_polynomials_cont_new(x1_big,x1_small,DF1_big,DF1_small,...
            F_new_big,F_new_small,[],Anew_big,Anew_small);
        
        use_intlab = 0;
        if flag==0
            error('When decreasing number of modes something went wrong')
        end
        
        DF1 = derivative_to_matrix(derivative(F_new,x1,0));
        Anew = inv(DF1);
        previous_iter = [];
    end
    
    
    if talkative>0
        % some talking to the user, when requested
        fprintf('\n ---- SUCCESS ----\n ')
        if talkative>1
            fprintf('The %i-th interval of validation is [%e, %e]\n',i,Imin,Imax);
        end
        if talkative>0 && mod(i,5) ==0
            fprintf('The %i-th step has been validated successfully\n',i);
            if talkative>1
                fprintf('The stepsize is %e\n',h)
            end
        end
        fprintf('\n')
    end
    
    %storage
    step_size(i) = h;
    number_modes(i+1) = x1.nodes;
    norm_x(:,i)    = vert(norm(x1));
    norm_x(1:x0.size_scalar,i) = x1.scalar;% to avoid absolute values here
    Interval(:,i)  = [Imin,Imax]';
    Z0_iter(:,i)   = vert(Z0vector);
    Z1_iter(:,i)   = vert(Z1vector);
    Z2_iter(:,i)   = vert(Z2vector);
    Y_iter(:,i)    = vert(Yvector);
    
    % iteration
    x_dot_0 = x_dot_1;
    Aold = Anew;
    F_old = F_new;
    h = h* new_step;
    DF0 = DF1;
    x0 = x1;
    
    if mod(i,250) == 0
        x_n = x1;
        x_dot_n = x_dot_1;
        s_temp = strcat(s,'temp');
        if ~isempty(saddle_index_var)
            save(s_temp, 'x_n', 'x_dot_n', 'step_size','number_modes', 'norm_x', 'Interval', 'Y_iter', 'Z0_iter', 'Z1_iter', 'Z2_iter', ...
                'saddle_x0_stored','saddle_x1_stored','saddle_index_var','saddle_iter');
        else
            save(s_temp, 'x_n', 'x_dot_n','step_size','number_modes', 'norm_x', 'Interval', 'Y_iter', 'Z0_iter', 'Z1_iter', 'Z2_iter');
        end
    end
end

x_n = x1;
x_dot_n = x_dot_1;
F_new_save = F_new;
F_new_save.scalar_equations = remove_lin_coef(F_new.scalar_equations,...
    F_new.scalar_equations.number_equations_lin);
F = F_new_save;

% save to location the elements that are interesting
if ~isempty(saddle_index_var)
    save(s, 'x_0','x_n', 'x_dot_n', 'step_size','number_modes', 'norm_x', 'Interval', 'Y_iter', 'Z0_iter', 'Z1_iter', 'Z2_iter', ...
        'saddle_x0_stored','saddle_x1_stored','saddle_index_var','saddle_iter','F');
else
    save(s, 'x_0','x_n', 'x_dot_n','step_size','number_modes', 'norm_x', 'Interval', 'Y_iter', 'Z0_iter', 'Z1_iter', 'Z2_iter','F');
end

use_intlab = temp_use_intval;

return
