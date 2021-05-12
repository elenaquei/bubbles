function [flag,Imin,x_end,Imax] = substepping(x_0,x_1,F0,F1,times_run, bool_refine_anyway)
% function [flag, Imin, x_end,Imax] = substepping(x_0,x_1,F0,F1,times_run,bool_refine_anyway)

global talkative
global use_intlab

if nargin <5 || isempty(times_run)
    times_run = 0;
end

if nargin<6 || isempty(bool_refine_anyway)
    bool_refine_anyway = 0;
end


% %%% debugging
% flag = 1;
% Imin=10^-6;
% return

temp_use_intval = use_intlab;

step_size = max(norm(x_0-x_1));

use_intlab = temp_use_intval;
% validate
if bool_refine_anyway
    flag = 0;
else
    [flag,Imin,Imax] = radii_polynomials_cont_new(x_0,x_1,[],[],...
        F0,F1);
end

if flag 
    return
    %flag = 1;
    %warning('DEBUGGING TIME')
    %error('Continuation from the Hopf to the standard continuation failed - first step')
end
use_intlab = 0;

n_steps = measure_of_ugliness(F0,F1,step_size);

dbstack
fprintf('%i extra steps needed for validation\n',n_steps)
talkative = talkative - 1;

x_iminus1 = x_0;
x_Delta = 1/n_steps * (x_1 - x_0);

coefs_0 = F0.scalar_equations.linear_coef;
coefs_1 = F1.scalar_equations.linear_coef;

coefs_iminus1 = coefs_0;
coefs_Delta = divide_cell(minus_cell(coefs_1,coefs_0),n_steps);

F_iminus1 = F0;
F_i = F_iminus1;
Imin = 0;
Imax = Inf;

for i = 1:n_steps
    x_i = x_iminus1 + x_Delta;
    coefs_i = add_cell(coefs_iminus1, coefs_Delta);
    
    
    F_i.scalar_equations.linear_coef= coefs_i;
    
    x_i = Newton_2(x_i,F_i);
    
    %[flag] = radii_polynomials_cont_new(x_iminus1,x_i,[],[],...
    %    F_iminus1, F_i);
    if i ==1 
        DH_i = derivative_to_matrix(derivative(F_i,x_i,0));
        A_i = inv(DH_i);
        use_intlab = temp_use_intval;
        [flag,Imin_loc,Imax_loc,info_iter] = ...
            radii_polynomials_cont_new(x_iminus1,x_i,[],DH_i,...
            F_iminus1, F_i,[],[],A_i);
        use_intlab = 0;
    else
        DH_i = derivative_to_matrix(derivative(F_i,x_i,0));
        A_i = inv(DH_i);
        use_intlab = temp_use_intval;
        [flag,Imin_loc,Imax_loc,info_iter] = ...
            radii_polynomials_cont_new(x_iminus1,x_i,DH_iminus1,DH_i,...
            F_iminus1, F_i,info_iter,A_iminus1,A_i);
        use_intlab = 0;
    end
    
    if ~flag
        fprintf('%i iteration out of %i FAILED (times run = %i)\n',i,n_steps,times_run)
        if times_run<40
            [flag,Imin_loc, x_i,Imax_loc] = substepping(x_iminus1,x_i,F_iminus1,F_i,times_run+1);
        else
            flag = 0;
            x_end = [];
            Imin =[];
            Imax =[];
            return
        end
    elseif mod(i,5) == 0 && talkative >1
        fprintf('%i iteration out of %i successful\n',i,n_steps)
    elseif mod(i,10) == 0 && talkative >0
        fprintf('%i iteration out of %i successful\n',i,n_steps)
    elseif mod(i,20) == 0
        fprintf('%i iteration out of %i successful at time %s\n',i,n_steps, datestr(now,'HH:MM:SS.FFF'))
    end
    
    x_iminus1 = x_i;
    coefs_iminus1 = coefs_i;
    F_iminus1 = F_i;
    DH_iminus1 = DH_i;
    A_iminus1 = A_i;
    Imin = max(Imin, Imin_loc);
    Imax = min(Imax, Imax_loc);
end

x_end = x_i;

if times_run == 0
    % check that it glues back together
    if norm(x_end - x_1) > Imin_loc
        [flag,Imin_end]=radii_polynomials(x_end,F_i,DH_i, A_i);
        if norm(x_end - x_1) > Imin_end
            fprintf('Norm of the distance between validated point and proposed point is %f',norm(x_end - x_1))
            error('Could not glue back together')
        end
    end
end

talkative = talkative + 1;
end

function cell_Delta = divide_cell(cell_3,N)

cell_Delta = cell(3,1);
for i = 1:3
    cell_Delta{i} = cell_3{i}/N;
end

end


function cell_Delta = add_cell(cell_3,cell_3_2)

cell_Delta = cell(3,1);
for i = 1:3
    cell_Delta{i} = cell_3{i}+cell_3_2{i};
end

end

function cell_Delta = minus_cell(cell_3,cell_3_2)

cell_Delta = cell(3,1);
for i = 1:3
    cell_Delta{i} = cell_3{i}-cell_3_2{i};
end

end

function n_step = measure_of_ugliness(F0,F1,~)
% function n_step = measure_of_ugliness(F0,F1,step_size)
% a very heuristic way to determine how many substeps we would need to
% avoid failure
% 
% C0 = F0.scalar_equations.linear_coef;
% C1 = F1.scalar_equations.linear_coef;
% n_step = 0*C0{3};
% for i = 1:size(C0{3})
%     c0 = Xi_vector(C0{1}(i,:),squeeze(C0{2}(i,:,:)));
%     c1 = Xi_vector(C1{1}(i,:),squeeze(C1{2}(i,:,:)));
% 
%     n_step(i)=222;%2*max(ceil(norm(c0-c1)/(20*step_size)));
% end
% 
% n_step = max(max(n_step),150);

% n_step = 222; % WORKS FOR HOPF GLUING LORENZE84

n_step = 50;

end
