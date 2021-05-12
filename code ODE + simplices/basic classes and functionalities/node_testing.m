function [best_n_nodes,x1,x_dot1,DH1,Imin,Imax,Y0,Z0,Z1,Z2]=node_testing(...
    x0_in,x_dot0_in,coefs_linear0_in, coefs_linear_in, alpha_coef,delta_step_in,...
    max_Newton_iter,res_Newton)
%function [best_n_nodes,x1,x_dot1,DH1,Imin,Imax,Y0,Z0,Z1,Z2]=node_testing(...
%    x0_in,x_dot0_in,coefs_linear0_in, coefs_linear_in, alpha_coef, delta_step_in...
%    max_Newton_iter,res_Newton)
%
% Testing the number of nodes in [n_nodes-2:2:n_nodes+2] and
% finding the best number of nodes. The best number of nodes, the next 
% validated solution with this number of nodes and the
% validation elements of this segment are returned
%
% INPUTS
% x0_in                Xi_vector, input solution
% x_dot0_in            vector, local tangent to the curve at x0_in
% coefs_linear0_in     cell, coefficients of the linear part associated to the x0_in solution
% coefs_linear_in      cell, coefficients of the linear part of the general problem
% alpha_coef           coefs, coefficient of the vector field
% delta_step_in        assumed step size for the next validation step
% max_Newton_iter      maximum number of iterations in Newton
% res_Newton           requested residual in Newton
%
%
% OUTPUT
% best_n_nodes         number of nodes for which the validation is faster
% x1                   validated next solution point with best_n_nodes
% x_dot1               local tangent of the curve in x1
% DH1                  derivative in x1
% Imin                 minimum validation radius for [x0_in, x1]
% Imax                 maximum validation radius for [x0_in, x1]
% Y0,Z0,Z1,Z2          bounds used in the radii polynomials in the
%                      validation of [x0_in,x1]


global use_intlab

temp_intlab=use_intlab;
n_nodes=x0_in.nodes;

% LONGER TEST: best_n_nodes = n_nodes+2
test_n_nodes = n_nodes+2;
previous_iter.Y=[];
previous_iter.Z1=[];
%n_nodes=test_n_nodes;
x0=reshape_Xi(x0_in,test_n_nodes);
coefs_linear=reshape_coefs(coefs_linear_in,test_n_nodes);

%x_dot0=Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot0_in,x0_in.size_scal,x0_in.size_vec,nodes),n_nodes));
x_dot0=Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot0_in,x0_in.size_scalar,x0_in.size_vector,...
    x0_in.nodes),x0.nodes));

% add one equation for archlength continuation
coefs_linear0=coefs_linear;
coefs_linear0{1}=[-x_dot0(1:2).';[0,0]];
coefs_linear0{2}(2,:,:)=coefs_linear{2}(1,:,:);
coefs_linear0{2}(1,:,:)=reshape(-x_dot0(3:end),2*x0.nodes+1,x0.size_vector).';
coefs_linear0{3}(1) = Xi_vec2vec(x0).'*x_dot0;
coefs_linear0{3}(2)=0;
coefs_linear0{3}(2)=real(-coefs_linear0{1}(2,:)*x0.scalar.'-...
    sum(sum(squeeze(coefs_linear0{2}(2,:,:)).*x0.vector,1)));

DH0=Function_derivative(x0,alpha_coef,coefs_linear0);

use_intlab=0;

[x1_longer,x_dot1_longer,DH1_longer,coefs_linear1] = compute_new_point(x0,delta_step_in,x_dot0,...
    alpha_coef,coefs_linear,max_Newton_iter,res_Newton);

use_intlab = temp_intlab;

% validate (x_iter-1 , x_iter)
tic;

[flag,Imin_longer,Imax_longer,~,Y0_longer,Z0_longer,Z1_longer,Z2_longer,...
    new_step_longer] = radii_polynomials_cont(x0,x1_longer,DH0,DH1_longer,...
    alpha_coef,coefs_linear0, coefs_linear1,previous_iter);
if flag>0
    delta_step_longer=delta_step_in*new_step_longer;
    speed_longer=delta_step_longer/toc;
else
    speed_longer=0;
end




%% SECOND POSSIBILITY : beat_n_nodes = n_nodes
x0=x0_in;
x_dot0=x_dot0_in;
use_intlab=0;
[x1_std,x_dot1_std,DH1_std,coefs_linear1] = compute_new_point(x0,delta_step_in,x_dot0,...
    alpha_coef,coefs_linear_in,max_Newton_iter,res_Newton);
DH0=Function_derivative(x0_in,alpha_coef,coefs_linear0_in);
use_intlab = temp_intlab;
% validate (x_iter-1 , x_iter)
tic;

[flag,Imin_std,Imax_std,~,Y0_std,Z0_std,Z1_std,Z2_std,new_step_std] = ...
    radii_polynomials_cont(x0,x1_std,DH0,DH1_std,alpha_coef,...
    coefs_linear0_in, coefs_linear1,previous_iter);
if flag>0
    delta_step_std=delta_step_in*new_step_std;
    speed_standard=delta_step_std/toc;
else
    speed_standard=0;
end

%% THIRD POSSIBILITY: best_n_nodes = n_nodes - 2
test_n_nodes = n_nodes - 2;
use_intlab=0;

% BIG dimension
x0_BIG=x0_in;
% x_dot0_BIG=x_dot0_in;
DH0_BIG=Function_derivative(x0_BIG,alpha_coef,coefs_linear0_in);
coefs_linear0_BIG=coefs_linear0_in;

% SMALL dimension
x0_shorter=reshape_Xi(x0_BIG,test_n_nodes);
x_dot0_shorter=Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot0_in,x0_BIG.size_scalar,x0_BIG.size_vector,...
    x0_BIG.nodes),x0_shorter.nodes));
%coefs_linear0_shorter=reshape_coefs(coefs_linear0_BIG,best_n_nodes);

% compute new step with less nodes
coefs_linear_shorter=reshape_coefs(coefs_linear_in,test_n_nodes);
[x1_shorter,x_dot1_shorter,DH1_shorter,coefs_linear1_shorter] = compute_new_point(x0_shorter,delta_step_in,x_dot0_shorter,...
    alpha_coef,coefs_linear_shorter,max_Newton_iter,res_Newton);

% padd with zeros
x1_BIG=reshape_Xi(x1_shorter,n_nodes);
coefs_linear1_BIG=reshape_coefs(coefs_linear1_shorter,n_nodes);
DH1_BIG =Function_derivative(x1_BIG,alpha_coef,coefs_linear1_BIG);

use_intlab=temp_intlab;

previous_iter.Y=[];
previous_iter.Z1=[];
[flag,Imin_shorter,Imax_shorter,~,Y0_shorter,Z0_shorter,Z1_shorter,Z2_shorter,new_step_shorter] = ...
    radii_polynomials_cont(x0_BIG,x1_BIG,DH0_BIG,DH1_BIG,alpha_coef,...
    coefs_linear0_BIG, coefs_linear1_BIG,previous_iter);

if flag>0
    delta_step_shorter=delta_step_in*new_step_shorter;
    speed_shorter=delta_step_shorter/toc;
else
    speed_shorter=0;
end


%% COMPARISON
if speed_shorter>speed_standard && speed_shorter>speed_longer
    best_n_nodes=n_nodes-2;
    
    fprintf('NUMBER OF NODES DECREASED TO %d\n\n',best_n_nodes);
    
    Imin=Imin_shorter;
    Imax=Imax_shorter;
    Y0=Y0_shorter;Z0=Z0_shorter;Z1=Z1_shorter;Z2=Z2_shorter;
    x1=x1_shorter;x_dot1=x_dot1_shorter;DH1=DH1_shorter;
    
elseif speed_longer>speed_standard && speed_longer>speed_shorter
    best_n_nodes=n_nodes+2;
    fprintf('NUMBER OF NODES INCREASED TO %d\n\n',best_n_nodes);
    
    Imin=Imin_longer;
    Imax=Imax_longer;
    Y0=Y0_longer;Z0=Z0_longer;Z1=Z1_longer;Z2=Z2_longer;
    x1=x1_longer;x_dot1=x_dot1_longer;DH1=DH1_longer;

elseif speed_standard>0
    best_n_nodes=n_nodes;

    Imin=Imin_std;
    Imax=Imax_std;
    Y0=Y0_std;Z0=Z0_std;Z1=Z1_std;Z2=Z2_std;
    x1=x1_std;x_dot1=x_dot1_std;DH1=DH1_std;
else
    error('something went wrong: it was impossible to validate the segment with a variety of numbers of nodes');

end

return
end
