function [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,simplex, list_of_nodes] = ...
    radii_polynomials_simplex(simplex, list_of_nodes,bool_Hopf)%xBar0,xBar1,xBar2,...
%alpha0,alpha1,alpha2,previous_iter0,previous_iter1)
% INPUT:
% simplex   instance of the class simplex - validation material
% list_of_nodes     cell of instances of the node class
%
% OUTPUT:
% flag      0  - program failed
%           1  - all good
%           2  - really good (consider decreasing delta)
% Imin      positive real value, left bound of the interval;
% Imax      positive real value, right bound of the interval;
% Yvector   value of the Y bound
% Z0vector  value of the Z0 bound
% Z1vector  value of the Z1 bound
% Z2vector  value of the Z2 bound
% simplex   updated with the validation stepsize
% list_of_nodes     updated with the validation results
%
% In this function, the radii polynomial refering to the given problem and
% solutions are computed.

global use_intlab
global Display
global talkative
global debug
% ERROR Spotting

if length(list_of_nodes) == 3
    indeces = [1,2,3];
else
    indeces = simplex.nodes_number;
end

xBar0 = list_of_nodes{indeces(1)}.solution;
xBar1 = list_of_nodes{indeces(2)}.solution;
xBar2 = list_of_nodes{indeces(3)}.solution;

    function alpha = square_it(alpha,x, old_alpha)
        if ~square(alpha)
            if bool_Hopf
                % adds two new scalar equations
                alpha.scalar_equations = F_update_Hopf(alpha.scalar_equations,x);
            else
                alpha.scalar_equations = fancy_scalar_condition(x, alpha.scalar_equations);
            end
            if nargin<3
                alpha = continuation_equation_simplex(alpha,x);
            else
                alpha = continuation_equation_simplex(alpha,x, old_alpha);
            end
        end
    end

alpha0 = list_of_nodes{indeces(1)}.problem;
alpha1 = list_of_nodes{indeces(2)}.problem;
alpha2 = list_of_nodes{indeces(3)}.problem;
temp_intlab = use_intlab;
use_intlab = 0;
alpha0 = square_it(alpha0,xBar0);
alpha1 = square_it(alpha1,xBar1, alpha0);
alpha2 = square_it(alpha2,xBar2, alpha0);
use_intlab = temp_intlab;

previous_iter0 = list_of_nodes{indeces(1)}.previous_validation;
previous_iter1 = list_of_nodes{indeces(2)}.previous_validation;
previous_iter2 = list_of_nodes{indeces(3)}.previous_validation;

if xBar0.size_scalar~=xBar1.size_scalar
    error('Scalar dimensions are not consistent.')
end
if xBar0.size_vector~=xBar1.size_vector
    error('Vector dimensions are not consistent.')
end
if  ~compatible_vec(alpha0,xBar0) || ~compatible_vec(alpha1,xBar0) || ~compatible_vec(alpha1,xBar1)
    error('Problems and solutions not compatible.')
end
if ~square(alpha1)
    error('Problem is not square, a unique solution does not exist.')
end

% default return in case of crash
flag=0;
Imin=[];
Imax=[];
Yvector=[];Z0vector=[];Z1vector=[];Z2vector=[];
new_iter=struct('Y',[],'Z1',[],'Z_norm',[]);
new_step=0.5; % in case of failure, start here

    function [A_struct, DF_small] = create_A_struct(alpha, x)
        [A_small, M, P, Q, R, phi, D3F2, DF_small] = A_delay_symplex(alpha, x);
        if use_intlab
            [A_small, M, P, Q, R, phi, D3F2, DF_small] = everything_intval(A_small, M, P, Q, R, phi, D3F2, DF_small);
        end
        A_struct.A = A_small;
        A_struct.M = M;
        A_struct.P = P;
        A_struct.Q = Q;
        A_struct.R = R;
        A_struct.phi = phi;
        A_struct.D3F2 = D3F2;
    end

% set up of the derivatives
if has_delay(alpha0)
    
    [A0_struct, DH0] = create_A_struct(alpha0, xBar0);
    [A1_struct, DH1] = create_A_struct(alpha1, xBar1);
    [A2_struct, DH2] = create_A_struct(alpha2, xBar2);
    
else
    DH0=derivative_to_matrix(derivative(alpha0,xBar0,0));
    DH1=derivative_to_matrix(derivative(alpha1,xBar1,0));
    DH2=derivative_to_matrix(derivative(alpha2,xBar2,0));
    
    A0 = inv(DH0);
    A1 = inv(DH1);
    A2 = inv(DH2);
    
    % symmetrise A0 and A1
    A0=symmetrise_A(A0,xBar0);
    A1=symmetrise_A(A1,xBar0);
    A2=symmetrise_A(A2,xBar0);
end

Adagger_delta1=DH1-DH0;
Adagger_delta2=DH2-DH0;

% change into intvals
if use_intlab
    [Adagger_delta1, Adagger_delta2]=everything_intval(Adagger_delta1, Adagger_delta2);
    %nu=intval(nu);
    [xBar0,xBar1,xBar2]=everything_intval(xBar0,xBar1,xBar2);
    for i=1:xBar0.size_vector
        alpha0.vector_field.value{i}=intval(alpha0.vector_field.value{i});
        alpha1.vector_field.value{i}=intval(alpha1.vector_field.value{i});
        alpha2.vector_field.value{i}=intval(alpha2.vector_field.value{i});
    end
    for i = 1:3
        alpha0.scalar_equations.linear_coef{i} = intval(alpha0.scalar_equations.linear_coef{i});
        alpha1.scalar_equations.linear_coef{i} = intval(alpha1.scalar_equations.linear_coef{i});
        alpha2.scalar_equations.linear_coef{i} = intval(alpha2.scalar_equations.linear_coef{i});
    end
    for i = 1:alpha0.scalar_equations.number_equations_pol
        alpha0.scalar_equations.polynomial_equations.value{i}=intval(alpha0.scalar_equations.polynomial_equations.value{i});
        alpha1.scalar_equations.polynomial_equations.value{i}=intval(alpha1.scalar_equations.polynomial_equations.value{i});
        alpha2.scalar_equations.polynomial_equations.value{i}=intval(alpha2.scalar_equations.polynomial_equations.value{i});
    end
end

if talkative>1
    fprintf('\nComputed validation set up, time %s\n',datestr(now,13));
end

% Z0 BOUND
if has_delay(alpha0)
    Z0vector = Z0_delay_simplex(alpha0, alpha1, alpha2, xBar0,xBar1,xBar2, A0_struct, A1_struct, A2_struct, DH0, DH1, DH2);
else
    [Z0vector,Z0s] = Z0_bound_simplex(DH0,DH1,DH2,A0,A1,A2,xBar0);
end

if talkative>1
    fprintf('\nComputed Z0, time %s\n',datestr(now,13));
end
if any(Z0vector>1)
    return
end

% Z1 BOUND
if has_delay(alpha0)
    Z1vector = Z1_delay_simplex(alpha0, alpha1, alpha2, xBar0,xBar1,xBar2, A0_struct, A1_struct, A2_struct);
else
    if ~isempty(previous_iter0) && ~isempty(previous_iter0.Z1) &&...
            ~isempty(previous_iter1) && ~isempty(previous_iter1.Z1)
        [Z1vector,new_iter_Z1,Z1s]=Z1_bound_simplex(A0,A1,A1,xBar0,xBar1,xBar1,alpha0,...
            alpha1,alpha1,Adagger_delta1,Adagger_delta2,previous_iter0.Z1,previous_iter1.Z1);
    else
        [Z1vector,new_iter_Z1,Z1s]=Z1_bound_simplex(A0,A1,A1,xBar0,xBar1,xBar2,alpha0,...
            alpha1,alpha2,Adagger_delta1,Adagger_delta2);
    end
end

if talkative>1
    fprintf('\nComputed Z1, time %s\n',datestr(now,13));
end

if any(Z1vector>1)
    return
end

% Z2 BOUND
if has_delay(alpha0)
    Z2vector = Z2_delay_simplex(alpha0, alpha1, alpha2,xBar0,xBar1,xBar2, A0_struct, A1_struct, A2_struct);
else
    [Z2vector,Z2s]=Z2_bound_simplex(A0,A1,A2,xBar0,xBar1,xBar2,alpha0);
end


if talkative>1
    fprintf('\nComputed Z2, time %s\n',datestr(now,13));
end


% Y BOUND
if has_delay(alpha0)
    % Yvector = Y_delay_simplex(alpha0, alpha1, alpha2, xBar0,xBar1,xBar2, A0_struct, A1_struct, A2_struct);
    upper_bound = (intval(Z1vector)+intval(Z0vector)-1).^2./(4*Z2vector);
    Yvector = Y_delay_refinement(alpha0, alpha1, alpha2, xBar0,xBar1,xBar2, A0_struct, A1_struct, A2_struct, upper_bound);
else
    if ~isempty(previous_iter0) && ~isempty(previous_iter0.Y) && ...
            ~isempty(previous_iter1) && ~isempty(previous_iter1.Y)
        [Yvector,new_iter_Y,Ys]=Y_bound_simplex(A0,A1,A2,xBar0,xBar1,xBar2,...
            alpha0,alpha1,alpha2,previous_iter0.Y,previous_iter1.Y);
    else
        [Yvector,new_iter_Y,Ys]=Y_bound_simplex(A0,A1,A2,xBar0,xBar1,xBar2,alpha0,alpha1,alpha2);
    end
end

if talkative>1
    fprintf('\nComputed Y, time %s\n',datestr(now,13));
end

if any(Yvector>1)
    return
end

if ~has_delay(alpha0)
    new_iter=struct('Y',new_iter_Y,'Z1',new_iter_Z1,'Z1_extrema',Z1s(:,end),...
        'Z0_extrema',Z0s(:,end),'Y_extrema',Ys(:,end),'Z2',Z2vector);
    list_of_nodes{indeces(3)}.previous_validation = new_iter;
end


% computation of the radius
try
    [bool,Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector);
catch
    return
end
if talkative>1
    fprintf('\n Computed validation interval, time %s\n\n',datestr(now,13));
end
if ~bool
    return
end

if any(size(Imin)<1) || min(Imin,Imax)<0
    return
end

% adapttive stepsize
% new_step = ideal_stepsize(Imin, Imax, Ys, Z0s, Z1s, Z2s);

% simplex.verification_coeff = new_step;

% plot
if Display
    p=@(r) Z0vector*r+Z1vector*r+Z2vector*r.^2-ones(size(Z0vector))*r+Yvector*ones(size(r));
    xplot=0:0.000001:Imax*1.5;
    figure;
    pplot=p(xplot);
    plot(xplot,pplot,xplot,0*xplot,'k');
    hold on
    plot(Imin,0,'k*',Imax,0,'k*')
end

% flag
if Imax-Imin>1e-4
    flag=2;
else
    flag=1;
end

return
end





function A_sim=symmetrise_A(A,x)
A_sim=A;

length_nodes=2*x.nodes+1;
sym_vec=@(v) ( v + conj(v(end:-1:1,:)))/2;
sym_mat=@(A) ( A + conj(A(end:-1:1,end:-1:1)))/2;

for j=1:x.size_vector
    vec=x.size_scalar+(j-1)*length_nodes+1:x.size_scalar+j*length_nodes;
    all_scal=1:x.size_scalar;
    A_sim(vec,all_scal)=sym_vec(A_sim(vec,all_scal));
end

for jj=1:x.size_vector
    for kk=1:x.size_vector
        vec1=x.size_scalar+(jj-1)*length_nodes+1:x.size_scalar+jj*length_nodes;
        vec2=x.size_scalar+(kk-1)*length_nodes+1:x.size_scalar+kk*length_nodes;
        A_sim(vec1,vec2)=sym_mat(A_sim(vec1,vec2));
    end
end

end