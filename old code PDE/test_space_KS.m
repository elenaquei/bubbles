function test_space_KS
global nu
global display
global nu_Hopf
global a_Hopf
display = 2;
nu = [1.1];

% new ordering of parameters: nu, lambda, a 

try
    intval(1);
catch
    addpath(genpath('/Users/queirolo/Dropbox/Hopf_Bifurcations/2D_examples/Intlab_V7.1'));
    startintlab
end

nu_Hopf = 0.132;

node_space = 10;
node_time = 3;

% results from JP
approx1D = [0
    0
   1.482795416033047
   0
  -0.453719432314703
   0
   0.058908791659034
   0
  -0.006357461422265
   0
   0.000592702920728].';

approx1D = -1i*approx1D(1:node_space+1); % if bigger change
approx = [conj(approx1D(end:-1:2)),approx1D];


vector1D = approx;
vector1D(node_space+1) = 0;
%x_old = 1i*[0,0,0,-3,-2,-1,0,1,2,3,0,0,0];
pyz0 = approx;
% zz = [0.000000000000000 + 0.006357461422265i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 - 0.058908791659034i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 + 0.453719432314703i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 - 1.482795416033047i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 + 1.482795416033047i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 - 0.453719432314703i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 + 0.058908791659034i
%   0.000000000000000 + 0.000000000000000i
%   0.000000000000000 - 0.006357461422265i].';

G2 = zeros(1,1+2*node_space);
K = -node_space:node_space;
G2(1:(2*node_space+1)) = -1i*pyz0.*K;

G = zeros(1, length(G2));
G(1,:) = G2;
const_G= [0]';

y = KS_loc(pyz0, G, const_G);
%y_test = KS(pyz0, G, const_G);

nu_Hopf = 0.131; %pos

x2 = Newton(pyz0, @(x) KS_loc(x,G,const_G), @(x) DKS_loc(x,G), @(x) fin_diff_KS(x,G));

% [l2,v2] = eigs(DKS(x2(:).',G));
% 
% E2_ = eig(DKS(x2(:).',G));

[V2, D2]= eig(DKS_loc(x2(:).',G));

nu_Hopf = 0.132; %neg
nu_Hopf_2 = nu_Hopf;

x2 = Newton(x2, @(x) KS_loc(x,G,const_G), @(x) DKS_loc(x,G), @(x) fin_diff_KS(x,G));

[V2, D2]= eig(DKS_loc(x2(:).',0*G));


z_hat = zeros(2*node_time+1, 2*node_space+1);

% D2, V2, x2 chosen solutions for nu_Hopf_2
% with DKS(x2(:).'

n = 4;
[temp_V, index] = sort(real(diag(D2)));
% small_eigns = V2(sort(index(end:-1:end-n)));

index_eig = index(end-2);

z_hat(node_time, : ) = V2(:,index_eig);
z_hat(node_time+2, : ) = conj(V2(end:-1:1,index_eig));
important_eig = D2(index_eig,index_eig);
if node_space~=10
    error('chosen column of V2 computed just for 10 spatial nodes')
end
z_hat = z_hat/sqrt(sum(sum(z_hat.*conj(z_hat))));

y_hat = x2;
a_hat = 10^-4;

y = x2;
z = z_hat;
nu_Hopf = nu_Hopf_2;
lambda = abs(imag(important_eig))^-1;
a = a_hat;

K = -node_space:node_space;
J = -node_time:node_time;

K2 = repmat(K,2*node_time+1,1);
J2 = repmat(J(:), 1, 2*node_space+1,1);

i_J_z_hat = 1i*J2.*z_hat(:,end:-1:1);
G1 = [0,0,0*y_hat(:).',i_J_z_hat(:).'];
const_G1 =0;
G3 = [0,0,G(:).',0*z_hat(:).'];
const_G3 = 0;
G2 = [0,0,0*y_hat(:).',conj(z_hat(:)).'];
const_G2 = -1;
y_2D = embedding(y,z);
y_plus_az = y_2D+a_hat*z_hat;
G4 = [0,0, 0*y_hat(:).', 1i*K2(:).'.*y_plus_az(:).'];
const_G4 = 0;

G_Hopf= [G1;G2;G3;G4];
const_G_Hopf = [const_G1,const_G2,const_G3,const_G4];

pyz0 = [lambda,a,y(:).',z(:).'];

Df = fin_diff_KS_Hopf(pyz0,G_Hopf, node_space);


% % section to comment: singular problem, not well defined over the 
% % amplitude
if 1 == 2%try
    if rcond(Df)<10^-6
        warning('Almost singular derivative, %e', rcond(Df))
    end
    
    KS_Hopf_vec_loc = @(x) KS_Hopf_vec(x,G_Hopf, const_G_Hopf, node_space);
    der_KS_Hopf = @(x) fin_diff_KS_Hopf(x,G_Hopf, node_space);
    
    
    pyz1 = Newton(pyz0, @(x) KS_Hopf_vec_loc(x),...
        @(x) der_KS_Hopf(x), @(x) der_KS_Hopf(x), @(x) norm_nu_vec(x,node_space), 10);
    
    pyz0_vector =PDE_vector([lambda, a], y, z.');
    pyz1_vector = PDE_vector(pyz1, 2, 1, node_space);
    
    pyz0_vector =PDE_vector([nu_Hopf,lambda, a], y.', shiftdim(z,-1));
    
    G_Hopf = cat(2,zeros(4,1),G_Hopf);
    y = evaluate_HopfKS(pyz0_vector,G_Hopf, const_G_Hopf,0*G_Hopf(1,:).',0);

end



% shift fixed parameter from nu to a

a_Hopf = a;
pyz0 = [nu_Hopf,lambda,y(:).',z(:).'];

Df = fin_diff_KS_Hopf_a(pyz0,G_Hopf, node_space);
if rcond(Df)<10^-7 % actually gets $worst$
    [v,d] = eig(Df);
    [~,index]= min(abs(diag(d)));
    G2 = v(:,index).';
    const_G2 = -sum(G2.*pyz0);
    G_Hopf= [G1;G2;G3;G4];
    const_G_Hopf = [const_G1,const_G2,const_G3,const_G4];
    warning('Almost singular derivative, %e', rcond(Df))
end



KS_Hopf_vec_loc = @(x) KS_Hopf_vec_a(x,G_Hopf, const_G_Hopf, node_space);
der_KS_Hopf = @(x) fin_diff_KS_Hopf_a(x,G_Hopf, node_space);

pyz1 = Newton(pyz0, @(x) KS_Hopf_vec_loc(x),...
    @(x) der_KS_Hopf(x), @(x) der_KS_Hopf(x), @(x) norm_nu_vec(x,node_space), 10);

long_pyz = [pyz1(1:2);a_Hopf;pyz1(3:end)];
long_G = [G_Hopf(:,1:2),zeros(size(G_Hopf,1),1),G_Hopf(:,3:end)];
Df = fin_diff_KS_Hopf_all_inputs(long_pyz,long_G, node_space);
Df_diff = der_KS_Hopf(pyz1);
Df_diff_big =[Df_diff(:,1),0*Df_diff(:,1),Df_diff(:,2:end)];
% spy(abs(Df_diff_big-Df)>10^-6)

% max(abs((KS_all_inputs(nu_Hopf,y,G,const_G) - KS(y,G,const_G))))
% works

% plot(abs(KS_Hopf_vec_all_inputs(long_pyz,long_G,const_G_Hopf,node)-...
%     KS_Hopf_vec_a(pyz1,G_Hopf,const_G_Hopf,node)))
% works

% continuation
step_size = 10^-6;
continuation_equation = kernel(Df);
approx_sol = long_pyz + step_size * continuation_equation;
big_G = [continuation_equation.';long_G];
const_cont_eq = - dot( approx_sol, continuation_equation);
big_const_G_Hopf = [const_cont_eq, const_G_Hopf];

KS_Hopf_vec_loc = @(x) KS_Hopf_vec_all_inputs(x,big_G, big_const_G_Hopf, node_space);
der_KS_Hopf = @(x) fin_diff_KS_Hopf_all_inputs(x,big_G, node_space);



precise_sol = Newton(approx_sol, @(x) KS_Hopf_vec_loc(x),...
    @(x) der_KS_Hopf(x), @(x) der_KS_Hopf(x), @(x) norm_nu_vec(x,node_space), 10);



X =PDE_vector(long_pyz(1:3), long_pyz(3+(1:node_space*2+1)).', ...
    shiftdim( reshape(long_pyz(node_space*2+5:end),[],2*node_space+1),-1));
[X1,Es1,const_Es1] = archlength(X, long_G, const_G_Hopf, continuation_equation, [], step_size);
max(abs(precise_sol - PDE_vec2vec(X1).'));

Df = derivative_KS_all_inputs(pyz1(2),y, G);
nu_Hopf = (pyz1(2));

Df_old = DKS_loc(y(:).',G);
% DF_test1 =  DKS(y(:).',G);
% DF_test2 =  DKS(nu_Hopf,y(:).',G);

nu_Hopf = precise_sol(1);
lambda = precise_sol(2);
a = precise_sol(3);
y = precise_sol(3+(1:2*node_space+1));
z = precise_sol(2*node_space+5:end);
z = reshape(z,[],2*node_space+1);

pyz1_vector =PDE_vector([nu_Hopf,lambda, a], y.', shiftdim(z,-1));

% WORKS
DHopf = derivative_Hopf(pyz1_vector, long_G, continuation_equation);
DHopf_loc = der_KS_Hopf(precise_sol);

pyz1_vector_standard_order =PDE_vector([nu_Hopf,lambda, a], y.', shiftdim(z,-1));
test = compute_linear_func(PDE_vec2vec(pyz1_vector_standard_order),big_G, big_const_G_Hopf);

Hopf_PDE = evaluate_HopfKS(pyz1_vector, long_G, const_G_Hopf, continuation_equation, const_cont_eq);
pyz1_better = Newton_HopfKS(pyz1_vector, long_G, const_G_Hopf, continuation_equation, const_cont_eq);

% test conv_operator
% node_loc = 3;
% w = reshape(1:9,3,3)+1i*reshape(0.5+(1:9),3,3);
% z = reshape(10:18,3,3)+1i*reshape(0.5+(10:18),3,3);
% f_x = @(x) conv2(x,w,'same');
% vec_to_mat = @(x) reshape(x,node_loc,node_loc);
% F_x = @(x) reshape( f_x(vec_to_mat(x)),[],1);
% DF = fin_diff(F_x,z(:));
% DF_test = conv_operator(w);

end


function y = KS_loc(x, G, const_G)
global nu_Hopf

node = (length(x) - 1)/2;

K = diag(-node:node);

x_conv_x = conv(x,x,'same');
y = - nu_Hopf * K^4* x(:) + K^2*x(:) + 1i*K*x_conv_x(:);

y(node+1) =sum(G(:).*x(:))+const_G;

end

function y_dot = KS_all_inputs(nu,y, G, const_G)
%works

node = (length(y) - 1)/2;

K = diag(-node:node);

x_conv_x = conv(y,y,'same');
y_dot = - nu * K^4* y(:) + K^2*y(:) + 1i*K*x_conv_x(:);

y_dot(node+1) =sum(G(:).*y(:))+const_G;

end

function Df = derivative_KS_all_inputs(nu,y, G)

node = (length(y) - 1)/2;

K = diag(-node:node);

Df = - nu * K^4 + K^2 + 2*1i*K*conv_operator(y,1);

Df(node+1,:) =G(:).';

end

function z_dot = KS_Hopf_expansion(lambda,a,y,z,G,const_G)
if size(z,2)~=length(y)
    error('Unmatching sizes')
end
global nu_Hopf

x = [lambda, a, y(:).', z(:).'];

par1 = sum(G(1,:).'.*x(:))+const_G(1);

node_space = (size(z,2)-1)/2;
node_time = (size(z,1)-1)/2;

K = -node_space:node_space;
J = -node_time:node_time;

K2 = repmat(K,2*node_time+1,1);
J2 = repmat(J(:), 1, 2*node_space+1,1);

y_2D = embedding(y,z);

z_dot = - 1i*J2.*z - nu_Hopf*lambda*K2.^4.*z + lambda*K2.^2.*z ...
    + lambda* 1i*K2.*(2*conv2(y_2D,z,'same')+a*conv2(z,z,'same'));

z_dot(node_time+1, node_space+1)=par1;
end

function z_dot = KS_Hopf_expansion_all_inputs(nu,lambda,a,y,z,G,const_G)
global one_int
if size(z,2)~=length(y)
    error('Unmatching sizes')
end

x = [nu, lambda, a,y(:).', z(:).'];

par1 = sum(G(1,:).'.*x(:))+const_G(1);

node_space = (size(z,2)-1)/2;
node_time = (size(z,1)-1)/2;

K = -node_space:node_space;
J = -node_time:node_time;

K2 = repmat(K,2*node_time+1,1);
J2 = repmat(J(:), 1, 2*node_space+1,1);

y_2D = embedding(y,z);
%%%% changed for debugging
if isempty(one_int) || one_int ==1 
z_dot = - 1i*J2.*z - nu*lambda*K2.^4.*z + lambda*K2.^2.*z ...
    + lambda* 1i*K2.*(2*conv2(y_2D,z,'same')+a*conv2(z,z,'same'));
else
    z_dot =conv2(z,z,'same');% - 1i*J2.*z - nu*lambda*K2.^4.*z + lambda*K2.^2.*z ...
    %+ lambda* 1i*K2.*(2*conv2(y_2D,z,'same')+a*conv2(z,z,'same'));
end

z_dot(node_time+1, node_space+1)=par1;
end

function z_dot = KS_Hopf_expansion_a(nu_Hopf,lambda,y,z,G,const_G)
global a_Hopf
a = a_Hopf;
if size(z,2)~=length(y)
    error('Unmatching sizes')
end
x = [nu_Hopf, lambda, y(:).', z(:).'];

par1 = sum(G(1,:).'.*x(:))+const_G(1);

node_space = (size(z,2)-1)/2;
node_time = (size(z,1)-1)/2;

K = -node_space:node_space;
J = -node_time:node_time;

K2 = repmat(K,2*node_time+1,1);
J2 = repmat(J(:), 1, 2*node_space+1,1);

y_2D = embedding(y,z);

z_dot = - 1i*J2.*z - nu_Hopf*lambda*K2.^4.*z + lambda*K2.^2.*z ...
    + lambda* 1i*K2.*(2*conv2(y_2D,z,'same')+a*conv2(z,z,'same'));

z_dot(node_time+1, node_space+1)=par1;
end

function [par, y_dot, z_dot] = KS_Hopf(lambda,a,y,z,G,const_G)
x = [lambda, a, y(:).', z(:).'];

par(1) = sum(G(1,:).'.*x(:))+const_G(1);
par(2) = sum(G(2,:).'.*x(:))+const_G(2);

y_dot = KS_loc(y,G(3,2+(1:length(y))), const_G(3));

z_dot = KS_Hopf_expansion(lambda,a,y,z,G(4,:),const_G(4));
end

function [par, y_dot, z_dot] = KS_Hopf_a(nu,lambda,y,z,G,const_G)
x = [nu, lambda, y(:).', z(:).'];

par(1) = sum(G(1,:).'.*x(:))+const_G(1);
par(2) = sum(G(2,:).'.*x(:))+const_G(2);

y_dot = KS_all_inputs(nu,y,G(3,2+(1:length(y))), const_G(3));

z_dot = KS_Hopf_expansion_a(nu,lambda,y,z,G(4,:),const_G(4));
end

function [par, y_dot, z_dot] = KS_Hopf_all_inputs(nu,lambda,a,y,z,G,const_G)
x = [nu, lambda, a, y(:).', z(:).'];
for i = 1:size(G,1)-2
    par(i) = sum(G(i,:).'.*x(:))+const_G(i);
end
y_dot = KS_all_inputs(nu,y,G(end-1,3+(1:length(y))), const_G(end - 1));

z_dot = KS_Hopf_expansion_all_inputs(nu,lambda,a,y,z,G(end,:),const_G(end));
end

function x_dot = KS_Hopf_vec(x,G,const_G,node)
lambda=x(1);
a= x(2);
y = x(2+(1:2*node+1));
z = x(2*node+4:end);
z = reshape(z,[],2*node+1);

[par, y_dot, z_dot] = KS_Hopf(lambda,a,y,z,G,const_G);
x_dot = [par(:).', y_dot(:).', z_dot(:).'];
end

function x_dot = KS_Hopf_vec_a(x,G,const_G,node)
nu= x(1);
lambda=x(2);
y = x(2+(1:2*node+1));
z = x(2*node+4:end);
z = reshape(z,[],2*node+1);

[par, y_dot, z_dot] = KS_Hopf_a(nu,lambda,y,z,G,const_G);
x_dot = [par(:).', y_dot(:).', z_dot(:).'];
end


function x_dot = KS_Hopf_vec_all_inputs(x,G,const_G,node)
nu= x(1);
lambda=x(2);
a = x(3);

y = x(3+(1:2*node+1));
z = x(2*node+5:end);
z = reshape(z,[],2*node+1);

[par, y_dot, z_dot] = KS_Hopf_all_inputs(nu,lambda,a,y,z,G,const_G);
x_dot = [par(:).', y_dot(:).', z_dot(:).'];
end


function Df = fin_diff_KS_Hopf(x,G, node)
global nu_Hopf

const_G = zeros(1,4);

sol_at_x = KS_Hopf_vec(x,G,const_G,node);

Df = zeros(length(x));
h = 10^-8;
for i = 1:length(x)
    x_i =x;
    x_i(i) = x(i) + h;
    diff_sols = KS_Hopf_vec(x_i, G, const_G,node)-sol_at_x;
    Df(:,i) = diff_sols/h;
end

lambda=x(1);
a= x(2);
y = x(2+(1:2*node+1));
z = x(2*node+4:end);
z = reshape(z,[],2*node+1);
y_2D = embedding(y,z);

node_space = (size(z,2)-1)/2;
node_time = (size(z,1)-1)/2;

K = -node_space:node_space;
J = -node_time:node_time;

K2 = repmat(K,2*node_time+1,1);

der_lambda = - nu_Hopf * K2.^4 .* z + K2.^2.* z + 1i*K2.*(2*conv2(z,y_2D,'same')+a*conv2(z,z,'same')); 
der_a = lambda * 1i * K2 .* conv2(z,z,'same');

index = 2 + length(y);
Df(index+1:end,2) = der_a(:);
Df(index+1:end,1) = der_lambda(:);
end

function Df = fin_diff_KS_Hopf_all_inputs(x,G, node)

const_G = zeros(1,4);

sol_at_x = KS_Hopf_vec_all_inputs(x,G,const_G,node);

Df = zeros(length(x)-5+size(G,1),length(x));
h = 10^-8;
for i = 1:length(x)
    x_i =x;
    x_i(i) = x(i) + h;
    diff_sols = KS_Hopf_vec_all_inputs(x_i, G, const_G,node)-sol_at_x;
    Df(:,i) = diff_sols/h;
end
end

function Df = fin_diff_KS_Hopf_a(x,G, node)
%global a_Hopf

const_G = zeros(1,4);

sol_at_x = KS_Hopf_vec_a(x,G,const_G,node);

Df = zeros(length(x));
h = 10^-8;
for i = 1:length(x)
    x_i =x;
    x_i(i) = x(i) + h;
    diff_sols = KS_Hopf_vec_a(x_i, G, const_G,node)-sol_at_x;
    Df(:,i) = diff_sols/h;
end

% lambda=x(1);
% nu_Hopf= x(2);
% a = a_Hopf;
% y = x(2+(1:2*node+1));
% z = x(2*node+4:end);
% z = reshape(z,[],2*node+1);
% y_2D = embedding(y,z);
% 
% node_space = (size(z,2)-1)/2;
% node_time = (size(z,1)-1)/2;
% 
% K = -node_space:node_space;
% J = -node_time:node_time;
% 
% K2 = repmat(K,2*node_time+1,1);
% 
% der_lambda = - nu_Hopf * K2.^4 .* z + K2.^2.* z + 1i*K2.*(2*conv2(z,y_2D,'same')+a*conv2(z,z,'same')); 
% der_nu = lambda*K2.^4 .* z;
% 
% index = 2 + length(y);
% Df(index+1:end,2) = der_nu(:);
% Df(index+1:end,1) = der_lambda(:);
end


function Df = DKS_loc(x,G)
global nu_Hopf

node = (length(x)-1)/2;

K = diag(-node:node);

x_first = x(node+1:end);
x_second = x(node+1:-1:1);

x_first = [x_first, 0*x_first(2:end)];
x_second = [x_second, 0*x_second(2:end)];

Df_y = - nu_Hopf * K^4 + K^2 + 1i * 2 * K * toeplitz(x_first, x_second);

Df = Df_y;
Df(1+node,:)= G(:);
end

function Df = fin_diff_KS(x,G)
Df = zeros(length(x));
h = 10^-6;
for i = 1:length(x)
    x_i =x;
    x_i(i) = x(i) + h;
    Df(:,i) = (KS_loc(x_i, G, 0)-KS_loc(x, G, 0))/h;
end
end

function x  = Newton(x0, f, df,fin_diff, norm_nu_loc, tries_max)
global display
if nargin<6
    tries_max = 100;
    if nargin<5
        norm_nu_loc = @(x) norm_nu(x);
    end
end
x = x0;
tries = 0;
while max(norm_nu_loc(f(x)))>10^-6 && tries<tries_max
    df_fin_dif = fin_diff(x);
    DF = df(x(:).');
    if display>1
        fprintf('Cond %e,    Norm %e\n',rcond(DF), norm_nu_loc(f(x)))
    end
    F = f(x);
    x = x(:) - DF\F(:);
    tries = tries +1;
end
if max(norm_nu_loc(f(x)))>10^-6
    error('Newton did not converge')
end
if any(isnan(x))
    error('Result contains NaN')
end
end

function x = vec_to_struct(v,node)
x.lambda = v(1);
x.a=v(2);
x.nu = v(3);
x.y = v(3+(1:2*node+1));
x.z = x(2*node+5:end);
x.z = reshape(x.z,[],2*node+1);
end

function v = struct_to_vec(x)
v = [x.nu, x.lambda, x.a, x.y(:).', x.z(:).'];
end

function Der = fin_diff(f,x)

fx = f(x);
Der = zeros(length(fx), length(x));

h = 10^-6;

for i = 1:length(x)
    dx = x;
    dx(i) = x(i)+h;
    fdx = f(dx);
    Der(:,i) = (fdx - fx)/h;
end


end


function n = norm_nu(x)
global nu
node = (length(x)-1)/2;
K = -node:node;
n = sum(nu.^abs(K).*abs(x(:).'));
end

function n = norm_nu_mat(z)
global nu
node_space = (size(z,2)-1)/2;
node_time = (size(z,1)-1)/2;
K = -node_space:node_space;
K = repmat(K,2*node_time+1, 1);

J = -node_time:node_time;
J = repmat(J(:),1,2*node_space+1);

n = sum(sum(nu.^(abs(K)+abs(J)).*abs(z)));
end

function n = norm_nu_vec(x,node)
num_par = mod(length(x) , 2*node+1);
for i = 1:num_par
    n(1) = abs(x(i));
end
y = x(num_par+(1:2*node+1));
z = x(2*node+num_par+2:end);
z = reshape(z,[],2*node+1);
n(3) = vector_norm(y);
n(4) = vector_norm(z);
n = max(n);
end