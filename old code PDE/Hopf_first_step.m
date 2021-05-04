function [X, G_Hopf, const_G_Hopf, continuation_equation, const_cont_eq]= Hopf_first_step...
    (node_space_out, node_time_out)
global display
global nu_Hopf
global a_Hopf
if isempty(display)
    display = 1;
end
% ordering of parameters: nu, lambda, a 
node_space = 10;
node_time = 3;


vec2PDE_vec_loc = @(x) PDE_vector(x(1:3), x(3+(1:node_space*2+1)).', ...
    shiftdim( reshape(x(node_space*2+5:end),[],2*node_space+1),-1));


try
    intval(1);
catch
    addpath(genpath('/Users/queirolo/Dropbox/Hopf_Bifurcations/2D_examples/Intlab_V7.1'));
    startintlab
end

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


G2 = zeros(1,1+2*node_space);
K = -node_space:node_space;
G2(1:(2*node_space+1)) = 1i*approx.*K;

G = zeros(1, length(G2));
G(1,:) = G2;
const_G= 0;

% better initial approx
nu_Hopf = 0.131; %pos

x2 = Newton(approx, @(x) KS(x,G,const_G), @(x) DKS(x,G));


nu_Hopf = 0.132; %neg

x2 = Newton(x2, @(x) KS(x,G,const_G), @(x) DKS(x,G));

x2 = (x2 + conj(x2(end:-1:1)))/2;

[V2, D2]= eig(DKS(x2(:).',0*G));

% D3 = D2;
% tries = 1;
% num_D2 = sum(real(diag(D2)>0));
% max_tries = 6300;
% while abs(num_D2 -sum((diag(D3))>0))~=2 && tries<max_tries
%     nu_Hopf = nu_Hopf+10^-3;
% 
%     x2 = Newton(x2, @(x) KS(x,G,const_G), @(x) DKS(x,G));
% 
%     [V2, D3]= eig(DKS(x2(:).',0*G));
%     num_D2 = sum(real(diag(D3)>0));
% 
%     tries = tries+1;
% end
% if tries< max_tries
%     D2 = D3;
% end

z_hat = zeros(2*node_time+1, 2*node_space+1);

% D2, V2, x2 chosen solutions for nu_Hopf_2
% with DKS(x2(:).'

[~, index] = sort(real(diag(D2)));

index_eig = index(end-2);

z_hat(node_time, : ) = V2(:,index_eig);
z_hat(node_time+2, : ) = conj(V2(end:-1:1,index_eig));
important_eig = D2(index_eig,index_eig);
if node_space~=10
    error('chosen column of V2 computed just for 10 spatial nodes')
end
z_hat = z_hat/sqrt(sum(sum(z_hat.*conj(z_hat))));

y_hat = x2;
a_hat = -10^-8;

y = x2;
z = z_hat;
%nu_Hopf = nu_Hopf;
lambda = abs(imag(important_eig))^-1;
a = a_hat;

K = -node_space:node_space;
J = -node_time:node_time;

K2 = repmat(K,2*node_time+1,1);
J2 = repmat(J(:), 1, 2*node_space+1,1);

i_J_z_hat = 1i*J2.*z_hat;
G1 = conj([0,0,0*y_hat(:).',i_J_z_hat(:).']);
const_G1 =0;
K_y_hat = 1i*y_hat.*K.';
G3 = conj([0,0,K_y_hat(:).',0*z_hat(:).']);
const_G3 = 0;
%G2 = [0,0,0*y_hat(:).',conj(z_hat(:)).']; %old
J2z = (J2.^2.*z);
G2 = conj([0,0,0*y_hat(:).',(J2z(:)).']); % NEW 
const_G2 = -1;
y_2D = embedding(y,z);
y_plus_az = y_2D+a_hat*z_hat;
G4 = conj([0,0, 0*y_hat(:).', 1i*K2(:).'.*y_plus_az(:).']);
const_G4 = 0;

G_Hopf= [G1;G2;G3;G4];
const_G_Hopf = [const_G1,const_G2,const_G3,const_G4];

% shift fixed parameter from nu to a

a_Hopf = a;
pyz0 = [nu_Hopf,lambda,y(:).',z(:).'];

KS_Hopf_vec_loc = @(x) KS_Hopf_vec_a(x,G_Hopf, const_G_Hopf, node_space);
der_KS_Hopf = @(x) fin_diff_KS_Hopf_a(x,G_Hopf, node_space);

pyz1 = Newton(pyz0, @(x) KS_Hopf_vec_loc(x),...
    @(x) der_KS_Hopf(x), @(x) norm_nu_vec(x,node_space), 10);

long_pyz = [pyz1(1:2);a_Hopf;pyz1(3:end)];
PDE_vec_pyz = vec2PDE_vec_loc(long_pyz);


% y = PDE_vec_pyz.vector1D;
% z = squeeze(PDE_vec_pyz.vector2D);
% a = a_Hopf;
% 
% K = -node_space:node_space;
% J = -node_time:node_time;
% 
% K2 = repmat(K,2*node_time+1,1);
% J2 = repmat(J(:), 1, 2*node_space+1,1);
% 
% i_J_z = 1i*J2.*z;
% G1 = conj([0,0,0*y(:).',i_J_z(:).']);
% const_G1 =0;
% K_y = 1i*y.*K;
% G3 = conj([0,0,K_y(:).',0*z(:).']);
% const_G3 = 0;
% %G2 = [0,0,0*y_hat(:).',conj(z_hat(:)).']; %old
% J2z = (J2.^2.*z);
% G2 = conj([0,0,0*y(:).',(J2z(:)).']); % NEW 
% const_G2 = -1;
% y_2D = embedding(y,z);
% y_plus_az = y_2D+a*z;
% G4 = conj([0,0, 0*y(:).', 1i*K2(:).'.*y_plus_az(:).']);
% const_G4 = 0;
% 
% G_Hopf= [G1;G2;G3;G4];
% const_G_Hopf = [const_G1,const_G2,const_G3,const_G4];

G_Hopf = [G_Hopf(:,1:2),zeros(size(G_Hopf,1),1),G_Hopf(:,3:end)];


Df = derivative_Hopf(PDE_vec_pyz,G_Hopf);

% continuation
continuation_equation = kernel(Df);
continuation_equation = 0 * continuation_equation; % NEW CONT EQ.
continuation_equation(3) = 1;
const_cont_eq = - PDE_vec_pyz.parameters(3);%dot( long_pyz, continuation_equation);

X =PDE_vector(long_pyz(1:3), long_pyz(3+(1:node_space*2+1)).', ...
    shiftdim( reshape(long_pyz(node_space*2+5:end),[],2*node_space+1),-1));

if node_space_out~=node_space || node_time_out~=node_time
    X = reshape(X, node_space_out, node_time_out);
    G_Hopf = reshape_linear_coef(G_Hopf,3, node_space,node_space_out, node_time_out);
    continuation_equation = reshape_linear_coef(continuation_equation.',3, node_space,node_space_out, node_time_out).';
    X = Newton_HopfKS(X,G_Hopf, const_G_Hopf, continuation_equation, const_cont_eq);
end
end

function y_dot = KS_all_inputs(nu,y, G, const_G)
%works

node = (length(y) - 1)/2;

K = diag(-node:node);

x_conv_x = conv(y,y,'same');
y_dot = - nu * K^4* y(:) + K^2*y(:) + 1i*K*x_conv_x(:);

y_dot(node+1) =sum(G(:).*y(:))+const_G;

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

function [par, y_dot, z_dot] = KS_Hopf_a(nu,lambda,y,z,G,const_G)
x = [nu, lambda, y(:).', z(:).'];

par(1) = sum(G(1,:).'.*x(:))+const_G(1);
par(2) = sum(G(2,:).'.*x(:))+const_G(2);

y_dot = KS_all_inputs(nu,y,G(3,2+(1:length(y))), const_G(3));

z_dot = KS_Hopf_expansion_a(nu,lambda,y,z,G(4,:),const_G(4));
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
end

function n = norm_nu(x)
global nu
node = (length(x)-1)/2;
K = -node:node;
n = sum(nu(end).^abs(K).*abs(x(:).'));
end

function n = norm_nu_mat(z)
global nu
node_space = (size(z,2)-1)/2;
node_time = (size(z,1)-1)/2;
K = -node_space:node_space;
K = repmat(K,2*node_time+1, 1);

J = -node_time:node_time;
J = repmat(J(:),1,2*node_space+1);

n = sum(sum(nu(end).^abs(K).*nu(1).^abs(J).*abs(z)));
end

function n = norm_nu_vec(x,node)
num_par = mod(length(x) , 2*node+1);
for i = 1:num_par
    n(1) = abs(x(i));
end
y = x(num_par+(1:2*node+1));
z = x(2*node+num_par+2:end);
z = reshape(z,[],2*node+1);
n(3) = norm_nu(y);
n(4) = norm_nu_mat(z);
n = max(n);
end