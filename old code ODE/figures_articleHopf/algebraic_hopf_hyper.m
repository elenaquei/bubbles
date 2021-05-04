function dummy = algebraic_hopf_hyper(sol_num)

% the code request you to specify fn_  derivatives_ second_der and
% third_der
random = 0;
if nargin<1
    sol_num = 1;
    random = 1;
end
if isempty(sol_num) || mod(sol_num,1)~=0 || sol_num<1 || sol_num>2
    sol_num =1;
    warning('wrong input number')
end


f = @fn_hyper; % We choose a map.
df = @derivatives_hyper;
ddf = @second_der_hyper;
dddf = @third_der_hyper;
N = 4;

X_sol = cell(2,1);
phi_sol = cell(2,1);

X_sol{1} =[-0.000000000000000
   1.414213562373095
   0.000000000000000
  -0.000000000000000
   0.000000000000000
   0.000000000000000
  -0.140851690016498
  -0.000000000000000
   0.311445987589876
  -0.234285486293461
   0.165664856091694
   0.000000000000000
   1.216197421320658
  -0.199194370304503];
phi_sol{1} =[0.7629
    0.4609
    0.8222
    0.6343];
% OR
X_sol{2} = [ -1.015513726331095
  1.464130774004092
  10.099019513592780
 -10.099019513592783
  -1.000000000000000
  30.610405390670862
  -0.607676658921156
   0.030310136744525
  -0.048061582141663
   3.964552188224825
   2.839481206874421
   0.062830225198066
   0.282991177919101
   0.788681629272493];
phi_sol{2} = [0.015584636790005
   0.609679882074390
   0.384647202347126
   0.867038269926880];
% OR
% X_sol{3} =[-1.015513725673144
%   -1.464130773588491
%   10.099019513592782
%  -10.099019513592785
%   -1.000000000000000
%   30.610405377381553
%   -0.920387258631270
%   -0.006032594052934
%   -0.087633536901453
%    1.098301280709103
%    0.775747371688678
%    0.028280980220684
%    0.080489086826049
%    1.265487850348961];
% phi_sol{3}=[0.512209857544997
%    0.949906183387119
%    0.143244539041576
%    0.445883892900819];

X = X_sol{sol_num};
phi = phi_sol{sol_num};

if random
    X = rand(size(X_sol{1}))*140-70;
    phi = rand(size(phi_sol{1}));
end

[X]=newton_Hopf(f,X,phi);

X(2) = abs(X(2));

if max(abs(X(1:6)-X_sol{1}(1:6)))>10^-6 && ....
        max(abs(X(1:6)-X_sol{2}(1:6)))>10^-6
    disp('maybe new solution')
end



alpha=X(1); % Parameter 
beta=X(2); % distance of the complex conjugate eigs from the real axis

x=X(2+(1:N)); % the variables from the model
v1=X(N+2+(1:N)); v2=X(2*N+2+(1:N)); % The real and imaginary parts of the eigenvector

%%%% validation
% for the validation, the exact derivative of f is necessary
DF = big_derivative(df,X,phi);
if rank(DF)-size(DF,1)<0
    error('Derivative not invertible')
end

A=inv(DF);

F=F_general_Hopf(f,X,phi,df);

% start intlab
Y = norm(A*F);

Z1 = norm( eye(length(X)) - A *DF);

R = 0.3; % first approximation of the maximum radius


XandR=infsup(X-R,X+R);
DDF= second_derivative_F(XandR,phi,df);

Z2 = norm(sup(abs(A*DDF)));

pol=@(r) Y+(Z1-1)*r+Z2*r.^2;

delta= (intval(Z1)-1)^2-4*intval(Z2)*intval(Y);
if delta<=0 
    error('Hopf not validated')
end
rmin=sup((-(Z1-1)-sqrt(delta))/(2*Z2));
rmax=inf((-(Z1-1)+sqrt(delta))/(2*Z2));

% check rmax<R 

if rmax>R
    warning('computation has to be rerun with higher R')
end

r=linspace(0,1.3*rmax,100);
plot(r,pol(r),'b',rmin,0,'ok',rmax,0,'or');


% compute the first Lyapunov coefficient as in Kuznetsov
q =infsup( v1 + 1i*v2-rmin, v1 + 1i*v2+rmin);

% p: DF(xH)^Tp=conj(beta)p
% p = null( DF(xH)^T-conj(beta)I)

df_mat=df(x,alpha);
% next line does not work
%p = null( df_mat.'-conj(1i*beta)*eye(size(df_mat)));

[all_p,all_beta]=eigs(df_mat.');
all_eigs=diag(all_beta);
[~,index_beta]=min(all_eigs-conj(1i*beta));
p=all_p(:,index_beta);
% verification of the eigenvalue
[l,p]=verifyeig(midrad(df_mat',rmin),conj(1i*beta),p);

complex_product=@(a,b) sum(conj(a).*b);


% look for k, the rescaling coefficient such that <kp,q>=1
% k = k1 + i k2

% k1/k2
ratio = real(complex_product(p,q))/imag(complex_product(p,q));
k2 = 1/( ratio* real(complex_product(p,q))+imag(complex_product(p,q)));
k1= k2*ratio;

p=(k1+1i*k2)*p;

beta_ver=infsup(beta-rmin,beta+rmin);

% compute the first lyapunov coefficient
l1 = 1/(2*beta_ver) * real(1i*complex_product(p,ddf(zeros(N,1),q,q)*complex_product(p,ddf(zeros(N,1),q,conj(q))))+...
    beta_ver*complex_product(p,dddf(zeros(N,1),q,q,conj(q))));

% check the FLC is nonzero
if intval(inf(l1))*intval(sup(l1))<0
    error('Could not validate Hopf, since the first Lyapunov coefficient is not verifyed to be non-zero')
end


% last check: all ther eigenvalues different then zero
[all_eigenvectors,all_eigenvalues]=eigs(df_mat.');
all_eigenvalues = diag(all_eigenvalues);
[~,index_beta1]=min(all_eigenvalues-(1i*beta));
[~,index_beta2]=min(all_eigenvalues+(1i*beta));
for i=1:length(all_eigs)
    if i~=index_beta1 && i~=index_beta2
        [L,~] = verifyeig(df_mat,all_eigenvalues(i),all_eigenvectors(:,i));
        if intval(inf(real(L)))*intval(sup(real(L)))<0
            error('Could not validate Hopf, since there is at least one eigenvalue that not verifyed to be non-zero')
        end
    end

end
fprintf('SUCCESS\n The Hopf bifurcation at [%1.3f,%1.3f], %1.3f has been verified\n',x(1),x(2),alpha)
fprintf('with a radius of %e to %e\n\n',rmin, rmax) 

x_star = X(3:6);
lambda_star = X(1);
eigenvec = v1+1i*v2;
eigenval = 1i*beta;
save('hopf_in_hyper', 'X','l1','x_star','lambda_star','eigenvec','eigenval');



return



