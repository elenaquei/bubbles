function algebraic_hopf_lorenz()

% the code request you to specify fn_  derivatives_ second_der and
% third_der

% ALL DATA FOR LORENZ 84
f = @fn_Lorenz84; % We choose a map.
df = @derivatives_Lorenz84;
ddf = @second_der_Lorenz84;
dddf = @third_der_Lorenz84;
N = 4;
DIM= N;
X = [0.02 0.52   0.005   0.6   0.91   0.2   0.8976   0.7768   0.72   0.155   0.538   0.2578   0.47   0.51]'; %lorenx84
phi = [0.347693611712440   0.362607042331847   0.999178390839916   0.649549106979193]'; %lorenx84

[X]=newton_Hopf(f,X,phi);

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
save('hopf_in_Rossler', 'X','l1','x_star','lambda_star','eigenvec','eigenval');



return



