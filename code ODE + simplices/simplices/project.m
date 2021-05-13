function x = project(x_approx, f)
% takes an approximation x_approx and projects it on the manifold f(x)=0

warning('not debugged')

tol = 10^-6;
DF_x = derivative_to_matrix(derivative(f,x_approx,0));

[Q1Q2,R0] = qr(DF_x.');

Q1 = Q1Q2(:,1:end-2);
% Q2 = Q1Q2(:,end-2:end);
R = R0(1:end-2,:);

xBar = x_approx;

for iter = 1:100
    fx = apply(f,xBar,0);
    
    fx=Xi_vec2vec(fx); % working
    x=Xi_vec2vec(xBar);
    
    z = R\fx;
    step = Q1 * z;
    x = x - step;
    if norm(step)<tol
        break
    end
    
    xBar=vec2Xi_vec(x,xBar.size_scalar,xBar.size_vector,xBar.nodes);
    xBar = symmetrise(xBar);
end
if norm(step)>tol
    warning('Did not converge')
end

end