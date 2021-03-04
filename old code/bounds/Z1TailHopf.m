function Z1_tail = Z1TailHopf(X)
% function Z1_tail = Z1TailHopf(X)

Z1_tail = zeros(5);

node_space = X.node_space;

x = norm(X);

y = x(4);
z = x(5);
nu = x(1);
lambda = x(2);
a = x(3);

B1 = tailBoundHopf(1,nu, lambda, node_space);%K2overLinear;

Z1_tail(4:5,:)= [ 0, 0, 0, 2*lambda*B1*y, 0;
    B1*(2*z*y+a*z^2),0, lambda*B1*z^2,B1*2*lambda*z,(y+a*z)*2*lambda*B1];
end