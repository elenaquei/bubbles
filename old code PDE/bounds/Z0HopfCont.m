function Z0cont = Z0HopfCont(Adag0,Adag1,X0,X1,G,Es0,Es1)

%not_big = 0;
A0_small =  Adag1^-1;%matrix_AHopf(X0, G, Es0, not_big);
Z0_extrema0 = abs(eye(size(A0_small,1)) - A0_small*Adag0);

A1_small = Adag0^-1;%matrix_AHopf(X1, G, Es1, not_big);
Z0_extrema1 = abs(eye(size(A1_small,1)) - A1_small*Adag1);

first_term = norm(compose_overoperator(max(Z0_extrema0,Z0_extrema1),...
    X0));

ADelta = A0_small-A1_small;
AdagDelta = Adag0-Adag1;

second_term_center = norm( compose_overoperator(ADelta*AdagDelta, X0));

nu0 = X0.parameters(1);
lambda0 = X0.parameters(2);
nu1 = X1.parameters(1);
lambda1 = X1.parameters(2);

K = X0.node_space+1;

nu_Delta = nu0 - nu1;
nu_lambda_Delta = nu0*lambda0 - nu1*lambda1;
lambda_Delta = lambda0 - lambda1;

tail_44_0 = bound_tail_AD_AdagD(nu_Delta,0,nu0,-1, X0.node_space);
tail_44_1 = bound_tail_AD_AdagD(nu_Delta,0,nu0,-1, X0.node_space);
tail_44 =tail_44_0 * tail_44_1;

lim_inf_0 = nu_lambda_Delta/(nu0*lambda0);
lim_inf_1 = nu_lambda_Delta/(nu1*lambda1);

lim_infty = lim_inf_0*lim_inf_1;

% f0=@(x) (-nu_lambda_Delta*x^2+lambda_Delta)/(-nu0*lambda0*x^2+lambda0);
% f1=@(x) (-nu_lambda_Delta*x^2+lambda_Delta)/(-nu1*lambda1*x^2+lambda1);
% f = @(x) f1(x)*f0(x);

nAtilde = max([abs(lim_infty), ...
    bound_tail_AD_AdagD(nu_lambda_Delta,-lambda_Delta,nu0*lambda0,-lambda0,  X0.node_space),...
    bound_tail_AD_AdagD(nu_lambda_Delta,-lambda_Delta,nu1*lambda1,-lambda1,  X0.node_space)]);%abs(f(X0.node_space+1))]);

second_term_tail = [0,0,0,tail_44,nAtilde].';

second_term = max(second_term_center,second_term_tail);

Z0cont = first_term + 1/4*second_term;
end

function nAtilde = bound_tail_AD_AdagD(a,b,c,d, nodes)

argmax = 0;
if -d/c>0
    argmax = max(argmax,sqrt(-d/c));
end
if argmax>nodes
    check_func = nodes:ceil(argmax);
else
    check_func = nodes;
end
f = @(k) (a*k.^2+b)./(c*k.^2+d);

nAtilde = max(abs(f(check_func)));

end