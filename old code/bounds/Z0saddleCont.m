function Z0cont = Z0saddleCont(Adag0,Adag1,X0,X1,G,Es0,Es1)

not_big = 0;
A0_small = matrix_Asaddle(X0, G, Es0, not_big);
Z0_extrema0 = abs(eye(size(A0_small,1)) - A0_small*Adag0);

A1_small = matrix_Asaddle(X1, G, Es1, not_big);
Z0_extrema1 = abs(eye(size(A1_small,1)) - A1_small*Adag1);

first_term = norm(compose_overoperator3(max(Z0_extrema0,Z0_extrema1),...
    X0));

ADelta = A0_small-A1_small;
AdagDelta = Adag0-Adag1;

second_term_center = norm( compose_overoperator3(ADelta*AdagDelta, X0));

l10 = X0.parameters(1);
l20 = X0.parameters(2);
l11 = X1.parameters(1);
l21 = X1.parameters(2);

a = -l21 *l11 + l10*l20;
c = l10*l20;

f= @(x) ( x.^2*a + l11 - l10)./(x.^2*c - l10);

nAtilde = max([abs(c/a), abs(f(X0.node_space+1))]);
second_term_tail = [zeros(1,9),nAtilde*ones(1,6)].';

second_term = max(second_term_center,second_term_tail);

Z0cont = first_term + 1/8*second_term;