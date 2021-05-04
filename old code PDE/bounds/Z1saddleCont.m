function Z1cont = Z1saddleCont(A0_in, A1_in, A_dagger0_in,A_dagger1_in,X0,X1,G,Es0,Es1)

if size(A_dagger0_in,1)<size(A0_in,1)
    A_dagger0_in = big_Adag_saddle(A_dagger0_in,X0);
end
if size(A_dagger1_in,1)<size(A0_in,1)
    A_dagger1_in = big_Adag_saddle(A_dagger1_in,X1);
end

A0 = fullopA_saddle(A0_in,X0);
A1 = fullopA_saddle(A1_in,X1);
A_dagger0 = fullopAdag_saddle(A_dagger0_in,X0);
A_dagger1 = fullopAdag_saddle(A_dagger1_in,X1);

ADelta = A0-A1;
As = intval(A0,A1);
A_daggerDelta = A_dagger0 - A_dagger1;
EsDelta = Es0 - Es1;
Xs = interval_PDE(X0,X1);
XDelta = X0 - X1;

DS0 = derivative_saddle_fullop(X0,G,Es0,EsDelta);
DS1 = derivative_saddle_fullop(X0,G,Es1,EsDelta);

Z1_extrema0 = A0*(A_dagger0-DS0);
Z1_extrema1 = A1*(A_dagger1-DS1);

central_bound = 2*ADelta*A_daggerDelta-2*ADelta*second_der_saddle_fullop(Xs,XDelta)-...
    As*third_der_saddle_fullop(Xs,XDelta,XDelta);

%Z1cont = norm(max(abs(Z1_extrema0),abs(Z1_extrema1)))+...
%    1/8*norm(central_bound);

extrema_operator = cnorm(max(abs(Z1_extrema0.center),abs(Z1_extrema1.center)));
extrema_nocenter0 = Z1_extrema0;
extrema_nocenter1 = Z1_extrema1;
extrema_nocenter0.center = 0*Z1_extrema0.center;
extrema_nocenter1.center = 0*Z1_extrema1.center;
extrema_tail = max(cnorm(extrema_nocenter0),cnorm(extrema_nocenter1));

extrema_norm = sum( max(extrema_operator, extrema_tail),2);

Z1cont = extrema_norm+...
    1/8*norm(central_bound);

end






% function Z1cont = Z1saddleCont(A0, A1, A_dagger0,A_dagger1,X0,X1,G,Es0,Es1)
% % function Z1cont = Z1HopfCont(A0, A1, A_dagger0,A_dagger1,X0,X1,G,Es0,Es1)
% 
% Z1_extrema = Z1extrema(A0,A1,A_dagger0,A_dagger1,X0,X1,G,Es0,Es1);
% 
% Xs = interval_PDE(X0,X1);
% 
% second_der_term = cont_termZ1(A0,A1,A_dagger0,A_dagger1,X0,X1,Xs);
% 
% Z1cont = Z1_extrema+1/8*second_der_term;
% end
% 
% function Z1_extrema = Z1extrema(A0,A1,A_dagger0,A_dagger1,X0,X1,G,Es0,Es1)
% % function Z1_extrema = Z1extrema(A0,A1,A_dagger0,A_dagger1,X0,X1,G,Es0,Es1)
% % compute Z1 at the extrema of the segment
% 
% big = 1 ;
% A_dagger_big0 = big_Adaggersaddle(A_dagger0,X0);
% DH0 =  overop2mat(derivative_saddle(X0,G,Es0,big));
% A_dagger_big1 = big_Adaggersaddle(A_dagger1,X1);
% DH1 =  overop2mat(derivative_saddle(X1,G,Es1,big));
% 
% max_op = max(abs(A0*(DH0-A_dagger_big0)),abs(A1*(DH1-A_dagger_big1)));
% 
% Z1_finite = cnorm(compose_overoperator3(max_op,reshape(X0,2*X0.node_space,2*X0.node_time)));
% Z1_tail = Z1Tailsaddle(max(abs(X0),abs(X1)));
% 
% Z1_extrema = sum(max(Z1_finite,Z1_tail),2); 
% end
% 
% function second_der_term = cont_termZ1(A0,A1,Adag0,Adag1,X0,X1,Xs)
% % function second_der_term = cont_termZ1(A0,A1,Adag0,Adag1,X0,X1,Xs)
% % term dealing with the second derivative
% 
% ADelta = A0-A1;
% A_dagger_big0 = big_Adagger(Adag0,X0);
% A_dagger_big1 = big_Adagger(Adag1,X1);
% AdagDelta = A_dagger_big0-A_dagger_big1;
% 
% first_term= norm_ADelta_AdagDelta(ADelta, AdagDelta,X0,X1);
% 
% second_term = sum(norm_ADelta_D2(A0,A1, X0, X1,Xs),2);
% 
% third_term =norm_As_D3(A0, A1,Xs, X0, X1);
% 
% second_der_term = first_term+second_term+third_term;
% 
% end
% 
% function second_term = norm_As_D3(A0, A1, Xs, X0, X1)
% % function second_term = norm_As_D3(A0, A1, Xs, X0, X1)
% %
% % second term: \| A_s D^3 H(xs) xDelta xDelta \|
% 
% X0_big = reshape(X0, X0.node_space*3, X0.node_time*3);
% X1_big = reshape(X1, X1.node_space*3, X1.node_time*3);
% Xs_big = reshape(Xs, Xs.node_space*3, Xs.node_time*3);
% 
% DDxHxD = overop2mat(third_derHopf(Xs_big,X0_big-X1_big)); 
% 
% A0 = big_A(A0,X0,X0_big);
% A1 = big_A(A1,X1,X1_big);
% 
% second_term_finite = cnorm( compose_overoperator((A0-A1)*DDxHxD, X0_big));
% second_term_infinite = tail_ADelta_D2(X0_big,X1_big);
% 
% second_term = sum(max(second_term_finite,second_term_infinite),2);
% end
% 
% function second_term = norm_ADelta_D2(A0, A1, X0, X1, Xs)
% % function second_term = norm_ADelta_D2(A0, A1, X0, X1, Xs)
% % \| ADelta D^2 H(xs) xDelta \|
% 
% X0_big = reshape(X0, X0.node_space*5, X0.node_time*5);
% X1_big = reshape(X1, X1.node_space*5, X1.node_time*5);
% Xs_big = reshape(Xs, Xs.node_space*5, Xs.node_time*5);
% 
% DDxHxD = second_derHopf(Xs_big,X0_big-X1_big); 
% 
% A0 = big_A(A0,X0,X0_big);
% A1 = big_A(A1,X1,X1_big);
% 
% second_term_finite = cnorm(compose_overoperator(A0-A1, X0_big)*DDxHxD);
% second_term_infinite = tail_ADelta_D2(X0,X1);
% 
% second_term = max(second_term_finite,second_term_infinite);
% end
% 
% function tail_ADelta_D2 = tail_ADelta_D2(X0,X1)
% % function tail_ADelta_D2 = tail_ADelta_D2(X0,X1)
% %
% % tail of \| ADelta D^2 H(xs) xDelta \|
% 
% tail_ADelta_D2 = zeros(5);
% 
% node_space = X0.node_space;
% 
% x0 = norm(X0);
% x1 = norm(X1);
% x = max(x0,x1);
% xD = norm(X0-X1);
% 
% y = x(4);
% z = x(5);
% l1 = x(1);
% l2 = x(2);
% %a = x(3);
% 
% yD = xD(4);
% zD = xD(5);
% l1D = xD(1);
% l2D = xD(2);
% aD = xD(3);
% 
% %y1 = x1(4);
% %z1 = x1(5);
% l11 = x1(1);
% l21 = x1(2);
% %a1 = x1(3);
% 
% B1 = tailBoundHopf(1,l1, l2, node_space);
% 
% a = l1D*l2+l1*l2D;
% b = l1D;
% c = l11*l21;
% d = l11;
% 
% f2 = @(k) (k^2*a+ b)/(k^2*c + d);
% B2 = max(abs(f2(node_space+1)),abs(a/c));
% 
% tail_ADelta_D2(4:5,:)= [ 0,0,0, B2+4*B1*(l1D*y+l1*yD), 0;
%     0,0,0,4*B1*(l1*zD+l1D*z),B2+4*B1*(l1D*y+2*l1D*a*z+2*l1*aD*z+l1*yD+2*l1*a*zD)];
% end
% 
% 
% 
% function first_term= norm_ADelta_AdagDelta(ADelta, AdagDelta,X0,X1)
% % function first_term= norm_ADelta_AdagDelta(ADelta, AdagDelta,X0,X1)
% % 
% % \| ADelta A^dagger Delta \|
% 
% center = norm( compose_overoperator(ADelta*AdagDelta, X0));
% 
% nAtilde = normAtilde(X0,X1);
% 
% tail = [0,0,0,nAtilde,nAtilde].';
% 
% first_term = max(center, tail);
% 
% end
% 
% % if holds
% % ADelta = ATilde * A1
% % we compute here the norm of ATilde
% 
% function nAtilde = normAtilde(X0,X1)
% % function nAtilde = normAtilde(X0,X1)
% %
% % \| Atilde \|
% 
% l10 = X0.parameters(1);
% l20 = X0.parameters(2);
% l11 = X1.parameters(1);
% l21 = X1.parameters(2);
% 
% 
% a = -l21 *l11 + l10*l20;
% c = l10*l20;
% 
% f= @(x) ( x.^2*a + l11 - l10)./(x.^2*c - l10);
% 
% nAtilde = max([abs(c/a), abs(f(X0.node_space+1))]);
% 
% end