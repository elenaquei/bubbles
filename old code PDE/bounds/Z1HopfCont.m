function Z1cont = Z1HopfCont(A0_in, A1_in, A_dagger0_in,A_dagger1_in,X0,X1,G,Es0,Es1)

length_big = length(reshape(X0, 2*X0.node_space, 2*X0.node_time));

if size(A_dagger0_in,1)<length_big
    A_dagger0_in = big_Adagger(A_dagger0_in,X0);
end
if size(A_dagger1_in,1)<length_big
    A_dagger1_in = big_Adagger(A_dagger1_in,X1);
end
if size(A0_in,1) < length_big
    A0_in = big_A(A0_in, X0);
end
if size(A1_in,1) < length_big
    A1_in = big_A(A1_in, X1);
end


A0 = fullopA_Hopf(A0_in,X0);
A1 = fullopA_Hopf(A1_in,X1);


As = intval(A0,A1);
Xs = interval_PDE(X0,X1);
XDelta = X0 - X1;


partial_sum = As*third_derHopf_fullop(Xs,XDelta,XDelta);
clear As

ADelta = A0-A1;
A_dagger0 = fullopAdag_Hopf(A_dagger0_in,X0);
A_dagger1 = fullopAdag_Hopf(A_dagger1_in,X1);
A_daggerDelta = A_dagger0 - A_dagger1;
central_bound = 2*ADelta*A_daggerDelta-2*ADelta*second_der_Hopf_fullop(Xs,XDelta)-...
    partial_sum;

clear ADelta A_daggerDelta Xs XDelta partial_sum

%Z1cont = norm(max(abs(Z1_extrema0),abs(Z1_extrema1)))+...
%    1/8*norm(central_bound);


DS0_minusAdag = derivative_Hopf_fullop_minus_Adag(X0);
DS1_minusAdag = derivative_Hopf_fullop_minus_Adag(X1);

Z1_extrema0 = A0*DS0_minusAdag;
Z1_extrema1 = A1*DS1_minusAdag;

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



% function Z1cont = Z1HopfCont(A0, A1, A_dagger0,A_dagger1,X0,X1,G,Es0,Es1)
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

function Z1_extrema = Z1extrema(A0,A1,A_dagger0,A_dagger1,X0,X1,G,Es0,Es1)
% function Z1_extrema = Z1extrema(A0,A1,A_dagger0,A_dagger1,X0,X1,G,Es0,Es1)
% compute Z1 at the extrema of the segment

big = 1 ;
A_dagger_big0 = big_Adagger(A_dagger0,X0);
DH0 =  overop2mat(Hopf_derivative(X0,G,Es0,big));
A_dagger_big1 = big_Adagger(A_dagger1,X1);
DH1 =  overop2mat(Hopf_derivative(X1,G,Es1,big));

max_op = max(abs(A0*(DH0-A_dagger_big0)),abs(A1*(DH1-A_dagger_big1)));

Z1_finite = cnorm(compose_overoperator(max_op,reshape(X0,2*X0.node_space,2*X0.node_time)));
Z1_tail = Z1TailHopf(max(abs(X0),abs(X1)));

Z1_extrema = sum(max(Z1_finite,Z1_tail),2); 
end

function second_der_term = cont_termZ1(A0,A1,Adag0,Adag1,X0,X1,Xs)
% function second_der_term = cont_termZ1(A0,A1,Adag0,Adag1,X0,X1,Xs)
% term dealing with the second derivative

ADelta = A0-A1;
A_dagger_big0 = big_Adagger(Adag0,X0);
A_dagger_big1 = big_Adagger(Adag1,X1);
AdagDelta = A_dagger_big0-A_dagger_big1;

first_term= norm_ADelta_AdagDelta(ADelta, AdagDelta,X0,X1);

second_term = sum(norm_ADelta_D2(A0,A1, X0, X1,Xs),2);

third_term =norm_As_D3(A0, A1,Xs, X0, X1);

second_der_term = first_term+second_term+third_term;

end

function second_term = norm_As_D3(A0, A1, Xs, X0, X1)
% function second_term = norm_As_D3(A0, A1, Xs, X0, X1)
%
% second term: \| A_s D^3 H(xs) xDelta xDelta \|

X0_big = reshape(X0, X0.node_space*3, X0.node_time*3);
X1_big = reshape(X1, X1.node_space*3, X1.node_time*3);
Xs_big = reshape(Xs, Xs.node_space*3, Xs.node_time*3);

DDxHxD = overop2mat(third_derHopf(Xs_big,X0_big-X1_big)); 

A0 = big_A(A0,X0,X0_big);
A1 = big_A(A1,X1,X1_big);

second_term_finite = cnorm( compose_overoperator((A0-A1)*DDxHxD, X0_big));
second_term_infinite = tail_ADelta_D2(X0_big,X1_big);

second_term = sum(max(second_term_finite,second_term_infinite),2);
end

function second_term = norm_ADelta_D2(A0, A1, X0, X1, Xs)
% function second_term = norm_ADelta_D2(A0, A1, X0, X1, Xs)
% \| ADelta D^2 H(xs) xDelta \|

X0_big = reshape(X0, X0.node_space*5, X0.node_time*5);
X1_big = reshape(X1, X1.node_space*5, X1.node_time*5);
Xs_big = reshape(Xs, Xs.node_space*5, Xs.node_time*5);

DDxHxD = second_derHopf(Xs_big,X0_big-X1_big); 

A0 = big_A(A0,X0,X0_big);
A1 = big_A(A1,X1,X1_big);

second_term_finite = cnorm(compose_overoperator(A0-A1, X0_big)*DDxHxD);
second_term_infinite = tail_ADelta_D2(X0,X1);

second_term = max(second_term_finite,second_term_infinite);
end

function tail_ADelta_D2 = tail_ADelta_D2(X0,X1)
% function tail_ADelta_D2 = tail_ADelta_D2(X0,X1)
%
% tail of \| ADelta D^2 H(xs) xDelta \|

tail_ADelta_D2 = zeros(5);

node_space = X0.node_space;

x0 = norm(X0);
x1 = norm(X1);
x = max(x0,x1);
xD = norm(X0-X1);

y = x(4);
z = x(5);
l1 = x(1);
l2 = x(2);
%a = x(3);

yD = xD(4);
zD = xD(5);
l1D = xD(1);
l2D = xD(2);
aD = xD(3);

%y1 = x1(4);
%z1 = x1(5);
l11 = x1(1);
l21 = x1(2);
%a1 = x1(3);

B1 = tailBoundHopf(1,l1, l2, node_space);

a = l1D*l2+l1*l2D;
b = l1D;
c = l11*l21;
d = l11;

f2 = @(k) (k^2*a+ b)/(k^2*c + d);
B2 = max(abs(f2(node_space+1)),abs(a/c));

tail_ADelta_D2(4:5,:)= [ 0,0,0, B2+4*B1*(l1D*y+l1*yD), 0;
    0,0,0,4*B1*(l1*zD+l1D*z),B2+4*B1*(l1D*y+2*l1D*a*z+2*l1*aD*z+l1*yD+2*l1*a*zD)];
end


% function DDDHxDxD = second_der_xs_xD(Xs, V)
% % INPUT
% % Xs    PDE_vector
% % V     PDE_vector
% % OUTPUT
% % DDDHxDxD  overoperator
% 
% K2_y = derivative_operator(2,Xs.vector1D);
% K2_z = derivative_operator(2,squeeze(Xs.vector2D));
% 
% % K2_yvec = -Xs.node_space:Xs.node_space;
% % K2_y = diag(K2_yvec);
% % K2_zvec = repmat(K2_yvec,1+2*Xs.node_time,1);
% % K2_z = diag(reshape(K2_zvec,1,[]));
% 
% ys = Xs.vector1D;
% zs = squeeze(Xs.vector2D);
% l1s = Xs.parameters(1);
% l2s = Xs.parameters(2);
% as = Xs.parameters(3);
% 
% y = V.vector1D;
% z = squeeze(V.vector2D);
% l1 = V.parameters(1);
% l2 = V.parameters(2);
% a = V.parameters(3);
% 
% D11 = zeros(3);
% D12 = zeros(3, 2*Xs.node_space+1);
% D13 = zeros(3, (2*Xs.node_space+1)*(2*Xs.node_time+1));
% 
% D21_1 = -K2_y^4*(ys.'*l2+l2s*y.')+4*1i*K2_y*Convo(ys,y).';
% D21_2 = -K2_y^4*(ys.'*l1+y.'*l1s);
% D21_3 = 0*D21_1;
% 
% D21 = [D21_1,D21_2, D21_3];
% 
% D22 = (-K2_y^4*l2s+K2_y^2)*l1-K2_y^4*l1s*l2+4*1i*K2_y*conv_operator(l1*ys+l1s*y,1);
% 
% zs_conv = zs;
% z_conv = z;
% y_conv = embedding(y, z_conv);%0*z_conv;
% %y_conv((size(y_conv,1)+1)/2,:) = y;
% ys_conv = embedding(ys, zs_conv);%0*zs_conv;
% %ys_conv((size(ys_conv,1)+1)/2,:) = ys;
% 
% zs = reshape(zs_conv,1,[]);
% z = reshape(z_conv,1,[]);
% ys = reshape(ys_conv,1,[]);
% y = reshape(y_conv,1,[]);
% 
% Convo_loc = @(x,y) reshape(Convo(x,y),1,[]).';
% 
% D31_1 = -K2_z*(l2*zs.'+l2s*z.')+4*1i*K2_z*(...
%     Convo_loc(ys_conv,z_conv)+Convo_loc(zs_conv,y_conv)+a*Convo_loc(zs_conv,zs_conv)+2*as*Convo_loc(z_conv,z_conv));
% D31_2 = -K2_z^4*(zs*l1+z*l1s).';
% D31_3 = 4*1i*K2_z*(l1*Convo_loc(zs_conv, zs_conv)+2*l1s*Convo_loc(z_conv,zs_conv));
% 
% D31 = [D31_1, D31_2, D31_3];
% 
% D32 = 4*1i*K2_z*conv_operator(l1s*z+l1*zs,1);
% 
% column_selection = Xs.node_time+1+(0:2*Xs.node_space)*(2*Xs.node_time+1);
% D32 = D32(:,column_selection);
% 
% vec_33 = l1*ys + l1s*y + 2*l1*as*zs+2*l1s*a*zs +2*l1s*as*z;
% D33 = (-K2_z^4*l2s+K2_z^2)*l1-K2_z^4*l1s*l2+4*1i*K2_z*conv_operator(vec_33,1);
% 
% D23 = 0*D32';
% 
% DDDHxDxD = [D11,D12, D13
%     D21, D22, D23
%     D31, D32, D33];
% end



function first_term= norm_ADelta_AdagDelta(ADelta, AdagDelta,X0,X1)
% function first_term= norm_ADelta_AdagDelta(ADelta, AdagDelta,X0,X1)
% 
% \| ADelta A^dagger Delta \|

center = norm( compose_overoperator(ADelta*AdagDelta, X0));

nAtilde = normAtilde(X0,X1);

tail = [0,0,0,nAtilde,nAtilde].';

first_term = max(center, tail);

end

% if holds
% ADelta = ATilde * A1
% we compute here the norm of ATilde

function nAtilde = normAtilde(X0,X1)
% function nAtilde = normAtilde(X0,X1)
%
% \| Atilde \|

l10 = X0.parameters(1);
l20 = X0.parameters(2);
l11 = X1.parameters(1);
l21 = X1.parameters(2);


a = -l21 *l11 + l10*l20;
c = l10*l20;

f= @(x) ( x.^2*a + l11 - l10)./(x.^2*c - l10);

nAtilde = max([abs(c/a), abs(f(X0.node_space+1))]);

end