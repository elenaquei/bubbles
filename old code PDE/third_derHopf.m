function DDDHxDxD = third_derHopf(X, V, W,flagFullSize,flagOp)
% INPUT
% X     PDE_vector
% V     PDE_vector, first direction
% w     PDE_vector, second direction
% flag  bool
% OUTPUT
% DDDHxDxD  overoperator

if nargin == 2
    W = V;
    flagFullSize = 0;
elseif nargin == 3
    if isfloat(W)
        flagFullSize = W;
        W = V;
    else
        flagFullSize = 0;
    end
end
if nargin<5
    flagOp = 0;
end

if flagFullSize
    X = reshape(X, X.node_space*2, X.node_time*2);
    V = reshape(V, V.node_space*2, V.node_time*2);
    W = reshape(W, W.node_space*2, W.node_time*2);
end

K2_y = derivative_operator(2,X.vector1D);
K2_z = derivative_operator(2,squeeze(X.vector2D));

%ys = Xs.vector1D;
z = squeeze(X.vector2D);
l1 = X.parameters(1);
%l2s = Xs.parameters(2);
a = X.parameters(3);

y1 = V.vector1D;
z1 = squeeze(V.vector2D);
l11 = V.parameters(1);
l21 = V.parameters(2);
a1 = V.parameters(3);

y2 = W.vector1D;
z2 = squeeze(W.vector2D);
l12 = W.parameters(1);
l22 = W.parameters(2);
a2 = W.parameters(3);

D11 = zeros(3);
D12 = zeros(3, 2*X.node_space+1);
D13 = zeros(3, (2*X.node_space+1)*(2*X.node_time+1));

D21_1 = 0 * y1.';
D21_2 = 0 * y1.';
D21_3 = 0 * D21_1;

D21 = [D21_1,D21_2, D21_3];

D22 = 0 * K2_y;

z_conv = z;
z1_conv = z1;
z2_conv = z2;
y1_conv = embedding(y1, z1_conv);
y2_conv = embedding(y2, z2_conv);

z = reshape(z_conv,1,[]);
z1 = reshape(z1_conv,1,[]);
z2 = reshape(z2_conv,1,[]);
%ys = reshape(ys_conv,1,[]);
y1 = reshape(y1_conv,1,[]);
y2 = reshape(y2_conv,1,[]);

Convo_loc = @(x,y) reshape(Convo(x,y),1,[]).';

D31_1 = - K2_z^4 *(l22 *z1.'+l21*z2.') + 2 *1i* K2_z*( Convo_loc(y2, z1) ... 
 + Convo_loc(y1,z2) +  a2 * Convo_loc(z, z1)...
 + a1 * Convo_loc(z,z2) + a * Convo_loc(z1,z2) );
%-2*K2_z*l21*z1.'+4*1i*K2_z*(4*a1*Convo_loc(z1_conv,z_conv)+...
%    2*Convo_loc(z1_conv,y_conv)+2*a*Convo_loc(z1_conv,z1_conv));
D31_2 = - K2_z^4*(l11*z2.'+l12*z1.');
D31_3 = 1i * K2_z * ( l12 * Convo_loc(z, z1) + ...
     l11 * Convo_loc( z,z2) + l1 * Convo_loc( z1, z2) );
%8*1i*K2_z*(2*l11*Convo_loc(z_conv, z1_conv)+l1*Convo_loc(z1_conv,z1_conv));

D31 = [D31_1, D31_2, D31_3];

D32 = 2 *1i* K2_z*(l12* conv_operator(z1,1) +...
    l11 *conv_operator(z2,1));

column_selection = X.node_time+1+(0:2*X.node_space)*(2*X.node_time+1);
D32 = D32(:,column_selection);

D33 = -K2_z^4 * (l11* l22 + l12 * l21) + 2* 1i *K2_z*( ... 
    l11* conv_operator(y2,1) + l12* conv_operator(y1,1) + ...
    a1 * l12 * conv_operator(z,1) + a  * l12 * conv_operator(z1,1) + ...
    a2 * l11 * conv_operator(z,1) + a2 * l1 *conv_operator(z1,1) + ...
    a1 * l1 * conv_operator(z2,1)+ a* l11 * conv_operator(z2,1) );
%-2*K2_z^4*l11*l21 + 4*1i*K2_z*conv_operator(l11*(2*y1+4*a*z1+4*a1*z)+ 4*l1*a1*z1,1);

D23 = 0*D32';

DDDHxDxD = [D11,D12, D13
    D21, D22, D23
    D31, D32, D33];
if size(DDDHxDxD,1)>10^3
    DDDHxDxD = sparse(DDDHxDxD);
end
if flagOp
    DDDHxDxD =compose_overoperator (DDDHxDxD, X.size_real, X.node_space, X.node_time);
end

end