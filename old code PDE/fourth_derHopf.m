function DDDHxDxD = fourth_derHopf(X, V, W, Z, flagFullSize)
% INPUT
% X    PDE_vector
% V     PDE_vector
% OUTPUT
% DDDHxDxD  overoperator

if nargin == 2
    W = V;
    Z = V;
    flagFullSize =0;
elseif nargin == 3
    if isfloat(W)
        flagFullSize = W;
        W = V;
        Z = V;
    else
        flagFullSize =0;
        Z=W;
    end
elseif nargin == 4
    if isfloat(Z)
        flagFullSize = Z;
        Z = W;
    else
        flagFullSize = 0;
    end
end

if flagFullSize
    X = reshape(X, X.node_space*2, X.node_time*2);
    V = reshape(V, V.node_space*2, V.node_time*2);
    W = reshape(W, W.node_space*2, W.node_time*2);
    Z = reshape(Z, Z.node_space*2, Z.node_time*2);
end

%K2_y = derivative_operator(2,X.vector1D);
K2_z = derivative_operator(2,squeeze(X.vector2D));

%ys = Xs.vector1D;
%z = squeeze(X.vector2D);
%l1 = X.parameters(1);
%l2 = X.parameters(2);
%a = X.parameters(3);

%y1 = V.vector1D;
z1 = squeeze(V.vector2D);
l11 = V.parameters(1);
%l21 = V.parameters(2);
a1 = V.parameters(3);

%y2 = V.vector1D;
z2 = squeeze(V.vector2D);
l12 = V.parameters(1);
%l22 = V.parameters(2);
a2 = V.parameters(3);

%y3 = V.vector1D;
z3 = squeeze(V.vector2D);
l13 = V.parameters(1);
%l23 = V.parameters(2);
a3 = V.parameters(3);

D11 = zeros(3);
D12 = zeros(3, 2*X.node_space+1);
D13 = zeros(3, (2*X.node_space+1)*(2*X.node_time+1));

D21 = zeros(2*X.node_space+1,3);
D22 = zeros(2*X.node_space+1);
D23 = zeros(2*X.node_space+1, (2*X.node_space+1)*(2*X.node_time+1));

%z_conv = z;
z1_conv = z1;
z2_conv = z2;
z3_conv = z3;
%y_conv = embedding(y, z1_conv);
%y1_conv = embedding(y1, z1_conv);
%y2_conv = embedding(y2, z1_conv);
%y3_conv = embedding(y3, z1_conv);

%z = reshape(z_conv,1,[]);
z1 = reshape(z1_conv,1,[]);
z2 = reshape(z2_conv,1,[]);
z3 = reshape(z3_conv,1,[]);
%ys = reshape(ys_conv,1,[]);
%y1 = reshape(y1_conv,1,[]);
%y2 = reshape(y2_conv,1,[]);
%y3 = reshape(y3_conv,1,[]);

Convo_loc = @(x,y) reshape(Convo(x,y),1,[]).';

D31_1 = 2 * a3 *1i* K2_z *Convo_loc(z1,z2) + ...
     2 * a2* 1i* K2_z* Convo_loc(z1,z3) + 2 * a1 *1i *K2_z* Convo_loc(z2, z3);
%24*a1*1i*K2_z*Convo_loc(z1_conv,z1_conv);
D31_2 = 0*D31_1;
D31_3 = 2* 1i *K2_z* l13* Convo_loc(z1,z2) + ...
    2 *1i* K2_z *l12* Convo_loc(z1,z3) + 2*1i* K2_z *l11 * Convo_loc(z2,z3);
%24*l11*1i*K2_z*Convo_loc(z1_conv,z1_conv);

D31 = [D31_1, D31_2, D31_3];

D32 = 0*K2_z;

column_selection = X.node_time+1+(0:2*X.node_space)*(2*X.node_time+1);
D32 = D32(:,column_selection);

z1 = conv_operator(z1,1);
z2 = conv_operator(z2,1);
z3 = conv_operator(z3,1);

D33 = 2*1i* K2_z *(l13 *( a2* z1 +  a1 *z2) +... 
 a3* ( l12 *z1 +l11 *z2) + ...
 ( a2 *l11 +a1* l12)* z3);
%40*1i*K2_z*conv_operator( l11*a1*z1,1);

if isa(D31,'intval')
    D1 = intval([D11,D12, D13]);
    D2 = intval([D21, D22, D23]);
else
    D1 = [D11,D12, D13];
    D2 = [D21, D22, D23];
end

DDDHxDxD = [D1
    D2
    D31, D32, D33];
end