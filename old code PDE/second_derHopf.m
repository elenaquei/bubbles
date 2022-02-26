function DDHxD = second_derHopf(Xs, V, flagFullSize,flagOp)
% INPUT
% Xs    PDE_vector
% V     PDE_vector
% OUTPUT
% DDDHxDxD  overoperator
if nargin>2 && flagFullSize
    Xs = reshape(Xs, Xs.node_space*2, Xs.node_time*2);
    V = reshape(V, V.node_space*2, V.node_time*2);
end
if nargin<4
    flagOp = 0;
end

K2_y = derivative_operator(2,Xs.vector1D);
K2_z = derivative_operator(2,squeeze(Xs.vector2D));

% K2_yvec = -Xs.node_space:Xs.node_space;
% K2_y = diag(K2_yvec);
% K2_zvec = repmat(K2_yvec,1+2*Xs.node_time,1);
% K2_z = diag(reshape(K2_zvec,1,[]));

ys = Xs.vector1D;
zs = squeeze(Xs.vector2D);
nus = Xs.parameters(1);
lambdas = Xs.parameters(2);
as = Xs.parameters(3);

y = V.vector1D;
z = squeeze(V.vector2D);
nu = V.parameters(1);
lambda = V.parameters(2);
a = V.parameters(3);

D11 = zeros(3);
D12 = zeros(3, 2*Xs.node_space+1);
D13 = zeros(3, (2*Xs.node_space+1)*(2*Xs.node_time+1));

D21_2 = -K2_y^4*y.';
D21_1 = 0*D21_2;
D21_3 = 0*D21_2;

D21 = [D21_1,D21_2, D21_3];

D22 = -K2_y^4*nu + 2 *1i*K2_y*conv_operator(y,1);

zs_conv = zs;
z_conv = z;
y_conv = embedding(y, z_conv);%0*z_conv;
%y_conv((size(y_conv,1)+1)/2,:) = y;
ys_conv = embedding(ys, zs_conv);%0*zs_conv;
%ys_conv((size(ys_conv,1)+1)/2,:) = ys;

zs = reshape(zs_conv,1,[]);
z = reshape(z_conv,1,[]);
ys = reshape(ys_conv,1,[]);
y = reshape(y_conv,1,[]);

Convo_loc = @(x,y) reshape(Convo(x,y),1,[]).';


% -K4 l21 z + 2 i k y1 z + 
% a1 i k z^2 + (K2 - K4 l2 + i k (2 y + 2 a z)) z1
D31_1 =  -K2_z^4*(zs*nu+z*nus).' + K2_z * z.'+ 1i*K2_z*(...
    2 * Convo_loc(ys_conv,z_conv)+ 2 * Convo_loc(zs_conv,y_conv)+...
    a*Convo_loc(zs_conv,zs_conv)+2*as*Convo_loc(zs_conv,z_conv));
D31_2 = -K2_z^4*(zs*nu+z*nus).';
D31_3 = 1i*K2_z*(lambda*Convo_loc(zs_conv, zs_conv)+2*lambdas*Convo_loc(z_conv,zs_conv));

D31 = [D31_1, D31_2, D31_3];

D32 = 2*1i*K2_z*conv_operator(lambdas*z+lambda*zs,1);

column_selection = Xs.node_time+1+(0:2*Xs.node_space)*(2*Xs.node_time+1);
D32 = D32(:,column_selection);

vec_33 = lambda*ys + lambdas*y + lambda*as*zs + lambdas*a*zs + lambdas*as*z;
D33 = -K2_z^4*lambdas*nu+K2_z^2*lambda-K2_z^4*nus*lambda+2*1i*K2_z*conv_operator(vec_33,1);

D23 = 0*D32';

DDHxD = [D11,D12, D13
    D21, D22, D23
    D31, D32, D33];

if size(DDHxD,1)>10^3
    DDHxD = sparse(DDHxD);
end

if flagOp
    DDHxD =compose_overoperator (DDHxD, Xs.size_real, Xs.node_space, Xs.node_time);
end
end