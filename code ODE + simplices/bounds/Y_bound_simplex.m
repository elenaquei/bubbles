function [Yvector,new_Y,Y_s]=Y_bound_simplex(A0,A1,A2,xBar0,xBar1,xBar2,alpha0,alpha1,alpha2,old_Y0,old_Y1)
% function Yvector=Y_bound_cont(A0,A1,xBar0,xBar1,alpha0,alpha1,old_Y)
%
% INPUT
% A0,A1,A2              complex matrices, approximations of the inverse of the
%                       function of the problem;
% xBar0,xBar1,xBar2     Xi_vector, numerical solutions of the ODE system;
% alpha0,alpha1,alpha2  full_problem
%
% OUTPUT
% Yvector   a standard vector of length = size_vector+size_scalar, the
%           values of this vector come from the Y bound in the associated pdf (we
%           refer to rigorous numerics for analytic solution of D.E..... for the
%           notation, mainly)
% new_Y     a standard vector of length = size_vector+size_scalar, Y_bound
%           (continuation excluded) useful for next run
% Y_s       Y0+ s^2 YDelta, dependency on the stepsize of Y (rough bound),
%           saved as standard vector of three elements [Y0, 0, YDelta]

% tested, working

global nu
global talkative
global use_intlab

% first, the three corners
if nargin>10
    fullY0=old_Y0;
else
    [~,fullY0]=Y_bound(A0,xBar0,alpha0);
end
if nargin>11
    fullY1=old_Y1;
else
    [~,fullY1]=Y_bound(A1,xBar1,alpha1);
end

[~,fullY2]=Y_bound(A2,xBar2,alpha2);
new_Y=fullY2;

Y0=cnorm_Xi_vector(max(max(abs(fullY0),abs(fullY1)),abs(fullY2)),nu);

if any(Y0>1)
    Yvector = Y0;
    Y_s = [];
    error('Bound at an edge is too big')
    return
end
%% 

xBarDelta1_short=xBar0-xBar1;
xBarDelta2_short=xBar0-xBar2;
%
% (!) requires use_intlab anyway (!)


coefs_linear0 = alpha0.scalar_equations.linear_coef;
coefs_linear1 = alpha1.scalar_equations.linear_coef;
coefs_linear2 = alpha2.scalar_equations.linear_coef;

xBarS1_short=interval_Xi(xBar0,xBar1);
xBarS2_short=interval_Xi(xBar0,xBar2);
xBarS_short = interval_Xi(xBarS1_short,xBarS2_short);

coefs_linearS1 = infsup_coefs(coefs_linear0,coefs_linear1);
coefs_linearS2 = infsup_coefs(coefs_linear0,coefs_linear2);
coefs_linearS = infsup_coefs(coefs_linearS1,coefs_linearS2);

xBarS1= reshape_Xi(xBarS1_short,xBarS1_short.nodes*(alpha0.vector_field.deg_vector-1));
xBarS2= reshape_Xi(xBarS2_short,xBarS2_short.nodes*(alpha0.vector_field.deg_vector-1));

xBarDelta1=reshape_Xi(xBarDelta1_short,xBarS1.nodes);
xBarDelta2=reshape_Xi(xBarDelta2_short,xBarS2.nodes);


if alpha0.vector_field.deg_vector>2
    pad=zeros(size(coefs_linearS1{2},1),size(coefs_linearS1{2},2),xBar0.nodes*(alpha0.vector_field.deg_vector-2));
    pad = intval(pad);
    coefs_linearS1{2}=intvalCAT(3,pad,intvalCAT(3,coefs_linearS1{2},pad));
    
    pad=zeros(size(coefs_linearS2{2},1),size(coefs_linearS2{2},2),xBar0.nodes*(alpha0.vector_field.deg_vector-2));
    pad = intval(pad);
    coefs_linearS2{2}=intvalCAT(3,pad,intvalCAT(3,coefs_linearS2{2},pad));
end

alphaS1 = alpha1;
alphaS1.scalar_equations.linear_coef = coefs_linearS1;

alphaS2 = alpha2;
alphaS2.scalar_equations.linear_coef = coefs_linearS2;

alphaS = alpha0;
alphaS.scalar_equations.linear_coef = coefs_linearS;

temp_intlab=use_intlab;
use_intlab=1;

%%
    function Y2 = Adelta_dxF_xDelta(A0, A1, xBarDelta, xBar0)
        
        new_nodes = ((size(DH_xs2,1)-xBar0.size_scalar)/xBar0.size_vector - 1)/2;
        xBarDelta_long = reshape_Xi(xBarDelta,new_nodes);
        
        full_DH_xs = DH_xs2*Xi_vec2vec(xBarDelta_long);
        
        Adiff=extend_approximate_inverse(A1,xBar0.size_scalar,...
            xBar0.size_vector,xBar0.nodes,new_nodes)-...
            extend_approximate_inverse(A0,xBar0.size_scalar,...
            xBar0.size_vector,xBar0.nodes,new_nodes);
        Y2=vec2Xi_vec(Adiff*full_DH_xs,...
            xBarDelta_long.size_scalar,xBarDelta_long.size_vector,xBarDelta_long.nodes);
    end
DH_xs2 = derivative_to_matrix(derivative(alphaS,xBarS_short));
Y2_d1d1 = Adelta_dxF_xDelta(A0, A1, xBarDelta1, xBar0);
Y2_d2d2 = Adelta_dxF_xDelta(A0, A2, xBarDelta2, xBar0);
% order doesn't matter, as long as we coer all options
Y2_d1d2 = Adelta_dxF_xDelta(A0, A1, xBarDelta2, xBar0); 
Y2_d2d1 = Adelta_dxF_xDelta(A0, A2, xBarDelta1, xBar0);

% Y2 = 2*cnorm_Xi_vector(max(max(Y2_d1d1,Y2_d2d2),max(Y2_d1d2,Y2_d2d1)),nu);

%%

As = interpolation(A0, A1, A2);
% 4 options for DDH_xdeltai_xDeltaj
DDH_xD1_xD1 = Function_directional_second_derivative(alpha1,xBarS_short,xBarDelta1_short,xBarDelta1_short);
DDH_xD1_xD2 = Function_directional_second_derivative(alpha1,xBarS_short,xBarDelta1_short,xBarDelta2_short);
% only 3 options for simmetry
DDH_xD2_xD2 = Function_directional_second_derivative(alpha1,xBarS_short,xBarDelta2_short,xBarDelta2_short);

% expand A0 to size of vec(DDH_xs)
nodes_new=DDH_xD2_xD2.nodes;
As_big=extend_approximate_inverse(As,xBar0.size_scalar,xBar0.size_vector,xBar0.nodes,nodes_new);

Y3_d1d1 = vec2Xi_vec(As_big*Xi_vec2vec(DDH_xD1_xD1),...
    xBarDelta1.size_scalar,xBarDelta1.size_vector,nodes_new);
Y3_d1d2 = vec2Xi_vec(As_big*Xi_vec2vec(DDH_xD1_xD2),...
    xBarDelta1.size_scalar,xBarDelta1.size_vector,nodes_new);
Y3_d2d1 = Y3_d1d2;
Y3_d2d2 = vec2Xi_vec(As_big*Xi_vec2vec(DDH_xD2_xD2),...
    xBarDelta1.size_scalar,xBarDelta1.size_vector,nodes_new);

% Y3_vec = max(sup(abs(Y3_d1d1)), max(abs(Y3_d1d2), abs(Y3_d2d2)));

% Y3=cnorm_Xi_vector(Y3_vec,nu);

%%

Hdiff1_xs = intval(0*Xi_vec2vec(xBar0));
Hdiff2_xs = intval(0*Xi_vec2vec(xBar0));
DxHdiff1_xsxD1 = 0*Xi_vec2vec(xBar0);
DxHdiff1_xsxD2 = 0*Xi_vec2vec(xBar0);
DxHdiff2_xsxD1 = 0*Xi_vec2vec(xBar0);
DxHdiff2_xsxD2 = 0*Xi_vec2vec(xBar0);
% addition for the new continuation equation and linear in s scalar eqs
for i = 1:alpha1.scalar_equations.number_equations_lin
    coef_1 = extract_lin_coef(alpha1.scalar_equations, i);
    coef_0 =  extract_lin_coef(alpha0.scalar_equations, i);
    c_1 = alpha1.scalar_equations.linear_coef{3}(i);
    c_0 = alpha0.scalar_equations.linear_coef{3}(i);
    coef_Delta1 = coef_1-coef_0;
    c_Delta1 = c_1-c_0;
    coef_2 = extract_lin_coef(alpha2.scalar_equations, i);
    c_2 = alpha2.scalar_equations.linear_coef{3}(i);
    coef_Delta2 = coef_2-coef_0;
    c_Delta2 = c_2-c_0;
    %DDH_xs.scalar(i) = 2*dot(Xi_vec2vec(xBarDelta_short),coef_Delta);
    Hdiff1_xs(i) = dot(Xi_vec2vec(xBarS_short),coef_Delta1)+c_Delta1;
    Hdiff2_xs(i) = dot(Xi_vec2vec(xBarS_short),coef_Delta2)+c_Delta2;
    
    DxHdiff1_xsxD1(i) = dot(Xi_vec2vec(xBarDelta1_short),coef_Delta1);
    DxHdiff1_xsxD2(i) = dot(Xi_vec2vec(xBarDelta1_short),coef_Delta2);
    DxHdiff2_xsxD1(i) = dot(Xi_vec2vec(xBarDelta2_short),coef_Delta1);
    DxHdiff2_xsxD2(i) = dot(Xi_vec2vec(xBarDelta2_short),coef_Delta2);
end

Adiff1_small = A1-A0;
Adiff2_small = A2-A0;

Y4_d1d1 = vec2Xi_vec(Adiff1_small*Hdiff1_xs,...
    xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);
Y4_d1d2 = vec2Xi_vec(Adiff1_small*Hdiff2_xs,...
    xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);
Y4_d2d1 = vec2Xi_vec(Adiff2_small*Hdiff1_xs,...
    xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);
Y4_d2d2 = vec2Xi_vec(Adiff2_small*Hdiff2_xs,...
    xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);

Y4_d1d1 = reshape_Xi(Y4_d1d1,Y2_d1d1.nodes);
Y4_d1d2 = reshape_Xi(Y4_d1d2,Y2_d1d1.nodes);
Y4_d2d1 = reshape_Xi(Y4_d2d1,Y2_d1d1.nodes);
Y4_d2d2 = reshape_Xi(Y4_d2d2,Y2_d1d1.nodes);

%Y4 = 2* cnorm_Xi_vector(max(max(abs(Y4_d1d1),abs(Y4_d1d2)),max(abs(Y4_d2d1),abs(Y4_d2d2))),nu);

Y5_d1d1 = vec2Xi_vec(As*DxHdiff1_xsxD1,...
    xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);
Y5_d1d2 = vec2Xi_vec(As*DxHdiff1_xsxD2,...
    xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);
Y5_d2d1 = vec2Xi_vec(As*DxHdiff2_xsxD1,...
    xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);
Y5_d2d2 = vec2Xi_vec(As*DxHdiff2_xsxD2,...
    xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);

Y5_d1d1 = reshape_Xi(Y5_d1d1,Y2_d1d1.nodes);
Y5_d1d2 = reshape_Xi(Y5_d1d2,Y2_d1d1.nodes);
Y5_d2d1 = reshape_Xi(Y5_d2d1,Y2_d1d1.nodes);
Y5_d2d2 = reshape_Xi(Y5_d2d2,Y2_d1d1.nodes);

% Y5 = cnorm_Xi_vector(max(max(abs(Y5_d1d1),abs(Y5_d1d2)),max(abs(Y5_d2d1),abs(Y5_d2d2))),nu);


% merging
Ys_d1d1 = Y2_d1d1 + Y3_d1d1 + 2*Y4_d1d1 + 2*Y5_d1d1;
Ys_d1d2 = Y2_d1d2 + Y3_d1d2 + 2*Y4_d1d2 + 2*Y5_d1d2;
Ys_d2d1 = Y2_d2d1 + Y3_d2d1 + 2*Y4_d2d1 + 2*Y5_d2d1;
Ys_d2d2 = Y2_d2d2 + Y3_d2d2 + 2*Y4_d2d2 + 2*Y5_d2d2;
YDelta = 1/8 * cnorm_Xi_vector(max(max(abs(Ys_d1d1),abs(Ys_d1d2)),max(abs(Ys_d2d1),abs(Ys_d2d2))),nu);


%%

use_intlab=temp_intlab;

%YDelta=(Y2+Y3+Y4+Y5)/8;
Y_s=[0*Y0, YDelta,0*Y0,Y0];

Yvector=Y0+YDelta;

%%
if any(Yvector>1)
    fprintf('Y_simplex computed, %d\n',max(Yvector));
    %error('Y_cont is bigger than 1, no interval found')
elseif any(isnan(Yvector))
    error('Y_simplex is NaN, big troubles!')
elseif talkative>0
%    fprintf('Y_cont computed, %d\n',Yvector);
%elseif talkative>0
    fprintf('Y_simplex computed, %d\n',max(Yvector));
end

end



%% help functions

function coefs_linearS = infsup_coefs(coefs_linear0,coefs_linear1)

coefs_linearS = cell(3,1);
for i = 1:3
    if ~isintval(coefs_linear1{i})
        min_coef1=min(real(coefs_linear1{i}),real(coefs_linear0{i}))+...
            1i*min(imag(coefs_linear1{i}),imag(coefs_linear0{i}));
        max_coef1=max(real(coefs_linear1{i}),real(coefs_linear0{i}))+...
            1i*max(imag(coefs_linear1{i}),imag(coefs_linear0{i}));
    else
        min_coef1=min(inf(real(coefs_linear1{i})),inf(real(coefs_linear0{i})))+...
            1i*min(inf(imag(coefs_linear1{i})),inf(imag(coefs_linear0{i})));
        max_coef1=max(sup(real(coefs_linear1{i})),sup(real(coefs_linear0{i})))+...
            1i*max(sup(imag(coefs_linear1{i})),sup(imag(coefs_linear0{i})));
    end
    coefs_linearS{i}=infsup(min_coef1,max_coef1);
end
coefs_linearS{1}(end)=0;
coefs_linearS{2}(end,:,:)=0;
coefs_linearS{3}(end)=0;
end
