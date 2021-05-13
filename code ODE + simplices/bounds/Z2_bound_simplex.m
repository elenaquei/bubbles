function [Z2vector,Z2_s]=Z2_bound_simplex(A0,A1,A2,xBar0,xBar1,xBar2,alpha)
% function Z2vector=Z2_bound_simplex(A0,A1,A2,xBar0,xBar1,xBar2,alpha)
% 
% INPUT:
% DHm     complex matrix, can also be called Adagger, is the approximate 
%         derivative in the numerical solution xBar (given);
% A       complex matrix, numerical inverse of DHm;
% xBar    Xi_vector, numerical solution. Theoretically, there is no need
%         of having the numerical solution in this function, the goal is to have
%         the data if the problem inherent the solution (the size, the number
%         of nodes, ecc.);
% alpha  coefs, exact data of the problem.
%
% OUTPUT:
% Z2_vector   a standard vector of length = size_vector+size_scalar, the
%             values of this vector come from the Z1 bound in the associated pdf (we
%             refer to rigorous numerics for analytic solution od D.E..... for the
%             notation, mainly)
% Z2_s        s^2 Z2, dependency on the stepsize of Z2 (rough bound),
%             saved as standard vector of three elements [0, Z2, 0, 0]
%
% In this function the coefficient of r.^2 in the radii polynomials is
% calculated.
% REMARK:
% A maximum radius need to be fixed in advance in order to
% simplify the calculations, it will be set as a global variable in the
% main, calling it RAD_MAX. If this is not done, the maximum radius is
% fixed by default to 10.^-3.
%
% It holds
%     Z2vector(i) = sum_{j=1}^{M+N} ||A_ij||  || [DDH(x+z)bc]_j ||
% and
%     || [DDH(x+z)bc]_j || <= sum_d sum_k sum_l | alpha_d,j d_k (d-e_k)_l |
%       ( |x|_c +r )^(d-e_k-e_l)

global RAD_MAX
global nu
%global use_intlab
global talkative

if isempty(RAD_MAX)
    RAD_MAX=1;
end
%RAD_MAX=10^-5;

xBar=max(max(cnorm_Xi_vector(xBar1,nu),cnorm_Xi_vector(xBar0,nu)),...
    cnorm_Xi_vector(xBar2,nu));

if isintval(A0)
    abs_A0=abs(A0);
    abs_A1=abs(A1);
    abs_A2=abs(A2);
    A_ij=max(component_matrix_norm(abs_A0.sup,xBar1),...
        max(component_matrix_norm(abs_A1.sup,xBar1),...
        component_matrix_norm(abs_A2.sup,xBar1)));    
else
    A_ij=max(component_matrix_norm(A2,xBar1),...
        max(component_matrix_norm(A0,xBar1),...
        component_matrix_norm(A1,xBar1))); 
end


DDH=Function_second_derivative(alpha,xBar1,RAD_MAX,xBar); 
% remark: the output WILL/SHOULD be a Xi_vector because only the convolution terms are
% left. BUT by now it's just a plain and normal vector. 


%if use_intlab
    Z2vector=A_ij*DDH;
%else
%    Z2vector=A_ij*DDH;
%end
Z2_s=[0*Z2vector,Z2vector,0*Z2vector,0*Z2vector];

if talkative>1
    fprintf('Z2_simplex computed, %d\n',max(Z2vector));
elseif talkative>0
    fprintf('Z2_simplex computed, %d\n',max(Z2vector));
end