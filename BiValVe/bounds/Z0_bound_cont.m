function [Z0_vector,Z0_s]=Z0_bound_cont(DH0,DH1,A0,A1,xBar0,~)
% function Z0_vector=Z0_bound(DH0,DH1,A0,A1,xBar0,xBar1)
% 
% INPUT:
% DH0,DH1      complex matrix, can also be called Adagger, is the approximate 
%              derivative in the numerical solution xBar (given);
% A0,A1        complex matrix, numerical inverse of DHm;
% xBar0,xBar1  Xi_vector, numerical solution. Theoretically, there is no need
%              of having the numerical solution in this function, the goal is to have
%              the data if the problem inherent the solution (the size, the number
%              of nodes, ecc.)
%
% OUTPUT:
% Z0_vector   a standard vector of length = size_vector+size_scalar, the
%             values of this vector come from the Z0 bound in the associated pdf (we
%             refer to rigorous numerics for analytic solution od D.E..... for the
%             notation, mainly)
% Z0_s        Z0+ s^2 Z0Delta, dependency on the stepsize of Z0 (rough bound),
%             saved as standard vector of three elements [0, Z0Delta,0,Z0]


global talkative
global use_intlab


if use_intlab
    B0=abs(intval(eye(size(A0)))-A0*DH0);
    B1=abs(intval(eye(size(A1)))-A1*DH1);
else
    B0=eye(size(A0))-A0*DH0;
    B1=eye(size(A1))-A1*DH1;
end

A_d=A1-A0;
DHd= DH1-DH0;

if any(size(A_d)~=size(DHd))
    disp('problem in sizes\n');
end

if use_intlab
    Z0_vector=norm_Ximat(max(sup(abs(B0)),sup(abs(B1)))+abs(1/8*A_d*DHd),xBar0);
    Z0=norm_Ximat(max(sup(abs(B0)),sup(abs(B1))),xBar0);
    Z0Delta=norm_Ximat(abs(1/8*A_d*DHd),xBar0);
else
    Z0_vector=norm_Ximat(max(abs(B0),abs(B1))+abs(1/8*A_d*DHd),xBar0);
    Z0=norm_Ximat(max(abs(B0),abs(B1)),xBar0);
    Z0Delta=norm_Ximat(abs(1/8*A_d*DHd),xBar0);
end

Z0_s=[0*Z0, Z0Delta,0*Z0,Z0];

% error handling
if any(Z0_vector>1)
    fprintf('Z0_cont computed, %d\n',max(Z0_vector));
    %error('Z0_cont is bigger than 1, no interval found')
elseif talkative>1
%    fprintf('Z0_cont computed, %d\n',Z0_vector);
%elseif talkative>0
    fprintf('Z0_cont computed, %d\n',max(Z0_vector));
end
