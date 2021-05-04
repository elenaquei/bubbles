function Norm = operator_norm(A)
% function Norm = operator_norm(A)
% 
% INPUT
% A         4D operator
% GLOBAL 
% nu        (nu1,n2)>(1,1)
% OUTPUT
% Norm      norm of the operator A
global nu

span_k1 = (size(A,1)-1)/2;
span_k2 = (size(A,2)-1)/2;
span_l1 = (size(A,3)-1)/2;
span_l2 = (size(A,4)-1)/2;

NU_K1 = nu(1).^(abs(-span_k1:span_k1));
NU_K2 = nu(2).^(abs(-span_k2:span_k2));

sup_mat = zeros(2*span_l1+1,2*span_l2+1);
for l1_tilde =-span_l1:span_l1
    for l2_tilde = -span_l2:span_l2
        A_loc = squeeze(A(:,:,l1_tilde,l2_tilde));
        sup_mat = 1/(nu(1)^(abs(l1_tilde))+nu(2)^(abs(l2_tilde)))* NU_K1 *abs(A_loc) *NU_K2';
    end
end

Norm = max(max(sup_mat));