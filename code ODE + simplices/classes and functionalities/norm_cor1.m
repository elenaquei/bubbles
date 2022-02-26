function norm_Aij=norm_cor1(A_ij,F,alpha_1)
% INPUT
% Aij        linear operator from l1_nu to l1_nu (considered as a square 
%            matrix of dimension 2F+1
% F          number of nodes
% alpha_1    period of the solution
% nu         coefficient of the nu-norm
% OUTPUT
% norm_Aij   application of corollary 1 of Hungria, Lessard, James
%
% based on the formula
%  ||A_ij|| \leq max (K_ij, 1/(F\alpha_1)
% 
% with
%  K_ij= max_{|n|<=F} 1/(nu^abs(n)) sum_{k=-F}^F |(A_ij)_{kn}| nu^{abs(k)}

global use_intlab

K_ij=InfNuNorm_mat(A_ij);

if ~use_intlab
    norm_Aij=max(K_ij, abs(1/(2*pi*F*alpha_1)) );
else
    
    norm_Aij=max(sup(K_ij), sup(abs(1/(intval(2*pi*F)*alpha_1))) );
end

return