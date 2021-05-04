function A_op = linear_part_A(l1,l2,node_space, node_time)
% function A_op = linear_part_A(l1,l2,node_space, node_time)
% 
% INPUT
% l1, l2                    parameters
% node_space,node_time      number of requested nodes
% OUTPUT 
% A_op                      overoperator with tail of A
A11 = zeros(3,3);
A12 = zeros(3,2*node_space+1);
A13 = zeros(3, (2*node_space+1)*(2*node_time+1));

A21 = zeros(2*node_space+1,3);
A31 = zeros((2*node_space+1)*(2*node_time+1),3);

A23 = zeros(2*node_space+1,(2*node_space+1)*(2*node_time+1));
A32 = A23.';

K_2 = -node_space:node_space;
A22 = diag((-K_2.^4*l2+K_2.^2).^-1);
A22(node_space+1,node_space+1) = 0;

J = -node_time:node_time;
J = J.';
J = repmat(J,1,1+2*node_space);
K_2 = repmat(K_2,1+2*node_time,1);

A33 = diag(reshape((-1i*J-K_2.^4*l1*l2+K_2.^2*l1).^-1,1,[]));
center_33 = ((2*node_space+1)*(2*node_time+1) +1)/2;
A33(center_33,center_33) = 0;

A_op = overoperator(A11, A12, A13,...
    A21, A22, A23,...
    A31, A32, A33);