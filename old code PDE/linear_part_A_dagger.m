function A_dagger_op = linear_part_A_dagger(l1,l2,node_space, node_time)

A_dagger_11 = zeros(3,3);
A_dagger_12 = zeros(3,2*node_space+1);
A_dagger_13 = zeros(3, (2*node_space+1)*(2*node_time+1));

A_dagger_21 = zeros(2*node_space+1,3);
A_dagger_31 = zeros((2*node_space+1)*(2*node_time+1),3);

A_dagger_23 = zeros(2*node_space+1,(2*node_space+1)*(2*node_time+1));
A_dagger_32 = A_dagger_23.';

K_2 = -node_space:node_space;
A_dagger_22 = diag((-K_2.^4*l2+K_2.^2));

K_1 = -node_time:node_time;
K_1 = K_1.';
K_1 = repmat(K_1,1,1+2*node_space);
K_2 = repmat(K_2,1+2*node_time,1);

A_dagger_33 = diag(reshape((-1i*K_1-K_2.^4*l1*l2+K_2.^2*l1),1,[]));

A_dagger_op = overoperator(A_dagger_11, A_dagger_12, A_dagger_13,...
    A_dagger_21, A_dagger_22, A_dagger_23,...
    A_dagger_31, A_dagger_32, A_dagger_33);