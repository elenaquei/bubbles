node_space = 2;
node_time = 3;

new_node_space = 3;
new_node_time = 4;

size_real = 1;

A = rand(size_real+2*node_space+1+(2*node_space+1)*(2*node_time+1));
A_op= compose_overoperator(A, size_real, node_space, node_time);

index_C = 1:size_real;
index_l1 = size_real+(1:2*node_space+1);
index_L1 = size_real+2*node_space+1 +(1:(2*node_space+1)*(2*node_time+1));

A_CC = A(index_C,index_C);
Aop_CC_small = operator(A_CC,'Cn','Cn');
A_Cl1 = A(index_l1,index_C);
Aop_Cl1_small = operator(A_Cl1,'Cn','l1');
A_CL1 = A(index_L1,index_C);
Aop_CL1_small = operator(A_CL1,'Cn','L1');
A_l1C = A(index_C,index_l1);
Aop_l1C_small = operator(A_l1C,'l1','Cn');
A_l1l1 = A(index_l1, index_l1);
Aop_l1l1_small = operator(A_l1l1,'l1','l1');
A_l1L1 = A(index_L1, index_l1);
Aop_l1L1_small = operator(A_l1L1,'l1','L1');
A_L1C = A(index_C, index_L1);
Aop_L1C_small = operator(A_L1C,'L1','Cn');
A_L1l1 = A(index_l1, index_L1);
Aop_L1l1_small = operator(A_L1l1,'L1','l1');
A_L1L1 = A(index_L1, index_L1);
Aop_L1L1_small = operator(A_L1L1,'L1','L1');

Aop_Cl1 = expand(Aop_Cl1_small,new_node_space, node_space, new_node_time);
Aop_CL1 = expand(Aop_CL1_small,new_node_space, node_space, new_node_time);
Aop_l1C = expand(Aop_l1C_small,new_node_space, node_space, new_node_time);
Aop_l1l1 = expand(Aop_l1l1_small,new_node_space, node_space, new_node_time);
Aop_l1L1 = expand(Aop_l1L1_small,new_node_space, node_space, new_node_time);
Aop_L1C = expand(Aop_L1C_small,new_node_space, node_space, new_node_time);
Aop_L1l1 = expand(Aop_L1l1_small,new_node_space, node_space, new_node_time);
Aop_L1L1 = expand(Aop_L1L1_small,new_node_space, node_space, new_node_time);

A_big_man = [Aop_CC_small.matrix, Aop_l1C.matrix, Aop_L1C.matrix;
    Aop_Cl1.matrix, Aop_l1l1.matrix, Aop_L1l1.matrix;
    Aop_CL1.matrix, Aop_l1L1.matrix, Aop_L1L1.matrix];

A_big_aut = reshape(A_op,new_node_space, new_node_time);

size(overop2mat(A_big_aut))
size(A_big_man)
spy(A_big_man);
spy(overop2mat(A_big_aut)-A_big_man)

% overoperator reshape works correctly