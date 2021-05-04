function Y_bound=YboundHopf(A,X,G,const_G, Es,const_Es)
% 
% INPUT
% A                 matrix
% X                 PDE_vector
% G, const_G        linear equation
% Es, const_Es      continuation equation
% PUTPUT
% Ybound            vector Y bound

%Y_vec = matrix_AHopf(X,G,Es,1)*PDE_vec2vec(evaluate_HopfKS(X,G,const_G,Es,const_Es,1)).';

y = evaluate_HopfKS(X,G,const_G,Es,const_Es,1);
op_A =  compose_overoperator(A,X.size_real, 2*X.node_space,2*X.node_time);

Y_arg = op_A * y ;

Y_bound = norm( Y_arg );
%Y_PDE = PDE_vector(Y_vec,X.size_real, X.size_vector, 2*X.node_space);
%Y_bound = norm(Y_PDE);
