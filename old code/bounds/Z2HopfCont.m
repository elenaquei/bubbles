function Z2cont = Z2HopfCont(A0, A1, X0, X1,Rstar)

Xs = interval_PDE(X0,X1);
Z2cont = max(Z2Hopf(A0,Xs,Rstar),Z2Hopf(A1,Xs,Rstar));

end
% bound_tail = zeros(5);
% 
% node_space = X0.node_space;
% 
% x0 = norm(X0);
% x1 = norm(X1);
% x = max(x0,x1);
% 
% y = x(4);
% z = x(5);
% nu = x(1);
% lambda = x(2);
% a = x(3);
% 
% % tail bound
% 
% % Bj bounds
% l1_sign = X0.parameters(1);
% l2_sign = X0.parameters(2);
% 
% % first line
% B1 = tailBoundHopf(1, 1, l2_sign, node_space);
% %B2 = tailBoundHopf(2, 1, l2_sign, node_space);
% B4 = tailBoundHopf(4, 1, l2_sign, node_space);
% 
% % fourth line of the global tail bounds
% bound_tail(4,:)= [ 0,B4,0, 4*B1+B4, 0];
% 
% % second line
% B1 = tailBoundHopf(1, l1_sign, l2_sign, node_space);
% B2 = tailBoundHopf(2, l1_sign, l2_sign, node_space);
% B4 = tailBoundHopf(4, l1_sign, l2_sign, node_space);
% 
% % fifth line of the global tail bounds
% bound_tail(5,:)= [B4*(z+lambda)+B2+2*B1*(z^2+2*z+2*y+2*a*z),B4*nu,2*B1*(z^2+2*nu*z),4*B1*(nu+z),B4*(nu+lambda)+B2+4*B1*(nu+y+nu*z+a*z+nu*a)];
% 
% 
% % center bound
% D2H_square_shape = cell(5);
% A_matching_size=cell(5,5);
% 
% A_abs = compose_overoperator(max(abs(A0),abs(A1)), X0);
% 
% for j = 1:3
%     for i = 1:3
%         D2H_square_shape{i,j} = 0;
%         A_matching_size{i,j} = A_abs.matrix_cell{1,1}.matrix(i,j);
%     end
% end
% 
% for j = 2:3
%     for i = 1:3
%         D2H_square_shape{i,j+2} = 0 * A_abs.matrix_cell{1,j}.matrix(1,:);
%         A_matching_size{i,j+2} = 0 * A_abs.matrix_cell{1,j}.matrix(1,:);
%         A_matching_size{j+2,i} = A_abs.matrix_cell{j,1}.matrix(:,i);
%     end
% end
% for i = 2:3
%     for j = 2:3
%          A_matching_size{i+2,j+2} = A_abs.matrix_cell{i,j}.matrix;
%     end
% end
% 
% K2_y = derivative_operator(2,X0.vector1D);
% K2_z = derivative_operator(2,squeeze(X0.vector2D));
% %K1_2D = derivative_operator(1,z2D);
% 
% % K2_yvec = -X0.node_space:X0.node_space;
% % K2_y = diag(K2_yvec);
% % K2_zvec = repmat(K2_yvec,1+2*X0.node_time,1);
% % K2_z = diag(reshape(K2_zvec,1,[]));
% 
% 
% D2H_square_shape{4,1} = K2_y^4*(y+nu)+4*K2_y*y;
% D2H_square_shape{4,2} =  K2_y^4 *(y+nu);
% % {4,3}
% D2H_square_shape{4,4} = K2_y^4*(lambda+nu)+K2_y^2+4*K2_y*(y+nu);
% % {4,5}
% D2H_square_shape{5,1} = K2_z^4*(z+lambda)+4*K2_z*(z^2+z+y+2*a*z);
% D2H_square_shape{5,2} = K2_z^4*(z+nu);
% D2H_square_shape{5,3} = 4 * K2_z *(z^2+2*nu*z);
% D2H_square_shape{5,4} = 4*K2_z*(lambda+z);
% 
% column_selection = X0.node_time+1+(0:2*node_space)*(2*X0.node_time+1);
% D2H_square_shape{5,4} = D2H_square_shape{5,4}(:,column_selection);
% 
% D2H_square_shape{5,5} = K2_z^4*(lambda+nu)+K2_z^2+4*K2_z*(nu+y+2*nu*z+2*a*z+2*nu*a);
% 
% D2H_square_shape{4,5} = zeros(size(D2H_square_shape{4,2},1),size(D2H_square_shape{5,5},2));
% 
% 
% ones_y = ones(size(D2H_square_shape{4,1},1),1);
% ones_z = ones(size(D2H_square_shape{5,1},1),1);
% 
% D2H_square_shape{4,1} = D2H_square_shape{4,1}*ones_y;
% D2H_square_shape{4,2} = D2H_square_shape{4,2}*ones_y;
% D2H_square_shape{4,3} = 0*D2H_square_shape{4,2};
% 
% D2H_square_shape{5,1} = D2H_square_shape{5,1}*ones_z;
% D2H_square_shape{5,2} = D2H_square_shape{5,2}*ones_z;
% D2H_square_shape{5,3} = D2H_square_shape{5,3}*ones_z;
% 
% 
% % product 
% prod_A_D2H = cell(5);
% 
% for i = 1:5
%     for j = 1:5
%         prod_A_D2H{i,j} = 0*A_matching_size{i,1}*D2H_square_shape{1,j};
%         for k = 1:5
%             prod_A_D2H{i,j} = prod_A_D2H{i,j} + A_matching_size{i,k}*D2H_square_shape{k,j};
%         end
%     end
% end
% 
% % reconstruct the product
% 
% prod_mat = [prod_A_D2H{1,1} prod_A_D2H{1,2} prod_A_D2H{1,3} prod_A_D2H{1,4} prod_A_D2H{1,5}
%     prod_A_D2H{2,1} prod_A_D2H{2,2} prod_A_D2H{2,3} prod_A_D2H{2,4} prod_A_D2H{2,5}
%     prod_A_D2H{3,1} prod_A_D2H{3,2} prod_A_D2H{3,3} prod_A_D2H{3,4} prod_A_D2H{3,5}
%     prod_A_D2H{4,1} prod_A_D2H{4,2} prod_A_D2H{4,3} prod_A_D2H{4,4} prod_A_D2H{4,5}
%     prod_A_D2H{5,1} prod_A_D2H{5,2} prod_A_D2H{5,3} prod_A_D2H{5,4} prod_A_D2H{5,5}];
% 
% prod_op = compose_overoperator(prod_mat, X0);
% bound_center = cnorm(prod_op);
% 
% % max of the two bounds
% Z2cont = max(bound_tail, bound_center);