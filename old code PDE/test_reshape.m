%test_reshape

node_space = 2;
node_time = 2;

new_node_space = 3;
new_node_time = 5;

diff_space = new_node_space - node_space;
diff_time = new_node_time - node_time;

a = rand(2*node_time+1, 2*node_space+1);
a_big = [zeros(diff_time, diff_space), zeros(diff_time, 2*node_space+1), zeros(diff_time, diff_space);
    zeros(2*node_time+1,diff_space), a, zeros(2*node_time+1,diff_space);
    zeros(diff_time, diff_space), zeros(diff_time, 2*node_space+1), zeros(diff_time, diff_space)];

op_small = conv_operator(a,1);
op_big = conv_operator(a_big,1);

Op_small = operator(op_small,'L1','L1');
Op_big = operator(op_big,'L1','L1');

Op_small = expand(Op_small,new_node_space, node_space, new_node_time);

Op_small_mat = Op_small.matrix;

% projection of big in small
op_big(1:2*new_node_time+1,:) = 0;
op_big(:,1:2*new_node_time+1) = 0;
op_big(end+1-(1:2*new_node_time+1),:) = 0;
op_big(:,end+1-(1:2*new_node_time+1)) = 0;

figure; spy(Op_small_mat); title('Reshaped')
figure; spy(op_big); title('Original')

if any(any(abs(op_big -Op_small_mat)>10^-10))
    figure; spy(op_big - Op_small_mat)
%    error('Reshape does not work')
end

% EXPAND FOR OPERATORS WORKS




