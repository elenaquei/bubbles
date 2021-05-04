function new_G = reshape_linear_coef(G,size_real, old_node_space,new_node_space, new_node_time)

size_G_new = size_real + 2*new_node_space+1 + (2*new_node_space+1)*(2*new_node_time+1);
new_G = zeros(size(G,1),size_G_new);

for i = 1:size(G,1)
    temp_X_PDE = PDE_vector(G(i,:),size_real, 1, old_node_space);
    temp_X_reshape = reshape(temp_X_PDE, new_node_space, new_node_time);
    new_G(i,:)= PDE_vec2vec(temp_X_reshape);

end