function par_vec = compute_linear_func(x_vec, G, const_G)
%function par_vec = compute_linear_func(x_vec, G, const_G)
%
% INPUT
% x_vec     vector of length m
% G         matrix of size n times m
% const_G   vector of length n
if isa(x_vec,'PDE_vector')
    x_vec = PDE_vec2vec(x_vec);
end

if length(x_vec)~= length(G(1,:))
    if length(x_vec)== length(G(:,1))
        G = G.';
    else
        error('Size input incompatible')
    end
end
   
par_vec = zeros(size(G,1),1);
for i = 1:size(G,1)
    par_vec(i) = x_vec * (G(i,:).')+const_G(i);
end