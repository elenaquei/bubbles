function A_norm = norm_matrix(A, size_real, size_vector, node_space, node_time)
% function N = norm_matrix(A, size_real, size_vector, node_space, node_time)
%
% INPUT 
% A             matrix of adequate dimension
% size_real     number of parameters
% size_vector   number of vector unknowns
% node_space    nodes in spacial direction
% node_time     nodes in time direction
%
% OR
%
% INPUT
% A
% X             PDE_vector

%global nu

if nargin ==2
    size_vector = size_real.size_vector;
    node_space = size_real.node_space;
    node_time = size_real.node_time;
    size_real = size_real.size_real;
end

if size(A,1)~=size(A,2)
    error('Matrix must be square')
end
if size(A,1)~= (size_real+size_vector*(2*node_space+1+(2*node_space+1)*(2*node_time+1)))
    error('Inputs do not match, dimensions incoherent')
end

A = compose_overoperator(A,size_real, node_space, node_time);
A_norm = norm(A);


% K2 = -node_space:node_space;
% K1 = -node_time:node_time;
% K1 = K1.';
% K_1and2 = repmat(K1,1,1+2*node_space).*repmat(K2,1+2*node_time,1);
% K12_vec = reshape(K_1and2,1,[]);
% 
% norm_nu = @(x) sum( nu.^(abs(K2).').*abs(x));
% norm_ellnu = @(x) sum( nu.^(abs(K12_vec).').*abs(x));
% norm_Bnu = @(x) max((nu.^abs(K2)*abs(x)).*nu.^(-abs(K2))); 
% norm_Bellnu = @(x) max((nu.^abs(K12_vec)*abs(x)).*nu.^(-abs(K12_vec))); 
% norm_nuInf = @(x) max(abs(x).*(nu.^-abs(K2)));
% norm_ellInf = @(x) max(abs(x).*(nu.^-abs(K12_vec)));
% norm_nu2ell = @(x) max((nu.^abs(K12_vec)*abs(x)).*nu.^(-abs(K2))); 
% norm_ell2nu = @(x) max((nu.^abs(K2)*abs(x)).*nu.^(-abs(K12_vec))); 
% % at least dimensionally meaningful
% 
% 
% node_real = 1:size_real;
% node_vec_space = node_real(end)+(1:(2*node_space+1)*size_vector);
% length_one_vec_space = (2*node_space+1);
% node_vec_time = node_vec_space(end) + (1:(2*node_space+1)*(2*node_time+1)*size_vector);
% length_one_vec_time = (2*node_space+1)*(2*node_time+1);
% 
% nodes_y_i_begin = size_real +((1:size_vector)-1)*length_one_vec_space+1;
% nodes_y_i_end = size_real +(1:size_vector)*length_one_vec_space;
% nodes_z_i_begin = size_real + length_one_vec_space*size_vector + ((1:size_vector)-1)*length_one_vec_time+1;
% nodes_z_i_end = size_real + length_one_vec_space*size_vector + (1:size_vector)*length_one_vec_time;
% 
% A_norm = zeros(size_real+size_vector*2);
% 
% A_norm(1:size_real,1:size_real) = abs(A(node_real,node_real));
% 
% for i=1:size_vector
%     for j = 1:size_real
%         A_norm(j, size_real+i) = norm_nuInf (A(j,nodes_y_i_begin(i):nodes_y_i_end(i)));
%         A_norm(j, size_real+size_vector+i) = norm_ellInf (A(j,nodes_z_i_begin(i):nodes_z_i_end(i)));
%         A_norm(size_real+i, j) = norm_nu (A(nodes_y_i_begin(i):nodes_y_i_end(i),j));
%         A_norm(size_real+size_vector+i, j) = norm_ellnu (A(nodes_z_i_begin(i):nodes_z_i_end(i),j));
%     end
%     for j = 1:size_vector
%         A_norm(size_real+j, size_real+i) = ...
%             norm_Bnu(A(nodes_y_i_begin(j):nodes_y_i_end(j),nodes_y_i_begin(i):nodes_y_i_end(i)));
%         A_norm(size_real+j, size_real+size_vector+i) = ...
%             norm_ell2nu(A(nodes_y_i_begin(j):nodes_y_i_end(j),nodes_z_i_begin(i):nodes_z_i_end(i)));
%         A_norm( size_real+size_vector+i,size_real+j)=...
%              norm_nu2ell(A(nodes_z_i_begin(i):nodes_z_i_end(i),nodes_y_i_begin(j):nodes_y_i_end(j)));
%         A_norm( size_real+size_vector+i,size_real+size_vector+j)=...
%              norm_Bellnu(A(nodes_z_i_begin(i):nodes_z_i_end(i),nodes_z_i_begin(j):nodes_z_i_end(j)));
%     end
% end



end