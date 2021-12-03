function n = norm_Xi(X , size_scalar, size_vector, nodes_col, nodes_row, nu_local)
% function n = norm_Xi(X , size_scalar, size_vector, nodes_col, nodes_row, nu_local)
%
% compute the norm of X in the \R^size_scalar \times (l^1_\nu
% norm)^size_vector norm, therefore return a vector of length
% size_scalar + size_vector
% Input can be a Xi_vector, a vector or a matrix. If one of the latter two,
% the inputs size_scalar, size_vector and nodes_col are required. Just if
% the matrix is not square, will the input nodes_row be requested. The last
% input nu_local is requested if the global value nu does not exist. The 
% current value nu overwrited the value of global nu.
% If X Xi_vector, use the norm(Xi_vector)!
% 
% INPUT
% X             being it a matrix, a vector or a Xi_vector
% size_scalar   integer
% size_vector   integer
% nodes_col     integer, such that size_scalar + (2*nodes_col)*size_vector
%                 is the length of the vector and the number of columns of
%                 the matrix
% nodes_row    integer, DEFAULT nodes_col, such that 
%                 size_scalar + (2*nodes_col)*size_vector
%                 is the number of rows of the matrix
% nu_local     value of nu, DEFAULT global nu
% 
% OUTPUT 
% n            norm(X)_(l^1_nu), vector of length size_scalar + size_vector
global nu

if nargin < 6
    nu_local = nu;
end

if nargin ==1 
    disp('Use norm(X) instead')
    n = norm(X);
    return
elseif isempty(size_scalar) && isempty(size_vector) && isempty(nodes_col) && isempty(nodes_row)
    n = cnorm_Xi_vector(X,nu_local);
    return
elseif min(size(X))==1
    if length(X) ~= size_scalar + (2*nodes_col)*size_vector
        error('Length of X and other inputs do not match')
    end
    n = cnorm_Xi_vector(vec2Xi_vec(X,size_scalar,size_vector, nodes_col),nu_local);
    return
end

if nargin==4
    nodes_row = nodes_col;
end

if (size(X,1)~=size_scalar + (2*nodes_row)*size_vector) || (size(X,2)~=size_scalar + (2*nodes_col)*size_vector)
    error('Dimension of matrix and other inputs do not match')
end

n = test_norm_mat(X,size_scalar,size_vector, nodes_col, nodes_row, nu_local);
n2 = norm_mat(X, size_scalar, size_vector, nodes_col, nodes_row, nu_local);
end 


function nor = norm_mat(X, size_scal, size_vec, nodes_col, nodes_row, nu)
% function nor = norm_mat(X, size_scal, size_vec, nodes_col, nodes_row, nu)
% 
% helper function

nor = zeros(size_scal + size_vec);

for i =1:size_scalar
    Inf_norm =0;
    index = 0;
    for j = 1:size_vector
        index = index(end) + (-nodes_col:nodes_col);
        Inf_norm = inf_nu(X(j,index),nodes_row,nu);
    end
    nor(i) = sum(abs(X(1:size_scalar,i))) + Inf_norm;
end

index_i = 0;

for i = size_scalar+(1:size_vector)
    Nu_norm = 0;
    Mat_norm = 0;
    index_i= index_i(end) + (-nodes_row:nodes_row);
    for j = 1:size_scalar
        Nu_norm = Nu_norm + nu_norm(X(index_i,j),nodes,nu);
    end
    index=0;
    for j = 1:size_vector
        index = index(end) + (-nodes_col:nodes_col);
        Mat_norm = Mat_norm + mat_nu(X(index_i, index),nodes_row, nodes_col, nu);
    end
    nor(i) = Nu_norm + Mat_norm;
end
end

function nor = inf_nu(v,nodes,nu)
nor = max ( abs(v) ./( nu.^(-nodes:nodes)));
end

function nor = nu_norm(v, nodes, nu)
nor = sum(abs(v) .* ( nu.^(-nodes:nodes)));
end

function nor = mat_nu(B,nodes_row, nodes_col, nu)
nor = max((nu.^abs(-nodes_row:nodes_row)*abs(B)).*nu.^(-abs(-nodes_col:nodes_col)));
end
