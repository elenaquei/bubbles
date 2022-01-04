function norm_mat=norm_Ximat(B,xBar_size_scal, xBar_size_vector, xBar_nodes_col, xBar_nodes_row)
% function norm_mat=norm_Ximat(B,xBar)

if nargin == 2
    xBar = xBar_size_scal;
    norm_mat= norm_mat_loc(B,xBar.size_scalar,xBar.size_vector,xBar.nodes);
else
    if nargin==4
        xBar_nodes_row=xBar_nodes_col;
    elseif nargin<4
        error('not enough input arguments');
    end
    norm_mat= norm_mat_loc(B,xBar_size_scal, xBar_size_vector,  xBar_nodes_col, xBar_nodes_row);
end

return

function [norm_vec]= norm_mat_loc(Z_mat,size_scalar,size_vector,nodes_col,nodes_row)
%function [norm_vec]=norm_mat_loc(Z_mat,size_scal,size_vec,nodes_col,nodes_row)
%
% given size_scal,size_vec, nodes_col, nodes_row
% generalization of the norm of the Xi matrices in case the vertical
% and horizontal number of components is not the same.
global nu
global use_intlab
global norm_weight

if isempty(norm_weight)
    norm_weight_loc = ones(size_scalar+size_vector,1);
elseif length(norm_weight) ~= size_scalar+size_vector
    if length(norm_weight) == 1
        norm_weight_loc = norm_weight * ones(size_scalar+size_vector,1);
    else
        error('The weight of the norm is incompatible with the size of the problem')
    end
else
    norm_weight_loc = norm_weight;
end

if nargin==4
    nodes_row=nodes_col;
elseif nargin<4
    error('not enough input arguments');
end

if use_intlab && ~isintval(Z_mat)
    Z_mat = intval(Z_mat);
end
scalar_weights = norm_weight_loc(1:size_scalar);
vector_weights = norm_weight_loc(1+size_scalar:end);
if ~use_intlab
    norm_scal=sum(diag(scalar_weights)*abs(Z_mat(1:size_scalar,1:size_scalar))*diag(scalar_weights.^-1),2);
else
    norm_scal=sum(sup(diag(scalar_weights)*abs(Z_mat(1:size_scalar,1:size_scalar))*diag(scalar_weights.^-1)),2);
end
for i=1:size_vector
    if ~use_intlab
        norm_scal=norm_scal+ sum(diag(scalar_weights)*abs(Z_mat(1:size_scalar,size_scalar+(i-1)*(2*nodes_col+1)+(1:2*nodes_col+1)))/vector_weights(i)...
            .*repmat((nu.^ abs(-nodes_col:nodes_col)),[size_scalar,1]),2);
    else
        
        norm_scal=norm_scal+ sum(sup(diag(scalar_weights)*abs(Z_mat(1:size_scalar,size_scalar+(i-1)*(2*nodes_col+1)+(1:2*nodes_col+1)))/vector_weights(i)...
            .*repmat(nu.^abs(-nodes_col:nodes_col),[size_scalar,1])),2);
    end
end
norm_vec=zeros(size_vector,1);
for j=1:size_vector
    for i=1:size_scalar
        Zji=Z_mat(size_scalar+(j-1)*(2*nodes_row+1)+(1:2*nodes_row+1),i);
        norm_vec(j)=norm_vec(j)+NuNorm(Zji,nu)*norm_weight_loc(j)/norm_weight_loc(i);
    end
    %disp(Z)
    for i=1:size_vector
        Zji=Z_mat(size_scalar+(j-1)*(2*nodes_row+1)+(1:2*nodes_row+1),size_scalar+(i-1)...
            *(2*nodes_col+1)+(1:2*nodes_col+1));
        norm_vec(j)=norm_vec(j)+InfNuNorm_mat(Zji)*norm_weight_loc(j)/norm_weight_loc(i);
    end
end
norm_vec=[norm_scal;norm_vec];
return