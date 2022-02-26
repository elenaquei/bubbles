function A_ij= component_matrix_norm(A,xBar)
% INPUT
% A       square matrix
% xBar    Xi vector
%
% size(A)=xBar.size_scalar + xBar.size_vector*(xBar.nodes*2+1)
%
% OUTPUT
% A_ij   square matrix, elementwise operator nu-norm of A
% size(A_ij) = xBar.size_scalar + xBar.size_vector
%
% nu must be a global variable

global nu

if size(A,1)~=xBar.size_scalar + xBar.size_vector*(xBar.nodes*2+1)
    error('Sizes do not coincide');
elseif size(A,1)~=size(A,2)
    error('The input matrix must be square');
end

A_ij= zeros(xBar.size_scalar+xBar.size_vector);
size_vec=(xBar.nodes*2+1);

%A_ij(1:xBar.size_scalar,1:xBar.size_scalar) = abs(A(1:xBar.size_scalar,1:xBar.size_scalar));

for i=1:xBar.size_vector
    % vector - scalar
    for j=1:xBar.size_scalar
        a_ij=A(xBar.size_scalar+1+ (i-1)*size_vec:xBar.size_scalar+i*size_vec,...
            j);
        A_ij(xBar.size_scalar+i,j)=NuNorm(a_ij,nu);
    end
    % scalar - vector
    for j=1:xBar.size_scalar
        a_ij=A(j,xBar.size_scalar+1+ (i-1)*size_vec:xBar.size_scalar+i*size_vec);
        A_ij(j,xBar.size_scalar+i)=NuNorm(a_ij,nu);
    end
    % vector - vector
    for j=1:xBar.size_vector
        a_ij=A(xBar.size_scalar+1+ (i-1)*size_vec:xBar.size_scalar+i*size_vec,...
            xBar.size_scalar+1+ (j-1)*size_vec:xBar.size_scalar+j*size_vec);
        A_ij(xBar.size_scalar+i,xBar.size_scalar+j)=norm_cor1(a_ij,xBar.nodes,xBar.scalar(1));
    end
end

return