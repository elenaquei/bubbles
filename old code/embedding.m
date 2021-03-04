function vector_extended = embedding(vector1D, new_dim)
% function vector_extended = embedding(small_vector, new_dim)
%
% embedding from 1D vectors into 2D vectors
%
% INPUT
% small_vector      1D vector
% new_dim           dimension of the output, can be EITHER 
%                       - an integer, number of nodes of the output, 
%                   OR  - a 2D vector, the output will have the same
%                   dimension
% OUTPUT
% vector_extended   2D vector of required dimension
if max(size(new_dim))~=1
    vector2D = new_dim;
    new_dim = (size(vector2D,1)-1)/2;
end
vector_extended = zeros(2*new_dim+1,max(size(vector1D)));

if isintval(vector1D)
    vector_extended = intval(vector_extended);
end

vector_extended(new_dim+1,:) = vector1D(:);