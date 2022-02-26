function Xi_vec = vec2Xi_vec(vec,size_scal,size_vec,nodes)
% function Xi_vec = vec2Xi_vec(vec,size_scal,size_vec,nodes)
%
% this function takes a vector and return a Xi_vector
% Inputs:
% vec          vector
% size_scal    number of scalar components
% size_vec     number of vector components
% nodes        number of nodes
% OR
% size_scalar  Xi_vector of appropriate shape

if nargin==2
    Xi_vec=size_scal;
    size_scal=Xi_vec.size_scalar;
    size_vec=Xi_vec.size_vector;
    nodes=Xi_vec.nodes;
end

k=size_scal+size_vec*(2*nodes+1);

if min(size(vec))~=1
    error('First input must be a vector')
end
if max(size(vec))~=k
    error('Dimensions do not agree')
end

scalars=vec(1:size_scal);
if size(scalars,1)>size(scalars,2)
    scalars=scalars.';
end
vectors=reshape(vec(size_scal+1:end),2*nodes+1,size_vec).';

Xi_vec=Xi_vector(scalars,vectors,size_scal,size_vec,nodes);