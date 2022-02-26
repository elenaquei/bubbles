function vector = Xi_vec2vec(Xi_v)
% function vector = Xi_vec2vec(Xi_vector)
% 
% This function takes as inout a Xi_vector and returns a vector, such that
% the first M components are the scalars, then the next 2m+1 are the
% Fourrier coefficients of the first vector and so on until the end

%global use_intlab

k=Xi_v.size_scalar+Xi_v.size_vector*(2*Xi_v.nodes+1);
if ~isintval(Xi_v.scalar)
    vector=zeros( k,1);
else
    vector=intval(zeros(k,1));
end

vector(1:Xi_v.size_scalar)=Xi_v.scalar;

for i=1:Xi_v.size_vector
    indeces=Xi_v.size_scalar+(i-1)*(2*Xi_v.nodes+1)+(1:(2*Xi_v.nodes+1));
    vector(indeces)=Xi_v.vector(i,:);
end

    
    