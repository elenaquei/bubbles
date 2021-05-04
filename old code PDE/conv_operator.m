function Z = conv_operator(z, flagFull)
% function Z = conv_operator(z, flagFull)
%
% takes as input a 1D or 2D vector and returns an operator compatible with
% it, equivalent to the convolution
%
% INPUT
% z         1D or 2D Fourier coefficients
% flagFull  (DEFAULT 0) complete comvolution or restricted one
% OUTPUT
% Z         corresponding convolution operator

if nargin<2
    flagFull = 0;
end
if length(size(z))>2
    error('Unable to treat this')
end
if min(size(z))~=1
    %z = z(:).';%reshape(z,1,[]);
    flag2D = 1;
else
    flag2D = 0;
end
if ~ flag2D % works - kind of obviously
    nodes = (length(z) -1)/2;
    z = z(nodes+1:-1:1);
    if flagFull
        z = z(:).';
        z = [z,z(2:end)*0];
    end
    Z = toeplitz(z);
else % diesn't work: why so surprised?
    z = z;
    Z = zeros(prod(size(z)));
    if isa(z,'intval')
        Z = intval(Z);
    end
    M = size(z,1);
    first_half = (size(z,1)+1)/2:-1:1;
    second_half = (size(z,1)+1)/2 :size(z,1);
    for i = 1: size(z,2)
        z1 = z(first_half,i);
        z2 = z(second_half,i);
        T = toeplitz([z1(:).',0*first_half(1:end-1)],[z2(:).',0*second_half(1:end-1)]);
        sign_dist_from_diag = (size(z,2)+1)/2 - i;
        distance_from_diag = abs( sign_dist_from_diag );
        %number_of_blocks = size(z,2) - distance_from_diag;
        if isa(z,'intval')
            add_Z = infsup(...
                kron(diag(ones(size(z,2)-distance_from_diag,1),-sign_dist_from_diag),inf(T)),...
                kron(diag(ones(size(z,2)-distance_from_diag,1),-sign_dist_from_diag),sup(T)));
        else
            add_Z = kron(diag(ones(size(z,2)-distance_from_diag,1),-sign_dist_from_diag),T);
        end
        Z = Z + add_Z;
        %         for j = 1: number_of_blocks
        %             k = distance_from_diag;
        %             if sign_dist_from_diag > 0
        %                 index_1_start = M * (k + j - 1)+1;
        %                 index_2_start = M * (j - 1) +1;
        %             else
        %                 index_1_start = M* (j - 1)+1;
        %                 index_2_start = M * (k + j -1)+1;
        %             end
        %             index_1 = index_1_start + (1:M) - 1;
        %             index_2 = index_2_start + (1:M) - 1;
        %             try
        %                 %Z(index_1, index_2) = T;
        %             catch
        %                 disp(k)
        %             end
        %         end
    end
    %Z2 = test(z);
    %if max(max(abs(Z2-Z)))>10^-10
    %    error('90')
    %end
    Z = Z.';
end
return
end

function Z = test(z)

Z = zeros(prod(size(z)));
if isa(z,'intval')
    Z = intval(Z);
end
M = size(z,1);
first_half = (size(z,1)+1)/2:-1:1;
second_half = (size(z,1)+1)/2 :size(z,1);
for i = 1: size(z,2)
    z1 = z(first_half,i);
    z2 = z(second_half,i);
    T = toeplitz([z1(:).',0*first_half(1:end-1)],[z2(:).',0*second_half(1:end-1)]);
    sign_dist_from_diag = (size(z,2)+1)/2 - i;
    distance_from_diag = abs( sign_dist_from_diag );
    number_of_blocks = size(z,2) - distance_from_diag;
    for j = 1: number_of_blocks
        k = distance_from_diag;
        if sign_dist_from_diag > 0
            index_1_start = M * (k + j - 1)+1;
            index_2_start = M * (j - 1) +1;
        else
            index_1_start = M* (j - 1)+1;
            index_2_start = M * (k + j -1)+1;
        end
        index_1 = index_1_start + (1:M) - 1;
        index_2 = index_2_start + (1:M) - 1;
        try
            Z(index_1, index_2) = T;
        catch
            disp(k)
        end
    end
end
end