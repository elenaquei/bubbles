function [conv_u_v] = Convo(u,v,flag)
% function [conv_u_v] = Convo(u,v,flag)
%
% any convolution depending on the inputs
% INPUT
%
% u,v       vectors, either 1D or 2D (can be different)
% flag      for dimension of output, if 1 maximal output, if 0 output equal input
%           (DEFAULT 0)
%
% OUTPUT
% conv_u_v  1D or 2D convolution
if nargin<3 || isempty(flag)
    flag =0;
end
if all(flag == 1) || (length(flag)==4 && all(flag == 'full'))
    option = 'full';
elseif all(flag == 0) || (length(flag)==4 && all(flag == 'same'))
    option = 'same';
else
    error('unrecognised input')
end
% if both 1D, use 1D
if min(size(u))==1 && min(size(v))==1
    if isintval(u) || isintval(v)
        u = intval(u);
        v = intval(v);
        conv_u_v = convINTVAL(u,v,option);
    else
        nodes_u = (length(u)-1)/2;
        nodes_v = (length(v)-1)/2;
        conv_u_v = ifft(fft(u,2*length(u)-1).*fft(v,2*length(v)-1));
        if all(option == 'same') && nodes_u == nodes_v
            conv_u_v = conv_u_v(nodes_u+(1:2*nodes_u+1));
        end
        %conv_u_v = conv(u,v,option);
    end
    return
end
% if both 2D, use 2D
if min(size(u))~=1 && min(size(v))~=1
    if isintval(u) || isintval(v)
        u = intval(u);
        v = intval(v);
        conv_u_v = conv2INTVAL(u,v,option);
    else
        %conv_u_v = conv2(u,v,option);
        fft2u = fft(fft(u,2*size(u,1)-1,1),2*size(u,2)-1,2);
        fft2v = fft(fft(v,2*size(v,1)-1,1),2*size(v,2)-1,2);
        conv_u_v = ifft(ifft(fft2u.*fft2v),[],2);
        index_1 = (size(conv_u_v,1)+1)/2 + (-(size(u,1)-1)/2:(size(u,1)-1)/2);
        index_2 = (size(conv_u_v,2)+1)/2 + (-(size(u,2)-1)/2:(size(u,2)-1)/2);
        if all(option == 'same') && all(size(u)==size(v))
            conv_u_v = conv_u_v(index_1,index_2);
        end
    end
    return
end
% if mixed sizes, make 1D into 2D
if min(size(u))==1
    vector2D = v;
    vector1D = u;
else
    vector2D = u;
    vector1D = v;
end
vector_extended=embedding(vector1D,vector2D);

% call 2D
if ~isintval(vector2D) && ~isintval(vector_extended)
    conv_u_v = Convo(vector2D,vector_extended,option);
else
    
    vector2D = intval(vector2D);
    vector_extended = intval(vector_extended);
    conv_u_v = Convo(vector2D, vector_extended,option);
    
end

