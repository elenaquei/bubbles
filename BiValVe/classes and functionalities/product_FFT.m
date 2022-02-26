function y_conv=product_FFT(x_conv,powers_vectors)
% function y_conv=product_FFT(x_conv,powers_vectors)
%
% returns x_conv.^powers_vector
% the shortest size of x_conv must be equal to the number of elements in
% powers_vectors

size_x=size(x_conv);
turn=0;

if length(size_x)>2
    error('here troubles');
end

if size_x(1)>size_x(2)
    x_conv=x_conv.';
    size_x=size(x_conv);
    turn=1;
end

if sum(powers_vectors)<=0
    y_conv=1/size_x(2)+x_conv(1,:)*0;
    return
elseif sum(powers_vectors)==1
    y_conv=x_conv(powers_vectors==1,:);
    return
end
y_conv=1;


for i=1:size_x(1)
    y_conv=y_conv.*(x_conv(i,:).^powers_vectors(i));
end

if turn
    y_conv=y_conv.';
end

l=length(y_conv);
y_conv=l^(sum(powers_vectors)-1)*y_conv;

return
end