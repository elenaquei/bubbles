function z=fft_in(z_in,sign)
% function z=fft_in(z_in,sign)
%
% this function computes the fft or the inverse fft of the vector,
% depending on the flag
%
% INPUT
% z_in      vector
% sign      if 1, fft(z_in) is returned, if -1, ifft(z_in) 
%           DEFAULT=1
%
% OUTPUT
% z         vector, either fft(z_in) or ifft(z_in)

if nargin==1
    sign=1;
end

if sign==1
    z=fft(z_in);
else
    z=ifft(z_in);
end
return