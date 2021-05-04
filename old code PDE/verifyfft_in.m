function Z = verifyfft_in(z,sign,flag_long)
% function Z = verifyfft_in(z,sign,flag_long)
%
% This functionis just a help function to be able to call verifyfft with
% any vector
% z is the vector on wich the FFT needs to be applied, written in the form
%    z=[a_{-m}.....a_0....a_{-m}]
% sign is 1 for FFt and -1 for IFFT DEFAULT=1
% Z is the FFT of z using the verifyfft
% if flag_long = 1, the solution will not be cut (! padded with zeros !).
% DEFAULT=0

%global use_intlab

if nargin==1
    sign=1;
    flag_long=0;
elseif nargin==2
    flag_long=0;
else
    %  error('Too many inputs');
end
if flag_long~=1 && flag_long~=0
    flag_long=0;
end
len=length(z);
cols=1;

%turn into column vector
if size(z,2)>1
    if size(z,1)==1
        z=z.';
    else
        cols=size(z,2);
    end
end

%determine the next power of 2
next_pow=ceil(log2(len));

if (2^next_pow-len) ==0
    Z = do_the_thingy(z, sign);
    Z=check_transpose(Z,z);

    return;
end

diff_len= 2^next_pow-len;


if mod(diff_len,2)==1
    m=(diff_len-1)/2;
    y=[zeros(m+1,cols);z;zeros(m,cols)];
    
    Z = do_the_thingy(y, sign);
    if ~flag_long
        Z=Z(m+2:end-m,:);
    end
    Z=check_transpose(Z,z);
    return
end

m=diff_len/2;

y=[zeros(m,cols);z;zeros(m,cols)];
Z = do_the_thingy(y, sign);
if ~flag_long
    Z=Z(m+1:end-m,:);
end
Z=check_transpose(Z,z);
return
end


function Z = do_the_thingy(y, sign)

if sign == -1
    y = ifftshift(y)*length(y);
end

if isa(y,'intval') 
    Z=verifyfft(y,sign);
else
    Z=fft_in(y,sign);
end

if sign == 1
    Z = fftshift(Z)/length(y);
end

end