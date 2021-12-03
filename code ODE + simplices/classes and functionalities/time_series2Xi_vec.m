function xXi = time_series2Xi_vec(t,y,nodes)
% function xXi = time_series2Xi_vec(t,y,nodes)
% 
% INPUT
% t      vector of times or double with approximate period
% y      matrix of the kind y(linspace(0,t(end),length(y))
% nodes  requested nodes in the output (DEFAULT: size of y /2)
% OUTPUT
% xXi   Xi_vector, such that ifft(xXi) = y

if any(abs(y(1,:)-y(end,:))>0.5)
%    error('The orbit should be closed in order to use Fourier')
end
size_vec = size(y,2);
size_scalar = 1;
if nargin<3 || isempty(nodes)
    n_nodes = floor(size(y,1)/2);
else
    n_nodes = nodes;
end
x=0*y;
for i=1:size_vec
    x(:,i)=1/(size(y,1))*fft(y(:,i));
end
xBar=cell(size_vec+size_scalar,1);
xBar{1}=t(end)/(2*pi);
m=n_nodes;

for i=1:size_vec
    xBar{i+size_scalar}=fftshift([x(1:m+1,i);x(end-m+1:end,i)]);
    xBar{i+size_scalar}(1:m)=conj(flip(xBar{i+size_scalar}(m+2:end)));
    xBar{i+size_scalar}(m+1)=real(xBar{i+size_scalar}(m+1));
end

scal=[t(end)/(2*pi)];
vec=[xBar{1+size_scalar}];
for i=2:size_vec
    vec=[vec,xBar{i+size_scalar}];
end
xXi=Xi_vector(scal,vec,1,size_vec,n_nodes);