function Mat_a = Convo_mat(a)
% INPUT 1D sequence 'a' in Fourier
% OUTPUT Mat_a matrix such that for compatible sequence 'b', the
% representation Convo(a,b) = Mat_a*b is valid. In other words, implements
% truncated convolution as a matrix operation.
N = (length(a)-1)/2;
a = flip(a);
M1 = triu(toeplitz(a,a.'));
M2 = toeplitz([flip(a(1:N+1));zeros(N,1)],[a(N+1:end);zeros(N,1)].');
Mat_a = [M2(1:N,:);M1(1:N+1,:)];
end