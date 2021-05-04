function Z1 = Z1Hopf(A,~,X,G,Es)
% function Z1 = Z1Hopf(A,A_dagger,X,G,Es)
flagBig = 1;
%A_dagger_big = big_Adagger(A_dagger,X); % STIL DOES NOT WORK!
DH =  derivative_Hopf(X,G,Es,flagBig);
DH_op = compose_overoperator(DH,X.size_real, 2*X.node_space, 2*X.node_time);
Dh_minus_Adag = overop2mat(tail(DH_op,X.node_space,X.node_time));
Dh_minus_Adag(1:1+size(Dh_minus_Adag,1):end) = 0;% cancels diagonal
Z1_finite = cnorm(compose_overoperator(A*(Dh_minus_Adag),...
    reshape(X,2*X.node_space,2*X.node_time)));
if any(abs(Z1_finite)>1)
    error('Something went wrong in the reshaping')
end
Z1_tail = Z1TailHopf(X);

Z1 = sum(max(Z1_finite,Z1_tail),2); 

if any(Z1>1)
    disp(Z1)
    warning('Z1 bound too big')
end

end

% compute non-linear derivative and delete center with reshape (reshape
% works)