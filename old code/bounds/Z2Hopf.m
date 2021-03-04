function Z2 = Z2Hopf(A,X,Rstar)
% function Z2 = Z2Hopf(A,X,Rstar)
% 

finite_part = Z2_finite_Hopf(A,X,Rstar);

infinite_part = Z2_infinite_Hopf(X,Rstar);

Z2 = max(finite_part,infinite_part);
end

function finite_part = Z2_finite_Hopf(A,X,Rstar)

X_norm = norm(X);
nu = X_norm(1)+Rstar;
lambda = X_norm(2)+Rstar;
a = X_norm(3)+Rstar;
y = X_norm(4)+Rstar;
z = X_norm(5)+Rstar;

k = -X.node_space:X.node_space;
K2_y = diag(k);
k_rep = repmat(k',1, 2*X.node_time+1);
K2_z = diag(k_rep(:));

A_split = cell(5,5);

nodes_4th = 3+(1:(2*X.node_space+1));
nodes_5th = nodes_4th(end)+(1:(2*X.node_space+1)*(2*X.node_time+1));

if nodes_5th~=size(K2_z,1)
    error('Problem with dimensions')
end

A_split{1,4} = A(1,nodes_4th);
A_split{2,4} = A(2,nodes_4th);
A_split{3,4} = A(3,nodes_4th);
A_split{4,4} = A(nodes_4th,nodes_4th);
A_split{5,4} = A(nodes_5th,nodes_4th);

A_split{1,5} = A(1,nodes_5th);
A_split{2,5} = A(2,nodes_5th);
A_split{3,5} = A(3,nodes_5th);
A_split{4,5} = A(nodes_4th,nodes_5th);
A_split{5,5} = A(nodes_5th,nodes_5th);

norm_full = zeros(5);

norm_ij = @(i,j,x) norm_IJ(X.node_space,X.node_time,i,j,x);

for i = 1:5
    AijKp = zeros(5,4);
    for j = 4:5
        if j ==4
            K_temp = K2_y;
        else
            K_temp = K2_z;
        end
        for p = [1,2,4]
            AijKp(j,p) = norm_ij(i,j,A_split{i,j}*K_temp.^p);
        end
    end
    norm_full(i,1) = AijKp(5,1)*(z^2+2*z+2*y+2*a*z)+...
        AijKp(5,2)+AijKp(5,4)*(z+nu);
    
    norm_full(i,2) = AijKp(4,4)+AijKp(5,4)*(z+lambda);
    
    norm_full(i,3) = AijKp(5,1)*(z^2+2*lambda*z);
    
    norm_full(i,4) = 2*AijKp(4,1)+AijKp(4,4)+...
        AijKp(5,1)*2*(lambda+z);
    
    norm_full(i,5) = AijKp(5,1)*(y+lambda+a*z+lambda*z+lambda*a)+AijKp(5,2)+...
        AijKp(5,4)*(nu+lambda);
end
finite_part = sum(norm_full,2);
end

function N = norm_IJ(node_space, node_time,i,j,x)
global nu

K2 = -node_space:node_space;
K1 = -node_time:node_time;
K1 = K1.';
%K_1and2 = repmat(K1,1,1+2*node_space).*repmat(K2,1+2*node_time,1);
%K12_vec = reshape(K_1and2,1,[]);

nu_K_1and2 = nu(1).^abs(repmat(K1,1,1+2*node_space)).*...
    nu(end).^abs(repmat(K2,1+2*node_time,1));
nu_K_1and2_vec = reshape(nu_K_1and2,1,[]);

nu_K2 = nu(end).^(abs(K2));

%warning('Norms in Z2 to check')

norm_nu = @(x) ( nu_K2*abs(x));
norm_ell = @(x) (nu_K_1and2_vec*abs(x));
norm_Bnu = @(x) max(nu_K2*abs(x).*(nu_K2.^-1)); 
norm_Bell = @(x) max((nu_K_1and2_vec*abs(x)).*(nu_K_1and2_vec.^-1)); 
norm_nuInf = @(x) max(abs(x).*(nu_K2.^-1));
norm_ellInf = @(x) max(abs(x).*(nu_K_1and2_vec.^-1));
norm_nu2ell = @(x) max((nu_K_1and2_vec*abs(x)).*(nu_K2.^-1)); 
norm_ell2nu = @(x) max((nu_K2*abs(x)).*(nu_K_1and2_vec.^-1)); 

switch i 
    case {1,2,3}
        switch j
            case {1,2,3}
                N = abs(x);
            case 4
                N = norm_nuInf(x);
            case 5 
                N = norm_ellInf(x);
        end
    case 4
        switch j
            case {1,2,3}
                N = norm_nu(x);
            case 4
                N = norm_Bnu(x);
            case 5 
                N = norm_ell2nu(x);
        end
    case 5
        switch j
            case {1,2,3}
                N = norm_ell(x);
            case 4
                N = norm_nu2ell(x);
            case 5 
                N = norm_Bell(x);
        end
end


end


function infinite_part = Z2_infinite_Hopf(X,Rstar)
lambda1 = X.parameters(1);
lambda2 = X.parameters(2);

B1 = tailBoundHopf(1,lambda1, lambda2, X.node_space);
B2 = tailBoundHopf(2,lambda1, lambda2, X.node_space);
B4 = tailBoundHopf(4,lambda1, lambda2, X.node_space);

X_norm = norm(X);
nu = X_norm(1)+Rstar;
lambda = X_norm(2)+Rstar;
a = X_norm(3)+Rstar;
y = X_norm(4)+Rstar;
z = X_norm(5)+Rstar;

term51 = B4*(z+nu)+B2+B1*(z^2+2*z+2*y+2*a*z);
term55 = B4*(nu+lambda)+B2+2*B1*(lambda+y+lambda*z+a*z+lambda*a);

mid_way = zeros(5);
mid_way(4:5,:) = ...
    [ 0,B4,0, 2*B1+B4, 0;
    term51,B4*(z+lambda),B1*(z^2+2*lambda*z),2*B1*(lambda+z),term55];
infinite_part = sum(mid_way,2);
end


