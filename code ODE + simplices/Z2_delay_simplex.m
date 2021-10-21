function [Z2vector,Z2_s]=Z2_delay_simplex(x0,x1,x2,alpha, RAD_MAX)
global nu

modes = x0.n_nodes;

[A_small0, M0, P0, Q0, R0, phi0, D3F20] = A_delay_symplex(alpha, x0);
[A_small1, M1, P1, Q1, R1, phi1, D3F21] = A_delay_symplex(alpha, x1);
[A_small2, M2, P2, Q2, R2, phi2, D3F22] = A_delay_symplex(alpha, x2);

A_small_int = three_intval(A_small0, A_small1, A_small2);

M_int = three_intval(M0, M1, M2);
P_int = three_intval(P0, P1, P2);
Q_int = three_intval(Q0, Q1, Q2);
R_int = three_intval(R0, R1, R2);
phi_int = three_intval(phi0, phi1, phi2);
D3F2_int = three_intval(D3F20, D3F21, D3F22);

block_norm_A = block_norm(M_int, P_int, Q_int, R_int, D3F2_int/(modes+1), phi_int, x0);

upper_bound_psi_DF = [1,2,3]; %TODO

Z22 = block_norm_A * upper_bound_psi_DF;

if nu^(2*modes+2) < exp(1)
    Lemma_bound = nu^(2*modes+2)/ (exp(1) * log(nu^(2*modes+2)));
else
    Lemma_bound = 1;
end

block_norm_AM = block_norm((modes+1)*M_int, P_int, (modes+1)*Q_int, R_int, D3F2_int, phi_int, x0);
upper_bound_DDF_without_dpsi = [1,2,3];
Z21 = block_norm_AM * upper_bound_DDF_without_dpsi;

Z2 = Z22 + Z21;
end

function upper_bound_psi_DF = dpsiDF(alpha, x, Rmax)

K = -x.nodes:x.nodes;
x_norm = norm(x) + Rmax;
lambda_norm = x_norm(1:x.size_scalar);
v_norm = x_norm(x.size_scalar:end);
Kx = x;
for j = 1:x.size_vector
    Kx.vector(j,:) = K.*x.vector;
end
Kx_norm = norm(Kx)+Rmax;
Kv_norm = Kx_norm(x.size_scalar+1:end);
KKx = x;
for j = 1:x.size_vector
    KKx.vector(j,:) = K.*Kx.vector;
end
KKx_norm = norm(KKx)+Rmax;
KKv_norm = KKx_norm(x.size_scalar+1:end);

DpsiDF = zeros(length(x),1);

for i=1:alpha.vector_field.n_equations % equation
    for j=1:alpha.vector_field.n_terms(i) % element of the equation
        d=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
        const=abs(alpha.vector_field.value{i}(j));
        N=length(x);
        e=eye(N);
        
        % DDH_i = sum_{k=1,d_k>0}^N sum_{j=1,(d-e_k)_n>0}^N
        % (d-e_k)_n d_k x^{d-e_k-e_n} * dir1^e_k * dir2^e_n
        
        % special case: second derivative w.r.t. psi
        for k = 1
            power = d - e(:,1) - e(:,1);
            if any(power<0)
                continue
            end
            has_delay = (alpha.delay{i}{j}~=0); 
            delay = alpha.delay{i}{j};
            delay_v_norm = norm_v.^has_delay;
            temp_dpsiDF_firsterm = const*d(1)*(d(1)-1) * x_norm(1)^(d(1)-2) * prod(x_norm(2:end).^power(2:end))*prod(v_norm.*has_delay);
            delay_derivative=0;
            for term = 1:x.size_vector
                delay_derivative = delay_derivative + delay(term)*Kv_norm(term)*prod(v_norm(1:term-1).*has_delay(1:term-1))*prod(v_norm(term+1:end).*has_delay(term+1:end));
            end
            temp_dpsiDF_secondterm =const*d(1)* x_norm(1)^(d(1)-1) * prod(x_norm(2:end).^power(2:end))*delay_derivative;
            temp_dpsiDF_thirdterm = const * prod(x_norm.^power) 
        end
        for k=2:length(x)
            %if d(k)<=0
            %    continue
            %end
            % second derivative, not needed here for n=1:length(x)
                %if d(n)-e(n,k)<=0
                %    continue
                %end
                power=d-e(k,:);
                power_scalar=power(1:size_scal);
                power_vector=power(size_scal+1:end);
                if any(power<0)
                    continue
                end
                if k<=size_scal
                    temp_dir1=dir1.scalar(k);%prod(dir1.scalar.^(e(k,1:size_scal)));
                else
                    temp_dir1=dir1_conv(k-size_scal,:);%product_FFT(dir1_conv,e(k,size_scal+1:end));
                end
                has_delay = (alpha.delay{i}{j}(k)~=0); 
                derivative_delay = prod(Kx.^has_delay);
                temp_DDF(i,:)= temp_DDF(i,:) + const*d(k)*...
                    prod(lambda_norm.^(power_scalar.'))*...
                    prod(v_norm.^(power_vector.')).*...
                    derivative_delay;
            end
        end
    end
end


% second derivative of G(x)
temp_DDG=zeros(alpha.scalar_equations.num_equations,length(x_conv));
if isintval(temp_DDF)
    temp_DDG = intval(temp_DDG);
elseif isintval(alpha.scalar_equations.polynomial_equations.value(1))
    temp_DDF = intval(temp_DDF);
    temp_DDG = intval(temp_DDG);
end


for i=1:alpha.scalar_equations.number_equations_pol % equation
    % alpha.scalar_equations.number_equations_lin+
    for j=1:alpha.scalar_equations.polynomial_equations.n_terms(i) % element of the equation
        d=[alpha.scalar_equations.polynomial_equations.power_scalar{i}(:,j).', alpha.scalar_equations.polynomial_equations.power_vector{i}{j}.'];
        const=alpha.scalar_equations.polynomial_equations.value{i}(j);
        N=length(d);
        e=eye(N);
        
        % DDH_i = sum_{k=1,d_k>0}^N sum_{j=1,(d-e_k)_n>0}^N
        % (d-e_k)_n d_k x^{d-e_k-e_n} * dir1^e_k * dir2^e_n
        
        for k=1:N
            %if d(k)<=0
            %    continue
            %end
            for n=1:N
                %if d(n)-e(n,k)<=0
                %    continue
                %end
                power=((d-e(k,:)-e(n,:)).');
                power_scalar=power(1:size_scal);
                power_vector=power(size_scal+1:end);
                if any(power<0)
                    continue
                end
                if k<=size_scal
                    temp_dir1=dir1.scalar(k);%prod(dir1.scalar.^(e(k,1:size_scal)));
                else
                    temp_dir1=dir1_conv(k-size_scal,:);%product_FFT(dir1_conv,e(k,size_scal+1:end));
                end
                if n<=size_scal
                    temp_dir2=dir2.scalar(n);%prod(dir2.scalar.^(e(n,1:size_scal)));
                else
                    temp_dir2=dir2_conv(n-size_scal,:);%product_FFT(dir2_conv,e(n,size_scal+1:end));
                end
                temp_DDG(alpha.scalar_equations.number_equations_lin+i,:)= ...
                    temp_DDG(alpha.scalar_equations.number_equations_lin+i,:) + const*d(k)*(d(n)-e(k,n))*...
                    prod(xBar.scalar.^(power_scalar.'))*...
                    product_FFT(x_conv,power_vector).*...
                    temp_dir1.*...
                    temp_dir2;
            end
        end
    end
end


end

function block = block_norm(M, P, Q, R, D3F2, psi, x)
size_scalar = x.size_scalar;
size_vector = x.size_vector;
nodes = x.nodes;

norm_P = norm_C_to_ell1(P, nodes, size_scalar, size_vector);
norm_Q = norm_ell1_to_C(Q, nodes, size_scalar, size_vector);
norm_M = norm_ell1_to_ell1(M, nodes, size_vector);
norm_R = abs(R);
norm_V1 = norm_P * abs(D3F2) / abs(psi);
norm_V2 = norm_R * abs(D3F2) / abs(psi);
norm_V3 = 1 / abs(psi);

block = [ norm_M + norm_V1 + norm_V3, norm_P
    norm_Q+norm_V2, n_R];
end

function upper_bound_psi_DF()
end

function upper_bound_DDF_without_psi()
end

function n = norm_C_to_ell1(P, nodes, size_scalar, size_vector)
global use_intval
global nu
n = zeros(size_vector,size_scalar);
if isintval(P) || use_intval
    n = intval(n);
end
K = (-nodes:nodes).';
K_big = repmat(K, 1, size_scalar);
for i = 1:size_vector
    indeces = (1:2*nodes+1) + (i-1) * (2*nodes+1);
    n(i,:) = sum(abs(P(indeces, :)).* nu.^(abs(K_big)),1);
end
end

function n = norm_ell1_to_C(Q, nodes, size_scalar, size_vector)
global use_intval
global nu
n = zeros(size_scalar, size_vector);
if isintval(Q) || use_intval
    n = intval(n);
end
K = (-nodes:nodes);
K_big = repmat(K, size_scalar, 1);
for i = 1:size_vector
    indeces = (1:2*nodes+1) + (i-1) * (2*nodes+1);
    n(:,i) = max(abs(Q(:, indeces)).* nu.^(-abs(K_big)), [], 2);
end
end

function n = norm_ell1_to_ell1(M, nodes, size_vector)
global use_intval
global nu
n = zeros(size_vector, size_vector);
if isintval(M) || use_intval
    n = intval(n);
end
K = (-nodes:nodes);
for i = 1:size_vector
    indeces_i = (1:2*nodes+1) + (i-1) * (2*nodes+1);
    for j = 1:size_vector
        indeces_j = (1:2*nodes+1) + (j-1) * (2*nodes+1);
        n(i,j) = max(nu.^(abs(K))*abs(M(indeces_i, indeces_j)).* nu.^(-abs(K)));
    end
end
end