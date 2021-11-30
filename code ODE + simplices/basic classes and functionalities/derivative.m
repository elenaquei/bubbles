classdef derivative
    properties %(SetAccess=private)
        size_scalar
        size_vector
        deg_vector
        nodes
        derivative_Glambda   % M x M
        derivative_Gx_v      % M1 x N*(2*nodes+1)    %% *NOT* a Xi_vector
        derivative_G_hat     % M2 x N
        derivative_Flambda   % M x N x (2*nodes+1)*deg
        derivative_Fx_diagonal % N x N
        derivative_Fx_toeplix  % N x N x (2*nodes+1)*deg
        derivative_Fx_delay   % N x N cell - each cell with as many cells 
                                    % as terms, each term with two elements: 
                                    % convolution and exp_coef
    end
    methods
        % DERIVATIVE
        function DF = derivative(alpha, x, bool)
            % function DF = derivative(alpha, x, bool)
            % 
            % INPUT
            % alpha     full_problem
            % x         Xi_vector
            % bool      boolean, if 0 output of same dimension of input
            %           DEFAULT 1
            %
            % OUTPUT
            % DF        instance of class derivative
            global talkative
            if ~compatible_vec(alpha,x)
                error('The inputs are incompatible')
            end
            if ~isa(alpha,'full_problem')
                error('wrong call of function, a full_problem instance is requested')
            end
            if alpha.scalar_equations.number_of_nodes ~= x.nodes
                if talkative>1
                disp('In the computation of the derivative,nodes of inputs differ, the linear coefficients are padded or truncated to fit')
                end
            end
            if nargin<3
                bool = 1;
            end
            DF.size_scalar=alpha.scalar_equations.size_scalar;
            DF.size_vector=alpha.scalar_equations.size_vector;
            DF.deg_vector= alpha.vector_field.deg_vector;
            DF.nodes = x.nodes * DF.deg_vector;
            alpha.scalar_equations = reshape(alpha.scalar_equations,DF.nodes);
            [scalar_der_G,vector_der_G, vector_der_G_hat]= compute_derivative(alpha.scalar_equations,x);
            DF.derivative_Glambda= scalar_der_G;
            DF.derivative_Gx_v = vector_der_G; % still in the form of (n_eqs, size_x, coefs)
            DF.derivative_G_hat = vector_der_G_hat;
            
            [scalar_der_F ,diagonal_der_F, toeplix_der_F, delay_der_F] = ...
                compute_derivative(alpha.vector_field,x);
%             DF.derivative_Flambda = 0*scalar_der_F; % SIZE xi_vec.size_scalar,a.n_equations,xi_vec.nodes*2*a.deg_vector+1
%             for i =1:DF.size_scalar
%                 for j = 1:alpha.vector_field.n_equations
%                     DF.derivative_Flambda(i,j,:) =  verifyfft_in(squeeze(scalar_der_F(i,j,:)),1);
%                 end
%             end
            DF.derivative_Flambda = scalar_der_F;
            
            DF.derivative_Fx_diagonal = diagonal_der_F;
            DF.derivative_Fx_delay = delay_der_F;
            
            DF.derivative_Fx_toeplix = toeplix_der_F; % SIZE xi_vec.size_vector,a.n_equations,xi_vec.nodes*2*a.deg_vector+1
            for i =1:DF.size_vector
                for j = 1:alpha.vector_field.n_equations
                    DF.derivative_Fx_toeplix(i,j,:) =  verifyfft_in(squeeze(toeplix_der_F(i,j,:)),1);
                end
            end
            
            if bool ==0 
                DF = reshape(DF, x.nodes);
            else
                DF = reshape(DF, x.nodes*(alpha.vector_field.deg_vector));
            end
        end
        % DERIVATIVE
        
        % DERIVATIVE_TO_MATRIX
        function A_mat = derivative_to_matrix(A, nodes)
            % function A_mat = derivative_to_matrix(A, nodes)
            % INPUT
            % A      derivative
            % nodes  number of nodes requested, DEFAULT A.n_nodes
            % OUTPUT
            % A_mat  complex matrix, of size M + N * (2*nodes+1)
            global use_intlab
            
            if nargin >1
                A = reshape(A,nodes);
            end
            
            M = size(A.derivative_Flambda,1);
            N= size(A.derivative_Flambda,2);
            
            One = ones(1,A.nodes*2+1);
            Diag = diag(-A.nodes:A.nodes);%repmat(diag(-A.nodes:A.nodes),N,N);
            
            DGx_v = reshape(permute(A.derivative_Gx_v,[1,3,2]),[],N*(A.nodes*2+1));
            
            if ~isintval(A.derivative_G_hat)
                DxG = [DGx_v;kron(A.derivative_G_hat,One)]; 
            else
                %DxG = [DGx_v;kron(A.derivative_G_hat,One)];
                sup_A = kron(sup(A.derivative_G_hat),One);
                inf_A = kron(inf(A.derivative_G_hat),One);
                int_A = infsup(inf_A,sup_A);
                DxG = [DGx_v;int_A];
            end
            DlF = reshape(permute(A.derivative_Flambda,[1,3,2]),[M, N*(2*A.nodes+1),1]).';
            DxF = zeros(N*(2*A.nodes+1));
            if use_intlab || isintval(A.derivative_Flambda)
                DxF = intval(DxF);
            end
            index = 1:(2*A.nodes+1);
            for i =1:N
                for j = 1:N
                    DxF((i-1)*length(index)+index,(j-1)*length(index) +index) = ...
                        toeplitz(...
                        [squeeze(A.derivative_Fx_toeplix(j,i,A.nodes+1:end));zeros(A.nodes,1)],...
                        [squeeze(A.derivative_Fx_toeplix(j,i,A.nodes+1:-1:1));zeros(A.nodes,1)]...
                        );
                    delay_part = zeros(size(Diag,1));
                    for term = 1:size(A.derivative_Fx_delay{i,j})
                        if ~isempty(A.derivative_Fx_delay{i,j}{term})
                        conv = A.derivative_Fx_delay{i,j}{term}.convolution;
                        exp_coef = A.derivative_Fx_delay{i,j}{term}.exp_coef;
                        conv_mat = toeplitz(...
                        [conv(A.nodes+1:end);zeros(A.nodes,1)],...
                        [conv(A.nodes+1:-1:1);zeros(A.nodes,1)]...
                        );
                        exp_mat_delay = diag(exp(exp_coef*(-A.nodes:A.nodes)));
                        delay_part = delay_part + conv_mat*exp_mat_delay;
                        end
                    end
                    DxF((i-1)*length(index)+index,(j-1)*length(index) +index) =...
                        DxF((i-1)*length(index)+index,(j-1)*length(index) +index) + delay_part;
                end
            end
            if ~isintval(A.derivative_Fx_diagonal)
                DxF = DxF+kron( 1i*A.derivative_Fx_diagonal, Diag);
            else
                %DxG = [DGx_v;kron(A.derivative_G_hat,One)];
                sup_A = kron(sup(1i*A.derivative_Fx_diagonal),Diag);
                inf_A = kron(inf(1i*A.derivative_Fx_diagonal),Diag);
                int_A = infsup(inf_A,sup_A);
                clear inf_A sup_A;
                DxF = DxF+int_A;
                clear int_A;
            end
            
            
            A_mat = [A.derivative_Glambda , DxG 
                DlF, DxF];
            
            if size(A_mat,1)~=size(A_mat,2)
                % debugging option
                %warning('system not square')
            end
            if size(A_mat,2)~= M+N*(2*A.nodes+1)
                warning('PROBABLY BUG')
            end
        end
        % end DERIVATIVE_TO_MATRIX
        
        % RESHAPE
        function B = reshape(A,n_nodes)
            % function B = reshape(A,n_nodes)
            %
            % INPUT
            % A          instance of derivative
            % n_nodes    integer
            % OUTPUT
            % B          derivative with number of nodes equal to n_nodes
            B = truncate(A,n_nodes);
            B = pad( B, n_nodes);
        end
        % end RESHAPE
        
        % TRUNCATE
        function B = truncate(A,n_node)
            % function B = truncate(A,n_node)
            
            % if A.nodes < n_node
            %     error('No truncation necessary, but padding')
            % elseif A.nodes == n_node
            %     B=A;
            %     return
            % end
            used_nodes =-n_node:n_node;
            B = A;
            
            if size(B.derivative_Gx_v,3) > 2*n_node+1
                center_node =floor(size(B.derivative_Gx_v,3)/2)+1;
                B.derivative_Gx_v = B.derivative_Gx_v(:,:,center_node+used_nodes);
            end
            
            if size(B.derivative_Flambda,3)> 2*n_node +1
                center_node = floor(size(B.derivative_Flambda,3)/2)+1;   % M x N x (2*nodes+1)*deg
                B.derivative_Flambda = B.derivative_Flambda(:,:,center_node+used_nodes);
            end
            
            if size(B.derivative_Fx_toeplix,3)> 2*n_node +1
                Size= size(B.derivative_Fx_toeplix,3);
                if mod(Size,2) ==0
                    B.derivative_Fx_toeplix = B.derivative_Fx_toeplix(:,:,2:end);
                end
                center_node =floor(size(B.derivative_Fx_toeplix,3)/2)+1;
                B.derivative_Fx_toeplix = B.derivative_Fx_toeplix(:,:,center_node+used_nodes); %derivative_Fx_toeplix   N x N x (2*nodes+1)*deg
            end
            
            % all delay derivatives also need truncation!
            for i = 1:size(A.derivative_Fx_delay,1)
                for j = 1:size(A.derivative_Fx_delay, 2)
                    for k = 1:size(A.derivative_Fx_delay{i,j})
                        if ~isempty(A.derivative_Fx_delay{i,j}{k})
                            if length(A.derivative_Fx_delay{i,j}{k}.convolution)>2*n_node +1
                                center_node =floor(length(A.derivative_Fx_delay{i,j}{k}.convolution)/2)+1;
                                B.derivative_Fx_delay{i,j}{k}.convolution = B.derivative_Fx_delay{i,j}{k}.convolution(center_node+used_nodes); 
                            end
                        end
                    end
                end
            end
            
            B.nodes = n_node;
        end
        % end TRUNCATE
        
        % PAD
        function B = pad(A, n_node)
            % function B = pad(A, n_node)
            
            % if A.nodes > n_node
            %     error('No padding necessary, but truncation')
            % elseif A.nodes == n_node
            %     B=A;
            %     return
            % end
            
            B = A;
            M1 = size(A.derivative_Gx_v,1);
            M = size(A.derivative_Flambda,1);
            N= size(A.derivative_Flambda,2);
            if size(B.derivative_Gx_v,3)< 2*n_node +1
                center_node =ceil(size(B.derivative_Gx_v,3)/2);
                padding = -center_node + 1 + n_node;
                if ~isintval(B.derivative_Gx_v)
                    B.derivative_Gx_v=padarray(B.derivative_Gx_v,[0,0,padding]); % derivative_Gx_v  M1 x N*(2*nodes+1)
                else
                    
                    inf_pad = padarray(inf(B.derivative_Gx_v),[0,0,padding]);
                    sup_pad = padarray(sup(B.derivative_Gx_v),[0,0,padding]);
                    B.derivative_Gx_v=infsup(inf_pad,sup_pad); % derivative_Gx_v  M1 x N*(2*nodes+1)
                end
            end
            
            if size(B.derivative_Flambda,3)< 2*n_node +1
                center_node = ceil(size(B.derivative_Flambda,3)/2);   % M x N x (2*nodes+1)*deg
                padding = -center_node + 1 + n_node;
                if ~isintval(B.derivative_Flambda)
                    B.derivative_Flambda = padarray(B.derivative_Flambda,[0,0,padding]);
                else
                    
                    inf_pad = padarray(inf(B.derivative_Flambda),[0,0,padding]);
                    sup_pad = padarray(sup(B.derivative_Flambda),[0,0,padding]);
                    B.derivative_Flambda = infsup(inf_pad,sup_pad);
                end
            end
            
            if size(B.derivative_Gx_v,3) < 2*n_node+1
                B.derivative_Gx_v = reshape_Xi(B.derivative_Gx_v,n_nodes);
            end
            
            if size(B.derivative_Fx_toeplix,3)< 2*n_node +1
                center_node =ceil(size(B.derivative_Fx_toeplix,3)/2);
                padding = -center_node + 1 + n_node;
                if ~isintval(B.derivative_Fx_toeplix)
                    B.derivative_Fx_toeplix = padarray(B.derivative_Fx_toeplix,[0,0,padding]); %derivative_Fx_toeplix   N x N x (2*nodes+1)*deg
                else
                    inf_pad = padarray(inf(B.derivative_Fx_toeplix),[0,0,padding]);
                    sup_pad = padarray(sup(B.derivative_Fx_toeplix),[0,0,padding]);
                    B.derivative_Fx_toeplix = infsup(inf_pad,sup_pad);
                end
            end
            
            % all delay derivatives also need padding!
            for i = 1:size(A.derivative_Fx_delay,1)
                for j = 1:size(A.derivative_Fx_delay, 2)
                    for k = 1:size(A.derivative_Fx_delay{i,j})
                        if ~isempty(A.derivative_Fx_delay{i,j}{k})
                            if length(A.derivative_Fx_delay{i,j}{k})>2*n_node +1
                                center_node =floor(length(A.derivative_Fx_delay{i,j}{k})/2)+1;
                                padding = -center_node + 1 + n_node;
                                if ~isintval(B.derivative_Fx_toeplix)
                                    B.derivative_Fx_delay{i,j}{k} = padarray(B.derivative_Fx_delay{i,j}{k},[0,0,padding]); 
                                else
                                    inf_pad = padarray(inf(B.derivative_Fx_delay{i,j}{k}),[0,0,padding]);
                                    sup_pad = padarray(sup(B.derivative_Fx_delay{i,j}{k}),[0,0,padding]);
                                    B.derivative_Fx_delay{i,j}{k} = infsup(inf_pad,sup_pad);
                                end
                                B.derivative_Fx_delay{i,j}{k} = B.derivative_Fx_delay{i,j}{k}(center_node+used_nodes);
                            end
                        end
                    end
                end
            end
            
            B.nodes = n_node;
        end
        % end PAD
        
        
        
        % % % ---- ----- TO CODE ----- --- % %
        % function norm
        % function product
        
        
    end
end