classdef operator
    properties
        domain
        codomain
        matrix
    end
    methods
        function A = operator(mat, dom, codom)
            if ~strcmp(dom,'Cn') && ~strcmp(dom,'l1') && ~strcmp(dom,'L1')
                error('Domain not supported, at the moment just Cn l1 and L1 supported')
            end
            if ~strcmp(codom,'Cn') && ~strcmp(codom,'l1') && ~strcmp(codom,'L1')
                error('Domain not supported, at the moment just C l1 and L1 supported')
            end
            if strcmp(dom,codom) && size(A,1)~=size(A,2)
                error('If end space equal start space, expected square matrix')
            end
            A.matrix = mat;
            A.domain = dom;
            A.codomain = codom;
        end
        
        function C = plus(A,B)
            if ~strcmp(A.domain,B.domain) ||~strcmp(A.codomain,B.codomain)
                error('To sum two operators they need to be compatible')
            end
            C = A;
            C.matrix = A.matrix + B.matrix;
        end
        
        function C = max(A,B)
            if ~strcmp(A.domain,B.domain) ||~strcmp(A.codomain,B.codomain)
                error('To sum two operators they need to be compatible')
            end
            C = A;
            C.matrix = max(A.matrix, B.matrix);
        end
        
        function C = minus(A,B)
            if ~strcmp(A.domain,B.domain) ||~strcmp(A.codomain,B.codomain)
                error('To substract two operators they need to be compatible')
            end
            C = A;
            C.matrix = A.matrix - B.matrix;
        end
        
        function C = intval(A,B)
            if ~strcmp(A.domain,B.domain) ||~strcmp(A.codomain,B.codomain)
                error('To substract two operators they need to be compatible')
            end
            C = A;
            %
            % matrix_inf_real = min(real(A.matrix), real(B.matrix));
            % matrix_sup_real = max(real(A.matrix), real(B.matrix));
            %
            % matrix_inf_imag = min(imag(A.matrix), imag(B.matrix));
            % matrix_sup_imag = max(imag(A.matrix), imag(B.matrix));
            %
            % matrix_inf = matrix_inf_real + 1i*matrix_inf_imag;
            % matrix_sup = matrix_sup_real + 1i*matrix_sup_imag;
            %
            C.matrix = infsup_help (A.matrix, B.matrix);
        end
        
        function C = min(A,B)
            if ~strcmp(A.domain,B.domain) ||~strcmp(A.codomain,B.codomain)
                error('To compare two operators they need to be compatible')
            end
            C = A;
            if any(size(A.matrix)~=size(B.matrix))
                error('At the moment not coded yet') %banana!
            end
            C.matrix = min(A.matrix, B.matrix);
        end
        
        
        function C = abs(A)
            C=A;
            C.matrix = abs(A.matrix);
        end
        
        function C = mtimes(A,B)
            if isa(A,'operator') && isa(B,'operator')
                if ~strcmp(A.domain,B.codomain)
                    error('Incompatible operators')
                end
                C = A;
                C.matrix = A.matrix * B.matrix;
                C.domain = B.domain;
            elseif isa(A,'operator') && isa(B,'double')
                if all(size(B)==1)
                    C = A;
                    C.matrix = B*A.matrix;
                else
                    try
                        C = A.matrix*B.';
                        C=C.';
                    catch
                        C = A.matrix*B;
                    end
                end
                
            elseif isa(B,'operator') && isa(A,'double')
                temp = A;
                A=B;
                B=temp;
                if all(size(B)==1)
                    C = A;
                    C.matrix = B*A.matrix;
                else
                    try
                        C = A.matrix*B.';
                    catch
                        C = A.matrix*B;
                    end
                end
            else
                error('Unrecognized input')
            end
        end
        
        function n = cnorm(A,node_space)
            global nu
            if ~strcmp(A.domain, 'Cn')
                error('not implemented')
            end
            if isintval(A.matrix)
                A.matrix = sup(abs(A.matrix));
            end
            if strcmp(A.codomain, 'Cn')
                n=abs(A.matrix);
                return
            end
            if strcmp(A.codomain, 'l1')
                node_space = (size(A.matrix,1)-1)/2;
                dim_C = size(A.matrix,2);
                K2 = -node_space:node_space;
                norm_nu = @(x) sum( nu(end).^(abs(K2).').*abs(x));
                n=zeros(1,dim_C);
                for i =1:dim_C
                    n(i) = norm_nu(A.matrix(:,i));
                end
                return
            end
            if strcmp(A.codomain, 'L1')
                dim_C = size(A.matrix,2);
                node_time = (size(A.matrix,1)/(2*node_space+1)-1)/2;
                K2 = -node_space:node_space;
                K1 = -node_time:node_time;
                K1 = K1.';
                nu_K_1and2 = nu(1).^abs(repmat(K1,1,1+2*node_space)).*...
                    nu(end).^abs(repmat(K2,1+2*node_time,1));
                nu_K_1and2_vec = reshape(nu_K_1and2,1,[]);
                norm_ell = @(x) sum( nu_K_1and2_vec.'.*abs(x));
                n=zeros(1,dim_C);
                for i =1:dim_C
                    n(i) = norm_ell(A.matrix(:,1));
                end
                return
            end
        end
        
        function n = norm(A,node_space)
            % function n = norm(A,node_space)
            %
            % applying the right norm for the requested operator
            global nu
            
            if isintval(A.matrix)
                A.matrix = sup(abs(A.matrix));
            end
            
            if strcmp(A.domain, 'Cn') && strcmp(A.codomain, 'Cn')
                n=sum(abs(A.matrix),2);
                return
            end
            if strcmp(A.domain, 'Cn') && strcmp(A.codomain, 'l1')
                node_space = (size(A.matrix,1)-1)/2;
                dim_C = size(A.matrix,2);
                K2 = -node_space:node_space;
                norm_nu = @(x) sum( nu(end).^(abs(K2).').*abs(x));
                n=0;
                for i =1:dim_C
                    n = n+ norm_nu(A.matrix(:,i));
                end
                return
            end
            if strcmp(A.domain, 'Cn') && strcmp(A.codomain, 'L1')
                dim_C = size(A.matrix,2);
                node_time = (size(A.matrix,1)/(2*node_space+1)-1)/2;
                K2 = -node_space:node_space;
                K1 = -node_time:node_time;
                K1 = K1.';
                nu_K_1and2 = nu(1).^abs(repmat(K1,1,1+2*node_space)).*...
                    nu(end).^abs(repmat(K2,1+2*node_time,1));
                nu_K_1and2_vec = reshape(nu_K_1and2,1,[]);
                norm_ell = @(x) nu_K_1and2_vec*abs(x);
                n=0;
                for i =1:dim_C
                    n = n+norm_ell(A.matrix(:,1));
                end
                return
            end
            if strcmp(A.domain, 'l1') && strcmp(A.codomain, 'Cn')
                dim_C = size(A.matrix,1);
                node_space = (size(A.matrix,2)-1)/2;
                K2 = -node_space:node_space;
                norm_nuInf = @(x) max(abs(x).*(nu(end).^-abs(K2)));
                n=zeros(dim_C,1);
                for i =1:dim_C
                    n(i) = norm_nuInf(A.matrix(i,:));
                end
                return
            end
            if strcmp(A.domain, 'l1') && strcmp(A.codomain, 'l1')
                node_space = (size(A.matrix,1)-1)/2;
                K2 = -node_space:node_space;
                norm_Bnu = @(x) max((nu(end).^abs(K2)*abs(x)).*nu(end).^(-abs(K2)));
                n = norm_Bnu(A.matrix);
                return
            end
            if strcmp(A.domain, 'l1') && strcmp(A.codomain, 'L1')
                node_space = (size(A.matrix,2)-1)/2;
                node_time =  (size(A.matrix,1)/(2*node_space+1)-1)/2;
                K2 = -node_space:node_space;
                K1 = -node_time:node_time;
                K1 = K1.';
                nu_K_1and2 = nu(1).^abs(repmat(K1,1,1+2*node_space)).*...
                    nu(end).^abs(repmat(K2,1+2*node_time,1));
                nu_K_1and2_vec = reshape(nu_K_1and2,1,[]);
                norm_nu2ell = @(x) max((nu_K_1and2_vec*abs(x)).*nu(end).^(-abs(K2)));
                n = norm_nu2ell(A.matrix);
                return
            end
            if strcmp(A.domain, 'L1') && strcmp(A.codomain, 'Cn')
                dim_C = size(A.matrix,1);
                node_time =  (size(A.matrix,2)/(2*node_space+1)-1)/2;
                K2 = -node_space:node_space;
                K1 = -node_time:node_time;
                K1 = K1.';
                nu_K_1and2 = nu(1).^abs(repmat(K1,1,1+2*node_space)).*...
                    nu(end).^abs(repmat(K2,1+2*node_time,1));
                nu_K_1and2_vec = reshape(nu_K_1and2,1,[]);
                norm_ellInf = @(x) max(abs(x).*(nu_K_1and2_vec.^-1));
                n=zeros(dim_C,1);
                for i =1:dim_C
                    n(i) = norm_ellInf(A.matrix(i,:));
                end
                return
            end
            if strcmp(A.domain, 'L1') && strcmp(A.codomain, 'l1')
                node_space = (size(A.matrix,1)-1)/2;
                node_time =  (size(A.matrix,2)/(2*node_space+1)-1)/2;
                K2 = -node_space:node_space;
                K1 = -node_time:node_time;
                K1 = K1.';
                nu_K_1and2 = nu(1).^abs(repmat(K1,1,1+2*node_space)).*...
                    nu(end).^abs(repmat(K2,1+2*node_time,1));
                nu_K_1and2_vec = reshape(nu_K_1and2,1,[]);
                norm_ell2nu = @(x) max((nu(end).^abs(K2)*abs(x)).*(nu_K_1and2_vec.^-1));
                n = norm_ell2nu(A.matrix);
                return
            end
            if strcmp(A.domain, 'L1') && strcmp(A.codomain, 'L1')
                node_time = (size(A.matrix,1)/(2*node_space+1)-1)/2;
                K2 = -node_space:node_space;
                K1 = -node_time:node_time;
                K1 = K1.';
                nu_K_1and2 = nu(1).^abs(repmat(K1,1,1+2*node_space)).*...
                    nu(end).^abs(repmat(K2,1+2*node_time,1));
                nu_K_1and2_vec = reshape(nu_K_1and2,1,[]);
                norm_Bell = @(x) max((nu_K_1and2_vec*abs(x)).*(nu_K_1and2_vec.^-1));
                n = norm_Bell(A.matrix);
                return
            end
            error('Unrecognized space')
        end
        
        function I = identity(A)
            I=A;
            I.matrix = eye(size(A.matrix));
        end
        
        function I = zeros(A)
            I = A;
            I.matrix = 0*I.matrix;
        end
        
        function A_big = expand(A,new_node_space, old_node_space, new_node_time)
            
            if nargin== 4
                [old_rows,old_cols,new_n_rows, new_n_cols] = expantion_indeces...
                    (A,new_node_space, old_node_space, new_node_time);
            elseif nargin== 3
                [old_rows,old_cols,new_n_rows, new_n_cols] = expantion_indeces...
                    (A,new_node_space, old_node_space);
            elseif nargin== 2
                [old_rows,old_cols,new_n_rows, new_n_cols] = expantion_indeces...
                    (A,new_node_space);
            else
                error('Not enough input arguments')
            end
            A_big = A;
            A_big.matrix = zeros(new_n_rows, new_n_cols);
            if isa(A.matrix,'intval')
                A_big.matrix = intval(A_big.matrix);
            end
            A_big.matrix(old_rows, old_cols) = A.matrix;
        end
        
        
        function [old_rows,old_cols,new_n_rows, new_n_cols] = expantion_indeces...
                (A,new_node_space, old_node_space, new_node_time)
            % function [old_rows,old_cols,new_n_rows, new_n_cols] = expantion_indeces...
            %    (A,new_node_space, old_node_space, new_node_time)
            caller = dbstack(1);
            caller_func = caller(1).name;
            switch A.codomain % rows
                case 'Cn'
                    old_n_rows = size(A.matrix,1);
                    new_n_rows = old_n_rows;
                    old_rows = 1:old_n_rows;
                case 'l1'
                    new_n_rows = 2*new_node_space+1;
                    old_node_space = (size(A.matrix,1)-1)/2;
                    if old_node_space>new_node_space && strcmp(caller_func, 'operator.expand')
                        error('Cannot expand to the new shape')
                    end
                    old_rows = new_node_space+1 +(-old_node_space:old_node_space);
                case 'L1'
                    old_node_time = ((size(A.matrix,1)/(2*old_node_space+1))-1)/2;
                    if old_node_space>new_node_space && strcmp(caller_func, 'operator.expand')
                        error('Cannot expand to the new shape')
                    end
                    if old_node_time>new_node_time && strcmp(caller_func, 'operator.expand')
                        error('Cannot expand to the new shape')
                    end
                    
                    new_n_rows = (2*new_node_space+1)*(2*new_node_time+1);
                    
                    vec1 = (new_node_time +1) +(-old_node_time:old_node_time);
                    vec2 = (2*new_node_time+1)*(new_node_space+(-old_node_space:old_node_space));
                    mat1 = repmat(vec1, length(vec2),1);
                    mat2 = repmat(vec2.',1, length(vec1));
                    mat = (mat1+mat2).';
                    old_rows = mat(:);
            end
            switch A.domain % cols
                case 'Cn'
                    old_n_cols = size(A.matrix,2);
                    new_n_cols = old_n_cols;
                    old_cols = 1:old_n_cols;
                case 'l1'
                    new_n_cols = 2*new_node_space+1;
                    old_node_space = (size(A.matrix,2)-1)/2;
                    if old_node_space>new_node_space && strcmp(caller_func, 'operator.expand')
                        error('Cannot expand to the new shape')
                    end
                    if old_node_space<new_node_space && strcmp(caller_func, 'operator.project')
                        error('Cannot project to the new shape')
                    end
                    old_cols = new_node_space+1 +(-old_node_space:old_node_space);
                case 'L1'
                    old_node_time = ((size(A.matrix,2)/(2*old_node_space+1))-1)/2;
                    if old_node_space>new_node_space && strcmp(caller_func, 'operator.expand')
                        error('Cannot expand to the new shape')
                    end
                    if old_node_time>new_node_time && strcmp(caller_func, 'operator.expand')
                        error('Cannot expand to the new shape')
                    end
                    if old_node_space<new_node_space && strcmp(caller_func, 'operator.project')
                        error('Cannot project to the new shape')
                    end
                    if old_node_time<new_node_time && strcmp(caller_func, 'operator.project')
                        error('Cannot project to the new shape')
                    end
                    
                    new_n_cols = (2*new_node_space+1)*(2*new_node_time+1);
                    
                    vec1 = (new_node_time +1) +(-old_node_time:old_node_time);
                    vec2 = (2*new_node_time+1)*(new_node_space+(-old_node_space:old_node_space));
                    mat1 = repmat(vec1, length(vec2),1);
                    mat2 = repmat(vec2.',1, length(vec1));
                    mat = (mat1+mat2).';
                    old_cols = mat(:);
            end
        end
        
        function [old_rows,old_cols,new_n_rows, new_n_cols] = projection_indeces...
                (A,new_node_space, old_node_space, new_node_time)
            % function [old_rows,old_cols,new_n_rows, new_n_cols] = projection_indeces...
            %    (A,new_node_space, old_node_space, new_node_time)
            caller = dbstack(1);
            caller_func = caller(1).name;
            switch A.codomain % rows
                case 'Cn'
                    old_n_rows = size(A.matrix,1);
                    new_n_rows = old_n_rows;
                    old_rows = 1:old_n_rows;
                case 'l1'
                    new_n_rows = 2*new_node_space+1;
                    old_node_space = (size(A.matrix,1)-1)/2;
                    if old_node_space<new_node_space && strcmp(caller_func, 'operator.project')
                        error('Cannot expand to the new shape')
                    end
                    old_rows = old_node_space+1 +(-new_node_space:new_node_space);
                case 'L1'
                    old_node_time = ((size(A.matrix,1)/(2*old_node_space+1))-1)/2;
                    if old_node_space<new_node_space && strcmp(caller_func, 'operator.project')
                        error('Cannot expand to the new shape')
                    end
                    if old_node_time<new_node_time && strcmp(caller_func, 'operator.project')
                        error('Cannot expand to the new shape')
                    end
                    
                    new_n_rows = (2*new_node_space+1)*(2*new_node_time+1);
                    
                    vec1 = (old_node_time +1) +(-new_node_time:new_node_time);
                    vec2 = (2*old_node_time+1)*(old_node_space+(-new_node_space:new_node_space));
                    mat1 = repmat(vec1, length(vec2),1);
                    mat2 = repmat(vec2.',1, length(vec1));
                    mat = (mat1+mat2).';
                    old_rows = mat(:);
            end
            switch A.domain % cols
                case 'Cn'
                    old_n_cols = size(A.matrix,2);
                    new_n_cols = old_n_cols;
                    old_cols = 1:old_n_cols;
                case 'l1'
                    new_n_cols = 2*new_node_space+1;
                    old_node_space = (size(A.matrix,2)-1)/2;
                    if old_node_space<new_node_space && strcmp(caller_func, 'operator.project')
                        error('Cannot expand to the new shape')
                    end
                    old_cols = old_node_space+1 +(-new_node_space:new_node_space);
                case 'L1'
                    old_node_time = ((size(A.matrix,2)/(2*old_node_space+1))-1)/2;
                    if old_node_space<new_node_space && strcmp(caller_func, 'operator.project')
                        error('Cannot expand to the new shape')
                    end
                    if old_node_time<new_node_time && strcmp(caller_func, 'operator.project')
                        error('Cannot expand to the new shape')
                    end
                    
                    new_n_cols = (2*new_node_space+1)*(2*new_node_time+1);
                    
                    vec1 = (old_node_time +1) +(-new_node_time:new_node_time);
                    vec2 = (2*old_node_time+1)*(old_node_space+(-new_node_space:new_node_space));
                    mat1 = repmat(vec1, length(vec2),1);
                    mat2 = repmat(vec2.',1, length(vec1));
                    mat = (mat1+mat2).';
                    old_cols = mat(:);
            end
        end
        
        
        function A_small = project(A,new_node_space, old_node_space, new_node_time)
            %             switch A.codomain % rows
            %                 case 'Cn'
            %                     old_n_rows = size(A.matrix,1);
            %                     new_n_rows = old_n_rows;
            %                     old_rows = 1:old_n_rows;
            %                 case 'l1'
            %                     new_n_rows = 2*new_node_space+1;
            %                     old_node_space = (size(A.matrix,1)-1)/2;
            %                     if old_node_space<new_node_space
            %                         error('Cannot project to the new shape')
            %                     end
            %                     old_rows = old_node_space+1 +(-new_node_space:new_node_space);
            %                 case 'L1'
            %                     old_node_time = ((size(A.matrix,1)/(2*old_node_space+1))-1)/2;
            %                     if old_node_space<new_node_space
            %                         error('Cannot project to the new shape')
            %                     end
            %                     if old_node_time<new_node_time
            %                         error('Cannot project to the new shape')
            %                     end
            %                     new_n_rows = (2*new_node_space+1)*(2*new_node_time+1);
            %                     span_old = 2*new_node_space*new_node_time + new_node_space+new_node_time;
            %                     old_node_center =  2*old_node_space*old_node_time + old_node_space+old_node_time;
            %                     old_rows = old_node_center+1 +...
            %                         (-span_old:span_old);
            %             end
            %             switch A.domain % cols
            %                 case 'Cn'
            %                     old_n_cols = size(A.matrix,2);
            %                     new_n_cols = old_n_cols;
            %                     old_cols = 1:old_n_cols;
            %                 case 'l1'
            %                     new_n_cols = 2*new_node_space+1;
            %                     old_node_space = (size(A.matrix,2)-1)/2;
            %                     if old_node_space<new_node_space
            %                         error('Cannot project to the new shape')
            %                     end
            %                     old_cols = old_node_space+1 +(-new_node_space:new_node_space);
            %                 case 'L1'
            %                     old_node_time = ((size(A.matrix,2)/(2*old_node_space+1))-1)/2;
            %                     if old_node_space<new_node_space
            %                         error('Cannot project to the new shape')
            %                     end
            %                     if old_node_time<new_node_time
            %                         error('Cannot project to the new shape')
            %                     end
            %                     new_n_cols = (2*new_node_space+1)*(2*new_node_time+1);
            %                     span_old = 2*new_node_space*new_node_time + new_node_space+new_node_time;
            %                     old_node_center =  2*old_node_space*old_node_time + old_node_space+old_node_time;
            %                     old_cols = old_node_center+1 +...
            %                         (-span_old:span_old);
            %             end
            if nargin== 4
                [old_rows,old_cols,new_n_rows, new_n_cols] = projection_indeces...
                    (A,new_node_space, old_node_space, new_node_time);
            elseif nargin== 3
                [old_rows,old_cols,new_n_rows, new_n_cols] = projection_indeces...
                    (A,new_node_space, old_node_space);
            elseif nargin== 2
                [old_rows,old_cols,new_n_rows, new_n_cols] = projection_indeces...
                    (A,new_node_space);
            else
                error('Not enough input arguments')
            end
            A_small = A;
            A_small.matrix = zeros(new_n_rows, new_n_cols);
            A_small.matrix = A.matrix(old_rows, old_cols);
        end
        
        function T = tail(A, center_node_space, all_node_space, center_node_time)
            if ~strcmp(A.domain,'L1') && ~strcmp(A.codomain, 'L1')
                T = A;
                B = project(A,center_node_space);
                B.matrix = 4+0*B.matrix;
                B = expand(B,all_node_space);
                B.matrix(B.matrix ==0) = 1;
                B.matrix(B.matrix == 4) = 0;
                T.matrix = A.matrix.*B.matrix;
            else
                if strcmp(A.codomain, 'L1')
                    size_A = size(A.matrix,1);
                else
                    size_A = size(A.matrix,2);
                end
                all_node_time = (size_A/(2*all_node_space+1)-1)/2;
                T = A;
                B = project(A,center_node_space, ...
                    all_node_space, center_node_time);
                B.matrix = 4+0*B.matrix;
                B = expand(B,all_node_space, ...
                    center_node_space, all_node_time);
                B.matrix(B.matrix ==0) = 1;
                B.matrix(B.matrix == 4) = 0;
                T.matrix = A.matrix.*B.matrix;
            end
        end
        
        function ND = non_diag(A)
            ND = A;
            ND.matrix(1:1+size(A.matrix+1,1):end) = 0; % put diagonal to zero
        end
        
        function M = op2mat(A)
            M = A.matrix;
        end
    end
end