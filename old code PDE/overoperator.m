classdef overoperator
    properties
        matrix_cell
        size_real
    end
    methods
        function B = overoperator(A11, A12, A13, A21, A22, A23, A31, A32, A33)
            if nargin <=1
                B.matrix_cell = cell(3,3);
                B.size_real = [];
                return
            end
            A = cell(3,3);
            A{1,1} = operator(A11,'Cn','Cn');
            A{1,2} = operator(A12, 'l1', 'Cn');
            A{1,3} = operator(A13, 'L1','Cn');
            
            A{2,1} = operator(A21,'Cn','l1');
            A{2,2} = operator(A22, 'l1', 'l1');
            A{2,3} = operator(A23, 'L1','l1');
            
            
            A{3,1} = operator(A31,'Cn','L1');
            A{3,2} = operator(A32, 'l1', 'L1');
            A{3,3} = operator(A33, 'L1','L1');
            
            B.matrix_cell = A;
            B.size_real = size(A11,1);
        end
        
        function C = plus(A,B)
            C=A;
            for i=1:3
                for j = 1:3
                    C.matrix_cell{i,j} = A.matrix_cell{i,j}+B.matrix_cell{i,j};
                end
            end
        end
        
        function C = minus(A,B)
            C=A;
            for i=1:3
                for j = 1:3
                    C.matrix_cell{i,j} = A.matrix_cell{i,j}-B.matrix_cell{i,j};
                end
            end
        end
        
        
        function C = intval(A,B)
            C=A;
            for i=1:3
                for j = 1:3
                    C.matrix_cell{i,j} = intval(A.matrix_cell{i,j},B.matrix_cell{i,j});
                end
            end
        end
        
        function node_space = node_space(A)
            node_space = (length(A.matrix_cell{1,2}.matrix(1,:))-1)/2;
        end
        
        function node_time = node_time(A)
            node_time = ((size(A.matrix_cell{3,3}.matrix,1) / node_space(A))-1)/2;
        end
        
        function s = size_center(A)
            s = zeros(1,2);
            s(1)= size(A.matrix_cell{2,2},1);
            s(2) = size(A.matrix_cell{3,3},1);
        end
        
        function n = norm(A)
            node_space = ( size(A.matrix_cell{2,2}.matrix,1)-1)/2;
            n=zeros(2+A.size_real,1);
            for j = 1:3
                n(1:A.size_real) = n(1:A.size_real)+norm(A.matrix_cell{1,j},node_space);
            end
            
            for i=2:3
                n_temp = zeros(1,3);
                for j = 1:3
                    n_temp(j) = norm(A.matrix_cell{i,j},node_space);
                end
                n(A.size_real+i-1) = sum(n_temp);
            end
        end
        
        function n = cnorm(A)
            node_space = ( size(A.matrix_cell{2,2}.matrix,1)-1)/2;
            n=zeros(2+A.size_real);
            n(1:A.size_real,1:A.size_real)= cnorm(A.matrix_cell{1,1},node_space);
            for j = 2:3
                n(1:A.size_real,A.size_real+j-1) = norm(A.matrix_cell{1,j},node_space);
                n(A.size_real+j-1,1:A.size_real) = cnorm(A.matrix_cell{j,1},node_space);
            end
            
            for i=2:3
                for j = 2:3
                    n(A.size_real+i-1,A.size_real+j-1) = norm(A.matrix_cell{i,j},node_space);
                end
            end
        end
        
        function B = abs(A)
            B=A;
            for i =1:3
                for j = 1:3
                    B.matrix_cell{i,j} = abs(A.matrix_cell{i,j});
                end
            end
        end
        
        function C = min(A,B)
            if ~isa(A,'overoperator') ||~isa(B,'overoperator')
                error('Comparison impossible')
            end
            C=A;
            for i =1:3
                for j = 1:3
                    C.matrix_cell{i,j} = min(A.matrix_cell{i,j},B.matrix_cell{i,j});
                end
            end
        end
        
        function C = max(A,B)
            if ~isa(A,'overoperator') ||~isa(B,'overoperator')
                error('Comparison impossible')
            end
            C=A;
            for i =1:3
                for j = 1:3
                    C.matrix_cell{i,j} = max(A.matrix_cell{i,j},B.matrix_cell{i,j});
                end
            end
        end
        
        function C = mtimes(A,B)
            if ~isa(A,'overoperator')
                temp = A;A = B; B= temp;
            end
            if isa(A,'overoperator')  && isa(B,'overoperator')
                C=A;
                for i=1:3
                    for j = 1:3
                        C.matrix_cell{i,j}.matrix = 0*C.matrix_cell{i,j}.matrix;
                        for k=1:3
                            C.matrix_cell{i,j} = C.matrix_cell{i,j}+A.matrix_cell{i,k}*B.matrix_cell{k,j};
                        end
                    end
                end
            elseif isa(A,'overoperator')  && isa(B,'PDE_vector')
                C=0*B;
                A = reshape(A,B.node_space,B.node_time);
                
                C.parameters = A.matrix_cell{1,1}*B.parameters + A.matrix_cell{1,2}*B.vector1D+...
                    A.matrix_cell{1,3}*reshape(B.vector2D,1,[]);
                
                C.vector1D = A.matrix_cell{2,1}*B.parameters + A.matrix_cell{2,2}*B.vector1D+...
                    A.matrix_cell{2,3}*reshape(B.vector2D,1,[]);
                
                vector2D = A.matrix_cell{3,1}*B.parameters + A.matrix_cell{3,2}*B.vector1D+...
                    A.matrix_cell{3,3}*reshape(B.vector2D,1,[]);
                
                C.vector2D = reshape(vector2D.', B.size_vector , [], ( 2*B.node_space + 1 ));
                
            elseif isa(A,'overoperator')  && isa(B,'double')
                C=A;
                for i=1:3
                    for j = 1:3
                        C.matrix_cell{i,j} = B*C.matrix_cell{i,j};
                    end
                end
            else
                error('unrecognized input type')
            end
        end
        
        function A_big = expand(A,new_node_space, new_node_time)
            A_big=A;
            old_node_space = (size(A.matrix_cell{2,2}.matrix,1)-1)/2;
            for i=1:3
                for j = 1:3
                    A_big.matrix_cell{i,j} = expand(A.matrix_cell{i,j},new_node_space,old_node_space, new_node_time);
                end
            end
        end
        
        function [old_rows,old_cols] = expantion_indeces(A,new_node_space, new_node_time)
            old_cols =0;
                for j = 1:3
                    [~,old_cols_temp] = expantion_indeces...
                        (A.matrix_cell{1,j},new_node_space,old_node_space, new_node_time);
                    old_cols = [old_cols, old_cols(end)+old_cols_temp];
                end
            old_cols = old_cols(2:end);
            old_rows = old_cols;
        end
        
        function A_small = project(A,new_node_space, new_node_time)
            A_small=A;
            old_node_space = (size(A.matrix_cell{2,2}.matrix,1)-1)/2;
            for i=1:3
                for j = 1:3
                    A_small.matrix_cell{i,j} = project(A.matrix_cell{i,j},new_node_space,old_node_space, new_node_time);
                end
            end
        end
        
        function A_tail = tail(A,center_node_space, center_node_time)
            A_tail=A;
            all_node_space = (size(A.matrix_cell{2,2}.matrix,1)-1)/2;
            for i=1:3
                for j = 1:3
                    A_tail.matrix_cell{i,j} = tail(A.matrix_cell{i,j},center_node_space, all_node_space, center_node_time);
                end
            end
        end
        
        function A_non_diag = non_diag(A)
            A_non_diag = A;
            for i = 1:3
                for j = 1:3
                   A_non_diag.matrix_cell{i,j} = non_diag(A.matrix_cell{i,j});
                end
            end
        end
        
        function B = reshape(A, new_node_space, new_node_time)
            old_node_space = (size(A.matrix_cell{2,2}.matrix,1)-1)/2;
            old_node_time = (size(A.matrix_cell{3,3}.matrix,1)/(2*old_node_space+1)-1)/2;
            if new_node_space==old_node_space && new_node_time==old_node_time
                B=A;
                return
            end
            if new_node_space>=old_node_space
                if new_node_time>=old_node_time
                    B = expand(A, new_node_space, new_node_time);
                else
                    B_mid = expand(A, new_node_space, old_node_time);
                    B = project(B_mid, new_node_space, new_node_time);
                end
            else
                if new_node_time<=old_node_time
                    B = project(A, new_node_space, new_node_time);
                else
                    B_mid = expand(A, old_node_space, new_node_time);
                    B = project(B_mid, new_node_space, new_node_time);
                end
            end
        end
        
        function I = identity(A)
            I = A;
            for i = 1:3
                for j =1:3
                    if i==j
                        I.matrix_cell{i,j} = identity(A.matrix_cell{i,j});
                    else
                        I.matrix_cell{i,j} = zeros(A.matrix_cell{i,j});
                    end
                end
            end
        end
        
        function M = overop2mat(A)
            M = [A.matrix_cell{1,1}.matrix A.matrix_cell{1,2}.matrix A.matrix_cell{1,3}.matrix
                A.matrix_cell{2,1}.matrix A.matrix_cell{2,2}.matrix A.matrix_cell{2,3}.matrix
                A.matrix_cell{3,1}.matrix A.matrix_cell{3,2}.matrix A.matrix_cell{3,3}.matrix];
        end
    end
end