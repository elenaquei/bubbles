classdef overoperator3
    properties
        op_cell
    end
    methods
        function B = overoperator3(A11, A12, A13, A21, A22, A23, A31, A32, A33)
            % inputs: 9 overoperators
            A = cell(3,3);
            A{1,1} = A11;
            A{1,2} = A12;
            A{1,3} = A13;
            
            A{2,1} = A21;
            A{2,2} = A22;
            A{2,3} = A23;
            
            
            A{3,1} = A31;
            A{3,2} = A32;
            A{3,3} = A33;
            
            B.op_cell = A;
            
        end
        
        function C = plus(A,B)
            C=A;
            for i=1:3
                for j = 1:3
                    C.op_cell{i,j} = A.op_cell{i,j}+B.op_cell{i,j};
                end
            end
        end
        
        function C = minus(A,B)
            C=A;
            for i=1:3
                for j = 1:3
                    C.op_cell{i,j} = A.op_cell{i,j}-B.op_cell{i,j};
                end
            end
        end
        
        function n = norm(A)
            n = sum([ norm(A.op_cell{1,1})  norm(A.op_cell{1,2})  norm(A.op_cell{1,3})
                norm(A.op_cell{2,1})  norm(A.op_cell{2,2})  norm(A.op_cell{2,3})
                norm(A.op_cell{3,1})  norm(A.op_cell{3,2})  norm(A.op_cell{3,3})],2);
        end
        
        function n = cnorm(A)
            n = [ cnorm(A.op_cell{1,1})  cnorm(A.op_cell{1,2})  cnorm(A.op_cell{1,3})
                cnorm(A.op_cell{2,1})  cnorm(A.op_cell{2,2})  cnorm(A.op_cell{2,3})
                cnorm(A.op_cell{3,1})  cnorm(A.op_cell{3,2})  cnorm(A.op_cell{3,3})];
        end
        
        function B = abs(A)
            B=A;
            for i =1:3
                for j = 1:3
                    B.op_cell{i,j} = abs(A.op_cell{i,j});
                end
            end
        end
        
        function s = size_center(A)
            s = zeros(6);
            for i = 1:3
                for j = 1:3
                    s((i-1)*2+1:(i*2-1),(j-1)*2+1:(j*2-1)) = size_center(A.op_cell{i,j});
                end
            end
        end
        
        
        function C = min(A,B)
            if ~isa(A,'overoperator3') ||~isa(B,'overoperator3')
                error('Comparison impossible')
            end
            C=A;
            for i =1:3
                for j = 1:3
                    C.op_cell{i,j} = min(A.op_cell{i,j},B.op_cell{i,j});
                end
            end
        end
        
        function C = max(A,B)
            if ~isa(A,'overoperator3') ||~isa(B,'overoperator3')
                error('Comparison impossible')
            end
            C=A;
            for i =1:3
                for j = 1:3
                    C.op_cell{i,j} = max(A.op_cell{i,j},B.op_cell{i,j});
                end
            end
        end
        
        
        function As = intval(A0,A1)
            As=A0;
            for i =1:3
                for j = 1:3
                    As.op_cell{i,j} = intval(A0.op_cell{i,j},A1.op_cell{i,j});
                end
            end
        end
        
        function C = mtimes(A,B)
            if ~isa(A,'overoperator3')
                temp = A;
                A = B;
                B = temp;
            end
            if isa(B,'overoperator3')
                C=A;
                for i=1:3
                    for j = 1:3
                        C.op_cell{i,j} = 0*C.op_cell{i,j};
                        for k=1:3
                            C.op_cell{i,j} = C.op_cell{i,j}+A.op_cell{i,k}*B.op_cell{k,j};
                        end
                    end
                end
            elseif isa(B,'PDE_vector')
                A = reshape(A,B.node_space,B.node_time);
                [B1,B2,B3] = split(B);
                
                C1 = A.op_cell{1,1}*B1 + A.op_cell{1,2}*B2 +...
                    A.op_cell{1,3}*B3;
                
                C2 = A.op_cell{2,1}*B1 + A.op_cell{2,2}*B2+...
                    A.op_cell{2,3}*B3;
                
                C3 = A.op_cell{3,1}*B1 + A.op_cell{3,2}*B2+...
                    A.op_cell{3,3}*B3;
                
                C=merge(C1,C2,C3);
                
            elseif isa(B,'double')
                C=A;
                for i=1:3
                    for j = 1:3
                        C.op_cell{i,j} = B*C.op_cell{i,j};
                    end
                end
            else
                error('unrecognized input type')
            end
        end
        
        function A_big = expand(A,new_node_space, new_node_time)
            A_big=A;
            for i=1:3
                for j = 1:3
                    A_big.op_cell{i,j} = expand(A.op_cell{i,j},new_node_space, new_node_time);
                end
            end
        end
        
        function A_small = project(A,new_node_space, new_node_time)
            A_small=A;
            for i=1:3
                for j = 1:3
                    A_small.op_cell{i,j} = project(A.op_cell{i,j},new_node_space, new_node_time);
                end
            end
        end
        
        function A_tail = tail(A,center_node_space, center_node_time)
            A_tail=A;
            for i=1:3
                for j = 1:3
                    A_tail.op_cell{i,j} = tail(A.op_cell{i,j},center_node_space, center_node_time);
                end
            end
        end
        
        function B = reshape(A, new_node_space, new_node_time)
            B=A;
            for i=1:3
                for j = 1:3
                    B.op_cell{i,j} = reshape(A.op_cell{i,j},new_node_space, new_node_time);
                end
            end
        end
        
        function I = identity(A)
            I = A;
            for i = 1:3
                for j =1:3
                    I.op_cell{i,j} = identity(A.op_cell{i,j});
                end
            end
        end
        
        function M = overop2mat(A)
            matrix_cell=cell(3);
            for i = 1:3
                for j = 1:3
                    matrix_cell{i,j} = overop2mat(A.op_cell{i,j});
                end
            end
            M = [matrix_cell{1,1} matrix_cell{1,2} matrix_cell{1,3}
                matrix_cell{2,1} matrix_cell{2,2} matrix_cell{2,3}
                matrix_cell{3,1} matrix_cell{3,2} matrix_cell{3,3}];
        end
        
        function nnode_space = node_space(A)
            nnode_space = (size(A.op_cell{2,2},1)-1)/2;
        end
        
        
%         function A_big_tail = add_tail(A, A_tail,new_node_space, new_node_time)
%             old_nspace = node_space(A);
%             old_ntime = node_time(A);
%             A_tail_op = tail_mat2op3(A_tail,old_nspace, old_ntime, new_node_space, new_node_time);
%             A_big = reshape(A, new_node_space, new_node_time);
%             A_big_tail = A_tail_op + A_big;
%         end
        
    end
end