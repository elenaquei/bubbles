
classdef tail_matrix
    properties
        mat % cell with polynomial or polynomial_fraction
    end
    methods
        function A = tail_matrix(cell_in)
            % function A = tail_matrix(cell_in)
            % INPUT
            % cell_in      square cell with polynomial or polynomial_fraction
            if ~isa(cell_in,'cell')
                error('uncexpected input')
            end
            for i =1:size(cell_in,1)
                for j = 1:size(cell_in,2)
                    if ~isa(cell_in{i,j}, 'polynomial') && ~isa(cell_in{i,j}, 'polynomial_fraction')
                        error('Unrecognized input type')
                    end
                end
            end
            A.mat = cell_in;
        end
        
        
        % PLUS
        function B = plus(A, M)
            if ~isa(M,'tail_matrix')
                temp = A;
                A = M;
                M = temp;
            end
            if isa(A,'operator')
                B = full_operator(A,M);
                return
            elseif ~isa(A,'tail_matrix')
                error('Unrecognised input')
            end
            if size(A.mat,1)~= size(M.mat,1) ||size(A.mat,2)~= size(M.mat,2)
                error('Incompatible sizes')
            end
            
            B = 0*M;
            
            for i =1:size(A.mat,1)
                for j = 1:size(M.mat,2)
                    B.mat{i,j} = A.mat{i,j} + M.mat{i,j};
                end
            end
            
        end
        % end PLUS
        
        
        % MINUS
        function B = minus(A, M)
            if ~isa(M,'tail_matrix')
                temp = A;
                A = M;
                M = temp;
            end
            if ~isa(A,'tail_matrix')
                error('Unrecognised input')
            end
            if size(A.mat,1)~= size(M.mat,1) ||size(A.mat,2)~= size(M.mat,2)
                error('Incompatible sizes')
            end
            
            B = 0*M;
            
            for i =1:size(A.mat,1)
                for j = 1:size(M.mat,2)
                    B.mat{i,j} = A.mat{i,j} - M.mat{i,j};
                end
            end
            
        end
        % end MINUS
        
        % ABS
        function C = abs(A)
            C = A;
            
            for i =1:size(A.mat,1)
                for j = 1:size(A.mat,2)
                    C.mat{i,j} = abs(A.mat{i,j});
                end
            end
        end
        % end ABS
        
        % INTVAL
        function As = intval(A0, A1)
            As = A0;
            
            for i =1:size(A0.mat,1)
                for j = 1:size(A1.mat,2)
                    As.mat{i,j} = intval(A0.mat{i,j}, A1.mat{i,j});
                end
            end
        end
        % end INTVAL
        
        % MTIMES
        function B = mtimes(A, M)
            if ~isa(M,'tail_matrix')
                temp = A;
                A = M;
                M = temp;
            end
            if isa(A,'float') || isa(A,'intval')
                if length(A)~=1
                    error('Not implemented')
                end
                B = M;
                
                for i =1:size(M.mat,1)
                    for j = 1:size(M.mat,2)
                        B.mat{i,j} = A * M.mat{i,j};
                    end
                end
            elseif isa(A,'tail_matrix')
                B = 0*M;
                
                for i =1:size(A.mat,1)
                    for j = 1:size(M.mat,2)
                        for k = 1:size(M.mat,1)
                            B.mat{i,j} = B.mat{i,j} + A.mat{i,k} * M.mat{k,j};
                        end
                    end
                end
            else
                error('Multiplication not possible')
            end
        end
        % end MTIMES
        
        % MAX
        function m = max(a,b)
            m = 0*a;
            
            for i =1:size(A.mat,1)
                for j = 1:size(M.mat,2)
                    m.mat{i,j} = max(a.mat{i,j},b.mat{i,j});
                    
                end
            end
        end
        % end MAX
        
        % EQ
        function m = eq(a,b)
            m = size(size(A.mat));
            for i =1:size(A.mat,1)
                for j = 1:size(M.mat,2)
                    m.mat{i,j} = eq(a.mat{i,j},b.mat{i,j});
                    
                end
            end
        end
        % end EQ
        
        % COMPOSE_TAIL
        function c = compose_tail(A11, A12, A13, A21, A22, A23, A31, A32, A33)
            c_cell = compose_cell(A11.mat, A12.mat, A13.mat,...
                A21.mat, A22.mat, A23.mat, A31.mat, A32.mat, A33.mat);
            c=tail_matrix(c_cell);
        end
        % end COMPOSE_TAIL
        
        % cnorm
        function b = cnorm(M,K)
            % function b = cbound(M,K)
            %
            % component wise bound of a tail matrix
            b = zeros(size(M.mat));
            for i = 1:size(M.mat,1)
                for j = 1:size(M.mat,2)
                    if nonzero(M.mat{i,j})
                        b(i,j) = upper_bound_fraction(M.mat{i,j},K);
                    else
                        b(i,j) = 0;
                    end
                end
            end
        end
        % end cnorm
        
    end
end