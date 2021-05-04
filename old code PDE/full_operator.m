classdef full_operator
    properties
        center % overoperator3
        tail_lin % tail_matrix
        tail_nonlin % tail_matrix
        coef_nonlin % cell(6,6) OR cell(3,3)
        dim % either 3 OR 6
    end
    methods
        function a = full_operator(center, tail_lin, tail_nonlin, coef_nonlin)
            if isa(center,'overoperator3')
                a.dim = 6;
            elseif isa(center,'overoperator')
                a.dim  = 2;
            else
                error('Unexpected input')
            end
            a.center = center;
            if nargin ==1
                tail_lin = zero_tail(a.dim);
            end
            a.tail_lin = tail_lin;
            if nargin>2
                if any(size(coef_nonlin)~=a.dim)
                    error('Are you sure? Revise class')
                end
                a.tail_nonlin = tail_nonlin;
                a.coef_nonlin = coef_nonlin;
            else
                a.tail_nonlin = 0*a.tail_lin;
                a.coef_nonlin = cell(6);
            end
        end
        
        function c = compose_fullop(A11, A12, A13, A21, A22, A23, A31, A32, A33)
            center_c = overoperator3(A11.center, A12.center, A13.center, ...
                A21.center, A22.center, A23.center,...
                A31.center, A32.center, A33.center);
            tail_lin_c = compose_tail(A11.tail_lin, A12.tail_lin, A13.tail_lin, ...
                A21.tail_lin, A22.tail_lin, A23.tail_lin,...
                A31.tail_lin, A32.tail_lin, A33.tail_lin);
            tail_nonlin_c = compose_tail(A11.tail_nonlin, A12.tail_nonlin, A13.tail_nonlin, ...
                A21.tail_nonlin, A22.tail_nonlin, A23.tail_nonlin,...
                A31.tail_nonlin, A32.tail_nonlin, A33.tail_nonlin);
            coef_nonlin_c = compose_cell(A11.coef_nonlin, A12.coef_nonlin, A13.coef_nonlin, ...
                A21.coef_nonlin, A22.coef_nonlin, A23.coef_nonlin,...
                A31.coef_nonlin, A32.coef_nonlin, A33.coef_nonlin);
            c = full_operator(center_c, tail_lin_c, tail_nonlin_c, coef_nonlin_c);
        end
        
        function c_tail = tail(A)
            c_tail = A;
            c_tail.center = tail(A.center);
        end
        
        function c = plus(a,b)
            if a.dim ~= b.dim
                error('incompatible sizes')
            end
            c = a;
            c.center = a.center+b.center;
            c.tail_lin = a.tail_lin + b.tail_lin;
            c.tail_nonlin = a.tail_nonlin;
            for i =1:a.dim
                for j = 1:a.dim
                    %if isa(a.tail_nonlin.mat{i,j},'polynomial_fraction')
                    a.tail_nonlin.mat{i,j} = simplify(a.tail_nonlin.mat{i,j});
                    %end
                    %if isa(a.tail_nonlin.mat{i,j},'polynomial_fraction')
                    b.tail_nonlin.mat{i,j} = simplify(b.tail_nonlin.mat{i,j});
                    %end
                    if nonzero(a.tail_nonlin.mat{i,j}) && ~nonzero(b.tail_nonlin.mat{i,j})
                        c.tail_nonlin.mat{i,j} = b.tail_nonlin.mat{i,j};
                        c.coef_nonlin{i,j} = b.coef_nonlin{i,j};
                    elseif ~nonzero(a.tail_nonlin.mat{i,j}) && ~nonzero(b.tail_nonlin.mat{i,j})
                        c.tail_nonlin.mat{i,j} = polynomial(0);
                        c.coef_nonlin{i,j} = [];
                    elseif ~nonzero(a.tail_nonlin.mat{i,j}) && nonzero(b.tail_nonlin.mat{i,j})
                        c.tail_nonlin.mat{i,j} = a.tail_nonlin.mat{i,j};
                        c.coef_nonlin{i,j} = a.coef_nonlin{i,j};
                    elseif a.tail_nonlin.mat{i,j} == b.tail_nonlin.mat{i,j}
                        c.tail_nonlin.mat{i,j} = a.tail_nonlin.mat{i,j};
                        c.coef_nonlin{i,j} = a.coef_nonlin{i,j}+b.coef_nonlin{i,j};
                    else
                        error('Incompatible difference')
                    end
                end
            end
        end
        
        
        function c = minus(a,b)
            if a.dim ~= b.dim
                error('incompatible sizes')
            end
            c = a;
            c.center = a.center-b.center;
            c.tail_lin = a.tail_lin - b.tail_lin;
            c.tail_nonlin = a.tail_nonlin;
            for i =1:a.dim
                for j = 1:a.dim
                    %if isa(a.tail_nonlin.mat{i,j},'polynomial_fraction')
                    a.tail_nonlin.mat{i,j} = simplify(a.tail_nonlin.mat{i,j});
                    %end
                    %if isa(a.tail_nonlin.mat{i,j},'polynomial_fraction')
                    b.tail_nonlin.mat{i,j} = simplify(b.tail_nonlin.mat{i,j});
                    %end
                    if nonzero(a.tail_nonlin.mat{i,j}) && ~nonzero(b.tail_nonlin.mat{i,j})
                        c.tail_nonlin.mat{i,j} = b.tail_nonlin.mat{i,j};
                        c.coef_nonlin{i,j} = -b.coef_nonlin{i,j};
                    elseif ~nonzero(a.tail_nonlin.mat{i,j}) && ~nonzero(b.tail_nonlin.mat{i,j})
                        c.tail_nonlin.mat{i,j} = polynomial(0);
                        c.coef_nonlin{i,j} = [];
                    elseif ~nonzero(a.tail_nonlin.mat{i,j}) && nonzero(b.tail_nonlin.mat{i,j})
                        c.tail_nonlin.mat{i,j} = a.tail_nonlin.mat{i,j};
                        c.coef_nonlin{i,j} = a.coef_nonlin{i,j};
                    elseif a.tail_nonlin.mat{i,j} == b.tail_nonlin.mat{i,j}
                        c.tail_nonlin.mat{i,j} = a.tail_nonlin.mat{i,j};
                        c.coef_nonlin{i,j} = a.coef_nonlin{i,j}-b.coef_nonlin{i,j};
                    else
                        error('Incompatible difference')
                    end
                end
            end
        end
        
        
        function c = abs(a)
            c = a;
            c.center = abs(a.center);
            c.tail_lin = abs(a.tail_lin);
            c.tail_nonlin = abs(a.tail_nonlin);
            for i =1:a.dim
                for j = 1:a.dim
                    a.tail_nonlin.mat{i,j} = simplify(a.tail_nonlin.mat{i,j});
                    
                    c.tail_nonlin.mat{i,j} = abs(a.tail_nonlin.mat{i,j});
                    c.coef_nonlin{i,j} = abs(a.coef_nonlin{i,j});
                end
            end
        end
        
        function c = mtimes(a,b)
            if ~isa(b,'full_operator')
                temp = b;
                b = a;
                a = temp;
            end
            if isa(a,'double')||isa(a,'intval')
                c = b;
                c.center = a*b.center;
                c.tail_lin = a*b.tail_lin;
                for i = 1:b.dim
                    for j = 1:b.dim
                        c.coef_nonlin{i,j} = a*b.coef_nonlin{i,j};
                    end
                end
                return
            end
            if ~isa(a,'full_operator')||~isa(b,'full_operator')
                error('Unexpected input typoe')
            end
            % ASSUMPTION : one of the two has zero non-linear tail
            if a.dim ~= b.dim
                error('incompatible sizes')
            end
            
            
            if ~iszerocoef(b)
                if ~iszerocoef(a)
                    error('Expected at most one non-zero non-linear tail')
                end
                temp = b;
                b = a;
                a = temp;
            end % a has the tail
            
            c = a;
            c.center = a.center*b.center;
            c.tail_lin = a.tail_lin * b.tail_lin;
            c.tail_nonlin = a.tail_nonlin * b.tail_lin;
            c.coef_nonlin = a.coef_nonlin;
        end
        
        function As = intval(A0,A1)
            if A0.dim ~= A1.dim
                error('Incompatible input')
            end
            As = A0;
            As.center=intval(A0.center,A1.center); % overoperator3
            As.tail_lin = intval(A0.tail_lin, A1.tail_lin); % tail_matrix
            for i =1:A0.dim
                for j = 1:A0.dim
                    if isempty(A0.coef_nonlin{i,j}) && isempty(A1.coef_nonlin{i,j})
                        As.coef_nonlin{i,j} = [];
                    else
                        As.coef_nonlin{i,j} = infsup_help(A0.coef_nonlin{i,j},A1.coef_nonlin{i,j});
                    end
                end
            end
        end
        
        
        
        function bool = iszerocoef(a)
            bool = 1;
            for i =1:a.dim
                for j = 1:a.dim
                    if ~isempty(a.coef_nonlin{i,j})
                        bool = 0;
                        return
                    end
                end
            end
        end
        
        function n = cnorm(a)
            n1 = cnorm(a.center);
            nnode_space = node_space(a.center);
            n2 = cnorm(a.tail_lin,nnode_space) + cnorm(a.tail_nonlin,nnode_space).*cnorm_coef(a.coef_nonlin,a.dim);
            n2_reshape = 0*n1;
            n2_reshape((end+1-size(n2,1)):end,(end+1-size(n2,1)):end) = n2; 
            
            n = max(n1,n2_reshape);
        end
        
        function n = norm(a)
            n = sum(cnorm(a),2);
        end
    end
end

function n = cnorm_coef(cell_dim, dim)
n = zeros(dim);
for i = 1:dim
    for j = 1:dim
        if isa(cell_dim{i,j},'intval')
            cell_dim{i,j} = sup(abs(cell_dim{i,j}));
        end
        if max(size(cell_dim{i,j}))==1
            n(i,j) = abs(cell_dim{i,j});
        else
            n(i,j) = vector_norm(cell_dim{i,j});
        end
    end
end
end

function  tail_lin = zero_tail(dim)
mat = cell(dim);
pol_zero = polynomial(0);
for i = 1:dim
    for j = 1:dim
        mat{i,j} = pol_zero;
    end
end
tail_lin = tail_matrix(mat);
end