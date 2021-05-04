classdef polynomial
    properties
        coefs
    end
    methods
        % CONSTRUCTOR
        function x=polynomial(coefs)
            %
            % polynomial constructor of the class polynomial
            %
            % [c_n, c_n-1, ... c_1, c0]
            x.coefs = coefs;
        end
        % end CONSTRUCTOR
        
        % MTIMES
        function c = mtimes(a,b)
            if ~isa(b,'polynomial')
                temp = a;
                a = b;
                b = temp;
            end
            if isa(a,'float') || isa(a,'intval')
                c = b;
                c.coefs = a * b.coefs;
            elseif isa(a,'polynomial')
                c = b;
                    if length(a.coefs) ==1 || length(b.coefs) ==1
                        c.coefs = a.coefs*b.coefs;
                        return
                    end
                if isa(a.coefs,'intval') || isa(b.coefs,'intval')
                    c.coefs = convINTVAL(a.coefs, b.coefs,'full');
                else
                    c.coefs = conv(a.coefs, b.coefs,'full');
                end
            elseif isa(a,'polynomial_fraction')
                b_frac = polynomial_fraction(b,polynomial(1));
                c = a*b_frac;
            else
                error('Incompatible product')
            end
        end
        % end MTIMES
        
        % PLUS
        function c = plus(a,b)
            if ~isa(b,'polynomial')
                temp = a;
                a = b;
                b = temp;
            end
            if isa(a,'float') || isa(a,'intval')
                c = b;
                c.coefs(end) = a + b.coefs(end);
            elseif isa(a,'polynomial')
                c = b;
                while length(a.coefs)<length(b.coefs)
                    a.coefs = [0,a.coefs];
                end
                while length(b.coefs)<length(a.coefs)
                    b.coefs = [0,b.coefs];
                end
                c.coefs = a.coefs + b.coefs;
            elseif isa(a,'polynomial_fraction')
                b_frac = polynomial_fraction(b,polynomial(1));
                c = a+b_frac;
            else
                error('Incompatible product')
            end
        end
        % end PLUS
        
        % MINUS
        function c = minus(a,b)
            c = a + (-1*b);
        end
        % end MINUS
        
        % INTVAL
        function cs = intval(a, b)
                while length(a.coefs)<length(b.coefs)
                    a.coefs = [0,a.coefs];
                end
                while length(b.coefs)<length(a.coefs)
                    b.coefs = [0,b.coefs];
                end
                cs = a;
                cs.coefs = infsup( min(a.coefs,b.coefs), max(a.coefs,b.coefs));
        end
        % end INTVAL
        
        
        % ABS
        function a_abs = abs(a)
            a_abs = a;
            a_abs.coefs = abs( a.coefs );
        end
        % end ABS
        
        % DERIVATIVE
        function b = derivative(a)
            b = a;
            if length(b.coefs) ==1
                b.coefs =0;
                return
            end
            coef = a.coefs;
            coef = ((length(coef)-1):-1:0).*coef;
            coef = coef(1:end-1);
            b.coefs = coef;
        end
        % end DERIVATIVE
        
        % UPPER_BOUND_INTEGER
        function K = upper_bound_integer(a)
            % function K = upper_bound_integer(a)
            %
            % returns an integer bound on the absolute value of the biggest zero
            if length(a.coefs)==1 
                if a.coefs == 0
                    K = 0;
                    return
                else
                    K =Inf;
                    return
                end
            end
            if isa(a.coefs,'intval')
                a.coefs = sup(a.coefs);
            end
            bound1 = max( [1, sum(abs(a.coefs))/abs(a.coefs(1))]);
            bound2 = 1 + max(abs(a.coefs(2:end)))/abs(a.coefs(1));
            
            K = min(bound1, bound2);
            K = ceil(K);
        end
        % end UPPER_BOUND_INTEGER
        
        % EVAL_POL
        function y = eval_pol(p, x)
            X = repmat(x(:),1,length(p.coefs));
            P = repmat(p.coefs(:).',length(x),1);
            Pow = repmat((length(p.coefs)-1):-1:0,length(x),1);
            y = sum( P.*(X.^Pow),2);
            if size(y,1)~=size(x,1)
                y = y.';
            end
        end
        % end EVAL_POL
        
        % MAX
        function m = max(a,b)
            m = 0*a;
            m.coefs= max(a.coefs, b.coefs);
        end
        % end MAX
        
        % EQ
        function bool = eq(a,b)
            try 
                bool = (a.coefs == b.coefs);
            catch
                bool = 0;
            end 
        end
        % end EQ
        
        % SIMPLIFY
        function a_simple = simplify(a)
            a_simple = a;
            while length(a_simple.coefs)>1 && a_simple.coefs(1) ==0
                a_simple.coefs(1) = [];
            end
        end
        % end SIMPLIFY
        
        % NONZERO
        function bool = nonzero(a)
            if isempty(a.coefs)
                bool = 0;
            elseif all(a.coefs==0)
                bool = 0;
            else
                bool = 1;
            end
        end
        % end NONZERO
    end
end