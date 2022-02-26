classdef polynomial_fraction
    properties
        nominator
        denominator
    end
    methods
        % CONSTRUCTOR
        function x=polynomial_fraction(p,q)
            %
            % polynomial constructor of the class polynomial
            if (isa(q,'float')  || isa(q,'intval') )&& length(q) ==1
                q = polynomial(q);
            end
            if ~isa(q,'polynomial') || ~isa(p,'polynomial')
                error('Inputs muust be two polynomials')
            end
            x.nominator = p;
            x.denominator = q;
            x = simplify(x);
        end
        % end CONSTRUCTOR
        
        % SIMPLIFY
        function pq_simple = simplify(pq) 
            
            pq.nominator = simplify(pq.nominator);
            pq.denominator = simplify(pq.denominator);
            
            p = pq.nominator.coefs;
            q = pq.denominator.coefs;
            if sum(abs(q))==0
                error('Dividing by zero')
            end
            if (isa(p,'double') && max(abs(p))==0) || (isa(p,'intval') && max(sup(abs(p)))==0)
                p=0;
                q=1;
            end
            
            while p(end) ==0 && q(end) == 0 && min(length(p), length(q))>1
                p(end) = [];
                q(end) = [];
            end
            if length(q) ==1 && q~=1
                p = p/q;
                q = 1;
            end
            pq_simple = pq;
            pq_simple.nominator.coefs = p;
            pq_simple.denominator.coefs = q;
        end
        % end SIMPLIFY
        
        % DERIVATIVE
        function pq_prime = derivative(pq)
            pq = simplify(pq);
            
            p = pq.nominator;
            q = pq.denominator;
            
            nominator_der = derivative(p)*q - derivative(q)*p;
            denominator_der = q * q;
            pq_prime = polynomial_fraction(nominator_der,denominator_der);
            
        end
        % end DERIVATIVE
        
        % MTIMES
        function pq3 = mtimes(pq1, pq2)
            if ~isa(pq2, 'polynomial_fraction')
                temp = pq1;
                pq1 = pq2;
                pq2 = temp;
            end
            if (isa(pq1, 'float')|| isa(pq1,'intval')) && pq1 == 0
                pq3 = pq2;
                pq3.nominator.coefs = 0;
                pq3.denominator.coefs = 1;
            elseif isa(pq1, 'float') || isa(pq1,'intval')|| isa(pq1, 'polynomial')
                pq3 = pq2;
                pq3.nominator = pq1*pq2.nominator;
            elseif isa(pq1, 'polynomial_fraction')
                pq3 = pq1;
                pq3.nominator = pq1.nominator * pq2.denominator;
            else
                error('Inputs not compatible')
            end 
        end
        % end MTIMES
        
        % PLUS
        function pq3 = plus(pq1, pq2)
            if ~isa(pq2, 'polynomial_fraction')
                temp = pq1;
                pq1 = pq2;
                pq2 = temp;
            end
            if (isa(pq1, 'float')|| isa(pq1,'intval')) || isa(pq1, 'polynomial')
                pq3 = pq2;
                pq3.nominator = pq1 * pq2.denominator + pq2.nominator;
            elseif isa(pq1, 'polynomial_fraction')
                pq3 = pq1;
                pq3.denominator = pq1.denominator * pq2.denominator...
                     + pq2.nominator * pq1.denominator;
            else
                error('Inputs not compatible')
            end 
        end
        % end PLUS
        
        % MINUS
        function pq3 = minus(pq1, pq2)
            pq3 = pq1 + (-1*pq2);
        end
        % end MINUS
        
        function pq = intval(pq1, pq2)
            pq1 = simplify(pq1);
            pq2 = simplify(pq2);
            if all(pq1.denominator.coefs == pq2.denominator.coefs)
                pq = pq1;
                pq.denominator = intval(pq1.nominator, pq2.nominator);
                return
            elseif pq1.nominator.coefs == 1 &&  pq2.nominator.coefs==1
                pq_den = pq1.denominator * pq2.denominator;
                pq_nom = intval(pq1.nominator, pq2.nominator);
                pq = polynomial_fraction(pq_nom, pq_den);
            else
                error('Not coded - was not needed')
            end
        end
        
        
        % EVAL_POL
        function y = eval_pol(pq,x)
            % function y = eval_pol(pq,x)
            %
            % polynomial_fraction evaluation at a given x 
            y_nom = eval_pol(pq.nominator, x);
            y_denom = eval_pol(pq.denominator, x);
            y = y_nom./y_denom;
        end
        % end EVAL_POL
        
        % EQ 
        function bool = eq(a,b)
            a = simplify(a);
            b = simplify(b);
            try 
                bool = (a.nominator == b.nominator)*(a.denominator == b.denominator);
            catch
                bool = 0;
            end
        end
        % end EQ
        
        % NONZERO
        function bool = nonzero(a)
            bool = 1;
            if isempty(a.nominator.coefs) && isempty(a.denominator.coefs)
                bool = 0;
            elseif ( isa(a.nominator.coefs,'double') && max(abs(a.nominator.coefs))==0) ||...
                    ( isa(a.nominator.coefs,'intval') && max(sup(abs(a.nominator.coefs)))==0)
                bool = 0;
            end
        end
        % end NONZERO
        
        % upper_bound_fraction
        function B = upper_bound_fraction(pq, K)
            % function B = upper_bound_fraction(pq, K)
            % 
            % INPUT
            % pq    polynomial_fraction
            % K     integer, bound on the non considered integers
            % OUTPUT
            % B > max( abs(pq(k)), k>=K)
            pq = simplify(pq);
            pq_prime = derivative(pq);
            pq_der_nom = pq_prime.nominator;
            k_star = upper_bound_integer(pq_der_nom);
            if k_star < K
                k = K;
            else
                try 
                    k = K:k_star;
                catch
                    disp(9)
                end
            end
            pq_k = eval_pol(pq,k);
            if isa(pq_k,'intval')
                pq_k = sup(pq_k);
            end
            B = max( abs( pq_k));
        end
        % end upper_bound_fraction
    end
end