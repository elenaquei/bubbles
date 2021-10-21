classdef polynomial_coefs
    properties
        size_scalar
        size_vector
        deg_scalar
        deg_vector
        n_equations
        n_terms % vector(n_equations)  -- number of terms per equation
        value % cell{n_equations}(n_terms) -- coefficient of the term
        dot % cell{n_equations}{n_terms}(variables,:) -- number of derivatives applied to the given variable 
        power_vector % cell{n_equations}{n_terms}(variables,:)
        delay % cell{n_equations}{n_terms}(variables)
        power_scalar % cell{n_equations}(scalar_variables,n_terms)
        monomial_per_term % cell{n_equations}(n_terms) -- number of monomials per each term considered
        %
        % 
        % EXAMPLE
        %
        % x_dot_dot(t-tau) [x(t-2)]^2x_dot^3 y(t-1)
        % 
        %  delay{n_equation, 1}{term} = [tau,2,0]
        %  power{n_equation, 1}{term} = [1,2,3]
        %  dot  {n_equation, 1}{term} = [2,0,1]
        %  
        %  delay{n_equation, 2}{term} = 1
        %  power{n_equation, 2}{term} = 1
        %  dot  {n_equation, 2}{term} = 0
    end
    methods
        % contructor POLYNOMIAL_COEFS
        function alpha = polynomial_coefs(size_Scal, size_Vec, number_equation, ...
                number_term, coefficients,power_Scal,power_Vec, derivatives, delays)
            % function alpha = polynomial_coefs(size_Scal, size_Vec, number_equation, ...
            %     number_term, coefficients, derivatives, power_Vec, delay, power_Scal)
            
            alpha.size_scalar=0;
            alpha.size_vector=0;
            alpha.deg_scalar=0;
            alpha.deg_vector=0;
            alpha.n_equations=0;
            alpha.n_terms=0;
            alpha.value=[]; % cell{n_equations}(n_terms) -- coefficient of the term
            alpha.dot=[]; % cell{n_equations}{n_terms}(variables,:) -- number of derivatives applied to the given variable
            alpha.power_vector=[]; % cell{n_equations}{n_terms}(variables,:)
            alpha.delay=[]; % cell{n_equations}{n_terms}(variables,:)
            alpha.power_scalar=[]; % cell{n_equations}(scalar_variables,n_terms)
            alpha.monomial_per_term=[];
            
            if nargin>0
                if ~(numel(size_Scal)==1 || mod(size_Scal,1)==0)
                    error('number of scalar variables must be an integer')
                end
                alpha.size_scalar = size_Scal;
            else
                return
            end
            if nargin>1
                if ~(numel(size_Vec)==1 || mod(size_Vec,1)==0)
                    error('number of vector variables must be an integer')
                end
                alpha.size_vector = size_Vec;
            end
            if nargin>2
                if ~(numel(number_equation)==1 || mod(number_equation,1)==0)
                    error('number of equations must be an integer')
                end
                alpha.n_equations = number_equation;
            end
            
            
            if nargin>3
                if length(number_term) ~= alpha.n_equations
                    error('The vector of the number of terms per equation need to be as long as the number of equations')
                end
                if any(mod(number_term,1)~=0)
                    error('Number of terms must be integers')
                end
                alpha.n_terms = number_term;
            end
            
            
            if nargin>4
                try
                    for i = 1: alpha.n_equations
                        if alpha.n_terms(i)>0 %july 5th
                            coefficients{i}(alpha.n_terms(i));
                        end
                    end
                catch
                    error('The coefficients must have the form cell{n_equations}(n_terms)')
                end
                % value % cell{n_equations}(n_terms)
                alpha.value = coefficients;
            end
            
            
            if nargin >5
                % power_Scal
                % power_scalar % cell{n_equations}(scal_variables,n_terms)
                try
                    for i = 1: alpha.n_equations
                        power_Scal{i}(1:alpha.size_scalar,1:alpha.n_terms(i));
                    end
                catch
                    error('The powers of scalars must have the form cell{n_equations}(scal_variables,n_terms)')
                end
                alpha.power_scalar = power_Scal;
            end
            
            if nargin >6
                % power_vector % cell{n_equations}{n_terms}(variables,:)
                try
                    for i = 1: alpha.n_equations
                        for j =1:alpha.n_terms(i)
                            power_Vec{i}{j}(1:alpha.size_vector,:);
                        end
                    end
                catch
                    error('The powers of vector must have the form cell{n_equations}{n_terms}(variables,:)')
                end
                
                alpha.power_vector = power_Vec;
            end
            
            
            alpha.monomial_per_term = cell(alpha.n_equations,1);
            for i = 1:alpha.n_equations
                alpha.monomial_per_term{i} = zeros(alpha.n_terms(i));
                for j = 1:alpha.n_terms(i)
                    alpha.monomial_per_term{i}(j) = size(alpha.power_vector{i}{j},2);
                end
            end
            
            if nargin >7
                % dot % cell{n_equations}{n_terms}(variables,:)
                try
                    for i = 1: alpha.n_equations
                        for j = 1:alpha.n_terms(i)
                            if any(size(derivatives{i}{j})~=size(alpha.power_vector{i}{j}))
                                error('Dimensions of derivatives and powers of vectors must be equal');
                            end
                        end
                    end
                catch
                    error('The derivatives must have the form cell{n_equations}{n_terms}(variables,:)')
                end
                % derivative % cell{n_equations}{n_terms}(variables,:)
                alpha.dot = derivatives;
            else
                
                alpha.dot = get_default_derivative(alpha);
            end
            
            if nargin >8
                % delay % cell{n_equations}{n_terms}(variables,:)
                try
                    for i = 1: alpha.n_equations
                        for j = 1:alpha.n_terms(i)
                            if any(size(delays{i}{j})~=size(alpha.power_vector{i}{j}))
                                error('Dimensions of delays and powers of vectors must be equal');
                            end
                        end
                    end
                catch
                    error('The delay must have the form cell{n_equations}{n_terms}(variables,:)')
                end
                % delay % cell{n_equations}{n_terms}(variables,:)
                alpha.delay = delays;
            else
                delays = cell(1,alpha.n_equations);
                for i = 1: alpha.n_equations
                    for j = 1: alpha.n_terms(i)
                        delays{i}{j}=0*(alpha.power_vector{i}{j});
                    end
                end
                alpha.delay = delays;
            end
            alpha = set_degrees(alpha);
        end
        % POLYNOMIAL_COEFS
        
        % COMPATIBLE
        function bool = compatible(alpha,beta)
            % function bool = compatible(alpha,beta)
            % this function checks if the two polynomial_coefs alpha and beta
            % have the same number of variables and define the same number of
            % equations
            bool = 0;
            try
                if(alpha.size_scalar~=beta.size_scalar)
                    return
                end
                if alpha.size_vector~= beta.size_vector
                    return
                end
                if alpha.n_equations~= beta.n_equations
                    return
                end
            catch
                return
            end
            bool = 1;
        end
        % end COMPATIBLE
        
        % COMPATIBLE_VEC
        function bool = compatible_vec(alpha,x_Xi_vec)
            % function bool = compatible_vec(alpha,x_Xi_vec)
            % this function checks if the two polynomial_coefs alpha and beta
            % have the same number of variables and define the same number of
            % equations
            bool = 0;
            try
                if(alpha.size_scalar~=x_Xi_vec.size_scalar)
                    return
                end
                if alpha.size_vector~= x_Xi_vec.size_vector
                    return
                end
            catch
                return
            end
            bool = 1;
        end
        % end COMPATIBLE_VEC
        
        % RESCALE
        function beta = rescale(alpha, index_lambda, rescaling)
            if index_lambda > alpha.size_scalar
                error('Not coded yet');
            end
            beta = alpha;
            for j = 1:alpha.n_equations
                T = find(alpha.power_scalar{j}(index_lambda,:));
                for t = T
                    beta.value{j}(t) = alpha.value{j}(t)/...
                        (rescaling^alpha.power_scalar{j}(index_lambda,t));
                end
            end
        end
        % end RESCALE
        
        
        % SET_DEGREES
        function beta = set_degrees(alpha)
            % function beta = set_degrees(alpha)
            % this function set the degrees of alpha
            beta = alpha;
            deg_scal = 0;
            deg_vec = 0;
            for i = 1:alpha.n_equations
                for j = 1:alpha.n_terms(i)
                    if sum(alpha.power_scalar{i}(:,j))>deg_scal
                        deg_scal = sum(alpha.power_scalar{i}(:,j));
                    end
                    deg_delay = sum(alpha.delay{i}{j}~=0);
                    deg_term_ij = sum(alpha.power_vector{i}{j}) + deg_delay;
                    if deg_term_ij>deg_vec
                        deg_vec = deg_term_ij;
                    end
                end
            end
            beta.deg_scalar = deg_scal;
            beta.deg_vector= deg_vec;
        end
        % end SET_DEGREES
        
        % get_default_monomial
        function monomial = get_default_monomial(alpha)
            % function monomial = get_default_monomial(alpha)
            %
            % set monomials all to 1
            monomial = cell(alpha.n_equations,1);
            for i = 1:alpha.n_equations
                monomial{i} = zeros(alpha.n_terms(i),1);
                for j = 1:alpha.n_terms(i)
                    monomial{i}(j) = 1;
                end
            end
            
        end
        % end get_default_monomial
        
        % get_default_delay
        function delay = get_default_delay(alpha)
            % function delay = get_default_delay(alpha)
            % get the default delay at 0 for everyone
            delay = cell(alpha.n_equations,1);
            for i = 1:alpha.n_equations
                delay{i} = cell(alpha.n_terms(i),1);
                for j = 1:alpha.n_terms(i)
                    delay{i}{j} = zeros(alpha.size_vector, alpha.monomial_per_term{i}{j});
                end
            end
        end
        % end get_default_delay
        
        % get_default_derivative
        function derivatives = get_default_derivative(alpha)
            % function derivatives = get_default_derivative(alpha)
            %
            % if alpha in the shape of a standard ODE, give a derivative to
            % all first terms of each equation, otherwise skip
            
            derivatives = cell(alpha.n_equations,1);
            for i = 1: alpha.n_equations
                derivatives{i}=cell(alpha.n_terms(i),1);
                for j = 1: alpha.n_terms(i)
                    derivatives{i}{j}=0*(alpha.power_vector{i}{j});
                end
                if size(derivatives{i}{1},2)==1 % if the first term of the equation is a monomial, derive it
                    derivatives{i}{1}(i,1)=1;
                else
                    warning('Not an ode in standard for, the first element of each equation should be a monomial');
                end
            end
        end
        % end get_default_derivative
        
        
        % set_zero_derivative
        function beta = set_zero_derivative(alpha)
            % function derivatives = set_zero_derivative(alpha)
            %
            % being alpha a polynomial coefficient, set all the derivatives
            % to zero
            beta = alpha;
            derivatives = cell(alpha.n_equations,1);
            for i = 1: alpha.n_equations
                derivatives{i}=cell(alpha.n_terms(i),1);
                terms_to_delete =[];
                for j = 1: alpha.n_terms(i)
                    
                    derivatives{i}{j}=0*(alpha.power_vector{i}{j});
                    if any(alpha.dot{i}{j}~=0) && all(alpha.power_vector{i}{j}==0)
                        terms_to_delete = [terms_to_delete, j];
                    end
                end
                if ~isempty(terms_to_delete)
                    beta.n_terms(i) = alpha.n_terms(i) - length(terms_to_delete);
                    beta.power_scalar{i}(:,terms_to_delete)=[];
                    beta.value{i}(terms_to_delete) =[];
                    
                    derivatives{i}{terms_to_delete} = [];
                    beta.power_vector{i}{terms_to_delete} = [];
                    beta.delay{i}{terms_to_delete} =[];
                    
                    
                    derivatives{i} = derivatives{i}(~cellfun('isempty',derivatives{i}));
                    beta.power_vector{i} = beta.power_vector{i}(~cellfun('isempty', beta.power_vector{i}));
                    beta.delay{i} = beta.delay{i}(~cellfun('isempty', beta.delay{i}));
                end
            end
            beta.dot = derivatives;
        end
        % end set_zero_derivative
        
        % APPLY
        function y = apply(a,xi_vec,bool)
            % function y = apply(a,xi_vec,bool)
            %
            % INPUTS
            % a         polynomial_coefs
            % xi_vec    Xi_vector
            % bool      boolean DEFAULT 0   -----   USELESS 
            %           if bool =0, output of the same size as input,
            %           otherwise output of size dim(x_xi)*deg(a)
            %
            % OUTPUT
            % y         complex vector size = a.n_equations, dim(x_xi)*deg(a),
            %           ifft of a applied to xi_vec
            global use_intlab
            if ~compatible_vec(a,xi_vec)
                error('Inputs not compatible')
            end
            
            if nargin<3
                bool = 1;
            end
            
            % if has_delay(a)
            %     y = apply_with_delay(a,xi_vec,bool);
            %     return
            % end
            
            xi_vec= set_ifft(xi_vec,a.deg_vector);
            %if nargin<3
            %    bool = 0;
            %end
            size_1 =a.n_equations;
            size_2 = size(xi_vec.ifft_vector,2);%xi_vec.nodes*2*a.deg_vector+1;
            K_nodes = size_2/2;
            K = -K_nodes:(K_nodes-1); %-(xi_vec.nodes*a.deg_vector):(xi_vec.nodes*a.deg_vector);
            help_xi_vec = reshape_Xi(xi_vec, K_nodes);%%xi_vec.nodes*a.deg_vector);
            
            y = zeros(size_1,size_2);
            if use_intlab
                y = intval(y);
            end
            for i =1 :size_1
                if any(any(a.power_scalar{i}<0))
                    error('No negative power is allowed')
                end
                for j = 1:a.n_terms(i)
                    if a.value{i}(j)==0
                        continue
                    end
                    term = a.value{i}(j)*prod(xi_vec.scalar.^(a.power_scalar{i}(:,j).'));
                    for k = 1:a.monomial_per_term{i}(j)
                        
                        if any(a.dot{i}{j}(:,k)>0)
                            if a.monomial_per_term{i}(j)==1 && sum(a.dot{i}{j}(:,k))==1
                                index = find(a.dot{i}{j}(:,k));
                                term = term * (verifyfft_in((1i*K).^a.dot{i}{j}(index,k).*help_xi_vec.vector(index,1:(end-1)),-1).');
                                continue
                            else
                                warning('General DERIVATIVES not implemented here')
                            end
                        end
                        if any(a.delay{i}{j}(:,k)~=0)
                            [~,ifft_prod] = power(xi_vec,a.power_vector{i}{j}(:,k));
                            delay_x = xi_vec.delay(a.delay{i}{j}(:,k));
                            delay_x= set_ifft(delay_x,a.deg_vector);
                            [~,ifft_prod_delay] = power(delay_x,ones(1:xi_vec.size_vector));
                            term = term .* ifft_prod .* ifft_prod_delay;
                        else
                            [~,ifft_prod] = power(xi_vec,a.power_vector{i}{j}(:,k));
                            term = term .* ifft_prod;
                        end
                        
                    end
                    y(i,:) = y(i,:) +term; % remember: it returns the ifft!
                end
                % y(i,:) = verifyfft_in(y(i,:),1);
            end
            
            if bool==0
                index = size(y,2)/2 + 1 + (-xi_vec.nodes:xi_vec.nodes);
                y = y(:,index);
            end
        end
        % end APPLY
        
        
        % APPLY_SUM
        function y = apply_sum(a,xi_vec)
            % function y = apply_sum(a,xi_vec)
            %
            % INPUTS
            % a         polynomial_coefs
            % xi_vec    Xi_vector
            %
            % OUTPUT
            % y         complex vector, a applied to sum(xi_vec), that is
            %           x(0)
            if ~compatible_vec(a,xi_vec)
                error('Inputs not compatible')
            end
            size_1 =a.n_equations;
            xi_vec_sum = sum(xi_vec.vector,2);
            y = zeros(size_1,1);
            if size_1>0
            if isintval(a.value{1}(1)) || isintval(xi_vec.scalar)
                y = intval(y);
            end
            end
            for i =1 :size_1
                for j = 1:a.n_terms(i)
                    term = a.value{i}(j)*prod(xi_vec.scalar.^(a.power_scalar{i}(:,j).'));
                    for k = 1:a.monomial_per_term{i}(j)
                        if any(a.delay{i}{j}(:,k)>0)
                            warning('DELAYS not implemented here')
                        end
                        if any(a.dot{i}{j}(:,k)>0)
                                warning('DERIVATIVES not implemented here')
                            
                        end
                        if any(a.power_vector{i}{j}(:,k)~=0)
                            term = term .* prod(xi_vec_sum.^a.power_vector{i}{j}(:,k));
                        end
                    end
                    y(i) = y(i) +term;
                end
            end
        end
        % end APPLY_SUM
        
        
        % COMPUTE_DERIVATIVE
        function [Dlambda, Dx_diag, Dx_vec, Dx_delay] = ...
                compute_derivative(a,xi_vec)
            % function [Dlambda, Dx] = compute_derivative(a,xi_vec,bool)
            % 
            % compute the derivative of a in xi_vec, of the same length of
            % xi_vec if bool==0
            %
            % INPUTS 
            % a           instance of polynomial_coefs
            % xi_vec      instance of Xi_vector
            % bool        boolean, if 0 the output has the same number of
            %             nodes as the input ---- NOT IMPLEMENTED YET 
            % OUTPUTS
            % Dlambda     derivative of the polynomial with respect to the
            %             scalar components
            % Dx_diag     coefficients in front of the diagonals (i.e. of
            %             the derivatives of x_dot)
            % Dx_vec      vectors that constitute the toeplix matrices
            global use_intlab
            size_nodes_time_deg = size_verifyfft(ones(xi_vec.nodes*2*a.deg_vector+1,1));
            Dlambda = zeros(xi_vec.size_scalar,a.n_equations,size_nodes_time_deg);
            Dx_vec = zeros(xi_vec.size_vector,a.n_equations,size_nodes_time_deg);
            Dx_diag = zeros(xi_vec.size_vector,a.n_equations);
            Dx_delay = cell(xi_vec.size_vector,a.n_equations);
            
            xi_vec= set_ifft(xi_vec,a.deg_vector);
            if use_intlab 
                Dlambda = intval(Dlambda);
                Dx_vec = intval(Dx_vec);
                Dx_diag = intval(Dx_diag);
            end
            for i = 1: xi_vec.size_scalar
                Dlambda(i,:,:)=apply(coef_derivative_scal(a,i),xi_vec,1);
            end
            
            for i =1 :xi_vec.size_vector
                Dx_vec(i,:,:) = apply(conv_coef_derivative_vec(a,i),xi_vec,1);
            end
            
            for k = 1:a.n_equations
                for j = 1:a.n_terms(k)
                    if any(a.dot{k}{j}(:,1)>0)
                        if any(abs(a.delay{k}{j})>0)
                            error('Cannot handle derivatives and delays together')
                        end
                        if a.monomial_per_term{k}(j)==1 && sum(a.dot{k}{j}(:,1))==1
                            index = find(a.dot{k}{j}(:,1));
                            Dx_diag(index,k) = Dx_diag(index,k)+ a.value{k}(j)*prod(xi_vec.scalar.^(a.power_scalar{k}(:,j).'));
                            continue
                        elseif sum(sum(a.dot{k}{j}(:,:)))~=1
                            warning('General DERIVATIVES not implemented here')
                        end
                    end
                end
            end
            psi = xi_vec.scalar(1); % the frequency
            for i = 1:xi_vec.size_vector % which unknown I'm taking derivatives of
                for j = 1:a.n_equations
                    Dx_delay_loc= cell(a.n_terms(j),1);
                    for k = 1:a.n_terms(j)
                        if a.delay{j}{k}(i) == 0 
                            % only consider the delay derivative if there
                            % is a delay in the right variable
                            Dx_delay_loc{k} = [];
                            continue
                        end
                        % convolution term
                        term = a.value{j}(k) * prod(xi_vec.scalar.^(a.power_scalar{j}(:,k).'));
                        [~,ifft_prod] = power(xi_vec,a.power_vector{j}{k}(:));
                        delete_i_delay = a.delay{j}{k}(:); % taking the derivative w.r.t. i
                        delete_i_delay(i) =0;
                        delay_x = xi_vec.delay(delete_i_delay);
                        delay_x= set_ifft(delay_x,a.deg_vector);
                        [~,ifft_prod_delay] = power(delay_x,ones(1:xi_vec.size_vector));
                        % [~,ifft_prod_delay] = power(xi_vec.delay(a.delay{j}{k}(:)),ones(1:xi_vec.size_vector));
                        Dx_delay_loc{k}.convolution =  verifyfft_in(term .* ifft_prod .* ifft_prod_delay, 1);
                        
                        % exponential term
                        tau = a.delay{j}{k}(i);
                        Dx_delay_loc{k}.exp_coef = 1i * psi * tau;
                    end
                    Dx_delay{j,i} = Dx_delay_loc;
                end
            end
        end
        % end COMPUTE_DERIVATIVE
        
        
        % COEF_DERIVATIVE_SCAL
        function beta = coef_derivative_scal(alpha, j)
            % function beta = coef_derivative_scal(alpha, j)
            % INPUT
            % alpha    instance of polynomial_coef
            % j        integer
            % OUTPUT
            % beta     instance of polynomial_coef such that 
            %       beta(x) = d_lambda_j alpha(x)
            
            if j > alpha.size_scalar
                error('Inputs are incompatible')
            end
            
            beta = alpha;
            
            for i = 1:alpha.n_equations
                for ii = 1:alpha.n_terms(i)
                    beta.value{i}(ii) = alpha.value{i}(ii) * alpha.power_scalar{i}(j,ii);
                    beta.power_scalar{i}(j,ii) = max(alpha.power_scalar{i}(j,ii)-1,0);
                end
            end
        end
        % end COEF_DERIVATIVE_SCAL
        
        % CONV_COEF_DERIVATIVE_VEC
        function beta = conv_coef_derivative_vec(alpha, j)
            if j > alpha.size_vector
                error('Inputs are incompatible')
            end
            
            beta = alpha;
            % delays get just copy-pasted
            
            for i = 1:alpha.n_equations
                for ii = 1:alpha.n_terms(i)
                    
                    this_value=alpha.value{i}(ii);
                    
                    if alpha.monomial_per_term{i}(ii) ==1 &&  sum(alpha.dot{i}{ii}(j,:))==1
                        this_value =0;
                        power_vec = 0;
                    elseif sum(alpha.dot{i}{ii}(j,:))>0
                        warning('code does not support derivatives in this format')
                    else%if sum(alpha.delay{i}{ii}(:,:))>0
                        % if there are delays, this term still can be
                        % computed without caring for them (other than they
                        % get copied)
                        %warning('Code does not support delays')
                    %else
                        power_vec = sum(alpha.power_vector{i}{ii}(:,:),2);
                        
                        this_value = this_value * power_vec(j);
                        power_vec(j) =max( power_vec(j) -1,0);
                        
                    end
                    beta.value{i}(ii) = this_value;
                    beta.power_vector{i}{ii} = power_vec;
                end
            end
        end
        % end COEF_DERIVATIVE_VEC
        
        % EQ
        function bool = eq(alpha, beta)
            % function bool = eq(alpha, beta)
            %
            % INPUT 
            % alpha     polynomial_coefs
            % beta      polynomial_coefs
            % OUTPUT
            % bool      1 if alpha == beta, 0 otherwise
            
            bool =0;
            try
            if alpha.size_scalar ~= beta.size_scalar
                return
            end
            if alpha.size_vector ~= beta.size_vector
                return
            end
            if alpha.deg_scalar ~= beta.deg_scalar
                return
            end
            if alpha.deg_vector ~= beta.deg_vector
                return
            end
            if alpha.n_equations ~= beta.n_equations
                return
            end
            if any(alpha.n_terms ~= beta.n_terms)
                return
            end
            for i=1:alpha.n_equations
                if any(alpha.value{i}~=beta.value{i})
                    return
                end
                if any(alpha.power_scalar{i}~=beta.power_scalar{i})
                    return
                end
                if any(alpha.monomial_per_term{i}~=beta.monomial_per_term{i})
                    return
                end
                for ii = 1:alpha.n_terms(i)
                    if any(alpha.dot{i}{ii}~=beta.dot{i}{ii})
                        return
                    end
                    if any(alpha.delay{i}{ii}~=beta.delay{i}{ii})
                        return
                    end
                    if any(alpha.power_vector{i}{ii}~=beta.power_vector{i}{ii})
                        return
                    end
                end
            end
            catch 
                return
            end
            bool =1;

        end
        % end EQ
        
        % NE
        function bool = ne(alpha,beta)
            % function bool = eq(alpha, beta)
            %
            % INPUT 
            % alpha     polynomial_coefs
            % beta      polynomial_coefs
            % OUTPUT
            % bool      1 if alpha ~= beta, 0 otherwise
            
            not_bool = eq(alpha,beta);
            if not_bool ==1
                bool=0;
            else 
                bool=1;
            end
        end
        % end NE
        
        % HAS_DELAY
        function bool = has_delay(alpha)
            % delay - cell{n_equations}{n_terms}(variables,:)
            bool = 0;
            for i=1:alpha.n_equations
                for j = 1:alpha.n_terms(i)
                    if any(alpha.delay{i}{j})
                        bool = 1;
                        return
                    end
                end
            end
        end
        % end HAS_DELAY
        
        
        
        
        % DISP
        function disp(alpha)
            % function disp(alpha)
            % a function to display in readable format the vector field
            
            if alpha.n_equations==0
                disp('       [Empty vector field,  0 equations]')
            end
            
            s = '     ';
            for i = 1:alpha.n_equations
                for j = 1:alpha.n_terms(i)
                    if alpha.value{i}(j)==0
                        continue
                    elseif alpha.value{i}(j) >0
                        s =  sprintf('%s + ',s);
                    elseif alpha.value{i}(j)<0
                        s =  sprintf('%s - ',s);
                    end
                    
                    if abs(alpha.value{i}(j))==1
                    else
                        s = strcat(s, best_double_print(abs(alpha.value{i}(j))));
                    end
                    term = '';
                    for k = 1:alpha.size_scalar
                        if alpha.power_scalar{i}(k,j)==0
                            continue
                        elseif alpha.power_scalar{i}(k,j)==1
                            term = sprintf('%s l%d',term ,k);
                        else
                            term = sprintf('%s l%d ^%d',term ,k,alpha.power_scalar{i}(k,j));
                        end
                    end
                    s = strcat(s,term);
                    term = '';
                    
                    if any(alpha.delay{i}{j}(:)>0)
                        index_delay = find(alpha.delay{i}{j}(:)>0);
                        for i_delay= 1:length(index_delay)
                            index = index_delay(i_delay);
                            delay_loc = alpha.delay{i_delay}{j}(index);
                            print_delay_loc = best_double_print(delay_loc);
                            term = sprintf('%s Delay( x%d, %s) ', term, index, print_delay_loc);
                        end
                    end
                    
                    for k = 1:alpha.monomial_per_term{i}(j)
                        if any(alpha.dot{i}{j}(:,k)>0)
                            if alpha.monomial_per_term{i}(j)==1 && sum(alpha.dot{i}{j}(:,k))==1
                                index = find(alpha.dot{i}{j}(:,k));
                                term = sprintf('%s dot x%d ',term,index);
                                continue
                            else
                                warning('General DERIVATIVES not implemented here')
                            end
                        end
                        for l = 1:alpha.size_vector
                            if alpha.power_vector{i}{j}(l,k) ==0
                                continue
                            elseif alpha.power_vector{i}{j}(l,k) ==1
                            term = sprintf('%s x%d ',term, l);
                            else
                            term = sprintf('%s x%d ^ %d',term, l,alpha.power_vector{i}{j}(l,k));
                            end
                        end
                    end
                    s = strcat(s,term);
                end
                s = sprintf('%s\n     ',s);
            end
            disp(s)
        end
        % end DISP
        
    end
end



function s = best_double_print(x_double)
if mod(x_double,1)==0
    s =  sprintf(' %d ', abs(x_double));
else
    if abs(x_double)>0.0001 && abs(x_double)< 10^3
        s =  sprintf(' %f ', x_double);
    else
        s =  sprintf(' %e ', x_double);
    end
end
end

    
    
    