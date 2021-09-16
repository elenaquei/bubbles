classdef full_problem
    properties %(SetAccess=private)
        scalar_equations % scalar_eq
        vector_field % polynomial_coef
    end
    methods
        % CONSTRUCTOR
        function h = full_problem(scalar_coef,vector_coef)
            % function h = full_problem(scalar_coef,vector_coef)
            %
            % constructing a full problem H, also for continuity, based on
            % the scalar coefficients and the vector coefficients
            %
            % INPUT
            % scalar_coef    instance of the class scalar_eq
            % vector_coef    instance of the class polynomial_coef
            % OUTPUT
            % h              instance of full_problem
            
            h.scalar_equations=scalar_eq();
            h.vector_field=polynomial_coefs();
            if nargin>0
                if ~isa(scalar_coef,'scalar_eq')
                    error('First input must be of the class scalar_eq');
                end
                h.scalar_equations = scalar_coef;
            end
            if nargin>1
                if ~isa(vector_coef, 'polynomial_coefs')
                    error('Second input must be of the class polynomial_coef');
                end
                h.vector_field = vector_coef;
            end
        end
        % end CONSTRUCTOR
        
        
        function y_xi_vector = feval(alpha, x_xi_vector, varagin)
            y_xi_vector = apply( alpha, x_xi_vector,varagin{:});
        end
        
        % APPLY
        function y_xi_vector = apply( alpha, x_xi_vector,bool)
            % function y_xi_vector = apply( alpha, x_xi_vector,bool)
            %
            % INPUTS
            % alpha          instance of full_problem
            % x_xi_vector    Xi_vector
            % bool           if bool ==0 output of same dimensions as input
            %
            % OUTPUT
            % y_xi_vector   Xi_vector
            
            if nargin<3
                bool=0;
            end
            
            %y_xi_vector = Xi_vector(x_xi_vector);
            %if bool
            old_nodes = x_xi_vector.nodes;
            new_nodes = x_xi_vector.nodes* alpha.vector_field.deg_vector;
            x_xi_vector = reshape_Xi(x_xi_vector,new_nodes);
            y_xi_vector = Xi_vector(alpha.scalar_equations.num_equations, alpha.vector_field.n_equations, x_xi_vector.nodes);
            alpha = reshape(alpha, new_nodes);
            %end
            y_xi_vector.scalar = apply(alpha.scalar_equations,x_xi_vector);
            % y_xi_vector.vector = apply(alpha.vector_field,x_xi_vector);
            y_xi_vector.ifft_vector = apply(alpha.vector_field,x_xi_vector);
            y_xi_vector.bool_ifft = 1;
            y_xi_vector = set_vector_from_ifft(y_xi_vector);
            y_xi_vector.vector = [y_xi_vector.vector,zeros(y_xi_vector.size_vector,1)];
            y_xi_vector.nodes = (size( y_xi_vector.vector,2)-1)/2;
            if bool==0
                y_xi_vector = reshape_Xi(y_xi_vector,old_nodes);
            else
                y_xi_vector = reshape_Xi(y_xi_vector,old_nodes*alpha.vector_field.deg_vector);
            end
        end
        % end APPLY
        
        % COMPATIBLE_VEC
        function bool = compatible_vec(alpha,x_Xi_vec)
            % function bool = compatible_vec(alpha,x_Xi_vec)
            % this function checks if the two polynomial_coefs alpha and beta
            % have the same number of variables and define the same number of
            % equations
            bool = 0;
            try
                if(alpha.scalar_equations.size_scalar~=x_Xi_vec.size_scalar)
                    return
                end
                if alpha.scalar_equations.size_vector~= x_Xi_vec.size_vector
                    return
                end
            catch
                return
            end
            bool = 1;
        end
        % end COMPATIBLE_VEC
        
        % SQUARE
        function bool = square(alpha)
            % function bool = square(alpha)
            % checking if the output has the same dimension of the input
            bool = 0;
            try
                if(alpha.scalar_equations.size_scalar~=alpha.scalar_equations.num_equations)
                    return
                end
                if alpha.scalar_equations.size_vector~= alpha.vector_field.n_equations
                    return
                end
            catch
                return
            end
            bool = 1;
        end
        % end SQUARE
        
        
        % RESHAPE
        function beta = reshape(alpha, new_nodes)
            % function beta = reshape(alpha, new_nodes)
            % change the scalar equations to thewanted number of nodes
            beta = alpha;
            beta.scalar_equations = reshape(alpha.scalar_equations,new_nodes);
        end
        % end RESHAPE
        
        
        % RESCALE
        function beta = rescale(alpha, index_lambda, rescaling)
            if index_lambda > alpha.scalar_equations.size_scalar + ...
                    alpha.scalar_equations.size_vector
                error('Index of rescaled variable too big')
            end
            beta = alpha;
            if index_lambda <= alpha.scalar_equations.size_scalar
                beta.scalar_equations = rescale(alpha.scalar_equations, index_lambda, rescaling);
                beta.vector_field = rescale(alpha.vector_field, index_lambda, rescaling);
            else
                error('not coded yet');
            end
            
        end
        % end RESCALE
        
        % FUNCTION_DIRECTIONAL_SECOND_DERIVATIVE
        function DDH=Function_directional_second_derivative(alpha,xBar,dir1,dir2)
            %  function DDH=Function_directional_second_derivative(alpha,xBar,dir1,dir2)
            %
            % INPUT:
            % xBar   Xi_vector
            % alpha  full_problem
            % dir1   Xi_vector
            % dir2   Xi_vector
            %
            % OUTPUT:
            % DDH     Xi_vector
            
            global use_intlab;
            
            DDH=Xi_vector();
            DDH.size_scalar=alpha.scalar_equations.num_equations;
            DDH.size_vector=alpha.vector_field.n_equations;
            DDH.scalar = zeros(1,DDH.size_scalar);
            if use_intlab
                DDH=intval(DDH);
            end
            
            
            size_scal=xBar.size_scalar;
            
            end_size_vec=xBar.nodes*alpha.vector_field.deg_vector;
            
            % zero padding starts here
            halfsize_vec=(xBar.nodes-1)*alpha.vector_field.deg_vector+2;
            DDH=reshape_Xi(DDH,halfsize_vec);
            
            if use_intlab
                DDH.vector=intval(DDH.vector);
            end
            
            
            % vector part
            % computing
            
            % way to go:
            % ifft(xvector)
            xBar=reshape_Xi(xBar,DDH.nodes);
            dir1=reshape_Xi(dir1,DDH.nodes);
            dir2=reshape_Xi(dir2,DDH.nodes);
            
            
            % preallocation for speed
            x_conv=get_ifft(xBar);
            dir1_conv=get_ifft(dir1);
            dir2_conv=get_ifft(dir2);
            
            % computation
            
            % second derivative of F(x)
            temp_DDF=zeros(alpha.vector_field.n_equations,length(x_conv)); %% JUST WORKS IF DIM < NODES
            
            if isintval(alpha.vector_field.value) || isintval(x_conv)
                temp_DDF = intval(temp_DDF);
            end
            
            for i=1:alpha.vector_field.n_equations % equation
                for j=1:alpha.vector_field.n_terms(i) % element of the equation
                    d=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
                    const=alpha.vector_field.value{i}(j);
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
                            temp_DDF(i,:)= temp_DDF(i,:) + const*d(k)*(d(n)-e(k,n))*...
                                prod(xBar.scalar.^(power_scalar.'))*...
                                product_FFT(x_conv,power_vector).*...
                                temp_dir1.*...
                                temp_dir2;
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
            
            % conversion back
            temp_DDH = [temp_DDG;temp_DDF];
            temp_y_vector=0*temp_DDH;
            for i=1:xBar.size_vector
                if all(temp_DDH(i,:)==0)
                    continue
                end
                temp_y_vector(i,:)=verifyfft_in(temp_DDH(i,:),1,1);
            end
            nodes_y_vector=(length(temp_y_vector(i,:))-2)/2;
            DDH.vector=temp_y_vector(:,2:end);
            DDH.nodes=nodes_y_vector;
            
            
            % reshape
            DDH=reshape_Xi(DDH,end_size_vec);
        end
        % end FUNCTION_DIRECTIONAL_SECOND_DERIVATIVE
        
        
        % FUNCTION_SECOND_DERIVATIVE
        function DDH=Function_second_derivative(alpha,xBar,RAD_MAX,xBar_norm)
            % function Function_second_derivative_new(alpha,xBar,RAD_MAX)  GOOD
            %
            % INPUT:
            % xBar     Xi_vector
            % alpha    full_problem
            % RAD_MAX  positive scalar
            %
            % OUTPUT:
            % DDH     vector
            
            global nu;
            global use_intlab;
            
            
            if ~isa(alpha, 'full_problem')
                error('input not suitable')
            end
            
            
            
            DDH=zeros(alpha.scalar_equations.num_equations+alpha.vector_field.n_equations,1);
            if use_intlab
                %    DDH=intval(DDH);
                % the output should be the maximum of the interval anyway
            end
            if nargin == 4
                XnormC = xBar_norm;
            else
                XnormC=cnorm_Xi_vector(xBar,nu);
            end
            
            xANDr=XnormC+RAD_MAX;
            
            size_scal=xBar.size_scalar;
            
            % here for the polynomial scalar equations
            
            % the derivative is null for hte linear part,
            % alpha.scalar_equation.number_equations_lin, but non-null for the scalar
            % polynomial part
            
            n_scalar_equations = alpha.scalar_equations.number_equations_lin;
            
            beta = alpha.scalar_equations.polynomial_equations;
            
            for i=1:alpha.scalar_equations.number_equations_pol % equation
                for j=1:beta.n_terms(i) % element of the equation
                    if any(size(beta.power_vector{i}{j},2)>1)
                        error('Second derivative not yet ready for this')
                    end
                    if beta.value{i}(j) == 0
                        continue
                    end
                    k_vec = j;
                    for kk = j+1:beta.n_terms(i)
                        if  all(beta.power_vector{i}{j} ==  beta.power_vector{i}{kk})
                            k_vec(end+1) = kk;
                        end
                    end
                    d=[beta.power_scalar{i}(:,j).', beta.power_vector{i}{j}.'];
                    N=length(d);
                    e=eye(N);
                    
                    scalar_part_x = xANDr(1:size_scal);
                    vector_part_x = xANDr(size_scal+1:end);
                            
                    % DDH_i = sum_{k=1,d_k>0}^N sum_{j=1,(d-e_k)_n>0}^N
                    % (d-e_k)_n d_k x^{d-e_k-e_n}
                    
                    for k=1:N
                        if d(k)<=0
                            continue
                        end
                        for n=1:N
                            if d(n)-e(n,k)<=0
                                continue
                            end
                            const_vec = 0*k_vec;
                            for kk = 1:length(k_vec)
                                power_scalar_kk=beta.power_scalar{i}(:,k_vec(kk))-e(k,1:size_scal).'-e(n,1:size_scal).';
                                d=[beta.power_scalar{i}(:,k_vec(kk)).', beta.power_vector{i}{j}.'];
                                if use_intlab
                                    const_vec(kk)=d(k)*(d(n)-e(k,n))*sup(abs(beta.value{i}(k_vec(kk))))*prod(scalar_part_x.^power_scalar_kk);
                                else
                                    const_vec(kk)=d(k)*(d(n)-e(k,n))*abs(beta.value{i}(k_vec(kk)))*prod(scalar_part_x.^power_scalar_kk);
                                end
                            end
                            const = sum(const_vec);
                            power_loc = (d-e(k,:)-e(n,:)).';
                            power_vector_x_k = power_loc(size_scal+1:end);
                            
                            DDH(i+n_scalar_equations)= DDH(i+n_scalar_equations) + abs(const*...
                                prod(vector_part_x.^power_vector_x_k));
                            
                        end
                    end
                    beta.value{i}(k_vec) = 0;
                    
                end
            end
            
            % here for the vector field part
            
            n_scalar_equations = alpha.scalar_equations.num_equations;
            alpha = alpha.vector_field;
            
            for i=1:xBar.size_vector % equation
                for j=1:alpha.n_terms(i) % element of the equation
            
                    if any(size(alpha.power_vector{i}{j},2)>1)
                        error('Second derivative not yet ready for this')
                    end
                    if alpha.value{i}(j) == 0
                        continue
                    end
                    k_vec = j;
                    for kk = j+1:alpha.n_terms(i)
                        if  all(alpha.power_vector{i}{j} ==  alpha.power_vector{i}{kk})
                            k_vec(end+1) = kk;
                        end
                    end
                    d=[alpha.power_scalar{i}(:,j).', alpha.power_vector{i}{j}.'];
                    N=length(d);
                    e=eye(N);
                    
                    scalar_part_x = xANDr(1:size_scal);
                    vector_part_x = xANDr(size_scal+1:end);
                            
                    % DDH_i = sum_{k=1,d_k>0}^N sum_{j=1,(d-e_k)_n>0}^N
                    % (d-e_k)_n d_k x^{d-e_k-e_n}
                    
                    for k=1:N
                        if d(k)<=0
                            continue
                        end
                        for n=1:N
                            if d(n)-e(n,k)<=0
                                continue
                            end
                            const_vec = 0*k_vec;
                            for kk = 1:length(k_vec)
                                power_scalar_kk=alpha.power_scalar{i}(:,k_vec(kk))-e(k,1:size_scal).'-e(n,1:size_scal).';
                                d=[alpha.power_scalar{i}(:,k_vec(kk)).', alpha.power_vector{i}{j}.'];
                                if use_intlab
                                    const_vec(kk)=d(k)*(d(n)-e(k,n))*sup(abs(alpha.value{i}(k_vec(kk))))*prod(scalar_part_x.^power_scalar_kk);
                                else
                                    const_vec(kk)=d(k)*(d(n)-e(k,n))*abs(alpha.value{i}(k_vec(kk)))*prod(scalar_part_x.^power_scalar_kk);
                                end
                            end
                            const = sum(const_vec);
                            power_loc = (d-e(k,:)-e(n,:)).';
                            power_vector_x_k = power_loc(size_scal+1:end);
                            
                            DDH(i+n_scalar_equations)= DDH(i+n_scalar_equations) + abs(const*...
                                prod(vector_part_x.^power_vector_x_k));
                            
                        end
                    end
                    alpha.value{i}(k_vec) = 0;
                end
                    
            end
        
        end
        % end FUNCTION_SECOND_DERIVATIVE
        
        % FUNCTION_THIRD_DERIVATIVE
        function DDDH=Function_third_derivative(alpha,xBar,RAD_MAX)
            % function Function_third_derivative(alpha,xBar,RAD_MAX)
            %
            % INPUTS
            % xBar     Xi_vector, numerical solution
            % alpha    full_problem, vector field and scalar equations
            % RAD_MAX  double, estimation of the maximum validation radius
            %
            % OUTPUT:
            % DDDH     vector
            
            global nu;
            global use_intlab;
            size_alpha = alpha.scalar_equations.num_equations + alpha.vector_field.n_equations;
            DDDH=zeros(size_alpha,1);
            
            XnormC=cnorm_Xi_vector(xBar,nu);
            
            xANDr=XnormC+RAD_MAX;
            
            % DDG
            for i=1:alpha.scalar_equations.number_equations_pol
                %1+alpha.scalar_equations.number_equations_lin : alpha.vector_field.n_equations% equation
                for j=1:alpha.scalar_equations.polynomial_equations.n_terms(i) % element of the equation
                    d=[alpha.scalar_equations.polynomial_equations.power_scalar{i}(:,j).', alpha.scalar_equations.polynomial_equations.power_vector{i}{j}.'];
                    if use_intlab
                        const=sup(abs(alpha.scalar_equations.polynomial_equations.value{i}(j)));
                    else
                        const=abs(alpha.scalar_equations.polynomial_equations.value{i}(j));
                    end
                    N=length(d);
                    e=eye(N);
                    
                    % DDH_i = sum_{k=1,d_k>0}^N sum_{n=1,(d-e_k)_n>0}^N sum_{m=1,(d-e_k-e_n)_m>0}^N
                    % (d-e_k-e_n)_m (d-e_k)_n d_k x^{d-e_k-e_n}
                    
                    for k=1:N
                        if d(k)<0
                            continue
                        end
                        for n=1:N
                            if d(n)-e(n,k)<=0
                                continue
                            end
                            
                            for m=1:N
                                if d(m)-e(m,k)-e(n,k)<=0
                                    continue
                                end
                                DDDH(i+alpha.scalar_equations.number_equations_lin)= DDDH(i+alpha.scalar_equations.number_equations_lin)...
                                    + abs(const*d(k)*(d(n)-e(k,n))*(d(m)-e(m,k)-e(n,k))*...
                                    prod(xANDr.^((d-e(k,:)-e(n,:)-e(m,:)).')));
                            end
                        end
                    end
                    
                    
                end
            end
            
            % DDF
            for i=1 : alpha.vector_field.n_equations% equation
                for j=1:alpha.vector_field.n_terms(i) % element of the equation
                    d=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
                    if use_intlab
                        const=sup(abs(alpha.vector_field.value{i}(j)));
                    else
                        const=abs(alpha.vector_field.value{i}(j));
                    end
                    N=length(d);
                    e=eye(N);
                    
                    % DDH_i = sum_{k=1,d_k>0}^N sum_{n=1,(d-e_k)_n>0}^N sum_{m=1,(d-e_k-e_n)_m>0}^N
                    % (d-e_k-e_n)_m (d-e_k)_n d_k x^{d-e_k-e_n}
                    
                    for k=1:N
                        if d(k)<0
                            continue
                        end
                        for n=1:N
                            if d(n)-e(n,k)<=0
                                continue
                            end
                            
                            for m=1:N
                                if d(m)-e(m,k)-e(n,k)<=0
                                    continue
                                end
                                DDDH(i+alpha.scalar_equations.num_equations)= DDDH(i+alpha.scalar_equations.num_equations)...
                                    + abs(const*d(k)*(d(n)-e(k,n))*(d(m)-e(m,k)-e(n,k))*...
                                    prod(xANDr.^((d-e(k,:)-e(n,:)-e(m,:)).')));
                            end
                        end
                    end
                    
                    
                end
            end
        end
        % end FUNCTION_THIRD_DERIVATIVE
        
        
        
        %function disp(alpha)
        %
        %end
        
    end
end