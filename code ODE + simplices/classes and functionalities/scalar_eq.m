classdef scalar_eq
    properties
        num_equations
        number_equations_lin
        number_equations_pol
        number_of_nodes
        size_scalar
        size_vector
        linear_coef  % cell(3,1), {1}(number_equations_lin,size_scalar), {2}(number_equations_lin,size_vector,1+2*number_of_nodes), {3}(number_equations_lin)
        polynomial_equations % polynomial_coef with number_equations = number_equations_pol
        number_equations_non_computable
        non_computable_function
    end
    methods
        %CONSTRUCTOR
        function alpha = scalar_eq(n_equations_lin, n_equations_vec, n_scalar,...
                n_vector, lin_coefficients, polynomials, n_noncomputable, function_handle)
            % function alpha = scalar_eq(n_equations_lin, n_equations_vec, n_scalar, n_vector, lin_coefficients, polynomials, n_noncomputable, function_handle)
            %
            % INPUT
            % n_equations_lin       number of linear equations
            % n_equations_vec       number of polynomial equations
            % n_scalar              number of scalar variables
            % n_vector              number of vector variables
            % lin_coefficients      cell(3,1) with linear coefficients
            % polynomials           polynomial_coef with size compatible to
            %                       the previous inputs
            % n_noncomputable       number of equations that shouldn't be
            %                       delt with directly (user input)
            % function_handle       handle to a computation of the
            %                       function and its derivative, test the 
            %                       output shape by calling
            %                       _test_noncomputable_output
            % OUTPUT
            % alpha                 scalar_eq
            %
            % lin_coefficients must be of the form cell(3,1), with
            % the first cell {1}(number_equations_lin,size_scalar),
            % the second {2}(number_equations_lin,size_vector_number_of_nodes)
            % and the third {3}(number_equations_lin)
            alpha.num_equations=0;
            alpha.number_equations_lin=0;
            alpha.number_equations_pol=0;
            alpha.size_scalar=0;
            alpha.size_vector=0;
            alpha.linear_coef=[];
            alpha.polynomial_equations=[];
            alpha.number_equations_non_computable = 0;
            alpha.non_computable_function = [];
            
            if nargin >0
                if ~isint(n_equations_lin)
                    error('number of linear eqautions must be integer')
                end
                alpha.number_equations_lin = n_equations_lin;
            end
            if nargin>1
                if ~isint(n_equations_vec)
                    error('Number of polynomial equations must be integer');
                end
                alpha.number_equations_pol=n_equations_vec;
            end
            alpha.num_equations = alpha.number_equations_pol+alpha.number_equations_lin;
            if nargin>2
                if ~isint(n_scalar)
                    error('Number of scalar variables must be integer');
                end
                alpha.size_scalar = n_scalar;
            end
            if nargin>3
                if ~isint(n_vector)
                    error('Number of vector variables must be integer');
                end
                alpha.size_vector =n_vector;
            end
            if nargin>4
                try
                    if alpha.number_equations_lin>0
                        lin_coefficients{1}(alpha.number_equations_lin,alpha.size_scalar);
                        lin_coefficients{2}(alpha.number_equations_lin,alpha.size_vector,1);
                        lin_coefficients{3}(alpha.number_equations_lin);
                    else
                        if prod(size(lin_coefficients{1}))~=0 || prod(size(lin_coefficients{2}))~=0 ...
                                || prod(size(lin_coefficients{3}))~=0
                            error('The coefficients of the linear equations are not of the right dimension');
                        end
                    end
                catch
                    error('The coefficients of the linear equations are not in the specified form');
                end
                alpha.linear_coef = lin_coefficients;
                alpha.number_of_nodes=(size(lin_coefficients{2},3)-1)/2;
                if alpha.number_equations_lin>0
                    if ~isint(alpha.number_of_nodes) || alpha.number_of_nodes ==0
                        error('The number of nodes defined by the linear equations is non-integer or not positive');
                    end
                end
                % lin_coef must be cell(3,1) with specific structure
            end
            if nargin>5
                if ~is_polynomial_coef(polynomials)
                    error('The polynomial equations must be instances of the class polynomial_coef');
                end
                alpha.polynomial_equations = polynomials;
            end
            if nargin>6
                alpha.number_equations_non_computable = n_noncomputable;
                alpha.non_computable_function = function_handle;
                alpha.num_equations = alpha.num_equations + n_noncomputable;
            end
        end
        % end CONSTRUCTOR
        
        % APPLY
        function y_scal = apply(alpha, xi_vec)
            % function y_scal = apply(alpha, xi_vec)
            %
            % INPUT
            % alpha   instance of scalar_eq
            % xi_vec  Xi_vector
            % OUTPUT
            % y_scal  complex vector
            global use_intlab
            y_scal = zeros(alpha.num_equations,1);
            if use_intlab || isintval(xi_vec) || isintval(alpha)
                y_scal = intval(zeros(alpha.num_equations,1));
            end
            
            % linear equations
            if alpha.number_equations_lin>0
                if (isempty(use_intlab) || ~use_intlab) && ~isintval(xi_vec) && ~isintval(alpha)
                    temp=zeros(size(alpha.linear_coef{3}));%zeros(alpha.size_scalar);
                    for i=1:alpha.number_equations_lin
                        temp(i)=sum(sum(squeeze(alpha.linear_coef{2}(i,:,:)).*xi_vec.vector,1));
                    end
                    y_scal(1:alpha.number_equations_lin)=alpha.linear_coef{3}+(alpha.linear_coef{1}*(xi_vec.scalar.'))+temp;
                else
                    temp=intval(zeros(size(alpha.linear_coef{3})));
                    for i=1:alpha.number_equations_lin
                        temp(i)=sum(sum(squeeze(intval(alpha.linear_coef{2}(i,:,:))).*xi_vec.vector,1));
                    end
                    y_scal(1:alpha.number_equations_lin)=intval(alpha.linear_coef{3})+(intval(alpha.linear_coef{1})*(xi_vec.scalar.'))+temp;
                end
            end
            
            % polynomial equations
            if alpha.number_equations_pol
                y_scal(alpha.number_equations_lin+(1:alpha.number_equations_pol))=apply_sum(alpha.polynomial_equations,xi_vec);
            end
            % non-computable equations
            if alpha.number_equations_non_computable>0
                y_temp = alpha.non_computable_function(xi_vec);
                
                y_scal(alpha.number_equations_lin+alpha.number_equations_pol+1:end) = y_temp;
            end
        end
        % end APPLY
        
        % DERIVATIVE
        function [DlambdaG, DxG_v, DxG_hat,DxG_non_comp] = compute_derivative(alpha,xi_vec)
            % function [DlambdaG, DxG_v, DxG_hat,DxG_non_comp] = compute_derivative(alpha,xi_vec)
            %
            % INPUTS
            % alpha    instance of scalar_eq
            % xi_vec   instance of Xi_vector
            %
            % OUTPUT
            % DlambdaG   finite x finite matrix
            % DxG_v      derivative on x of vec(lambda,x) * vec(v)
            % DxG_hat    infinite derivative, such that
            % tensor_product(DxG_hat, infinite_sequence_of_1) give the
            % derivative with respect to x of p(lambda,x(0))
            global use_intlab
            if ~compatible(alpha,xi_vec)
                error('Inputs not compatible')
            end
            DlambdaG = zeros(alpha.num_equations,xi_vec.size_scalar);
            DxG_hat = zeros(alpha.number_equations_pol+alpha.number_equations_non_computable,xi_vec.size_vector);
            if isintval(alpha.linear_coef{1}) || isintval(xi_vec.scalar) || use_intlab
                DlambdaG = intval(DlambdaG);
                DxG_hat = intval(DxG_hat);
            end
            DlambdaG(1:alpha.number_equations_lin,:) = alpha.linear_coef{1};
            DxG_v = alpha.linear_coef{2};%, alpha.number_equations_lin,[]);
            
            if alpha.number_equations_pol>0
                for i = 1: xi_vec.size_scalar
                    DlambdaG(alpha.number_equations_lin+(1:alpha.number_equations_pol),i)=...
                        apply_sum(coef_derivative_scal(alpha.polynomial_equations,i),xi_vec);
                end
                for i =1 :xi_vec.size_vector
                    DxG_hat(1:alpha.number_equations_pol,i) = apply_sum(conv_coef_derivative_vec(alpha.polynomial_equations,i),xi_vec);
                end
            end
            if alpha.number_equations_non_computable>0
                [~, DF] = alpha.non_computable_function(xi_vec);
                DlambdaG(alpha.number_equations_lin+alpha.number_equations_pol+1:end,:)=...
                    DF(:,1:alpha.size_scalar);
                DxG_hat(alpha.number_equations_pol+1:end,:)= DF(:,alpha.size_scalar+1:end);
            end
            
        end
        % end DERIVATIVE
        
        % COMPATIBLE
        function bool = compatible(alpha, xi_vec)
            % function bool = compatible(alpha, xi_vec)
            %
            % INPUT
            % alpha   scalar_eq
            % xi_vec  Xi_vector
            % OUTPUT
            % bool    1 if compatible
            %         0 otherwise
            
            bool = 0;
            try
                if alpha.size_scalar ~= xi_vec.size_scalar
                    return
                elseif alpha.size_vector ~=xi_vec.size_vector
                    return
                end
            catch
                return
            end
            bool =1;
        end
        % end COMPATIBLE
        
        % INTVAL
        function alpha_intval = intval(alpha)
            alpha_intval = alpha;
            alpha_intval.polynomial_equations = intval(alpha.polynomial_equations);
            for i = 1:3
                alpha_intval.linear_coef{i} = intval(alpha.linear_coef{i});
            end
        end
        % end INTVAL
        
        % ISINTVAL
        function bool = isintval(alpha)
            bool = isintval(alpha.linear_coef{1});
        end
        % end ISINTVAL
        
        % RESHAPE
        function beta = reshape(alpha, new_nodes)
            % function beta = reshape(alpha, n_nodes)
            beta = alpha;
            if alpha.number_of_nodes ==new_nodes
                return
            end
            old_nodes=alpha.number_of_nodes;%(length(alpha.linear_coef{2})-1)/2;
            
            if old_nodes<new_nodes
                if isintval(beta.linear_coef{2})
                    pad=zeros(size(sup(alpha.linear_coef{2}),1),size(sup(alpha.linear_coef{2}),2),new_nodes-old_nodes);
                    beta.linear_coef{2}=intvalCAT(3,intval(pad),intvalCAT(3,beta.linear_coef{2},intval(pad)));
                    if size(sup(beta.linear_coef{2}),3)~= beta.number_of_nodes*2+1
                        problem = 1;
                        error('')
                    end
                    
                else
                    pad=zeros(size(alpha.linear_coef{2},1),size(alpha.linear_coef{2},2),new_nodes-old_nodes);
                    beta.linear_coef{2}=cat(3,pad,cat(3,beta.linear_coef{2},pad));
                    if size(beta.linear_coef{2},3)~= beta.number_of_nodes*2+1
                        problem = 1;
                        error('')
                    end
                end
                beta.number_of_nodes= new_nodes;
                
            elseif old_nodes==new_nodes
                return
            else
                diff=old_nodes-new_nodes;
                beta.linear_coef{2}=alpha.linear_coef{2}(:,:,diff+1:end-diff);
                beta.number_of_nodes= new_nodes;
            end
        end
        % end RESHAPE
        
        
        % RESCALE
        function beta = rescale(alpha, index_lambda, rescaling)
            beta = alpha;
            beta.linear_coef{1}(:,index_lambda) = alpha.linear_coef{1}(:,index_lambda) / rescaling;
            beta.polynomial_equations = rescale(alpha.polynomial_equations, index_lambda, rescaling);
        end
        % end RESCALE
        
        
        % CHANGE_LIN_COEF
        function beta = change_lin_coef(alpha,lin_coef,n_lin_coef)
            % function beta = change_lin_coef(alpha,lin_coef,n_lin_coef)
            % INPUT
            % alpha       scalar_eq
            % lin_coef    cell{3}, compatible with linear_coef, that is
            %             {1}(1,size_scalar), {2}(1,size_vector,1+2*number_of_nodes), {3}(1)
            % n_lin_coef  number of hte equation to be rewritten
            % OUTPUT
            % beta        scalar_eq with the nw line
            
            if ~iscell(lin_coef) ||...
                    size(lin_coef{1},1) ~=1 || size(lin_coef{1},2) ~= alpha.size_scalar ||...
                    size(lin_coef{2},1)~=1 || size(lin_coef{2},2)~=alpha.size_vector || size(lin_coef{2},3)~= 1+ 2*alpha.number_of_nodes||...
                    max(size(lin_coef{3}))~=1
                error('Inputs not coeherent')
            elseif n_lin_coef > alpha.number_equations_lin
                if  n_lin_coef == 1+ alpha.number_equations_lin
                    alpha.number_equations_lin =1+ alpha.number_equations_lin;
                    alpha.num_equations =1+ alpha.num_equations;
                elseif  n_lin_coef == 1+ alpha.num_equations
                    n_lin_coef = 1+ alpha.number_equations_lin;
                    alpha.number_equations_lin =1+ alpha.number_equations_lin;
                    alpha.num_equations =1+ alpha.num_equations;
                else
                    error('It is not possible to rewrite or append such equation, because the refered number is too big')
                end
            end
            beta = alpha;
            beta.linear_coef{1}(n_lin_coef,:) = lin_coef{1};
            beta.linear_coef{2}(n_lin_coef,:,:) = lin_coef{2};
            beta.linear_coef{3}(n_lin_coef) = lin_coef{3};
            if size(beta.linear_coef{3},1)==1 && size(beta.linear_coef{3},2)>1
                beta.linear_coef{3}=beta.linear_coef{3}.';
            end
        end
        % end CHANGE_LIN_COEF
        
        % CHANGE_LIN_COEF_VECTOR
        function beta = change_lin_coef_vector(alpha,lin_coef_vec,n_lin_coef)
            % function beta = change_lin_coef(alpha,lin_coef,n_lin_coef)
            % INPUT
            % alpha         scalar_eq
            % lin_coef_vec  vector of length
            %               size_scalar+size_vector*(2*nodes+1) +1
            %               the last element being the constant
            % n_lin_coef    number of hte equation to be rewritten
            % OUTPUT
            % beta          scalar_eq with the new line
            
            if ~isa(lin_coef_vec,'double')
                error('New linear coefficients must be stored in a vector format')
            end
            if min(size(lin_coef_vec))~=1
                error('New linear coefficients must be stored in a vector format')
            end
            if length(lin_coef_vec)~= alpha.size_scalar+ alpha.size_vector*(2*alpha.number_of_nodes+1)+1
                error('Lenght of linear coefficients incompatible')
            end
            if nargin<3
                n_lin_coef = alpha.number_equations_lin+1;
            end
            coefs_linear_lor = cell(1,3);
            coefs_linear_lor{1} = horiz(lin_coef_vec(1:alpha.size_scalar));
            coefs_linear_lor{2} = shiftdim(reshape(lin_coef_vec(alpha.size_scalar+1:end-1),2*alpha.number_of_nodes+1,alpha.size_vector).',-1);
            coefs_linear_lor{3} = lin_coef_vec(end);
            beta = change_lin_coef(alpha,coefs_linear_lor,n_lin_coef);
            
        end
        % end CHANGE_LIN_COEF_VECTOR
        
        % CHANGE_LIN_COEF_XI
        function beta = change_lin_coef_Xi(alpha, lin_coef_Xi_vec, const, n_lin_coef)
            % function beta = change_lin_coef_Xi(alpha, lin_coef_Xi_vec, const, n_lin_coef)
            %
            % INPUT
            % alpha             scalar_eq
            % lin_coef_Xi_vec   Xi_vector, compatible with alpha
            % const             scalar
            % n_lin_coef        number of the equation to be changed
            % OUTPUT
            % beta              scalar_eq, with the requested equation
            % changed into x.*lin_coef_Xi_vec - const = 0
            if any(size(const)~=1)
                error('The constant must be a double, not a vector');
            end
            if length(lin_coef_Xi_vec)~= alpha.size_scalar+ alpha.size_vector*(2*alpha.number_of_nodes+1)
                error('Lenght of linear coefficients incompatible')
            end
            if nargin<3
                n_lin_coef = alpha.number_equations_lin+1;
            end
            lin_coef_vec= cat(Xi_vec2vec(lin_coef_Xi_vec),const);
            beta = change_lin_coef_vector(alpha,lin_coef_vec,n_lin_coef);
        end
        % end CHANGE_LIN_COEF_XI
        
        % REMOVE_LIN_COEF
        function beta = remove_lin_coef(alpha,n_lin_coef)
            % function beta = remove_lin_coef(alpha,n_lin_coef)
            % INPUT
            % alpha       scalar_eq
            % n_lin_coef  number of the equation to be removed
            % OUTPUT
            % beta        scalar_eq with the nw line
            
            if  n_lin_coef > alpha.number_equations_lin
                error('It is not possible to remove such equation, because the refered number is too big')
            end
            beta = alpha;
            beta.linear_coef{1}(n_lin_coef,:) = [];
            beta.linear_coef{2}(n_lin_coef,:,:) = [];
            beta.linear_coef{3}(n_lin_coef) = [];
            beta.num_equations= beta.num_equations -1;
            beta.number_equations_lin = beta.number_equations_lin -1;
        end
        % end REMOVE_LIN_COEF
        
        
        % EXTRACT_LIN_COEF
        function v = extract_lin_coef(alpha, n_lin_coef)
            % function v = extract_lin_coef(alpha, n_lin_coef)
            %
            % this function returns the vector v that defines the
            % n_lin_coef linear equation, that is, the linear equation
            % n_lin_coef of alpha is <v,x> + c
            % INPUT
            % alpha         scalar_eq
            % n_lin_coef    integer
            % OUTPUT
            % v             complex vector
            
            if n_lin_coef > alpha.number_equations_lin
                error('Inputs are incompatible, you must request an equation that exists already')
            end
            
            v_1 = alpha.linear_coef{1}(n_lin_coef,:).';
            v_2 = reshape(squeeze(alpha.linear_coef{2}(n_lin_coef,:,:)).',[],1);
            
            v = [v_1;v_2];
            
        end
        % end EXTRACT_LIN_COEF
        
        % EXTRACT_ALL_LIN_COEF
        function v = extract_all_lin_coef(alpha)
            % function v = extract_all_lin_coef(alpha)
            %
            % INPUT
            % alpha     scalar_eq
            % OUTPUT
            % v         matrix alpha.number_equations_lin x (alpha.nodes+1)*2
            
            v = zeros( alpha.number_equations_lin , alpha.size_scalar+alpha.size_vector*(alpha.number_of_nodes*2+1));
            
            for i = 1:alpha.number_equations_lin
                coef = extract_lin_coef(alpha, i);
                v(i,:) =coef;
            end
        end
        % end EXTRACT_ALL_LIN_COEF
        
        % DISP
        function disp(alpha)
            fprintf('   Number of scalar equations: %i,   of which linear: %i,   polynomial: %i,   others: %i\n',alpha.num_equations, alpha.number_equations_lin, alpha.number_equations_pol, alpha.number_equations_non_computable)
            fprintf('   Expects a Xi_vector of scalar size: %i,   vector size: %i,   number of nodes: %i\n', alpha.size_scalar, alpha.size_vector, alpha.number_of_nodes)
            if alpha.number_equations_pol>0
                disp('   Scalar nonlinear equations:'); disp(alpha.polynomial_equations);
            end
        end
        % end DISP
    end
end
