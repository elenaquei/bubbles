classdef Xi_vector
    properties
        size_scalar
        size_vector
        nodes
        scalar
        vector % N x (nodes*2+1)
        bool_ifft
        %end
        %properties(SetAccess = private)
        ifft_vector
    end
    methods
        % CONSTRUCTOR
        function x=Xi_vector(scal,vec,size_scal,size_vec,nodes)
            %
            % Xi_vector   -  constructor of the class Xi_vector, can be
            % called in the following formats
            %    Xi_vector()  or Xi_vector,  created an empty Xi_vector
            %    Xi_vector(Xi)    - Xi is a Xi_vector, the result will be a
            %         null Xi_vector with same dimensions as Xi
            %    Xi_vector(double[] x1,double[][] x2)  - it will create a
            %         Xi_vector with scalar part x1 and vector part x2
            %    Xi_vector(double[] x1, double[][] x2, int dim1, int dim2) -
            %         creates a Xi_vector with scalar part x1 and vector
            %         part x2 testing that the dimensions coincide with the
            %         given, therefore dim(x1)=dim1, dim(x2,1)=dim2
            %    Xi_vector(double[] x1, double[][] x2, int dim1, int dim2, int dim3) -
            %         creates a Xi_vector with scalar part x1 and vector
            %         part x2 testing that the dimensions coincide with the
            %         given, therefore dim(x1)=dim1, dim(x2,1)=dim2,
            %         dim(x2,2)=dim3
            %    Xi_vector(int, int, int) - creates a null Xi_vector with
            %         the given dimensions, in order number of scalars, of
            %         vectors and of nodes
            %
            if nargin>0
                if nargin ==3
                    x.scalar=zeros(scal,1);
                    x.vector=zeros(vec,size_scal*2+1);
                    x.size_scalar=scal;
                    x.size_vector=vec;
                    x.nodes=size_scal;
                    return
                end
                if nargin==1
                    x.size_scalar=scal.size_scalar;
                    x.size_vector=scal.size_vector;
                    x.nodes=scal.nodes;
                    x.scalar=0*scal.scalar;
                    x.vector=0*scal.vector;
                else
                    if min(size(scal))>1
                        error('need a vector inout, not a matrix')
                    end
                    x.scalar=scal;
                    x.size_scalar=length(scal);
                end
                if nargin>1
                    x.vector=vec;
                    x.nodes=(size(vec,2)-1)/2;
                    x.size_vector=size(vec,1);
                    if nargin>2
                        if ~((isequal(size(scal),[size_scal,1])) || isequal(size(scal),[1,size_scal]))
                            error('dimensions do not agree')
                        end
                        x.size_scalar=size_scal;
                        if nargin>3
                            if (size(vec,1) ~= size_vec)
                                if (size(vec,2) ~= size_vec)
                                    error('dimensions do not agree')
                                else
                                    x.vector=vec.';
                                end
                            end
                            x.size_vector=size_vec;
                            if nargin>4
                                if (size(x.vector,2) ~= nodes*2+1)
                                    error('dimensions do not agree')
                                end
                                x.nodes=nodes;
                            end
                        end
                    end
                end
                x.bool_ifft = 0;
                %x=x.set_ifft();
            end
        end %CONSTRUCTOR
        
        % SET_IFFT
        function x = set_ifft(x,deg)
            % function x = set_fft(x,deg)
            % compute fft(x.vector) and store it - if not computed yet
            global use_intlab
            if nargin==1
                deg=1;
            end
            
            dim_ifft = size_verifyfft([ zeros(1, (deg-1)*x.nodes), x.vector(1,:), zeros(1, (deg-1)*x.nodes)]);
            
            if all(size(x.bool_ifft)==[1,1]) && x.bool_ifft ==1 && max(size(x.ifft_vector)) > dim_ifft
                return
            else
                x_vector = [ zeros(x.size_vector, (deg-1)*x.nodes),x.vector, zeros(x.size_vector, (deg-1)*x.nodes)]; % appending zeros
                x.ifft_vector = zeros(x.size_vector, dim_ifft);
                if use_intlab
                    x.ifft_vector = intval(x.ifft_vector);
                end
                for i=1:x.size_vector
                    x.ifft_vector(i,:) = verifyfft_in(x_vector(i,:),-1,1);
                end
                x.bool_ifft =1;
            end
        end
        % end SET_IFFT
        
        % SET_VECTOR_FROM_IFFT
        function x = set_vector_from_ifft(x)
            % function x = set_vector_from_ifft(x)
            % compute fft(x.ifft_vector) and store it
            if x.bool_ifft ==0
                error('This function takes the fft and returns the initial vector, not possible without the fft');
            else
                x.vector = 0*x.ifft_vector;
                for i=1:x.size_vector
                    x.vector(i,:) = verifyfft_in(x.ifft_vector(i,:),1);
                end
                x.bool_ifft =1;
            end
        end
        % end SET_VECTOR_FROM_IFFT
        
        % GET_IFFT
        function v = get_ifft(x)
            %function v = get_ifft(x)
            % get the fft(x.vector) , by computing it if necessary
            if x.bool_ifft==0 || isempty(x.ifft_vector)
                x=x.set_ifft();
            end
            v=x.ifft_vector;
        end
        % end GET_IFFT
        
        
        % POWER
        function [y,ifft_y] = power(x,p)
            % function y = power(x,p)
            % returns a vector such that fft(y) is equal to the
            % convolutions of x taken p times
            % ex: x = (x1, x2, x3), p = [1,2,0]
            % y is such that fft(y) = conv( x1,x2,x2)
            %
            % Remark: using x.ifft, it means that y= x.ifft_vector.^p
            if x.bool_ifft==0
                x=x.set_ifft();
            end
            z = 5;
            if length(p)~=x.size_vector
                error('Powers are not compatible')
            end
            
            y = 1+0*x.ifft_vector(1,:);
            if sum(p)==0
                ifft_y = y;
                y = verifyfft_in(y,1);
                return
            end
            for i = 1:length(p)
                y = y.*(x.ifft_vector(i,:).^p(i));
            end
            
            ifft_y = y;
            if nargout <= 1
                y = verifyfft_in(ifft_y,1);
            else
                y = NaN;
            end
        end
        % end POWER
        
        
        
        % EQUAL_XI
        function y=equal_Xi(x)
            % x Xi_vector, y identical Xi_vector
            y=Xi_vector();
            y.size_scalar=x.size_scalar;
            y.size_vector=x.size_vector;
            y.nodes=x.nodes;
            y.scalar=x.scalar;
            y.vector=x.vector;
        end
        % EQUAL_XI
        
        % NORM
        function nor= norm(xi)
            % function nor= norm(xi)
            global nu
            nor=cnorm_Xi_vector(xi,nu);
        end
        % NORM
        
        % CNORM_XI_VECTOR
        function nor=cnorm_Xi_vector(xi,nu)
            % cnorm_Xi_vector(Xi_vector,double)
            % in this function, the component-wise norm of a Xi_vector is
            % taken, that means that the output will be a positive double array
            % of size size_scalar+size_vector, storing in the first part the
            % absolute value of the Xi_vector, in the second part the l1-norm
            % of the vectors of the Fourier coefficients. This routine is
            % called by norm_Xi_vector, that returns the maximum of
            % cnorm_Xi_vector
            
            global use_intlab
            
            nor=zeros(xi.size_scalar+xi.size_vector,1);
            if ~use_intlab
                nor(1:xi.size_scalar)=abs(xi.scalar);
            else
                nor(1:xi.size_scalar)=sup(abs(xi.scalar));
            end
            %nor(xi.size_scalar+1:end)=abs(xi.vector(:,xi.nodes+1));
            K=-xi.nodes:xi.nodes;
            if use_intlab
                K=intval(K);
                if ~isintval(nu)
                    nu=intval(nu);
                end
            else
                if ~isfloat(nu)
                    nu=sup(nu);
                end
            end
            for i=1:xi.size_vector
                if ~use_intlab
                    nor(xi.size_scalar+i)=sum(abs(xi.vector(i,:)).*(nu.^abs(K)));
                else
                    nor(xi.size_scalar+i)=sup(sum(abs(xi.vector(i,:)).*(nu.^abs(K))));
                end
            end
        end %CNORM
        
        
        % NORM_XI_VECTOR
        function nor=norm_Xi_vector(xi,nu)
            % norm_Xi_vector(Xi_vector,double)
            % return the maximum-norm of the Xi_vector xi
            nor=max(cnorm_Xi_vector(xi,nu));
        end
        % NORM_XI_VECTOR
        
        % C_DIVISION
        function xi_div=cdivision_Xi_vector(xi,yi,flag)
            % cdivision_Xi_vector(Xi_vector,Xi_vector,flag)
            % the first vector is divided component-wise by the second vector,
            % with the convention that division by 0 is 0. The two Xi_vectors
            % must have the same size_scalar and size_vector. If they have
            % different number of nodes, the smallest is chosen for flag=1,
            % otherwie the biggest, padding with 0s the tails of the
            % smallest vector.
            if xi.size_scalar ~= yi.size_scalar
                error('Scalar dimensions must agree')
            elseif xi.size_vector ~=yi.size_vector
                error('Vector dimensions must agree')
            elseif xi.nodes ~= yi.nodes
                warning('Number of nodes must agree')
                if flag
                    if xi.nodes<yi.nodes
                        yi=reshape_Xi(yi,xi.nodes);
                    else
                        xi=reshape_Xi(xi,yi.nodes);
                    end
                else
                    if xi.nodes<yi.nodes
                        xi=reshape_Xi(xi,yi.nodes);
                    else
                        yi=reshape_Xi(yi,xi.nodes);
                    end
                end
            end
            xi_div=Xi_vector(xi);
            I=find(yi.scalar);
            if size(xi.scalar)==size(yi.scalar)
                xi_div.scalar(I)=xi.scalar(I)./yi.scalar(I);
            else
                xi_div.scalar(I)=(xi.scalar(I).')./yi.scalar(I);
            end
            for i=1:xi.size_vector
                I=find(yi.vector(i,:));
                xi_div.vector(i,I)=xi.vector(i,I)./yi.vector(i,I);
            end
            xi_div.bool_ifft = 0;
            %[I,J]=find(yi.vector);
            %xi_div.vector(I,J)=xi.vector(I,J)./yi.vector(I,J);% this notation do NOT work
        end
        %C_DIVISION
        
        % symmetrise
        function y_vec=symmetrise(xi_vec)
            %function y_vec=symmetrised(xi_vec)
            y_vec=xi_vec;
            y_vec.scalar=(xi_vec.scalar+conj(xi_vec.scalar))/2;
            %m=xi_vec.nodes;
            for i=1:xi_vec.size_vector
                %y_vec.vector(i,1:m)=conj(flip(xi_vec.vector(i,m+2:end)));
                %y_vec.vector(i,m+1)=real(xi_vec.vector(i,m+1));
                y_vec.vector(i,:)=(xi_vec.vector(i,:)+conj(flip(xi_vec.vector(i,:),2)))/2;
            end
            y_vec.bool_ifft = 0;
        end
        % symmetrise
        
        
        
        % CREAT_XI_VECTOR (CREATION from FFT)
        function xi=create_Xi_vector(lambda,x)
            % creation of a Xi_vector starting from a scalar vector lambda consisting
            % of all the scalar paramenters and a matrix x of the values of the
            % solution depending on time t
            % x such be such that
            %   x(t,:) is the solution of the ODE at time t
            %
            % we assume than min(size(x)) is the number of vector unknowns, this means
            % that the number of vetor unknows is less than the size of the simulation
            
            
            xi=Xi_vector;
            
            if min(size(lambda))~=1
                error('Dimension of lambda is wrong')
            end
            xi.size_scalar=max(size(lambda));
            xi.scalar=lambda;
            
            xi.size_vector=size(min(x));
            
            xi.vector=verifyfft(x); % in order to apply FFT here, the size of x must be
            xi.bool_ifft = 0;
            % an exact power of 2
        end
        % CREAT_XI_VECTOR (CREATION from FFT)
        
        % SUM_XI_VECTOR
        function z_vec=sum_Xi_vector(x_vec,y_vec,flag)
            % z_vec=x_vec+y_vec
            % flag = 1 (DEFAULT)  the output is of the size of the smallest
            % vector, othewise maximum size
            % NOTE: the two vectors need to have the same size_scalar and
            % size_vector
            
            if x_vec.size_vector~=y_vec.size_vector
                error('sizes of inputs do not match')
            end
            if x_vec.size_scalar~=y_vec.size_scalar
                error('sizes of inputs do not match')
            end
            if nargin<3
                flag=1;
            end
            
            if x_vec.nodes>y_vec.nodes
                max_vec=equal_Xi(x_vec);
                min_vec=equal_Xi(y_vec);
            else
                max_vec=equal_Xi(y_vec);
                min_vec=equal_Xi(x_vec);
            end
            if flag
                z_vec=Xi_vector(min_vec);
                max_vec=reshape_Xi(max_vec,min_vec.nodes);
            else
                z_vec=Xi_vector(max_vec);
                min_vec=reshape_Xi(min_vec,max_vec.nodes);
            end
            if size(min_vec.scalar)~=size(max_vec.scalar)
                if size(min_vec.scalar.')==size(max_vec.scalar)
                    min_vec.scalar=min_vec.scalar.';
                else
                    error('Scalar sizes not compatible')
                end
            end
            z_vec.scalar=min_vec.scalar+max_vec.scalar;
            z_vec.vector=min_vec.vector+max_vec.vector;
            z_vec.bool_ifft = 0;
        end
        % SUM_XI_VECTOR
        
        % PLUS
        function z_vec=plus(x_vec,y_vec)
            % z_vec=x_vec-y_vec
            % flag = 1 (DEFAULT)  the output is of the size of the smallest
            % vector, othewise maximum size
            % NOTE: the two vectors need to have the same size_scalar and
            % size_vector
            
            z_vec=sum_Xi_vector(x_vec,y_vec);
            if x_vec.nodes~=y_vec.nodes
                warning('Number of nodes do not coincide, chosen the smallest');
            end
            z_vec.bool_ifft = 0;
        end
        % PLUS
        
        % EQ
        function bool=eq(x_vec,y_vec)
            % bool= (x_vec == y_vec)
            bool=0;
            if x_vec.nodes~=y_vec.nodes
                return
            end
            if x_vec.size_scalar~=y_vec.size_scalar
                return
            end
            if x_vec.size_vector ~= y_vec.size_vector
                return
            end
            if x_vec.scalar ~= y_vec.scalar
                return
            end
            if x_vec.vector ~= y_vec.vector
                return
            end
            bool=1;
            return
        end
        % EQ
        
        % NE
        function bool=ne(x_vec,y_vec)
            equality=eq(x_vec,y_vec);
            if equality==0
                bool=1;
            elseif equality==1
                bool=0;
            else
                error('wrong return eq');
            end
        end
        % NE
        
        % DIF_XI_VECTOR
        function z_vec=dif_Xi_vector(x_vec,y_vec,flag)
            % z_vec=x_vec-y_vec
            % flag = 1 (DEFAULT)  the output is of the size of the smallest
            % vector, othewise maximum size
            % NOTE: the two vectors need to have the same size_scalar and
            % size_vector
            
            if x_vec.size_vector~=y_vec.size_vector
                error('sizes of inputs do not match')
            end
            if x_vec.size_scalar~=y_vec.size_scalar
                error('sizes of inputs do not match')
            end
            if nargin<3
                flag=1;
            end
            
            if x_vec.nodes>y_vec.nodes
                max_vec=equal_Xi(x_vec);
                min_vec=equal_Xi(y_vec);
            else
                max_vec=equal_Xi(y_vec);
                min_vec=equal_Xi(x_vec);
            end
            if flag
                z_vec=Xi_vector(min_vec);
                max_vec=reshape_Xi(max_vec,min_vec.nodes);
            else
                z_vec=Xi_vector(max_vec);
                min_vec=reshape_Xi(min_vec,max_vec.nodes);
            end
            
            z_vec.scalar=min_vec.scalar-max_vec.scalar;
            z_vec.vector=min_vec.vector-max_vec.vector;
            z_vec.bool_ifft = 0;
        end
        % DIF_XI_VECTOR
        
        % MINUS
        function z_vec=minus(x_vec,y_vec)
            % z_vec=x_vec-y_vec
            % flag = 1 (DEFAULT)  the output is of the size of the smallest
            % vector, othewise maximum size
            % NOTE: the two vectors need to have the same size_scalar and
            % size_vector
            
            z_vec=dif_Xi_vector(x_vec,y_vec);
            if x_vec.nodes~=y_vec.nodes
                warning('Number of nodes do not coincide, chosen the smallest');
            end
            z_vec.bool_ifft = 0;
        end
        % MINUS
        
        % MTIMES
        function z_vec=mtimes(scal,x_vec)
            % function z_vec=mtimes(scal,x_vec)
            if any((size(scal))>1) || length(size(scal))>2
                error('just scalar multiplication is possible');
            end
            z_vec=Xi_vector(x_vec);
            z_vec.scalar=scal*x_vec.scalar;
            z_vec.vector=scal*x_vec.vector;
            z_vec.bool_ifft = 0;
        end
        % MTIMES
        
        % MRDIVIDE
        function z_vec = mrdivide(x_vec, scal)
            % function z_vec=mrdivide(x_vec, scal)
            if any((size(scal))>1) || length(size(scal))>2
                error('just scalar multiplication is possible');
            end
            z_vec=Xi_vector(x_vec);
            z_vec.scalar=scal/x_vec.scalar;
            z_vec.vector=scal/x_vec.vector;
            z_vec.bool_ifft = 0;
        end
        % MRDIVIDE
        
        function Xi=reshape(x,new_nodes)
            Xi=reshape_Xi(x,new_nodes);
        end
            
        % RESHAPE_XI
        function Xi=reshape_Xi(x,new_nodes)
            % x          Xi_vector
            % new_nodes  number of nodes used in the output
            % Xi         Xi_vector with new_nodes number of nodes
            %            if x.nodes>new_nodes, the vectors are cutted,
            %            otherwise they are padded with zeros
            
            %global use_intlab
            
            if x.nodes==new_nodes
                Xi=x;
                return
            elseif x.nodes<new_nodes
                dif_nodes=new_nodes-x.nodes;
                Xi=Xi_vector(x);
                Xi.scalar=x.scalar;
                Xi.nodes=new_nodes;
                if isa(x.vector,'intval')
                    Xi.vector=[intval(zeros(x.size_vector,dif_nodes)),x.vector,intval(zeros(x.size_vector,dif_nodes))];
                else
                    Xi.vector=[zeros(x.size_vector,dif_nodes),x.vector,zeros(x.size_vector,dif_nodes)];
                end
                Xi.bool_ifft = 0;
            else
                dif_nodes=x.nodes-new_nodes;
                Xi=Xi_vector(x);
                Xi.scalar=x.scalar;
                Xi.nodes=new_nodes;
                Xi.vector=x.vector(:,dif_nodes+1:end-dif_nodes);
                Xi.bool_ifft = 0;
            end
            
        end
        % RESHAPE_XI
        
        
        % XI_VEC2VEC_COMPATIBLE
        function bool = Xi_vec2vec_compatible(x_Xi, x_vec)
            % function bool = Xi_vec2vec_compatible(x_Xi, x_vec)
            %
            % This function tests the compatibility between a vector and a Xi_vector
            % This function returns a bool that is true if the two elements are
            % compatible and false otherwise. The two elements are compatible if the
            % length of the vector is equal to the lenght of Xi_vec2vec(x_Xi), or
            % considered otherwise x_Xi.size_scalar+x_Xi.size_vector*(x_Xi.nodes*2+1)
            %
            % INPUTS
            % x_Xi          Xi_vector
            % x_vec         standard vector
            %
            % OUTPUT
            % bool          check if coherent
            bool=false;
            if  x_Xi.size_scalar+x_Xi.size_vector*(x_Xi.nodes*2+1) == max(size(x_vec)) && min(size(x_vec))==1
                bool=true;
            end
        end
        % XI_VEC2VEC_COMPATIBLE
        
        
        % FINITE
        function Xi=finite(x,F_nodes)
            % x         Xi_vector
            % F_nodes   integer, F_nodes<x.nodes
            % Xi        Xi_vector, finite part of x.
            %   Xi_k=0  if k>F_nodes
            Xi=equal_Xi(x);
            %Xi.scalar(:)=0;
            if F_nodes<x.nodes
                Xi.vector(:,1:x.nodes-F_nodes)=0;
                Xi.vector(:,2*x.nodes-F_nodes+3:end)=0;
            end
            Xi.bool_ifft = 0;
        end
        % FINITE
        
        % TAIL
        function Xi=tail(x,F_nodes)
            % x         Xi_vector
            % F_nodes   integer, F_nodes<x.nodes
            % Xi        Xi_vector, tail part of x.
            %   Xi_k=0  if k>F_nodes
            global use_intlab
            Xi=equal_Xi(x);
            Xi.scalar(:)=0;
            if F_nodes<x.nodes
                if use_intlab
                    Xi.vector(:,x.nodes-F_nodes+1:x.nodes+F_nodes+1)=intval(0);
                else
                    Xi.vector(:,x.nodes-F_nodes+1:x.nodes+F_nodes+1)=0;
                end
            end
            Xi.bool_ifft = 0;
        end
        % TAIL
        
        % TILDE_XI
        function Xi=tilde_Xi(x)
            % x   Xi_vector
            % Xi  Xi_vector, such that in the vector part
            %   Xi_vector_k=k*x_k
            global use_intlab
            Xi=equal_Xi(x);
            K=-x.nodes:x.nodes;
            if use_intlab
                K=intval(K);
            end
            for i=1:x.size_vector
                Xi.vector(i,:)=x.vector(i,:).*K;
            end
            Xi.bool_ifft = 0;
        end
        % TILDE_XI
        
        
        % SPECIAL FUNCTIONS FOR DELAYS
        
        % KNORM
        function nor = Knorm(xi)
            nor = norm(tilde_Xi(xi));
        end
        % KNORM
        
        % DELAY
        function Xi = delay(x, tau)
            global use_intlab
            if length(tau)~=x.size_vector
                error('Dimensions not matching')
            end
            omega = x.scalar(1);
            K=-x.nodes:x.nodes;
            if use_intlab
                K=intval(K);
            end
            Xi = x;
            Xi.scalar = 0*x.scalar;
            const_func = 0*K;
            const_func(x.nodes+1) = 1;
            for i=1:x.size_vector
                if tau(i) == 0 
                    Xi.vector(i,:) = const_func;
                else
                    Xi.vector(i,:)=(exp(1i*K*tau(i)*omega) .* x.vector(i,:));
                end
            end
            Xi.bool_ifft = 0;
            Xi.ifft_vector = [];
        end
        % DELAY
        % end special functions for delay
        
        %SPROD_XI
        function Xi=sprod_Xi(x,alpha)
            %function Xi=sprod_Xi(x,alpha)
            % scalar product of a Xi_vector
            % INPUT
            % x     Xi_vector
            % alpha double
            % OUTPUT
            % Xi    Xi_vector  Xi=x*alpha
            Xi=equal_Xi(x);
            if ~isequal(size(alpha),[1,1])
                error('the second input must be a scalar')
            end
            Xi.vector=x.vector.*alpha;
            Xi.scalar=x.scalar.*alpha;
            Xi.bool_ifft = 0;
        end % SPROD_XI
        
        %INTVAL
        function y_Xi=intval(x_Xi)
            %turn x_Xi into an intval Xi_vector
            y_Xi=Xi_vector(x_Xi);
            y_Xi.scalar=intval(x_Xi.scalar);
            y_Xi.vector=intval(x_Xi.vector);
            y_Xi.bool_ifft = 0;
        end
        %INTVAL
        
        % INTERVAL_XI
        function xBarS=interval_Xi(xBar0,xBar1)
            % function xBarS=interval_Xi(xBar0,xBar1)
            %
            % taking as input two Xi_vectors, it construct an use_intlab Xi_vector having
            % the two inputs as borders
            
            %global use_intlab
            
            if isintval(xBar0.scalar)
                xBarS=Xi_vector(xBar0);
                scal_max=max(sup(real(xBar0.scalar)),sup(real(xBar1.scalar)))+...
                    1i*max(sup(imag(xBar0.scalar)),sup(imag(xBar1.scalar)));
                scal_min=min(inf(real(xBar0.scalar)),inf(real(xBar1.scalar)))+...
                    1i*min(inf(imag(xBar0.scalar)),inf(imag(xBar1.scalar)));
                xBarS.scalar=infsup(scal_min,scal_max);
                
                vec_max=max(sup(real(xBar0.vector)),sup(real(xBar1.vector)))+...
                    1i*max(sup(imag(xBar0.vector)),sup(imag(xBar1.vector)));
                vec_min=min(inf(real(xBar0.vector)),inf(real(xBar1.vector)))+...
                    1i*min(inf(imag(xBar0.vector)),inf(imag(xBar1.vector)));
                xBarS.vector=infsup(vec_min,vec_max);
            else
                xBarS=Xi_vector(xBar0);
                scal_max=max((real(xBar0.scalar)),(real(xBar1.scalar)))+...
                    1i*max((imag(xBar0.scalar)),(imag(xBar1.scalar)));
                scal_min=min((real(xBar0.scalar)),(real(xBar1.scalar)))+...
                    1i*min((imag(xBar0.scalar)),(imag(xBar1.scalar)));
                xBarS.scalar=infsup(scal_min,scal_max);
                
                vec_max=max((real(xBar0.vector)),(real(xBar1.vector)))+...
                    1i*max((imag(xBar0.vector)),(imag(xBar1.vector)));
                vec_min=min((real(xBar0.vector)),(real(xBar1.vector)))+...
                    1i*min((imag(xBar0.vector)),(imag(xBar1.vector)));
                xBarS.vector=infsup(vec_min,vec_max);
                %                 xBarS=Xi_vector(xBar0);
                %                 scal_max=max(xBar0.scalar,xBar1.scalar);
                %                 scal_min=min(xBar0.scalar,xBar1.scalar);
                %                 xBarS.scalar=infsup(scal_min,scal_max);
                %
                %                 vec_max=max(xBar0.vector,xBar1.vector);
                %                 vec_min=min(xBar0.vector,xBar1.vector);
                %                 xBarS.vector=infsup(vec_min,vec_max);
                
            end
            xBarS.bool_ifft = 0;
        end %INTERVAL_XI
        
        %ABS
        function y_Xi=abs(x_Xi)
            %turn x_Xi into an abs Xi_vector
            y_Xi=Xi_vector(x_Xi);
            y_Xi.scalar=abs(x_Xi.scalar);
            y_Xi.vector=abs(x_Xi.vector);
            y_Xi.bool_ifft = 0;
        end
        %ABS
        
        %LENGTH
        function l = length(x_Xi)
            l = x_Xi.size_scalar+x_Xi.size_vector*(2*x_Xi.nodes+1);
        end
        %LENGTH
        
        %RDIVIDE
        function y_Xi=rdivide(x_Xi,a)
            %turn x_Xi into an abs Xi_vector
            y_Xi=Xi_vector(x_Xi);
            y_Xi.scalar=rdivide(x_Xi.scalar,a);
            y_Xi.vector=rdivide(x_Xi.vector,a);
            y_Xi.bool_ifft = 0;
        end
        %RDIVIDE
        
        %MAX
        function z_Xi=max(x_Xi,y_Xi)
            %turn x_Xi into an intval Xi_vector
            global use_intlab
            
            if isfloat(y_Xi)
                z_Xi=max(max(norm(x_Xi)),y_Xi);
                return
            end
            
            z_Xi=Xi_vector(x_Xi);
            if x_Xi.size_scalar~=y_Xi.size_scalar || x_Xi.size_vector~=y_Xi.size_vector
                error('Sizes do not match')
            elseif x_Xi.nodes~=y_Xi.nodes
                if x_Xi.nodes>y_Xi.nodes
                    y_Xi=reshape_Xi(y_Xi,x_Xi.nodes);
                else
                    x_Xi=reshape_Xi(x_Xi,y_Xi.nodes);
                end
            end
            if isintval(x_Xi.scalar)
                if ~isintval(y_Xi.scalar)
                    y_Xi = intval(y_Xi);
                end
                z_Xi.scalar=max(x_Xi.scalar.sup,y_Xi.scalar.sup);
                z_Xi.vector=max(x_Xi.vector.sup,y_Xi.vector.sup);
            else
                if isintval(y_Xi.scalar)
                    x_Xi = intval(x_Xi);
                    z_Xi.scalar=max(x_Xi.scalar.sup,y_Xi.scalar.sup);
                    z_Xi.vector=max(x_Xi.vector.sup,y_Xi.vector.sup);
                else
                    z_Xi.scalar=max(x_Xi.scalar,y_Xi.scalar);
                    z_Xi.vector=max(x_Xi.vector,y_Xi.vector);
                end
            end
            z_Xi.bool_ifft = 0;
        end
        %MAX
        
        % TIME_SERIES
        function xt = time_series(X)
            if abs(X.scalar)>0.1
                npoints=200*ceil(1/real(X.scalar(1)));
                t1=linspace(0,1/X.scalar(1),npoints);
            else
                npoints=200*ceil(1/0.1);
                t1=linspace(0,1/0.1,npoints);
            end
            xt=zeros(npoints,X.size_vector);
            m=X.nodes;
            for i=1:npoints
                t=t1(i);
                for j=1:X.size_vector
                    xt(i,j)=sum(X.vector(j,:).*exp(1i*[-m:m]*2*pi*t*X.scalar(1)));
                end
            end
        end
        % TIME_SERIES
        
        
        % REDUCE_DIM
        function y = reduce_dim(x,new_dim,other_dim)
            if nargin <3
                first_dim = 1;
                last_dim = new_dim;
            else
                first_dim = min(new_dim, other_dim);
                last_dim = max(new_dim, other_dim);
            end
            if new_dim>x.size_vector
                error('Dimensions are reducable just to smaller phase planes')
            end
            y_vector = x.vector(first_dim:last_dim, :);
            y = Xi_vector(x.scalar,y_vector);
        end
        % end REDUCE_DIM
        
        
        % DISP
        function disp(x)
            fprintf('Scalar size (M): %d,',x.size_scalar)
            fprintf('Vector size (N): %d,',x.size_vector)
            fprintf('Nodes (F): %d.\n',x.nodes)
            fprintf('Scalar:\n')
            disp(x.scalar)
            disp('')
            fprintf('Vector:\n')
            disp(x.vector.')
            %fprintf('Scalar size (M): %d,\n Vector size (N): %d,\n Scalar: %f,\n Vector: %f\n',...
            %    x.size_scalar,x.size_vector,x.scalar,x.vector)
        end % DISP
        
        %PLOT
        function plot(X,varargin)
            % default plotting option for Xi_vector
            if abs(X.scalar(1))>0.1
                npoints=200*ceil(abs(1/real(X.scalar(1))));
                t1=linspace(0,1/abs(X.scalar(1)),npoints);
            else
                npoints=200*ceil(1/0.1);
                t1=linspace(0,1/0.1,npoints);
            end
            xt1=zeros(npoints,X.size_vector);
            %xt2=zeros(200,1);
            m=X.nodes;
            for i=1:npoints
                t=t1(i);
                for j=1:X.size_vector
                    xt1(i,j)=sum(X.vector(j,:).*exp(1i*[-m:m]*2*pi*t*X.scalar(1)));
                end
                %xt2(i)=sum(xBar{3}.*exp(1i*[-m:m]'*2*pi*t*xBar{1}));
            end
            
            t_plot = linspace(0, 2*pi * abs(X.scalar(1)),npoints);
            setround(0)
            
            %scrsz = get(0,'ScreenSize');
            %figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
            for j=1:X.size_vector
                subplot(1,X.size_vector,j);
                plot(t_plot,real(xt1(:,j)),varargin{:})
                xlabel('t')
                ylabel(sprintf('component %d',j));
                %xlim([0,2*pi * abs(X.scalar(1))])
                %ylim([1.2*min(real(xt1(:,j))),1.2*max(real(xt1(:,j)))]);
                
                %xTick_float = linspace(0,2*pi * abs(X.scalar(1)),...
                %    max(2,ceil(8/X.size_vector)));
                %yTick_float = sort([linspace(min(real(xt1(:,j))),...
                %    max(real(xt1(:,j))),4),0]);
                %set(gca,'XTick',xTick_float)
                %set(gca,'XTickLabel',sprintf('%5.3f|',xTick_float))
                %set(gca,'YTick',yTick_float)
                %set(gca,'YTickLabel',sprintf('%1.2f|',yTick_float))
            end
        end
        %PLOT
        
        %PLOT2
        function plot2(X,varargin)
            % 2D plotting option for Xi_vector
            if X.size_vector ~= 2
                error('This function is suppose to bejust for 2D vectors');
            end
            if abs(X.scalar(1))>0.1
                npoints=200*ceil(1/real(X.scalar(1)));
                t1=linspace(0,1/X.scalar(1),npoints);
            else
                npoints=200*ceil(1/0.1);
                t1=linspace(0,1/0.1,npoints);
            end
            xt1=zeros(npoints,X.size_vector);
            %xt2=zeros(200,1);
            m=X.nodes;
            for i=1:npoints
                t=t1(i);
                for j=1:X.size_vector
                    xt1(i,j)=sum(X.vector(j,:).*exp(1i*[-m:m]*2*pi*t*X.scalar(1)));
                end
                %xt2(i)=sum(xBar{3}.*exp(1i*[-m:m]'*2*pi*t*xBar{1}));
            end
            
            %scrsz = get(0,'ScreenSize');
            %figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
            
            plot(real(xt1(:,1)),real(xt1(:,2)),varargin{:})
            xlabel('x','Interpreter','Latex','FontSize',20)
            ylabel('y','Interpreter','Latex','FontSize',20)
            xlim([-0.6/abs(X.scalar(1)),0.6/abs(X.scalar(1))])
            axis equal
            if max(abs(imag(xt1)))>max(abs(real(xt1)))
                warning('The solution plotted has an iportant imaginary part')
            end
        end
        %PLOT2
        
        %PLOT3
        function plot3(X,varargin)
            % default plotting option for Xi_vector
            
            if X.size_vector~=3
                if X.size_vector<=2
                    warning('the vector has just two components, the standard plot is used');
                    plot(X)
                else
                    error('Wrong number of components')
                end
                return
            end
            npoints=2000*ceil(1/real(X.scalar(1)));
            t1=linspace(0,1/X.scalar(1),npoints);
            xt1=zeros(npoints,X.size_vector);
            m=X.nodes;
            for i=1:npoints
                t=t1(i);
                for j=1:X.size_vector
                    xt1(i,j)=sum(X.vector(j,:).*exp(1i*[-m:m]*2*pi*t*X.scalar(1)));
                end
            end
            
            %scrsz = get(0,'ScreenSize');
            %figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
            plot3(real(xt1(:,1)),real(xt1(:,2)),real(xt1(:,3)),varargin{:})
        end
        %PLOT3
        
        
    end
end