classdef PDE_vector
    properties
        size_real
        size_vector
        node_time
        node_space
        
        parameters % size_real
        vector1D % size_vector x ( 2*node_space + 1 )
        vector2D % size_vector x ( 2*node_time + 1 ) x ( 2*node_space + 1 )
    end
    methods
        % CONSTRUCTOR
        function x=PDE_vector(scalars, vector_1D, vector_2D, nodes_space)
            %
            % PDE_vector constructor of the class PDE_vector
            %
            % x=PDE_vector(scalars, vector_1D, vector_2D)
            % INPUT
            % scalars       vector or scalar parameters
            % vector_1D     vector of 1D Fourier coefficients
            % vector_2D     vector of 2D Fourier coefficients
            %
            % x = PDE_vector(vector, size_parameters, size_vector, nodes_space)
            % INPUT
            % vector        long vector with real parameters, Fourier
            %                   coefficients 1D and 2D
            % size_parameters size of the real paramters
            % size_vector   number of Fourier series in 1D ( and 2D)
            % nodes_space   number of nodes in space direction
            % x = PDE_vector(vector, x_PDE_vector)
            % INPUT
            % vector        huge vector with real parameters, Fourier
            %                   coefficients 1D and 2D
            % x_PDE_vector  PDE_vector of corresponding dimension
            
            if nargin == 2
                if ~isa(vector_1D,'PDE_vector')
                    error('Wrong inputs')
                else
                    vector = scalars;
                    PDE_vec = vector_1D;
                end
                x = PDE_vector(vector,PDE_vec.size_real, PDE_vec.size_vector,PDE_vec.node_space);
                return
            end
            if nargin==3
                if min(size(scalars))~=1
                    error('The parameters must be in a vector notation')
                end
                size_Vec1 = size(vector_1D);
                size_Vec2 = size(vector_2D);
                if size_Vec1(1)~= size_Vec2(1) || (length(size_Vec2)>2 &&size_Vec1(2)~= size_Vec2(3))
                    error('The space vector and time vector must agree')
                end
                %                 if size(vector_1D,2)==1
                %                     vector_1D=vector_1D.';
                %                 end
                %                 if min(size(vector_2D))==1
                %                     error('The time and space vector must be in a matrix notation')
                %                 end
                if size(vector_1D,2)~=size(vector_2D,3) || size(vector_1D,1)~=size(vector_2D,1)
                    error('Vector sizes incompatible')
                end
                x.size_vector = size(vector_1D,1);
                x.parameters = scalars;
                x.size_real = length(scalars);
                x.vector1D = vector_1D;
                x.node_space = (length(vector_1D)-1)/2;
                x.vector2D = vector_2D;
                x.node_time = (size(vector_2D,2)-1)/2;
                
            elseif nargin ==4
                
                x.size_real = vector_1D;
                x.size_vector = vector_2D;
                x.node_space = nodes_space;
                vector = scalars;
                x.parameters = vector(1:x.size_real);
                vector(1:x.size_real)=[];
                x.vector1D = vector(1:x.size_vector*(2*x.node_space+1));
                x.vector1D = reshape(x.vector1D,x.size_vector,[]);
                vector(1:x.size_vector*(2*x.node_space+1))=[];
                x.vector2D = vector;
                if isa(x.vector2D,'intval')
                    x_vector2D_inf = inf(x.vector2D);
                    x_vector2D_sup = sup(x.vector2D);
                    x_vector2D_inf = reshape(x_vector2D_inf, x.size_vector , [], ( 2*x.node_space + 1 ));
                    x_vector2D_sup = reshape(x_vector2D_sup, x.size_vector , [], ( 2*x.node_space + 1 ));
                    x.vector2D = infsup(x_vector2D_inf,x_vector2D_sup);
                else
                    x.vector2D = reshape(x.vector2D, x.size_vector , [], ( 2*x.node_space + 1 ));
                end
                x.node_time = (size(x.vector2D, 2) -1)/2;
                
            else
                error('This input number is not recognised')
            end
        end
        %CONSTRUCTOR
        
        % NORM
        function nor= norm(xi)
            % function nor= norm(xi)
            global nu
            nor=cnorm_PDE(xi,nu);
        end
        % NORM
        
        % ABS
        function abs_x = abs(x)
            abs_x = x;
            abs_x.parameters = abs(x.parameters);
            abs_x.vector1D = abs(x.vector1D);
            abs_x.vector2D = abs(x.vector2D);
            
        end
        % ABS
        
        % CNORM_PDE
        function nor=cnorm_PDE(xi,nu)
            % cnorm_Xi_vector(Xi_vector,double)
            % in this function, the component-wise norm of a Xi_vector is
            % taken, that means that the output will be a positive double array
            % of size size_scalar+size_vector, storing in the first part the
            % absolute value of the Xi_vector, in the second part the l1-norm
            % of the vectors of the Fourier coefficients. This routine is
            % called by norm_Xi_vector, that returns the maximum of
            % cnorm_Xi_vector
            
            global use_intlab
            
            if use_intlab
                if ~isintval(nu)
                    nu=intval(nu);
                end
            else
                if ~isfloat(nu)
                    nu=sup(nu);
                end
            end
            
            if length(nu) ==1
                nu1 = nu(1);
                nu2 = nu(1);
            elseif length(nu) ==2
                nu1 = nu(1);
                nu2 = nu(2);
            else
                error('Length of weight vector nu must be 1 or 2')
            end
            
            nor=zeros(xi.size_real+2*xi.size_vector,1);
            if ~isa(xi.parameters(1),'intval')
                nor(1:xi.size_real)=abs(xi.parameters);
            else
                nor(1:xi.size_real)=sup(abs(xi.parameters));
            end
            %nor(xi.size_scalar+1:end)=abs(xi.vector(:,xi.nodes+1));
            K2=-xi.node_space:xi.node_space;
            K1=-xi.node_time:xi.node_time;
            K1 = K1.';
            
            nuK1 = nu1.^ abs(K1);
            nuK2 = nu2.^ abs(K2);
            
            nuK_1and2_temp = repmat(nuK1,1,1+2*xi.node_space).*repmat(nuK2,1+2*xi.node_time,1);
            nuK_1and2_temp = shiftdim(nuK_1and2_temp,-1);
            nuK_1and2 = repmat(nuK_1and2_temp, [1, 1,1]);
            if use_intlab
                nuK2=intval(nuK2);
                nuK_1and2=intval(nuK_1and2);
                
            end
            
            
            for i=1:xi.size_vector
                if ~isa(xi.vector1D(1,1),'intval')
                    if xi.size_vector ==1
                        nor(xi.size_real+i)=sum(abs(xi.vector1D(:).').*(nuK2));
                        nor(xi.size_real+xi.size_vector+i)=sum(sum(abs(xi.vector2D).*(nuK_1and2)));
                    else
                        nor(xi.size_real+i)=sum(abs(xi.vector1D(i,:)).*(nuK2));
                        nor(xi.size_real+xi.size_vector+i)=sum(sum(abs(xi.vector2D(i,:,:)).*(nuK_1and2)));
                    end
                else
                    if xi.size_vector ==1
                        nor(xi.size_real+i)=sup(sum(abs(xi.vector1D(:).').*(nuK2)));
                        nor(xi.size_real+xi.size_vector+i)=sup(sum(sum(squeeze(abs(xi.vector2D).*(nuK_1and2)))));
                    else
                        nor(xi.size_real+i)=sup(sum(abs(xi.vector1D(i,:)).*(nuK2)));
                        nor(xi.size_real+xi.size_vector+i)=sum(sum(sup(abs(xi.vector2D(i,:,:)).*(nuK_1and2))));
                    end
                end
            end
        end
        %CNORM_PDE
        
        
        % NORM_PDE
        function nor=norm_PDE(xi)
            % norm_Xi_vector(Xi_vector,double)
            % return the maximum-norm of the Xi_vector xi
            global nu
            nor=max(cnorm_PDE(xi,nu));
        end
        % NORM_PDE
        
        % symmetrise
        function x_sym = symmetrise(x)
            x_sym = x;
            
            x_sym.parameters = (x.parameters + conj(x.parameters))/2;
            
            y_sym = (x.vector1D(:,:) + conj(x.vector1D(:,end:-1:1)))/2;
            
            z_sym = (x.vector2D(:,:,:) + conj(x.vector2D(:,end:-1:1,end:-1:1)))/2;
            
            x_sym.vector1D = y_sym;
            
            x_sym.vector2D = z_sym;
            % %function y_vec=symmetrised(xi_vec)
            % y_vec=xi_vec;
            % y_vec.parameters=(xi_vec.parameters+conj(xi_vec.parameters))/2;
            % %m=xi_vec.nodes;
            % for i=1:xi_vec.size_vector
            %     y_vec.vector1D(i,:)=(xi_vec.vector1D(i,:)+conj(flip(xi_vec.vector1D(i,:),2)))/2;
            %     y_vec.vector2D(i,:,:)=(xi_vec.vector2D(i,:,:)+conj(flip(xi_vec.vector2D(i,:,:),3))+...
            %         conj(flip(xi_vec.vector2D(i,:,:),2))+flip(flip(xi_vec.vector2D(i,:,:),3),2))/4;
            % end
            % %y_vec.bool_ifft = 0;
        end
        % symmetrise
        
        
        % PLUS
        function z_vec=plus(x_vec,y_vec)
            % z_vec=x_vec+y_vec
            % flag = 1 (DEFAULT)  the output is of the size of the smallest
            % vector, othewise maximum size
            % NOTE: the two vectors need to have the same size_scalar and
            % size_vector
            
            if x_vec.size_vector~=y_vec.size_vector
                error('sizes of inputs do not match')
            end
            if x_vec.size_real~=y_vec.size_real
                error('sizes of inputs do not match')
            end
            
            if x_vec.node_space>y_vec.node_space
                max_node_space=x_vec.node_space;
            else
                max_node_space=x_vec.node_space;
            end
            if x_vec.node_time>y_vec.node_time
                max_node_time=x_vec.node_time;
            else
                max_node_time=x_vec.node_time;
            end
            x_vec = reshape(x_vec,max_node_space,max_node_time);
            y_vec = reshape(y_vec,max_node_space,max_node_time);
            z_vec = x_vec;
            z_vec.parameters=x_vec.parameters+y_vec.parameters;
            z_vec.vector1D=x_vec.vector1D+y_vec.vector1D;
            z_vec.vector2D=x_vec.vector2D+y_vec.vector2D;
        end
        % PLUS
        
        % EQ
        function bool=eq(x_vec,y_vec)
            % bool= (x_vec == y_vec)
            bool=0;
            if x_vec.size_real~=y_vec.size_real
                return
            end
            if x_vec.size_vector~=y_vec.size_vector
                return
            end
            if x_vec.node_time ~= y_vec.node_time
                return
            end
            if x_vec.node_space ~= y_vec.node_space
                return
            end
            if x_vec.parameters ~= y_vec.parameters
                return
            end
            if x_vec.vector1D ~= y_vec.vector1D
                return
            end
            if x_vec.vector2D ~= y_vec.vector2D
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
        
        % MINUS
        function z_vec=minus(x_vec,y_vec)
            % z_vec=x_vec+y_vec
            % flag = 1 (DEFAULT)  the output is of the size of the smallest
            % vector, othewise maximum size
            % NOTE: the two vectors need to have the same size_scalar and
            % size_vector
            
            if x_vec.size_vector~=y_vec.size_vector
                error('sizes of inputs do not match')
            end
            if x_vec.size_real~=y_vec.size_real
                error('sizes of inputs do not match')
            end
            
            if x_vec.node_space>y_vec.node_space
                max_node_space=x_vec.node_space;
            else
                max_node_space=x_vec.node_space;
            end
            if x_vec.node_time>y_vec.node_time
                max_node_time=x_vec.node_time;
            else
                max_node_time=x_vec.node_time;
            end
            x_vec = reshape(x_vec,max_node_space,max_node_time);
            y_vec = reshape(y_vec,max_node_space,max_node_time);
            z_vec = x_vec;
            z_vec.parameters=x_vec.parameters-y_vec.parameters;
            z_vec.vector1D=x_vec.vector1D-y_vec.vector1D;
            z_vec.vector2D=x_vec.vector2D-y_vec.vector2D;
        end
        % MINUS
        
        % MAX
        function z=max(x,y)
            % function z=max(x,y)
            %
            % elementwise maximum for PDE
            global use_intlab
            
            if ~isa(y,'PDE_vector')
                error('imposssible comparison')
            end
            if ~PDE2PDE_compatible(x,y)
                error('Incompatible PDE vectors to maximize')
            end
            if PDE2PDE_compatible(x,y)==1
                warning('PDE vectors with different nodes, biggest chosen')
                nodes_space = max(x.node_space,y.node_space);
                nodes_time = max(x.node_time, y.node_time);
                x = reshape(x, nodes_space, nodes_time);
                y = reshape(y, nodes_space, nodes_time);
            end
            z=x;
            
            if use_intlab
                xscal=abs(x.parameters);
                yscal=abs(y.parameters);
                z.parameters=max(xscal.sup,yscal.sup);
                xvec=abs(y.vector1D);
                yvec=abs(x.vector1D);
                z.vector1D=max(xvec.sup,yvec.sup);
                xvec=abs(y.vector2D);
                yvec=abs(x.vector2D);
                z.vector2D=max(xvec.sup,yvec.sup);
            else
                z.parameters = max(x.parameters,y.parameters);
                z.vector1D = max(x.vector1D , y.vector1D);
                z.vector2D = max(x.vector2D , y.vector2D);
            end
        end
        % MAX
        
        % MIN
        function z=min(x,y)
            % function z=max(x,y)
            %
            % elementwise minimum for PDE
            global use_intlab
            
            if ~isPDEvec(y)
                error('imposssible comparison')
            end
            if ~PDE2PDE_compatible(x,y)
                error('Incompatible PDE vectors to maximize')
            end
            if PDE2PDE_compatible(x,y)==1
                warning('PDE vectors with different nodes, biggest chosen')
                nodes_space = min(x.node_space,y.node_space);
                nodes_time = min(x.node_time, y.node_time);
                x = reshape(x, nodes_space, nodes_time);
                y = reshape(y, nodes_space, nodes_time);
            end
            z=x;
            
            if use_intlab
                xscal=abs(x.parameters);
                yscal=abs(y.parameters);
                z.parameters=min(xscal.sup,yscal.sup);
                xvec=abs(y.vector1D);
                yvec=abs(x.vector1D);
                z.vector1D=min(xvec.sup,yvec.sup);
                xvec=abs(y.vector2D);
                yvec=abs(x.vector2D);
                z.vector2D=min(xvec.sup,yvec.sup);
            else
                z.parameters = min(x_vec.parameters,y_vec.parameters);
                z.vector1D = min(x_vec.vector1D , y_vec.vector1D);
                z.vector2D = min(x_vec.vector2D , y_vec.vector2D);
            end
        end
        % MIN
        
        
        % MTIMES (both scalar and PDE_vector)
        function z_vec=mtimes(x_vec,scal)
            % function z_vec=mtimes(scal,x_vec)
            if ~isa(x_vec, 'PDE_vector')
                temp = x_vec;
                x_vec = scal;
                scal = temp;
            end
            if isa(scal,'double') || isa(scal,'intval')
                if length(scal) ==1
                z_vec=x_vec;
                z_vec.parameters=scal*x_vec.parameters;
                z_vec.vector1D=scal*x_vec.vector1D;
                z_vec.vector2D=scal*x_vec.vector2D;
                else
                    z_vec = PDE_vector(scal*PDE_vec2vec(x_vec).',x_vec);
                end
            elseif isPDEvec(scal)
                y_vec=scal;
                if ~PDE2PDE_compatible(y_vec,x_vec)
                    error('Incompatible PDE vectors to multiply')
                end
                if PDE2PDE_compatible(y_vec,x_vec)==1
                    warning('PDE vectors with different nodes, biggest chosen')
                    nodes_space = max(y_vec.node_space,x_vec.node_space);
                    nodes_time = max(y_vec.node_time, x_vec.node_time);
                    x_vec = reshape(x_vec, nodes_space, nodes_time);
                    y_vec = reshape(y_vec, nodes_space, nodes_time);
                end
                z_vec = x_vec;
                z_vec.parameters = x_vec.parameters.*y_vec.parameters;
                z_vec.vector1D = x_vec.vector1D .* y_vec.vector1D;
                z_vec.vector2D = x_vec.vector2D .* y_vec.vector2D;
            else
                error('Unknown multiplication type')
            end
        end
        % MTIMES
        
        % % DOT
        % function y = dot(x_PDE, vec)
        %     if isa(vec,'PDE_vector')
        %         vec = PDE_vec2vec(vec);
        %     end
        %     x = PDE_vec2vec(x_PDE);
        %     y = x(:).'*vec(:);
        % end
        % % end DOT
        
        % isintval
        function bool = isintval(x_PDE)
            bool = (isintval(x_PDE.parameters)) ||...
                (isintval(x_PDE.vector1D)) || ...
                (isintval(x_PDE.vector2D));
        end
        % isintval
        
        % RESHAPE
        function Xi=reshape(x,node_space, node_time)% FINISH
            % INPUT
            % x             PDEvector
            % node_space    number of nodes used in the output
            % node_time
            % OR
            % INPUT
            % x
            % y             PDEvector, output should fit it
            % OUTPUT
            % Xi         PDEvector with requested number of nodes
            %            if x.nodes>new_nodes, the vectors are cutted,
            %            otherwise they are padded with zeros
            if nargin ==2
                y_vec = node_space;
                node_space = y_vec.node_space;
                node_time = y_vec.node_time;
            elseif isempty(node_space)
                node_space = x.node_space;
            elseif isempty(node_time)
                node_time = x.node_time;
            end
            
            if x.node_space==node_space && x.node_time == node_time
                Xi=x;
                return
            end
            Xi=x;
            if x.node_space<node_space
                dif_nodes=node_space-x.node_space;
                Xi.node_space=node_space;
                if isa(x.vector1D,'intval')
                    Xi.vector1D=[intval(zeros(x.size_vector,dif_nodes)),x.vector1D,intval(zeros(x.size_vector,dif_nodes))];
                    Xi.vector2D=intvalCAT(3,intval(zeros(x.size_vector,size(x.vector2D,2),dif_nodes)),x.vector2D,intval(zeros(x.size_vector,size(x.vector2D,2),dif_nodes)));
                else
                    Xi.vector1D=[zeros(x.size_vector,dif_nodes),Xi.vector1D,zeros(x.size_vector,dif_nodes)];
                    Xi.vector2D=cat(3,zeros(x.size_vector,size(x.vector2D,2),dif_nodes),Xi.vector2D,zeros(x.size_vector,size(x.vector2D,2),dif_nodes));
                end
            elseif x.node_space>node_space
                dif_nodes=x.node_space-node_space;
                Xi.node_space=node_space;
                Xi.vector1D=Xi.vector1D(:,dif_nodes+1:end-dif_nodes);
                Xi.vector2D=Xi.vector2D(:,:,dif_nodes+1:end-dif_nodes);
            end
            
            if x.node_time<node_time
                dif_nodes=node_time-x.node_time;
                Xi.node_time=node_time;
                if isa(x.vector1D,'intval')
                    Xi.vector2D=intvalCAT(2,intval(zeros(x.size_vector,dif_nodes,size(Xi.vector2D,3))),Xi.vector2D,intval(zeros(x.size_vector,dif_nodes,size(Xi.vector2D,3))));
                else
                    Xi.vector2D=cat(2,zeros(x.size_vector,dif_nodes,size(Xi.vector2D,3)),Xi.vector2D,zeros(x.size_vector,dif_nodes,size(Xi.vector2D,3)));
                end
            elseif x.node_time>node_time
                dif_nodes=x.nodes-node_time;
                Xi=Xi_vector(x);
                Xi.node_time=node_time;
                Xi.vector2D=Xi.vector2D(:,dif_nodes+1:end-dif_nodes,:);
            end
            
        end
        % RESHAPE
        
        
        % PDE_VEC2VEC_COMPATIBLE
        function bool = PDE2vec_compatible(x_Xi, x_vec)
            % function bool = Xi_vec2vec_compatible(x_Xi, x_vec)
            %
            % This function tests the compatibility between a vector and a Xi_vector
            % This function returns a bool that is true if the two elements are
            % compatible and false otherwise. The two elements are compatible if the
            % length of the vector is equal to the lenght of Xi_vec2vec(x_Xi), or
            % considered otherwise x_Xi.size_scalar+x_Xi.size_vector*(x_Xi.nodes*2+1)
            %
            % INPUTS
            % x_Xi          PDE_vector
            % x_vec         standard vector
            %
            % OUTPUT
            % bool          check if coherent
            bool=false;
            if  x_Xi.size_real+x_Xi.size_vector*(x_Xi.node_space*2+1)+...
                    x_Xi.size_vector*(x_Xi.node_space*2+1)*(x_Xi.node_time*2+1)...
                    == max(size(x_vec)) && min(size(x_vec))==1
                bool=true;
            end
        end
        % PDE_VEC2VEC_COMPATIBLE
        
        % PDE_VEC2VEC
        function vec = PDE_vec2vec(x_PDE)
            % function vec = PDE_vec2vec(x_PDE)
            %
            % INPUT
            % x_PDE     PDE_vector
            % OUTPUT
            % vec       standard vector, equivalent to x_PDE
            vec = [x_PDE.parameters(:).', x_PDE.vector1D(:).', x_PDE.vector2D(:).'];
        end
        % PDE_VEC2VEC
        
        % LENGTH
        function l = length(x_PDE)
            l = x_PDE.size_real + x_PDE.size_vector*(...
                2*x_PDE.node_space +1 + (2*x_PDE.node_space +1)*(2*x_PDE.node_time +1));
        end
        % end LENGTH
        
        % PDE2PDE_compatible
        function bool = PDE2PDE_compatible(x_PDE, y_PDE)
            % function bool = PDE_vecs_compatible(x_PDE, y_PDE)
            %
            % check if they have the same sizes
            % OUTPUT
            % bool == 0 if they are incompatible
            % bool == 1 if the sizes are compatible but the number of nodes
            % are different
            % bool == 2 if the sizes and the number of nodes are equal
            bool =0;
            if ~isa(x_PDE,'PDE_vector') || ~isa(y_PDE,'PDE_vector')
                return
            elseif  x_PDE.size_real ~= y_PDE.size_real
                return
            elseif x_PDE.size_vector ~= y_PDE.size_vector
                return
            end
            bool = 1;
            if x_PDE.node_space == y_PDE.node_space && x_PDE.node_time == y_PDE.node_time
                bool = 2;
            end
        end
        % PDE2PDE_compatible
        
        %INTVAL
        function y_PDE=intval(x0_PDE, x1_PDE)
            %turn x_Xi into an intval PDE_vector
            if nargin == 1
                y_PDE = x0_PDE;
                y_PDE.parameters=intval(x0_PDE.parameters);
                y_PDE.vector1D=intval(x0_PDE.vector1D);
                y_PDE.vector2D=intval(x0_PDE.vector2D);
                return
            end
            y_PDE= interval_PDE(x0_PDE, x1_PDE);
        end
        %INTVAL
        
        % INTERVAL_PDE
        function xS=interval_PDE(x_PDE,y_PDE)
            % function xBarS=interval_Xi(xBar0,xBar1)
            %
            % taking as input two Xi_vectors, it construct an use_intlab Xi_vector having
            % the two inputs as borders
            
            %global use_intlab
            if ~PDE2PDE_compatible(x_PDE, y_PDE)
                error('Incompatibel inputs to create an interval')
            elseif PDE2PDE_compatible(y_PDE,x_PDE)==1
                warning('PDE vectors with different nodes, biggest chosen')
                nodes_space = max(y_vec.node_space,x_PDE.node_space);
                nodes_time = max(y_vec.node_time, x_PDE.node_time);
                x_PDE = reshape(x_PDE, nodes_space, nodes_time);
                y_PDE = reshape(y_PDE, nodes_space, nodes_time);
                
            end
            xS=x_PDE;
            
            if isintval(x_PDE.parameters)
                
                scal_max=max(sup(real(x_PDE.parameters)),sup(real(y_PDE.parameters)))+...
                    1i*max(sup(imag(x_PDE.parameters)),sup(imag(y_PDE.parameters)));
                scal_min=min(inf(real(x_PDE.parameters)),inf(real(y_PDE.parameters)))+...
                    1i*min(inf(imag(x_PDE.parameters)),inf(imag(y_PDE.parameters)));
                xS.parameters=infsup(scal_min,scal_max);
                
                vec_max=max(sup(real(x_PDE.vector1D)),sup(real(y_PDE.vector1D)))+...
                    1i*max(sup(imag(x_PDE.vector1D)),sup(imag(y_PDE.vector1D)));
                vec_min=min(inf(real(x_PDE.vector1D)),inf(real(y_PDE.vector1D)))+...
                    1i*min(inf(imag(x_PDE.vector1D)),inf(imag(y_PDE.vector1D)));
                xS.vector1D=infsup(vec_min,vec_max);
                
                
                vec_max=max(sup(real(x_PDE.vector2D)),sup(real(y_PDE.vector2D)))+...
                    1i*max(sup(imag(x_PDE.vector2D)),sup(imag(y_PDE.vector2D)));
                vec_min=min(inf(real(x_PDE.vector2D)),inf(real(y_PDE.vector2D)))+...
                    1i*min(inf(imag(x_PDE.vector2D)),inf(imag(y_PDE.vector2D)));
                xS.vector2D=infsup(vec_min,vec_max);
            else
                scal_max=max((real(x_PDE.parameters)),(real(y_PDE.parameters)))+...
                    1i*max((imag(x_PDE.parameters)),(imag(y_PDE.parameters)));
                scal_min=min((real(x_PDE.parameters)),(real(y_PDE.parameters)))+...
                    1i*min((imag(x_PDE.parameters)),(imag(y_PDE.parameters)));
                xS.parameters=infsup(scal_min,scal_max);
                
                vec_max=max((real(x_PDE.vector1D)),(real(y_PDE.vector1D)))+...
                    1i*max((imag(x_PDE.vector1D)),(imag(y_PDE.vector1D)));
                vec_min=min((real(x_PDE.vector1D)),(real(y_PDE.vector1D)))+...
                    1i*min((imag(x_PDE.vector1D)),(imag(y_PDE.vector1D)));
                xS.vector1D=infsup(vec_min,vec_max);
                
                
                vec_max=max((real(x_PDE.vector2D)),(real(y_PDE.vector2D)))+...
                    1i*max((imag(x_PDE.vector2D)),(imag(y_PDE.vector2D)));
                vec_min=min((real(x_PDE.vector2D)),(real(y_PDE.vector2D)))+...
                    1i*min((imag(x_PDE.vector2D)),(imag(y_PDE.vector2D)));
                xS.vector2D=infsup(vec_min,vec_max);
                
            end
        end
        %INTERVAL_PDE
        
        
        %RDIVIDE
        function z_vec=rdivide(x_vec,scal)
            % divide a PDE_vector by a scalar OR another PDE vector
            if any((size(scal))>1) || length(size(scal))>2
                z_vec=x_vec;
                z_vec.parameters=x_vec.parameters/scal;
                z_vec.vector1D=x_vec.vector2D/scal;
                z_vec.vector1D=x_vec.vector2D/scal;
            elseif isPDEvec(scal)
                y_vec=scal;
                if ~PDE2PDE_compatible(y_vec,x_vec)
                    error('Incompatible PDE vectors to multiply')
                end
                if PDE2PDE_compatible(y_vec,x_vec)==1
                    warning('PDE vectors with different nodes, biggest chosen')
                    nodes_space = max(y_vec.node_space,x_vec.node_space);
                    nodes_time = max(y_vec.node_time, x_vec.node_time);
                    x_vec = reshape(x_vec, nodes_space, nodes_time);
                    y_vec = reshape(y_vec, nodes_space, nodes_time);
                end
                z_vec = x_vec;
                z_vec.parameters = x_vec.parameters./y_vec.parameters;
                z_vec.vector1D = x_vec.vector1D ./ y_vec.vector1D;
                z_vec.vector2D = x_vec.vector2D ./ y_vec.vector2D;
            else
                error('Unknown multiplication type')
            end
        end
        %RDIVIDE
        
        % DISP
        function disp(x)
            fprintf('Scalar size : %d,',x.size_real)
            fprintf('Vector size : %d,',x.size_vector)
            fprintf('Nodes in time and space : %d  %d.\n',x.node_time,x.node_space)
            fprintf('Scalar:\n')
            disp(x.parameters)
            disp('')
            fprintf('Vector1D:\n')
            if x.node_space>5
                disp(x.vector1D(:,x.node_space+1+[-5:5]))
            else
                disp(x.vector1D)
            end
            fprintf('Vector2D:\n')
            if x.node_space>5
                disp_vec2 = x.vector2D(:,:,x.node_space+[-5:5]);
            else
                disp_vec2 = x.vector2D;
            end
            if x.node_time>5
                disp_vec2 = disp_vec2(:,x.node_time+[-5:5],:);
            end
            disp_vec2 = squeeze(disp_vec2);
            disp(disp_vec2)
        end
        % DISP
        
        
        % MERGE
        function w = merge(v1, v2, varargin)
            % function w = merge(v1, v2, varargin)
            %
            % takes as input a list of PDE_vector and return a single
            % merged PDE_vector
            if nargin>2
                w = merge(merge(v1,v2), varargin{:});
                return
            end
            if v1.node_space~=v2.node_space
                v1 = reshape(v1, max(v1.node_space,v2.node_space), max(v1.node_time,v2.node_time));
                v2 = reshape(v2, max(v1.node_space,v2.node_space), max(v1.node_time,v2.node_time));
            end
            new_par = [v1.parameters(:).',v2.parameters(:).'];
            if ~isa(v1.vector1D,'intval')
                new_vec1D = cat(1,v1.vector1D,v2.vector1D);
                new_vec2D = cat(1,v1.vector2D,v2.vector2D);
            else
                new_vec1D = intvalCAT(1,v1.vector1D,v2.vector1D);
                new_vec2D = intvalCAT(1,v1.vector2D,v2.vector2D);
            end
            w = PDE_vector(new_par,new_vec1D,new_vec2D);
        end
        % MERGE
        
        % SPLIT
        function [varargout] = split(v)
            % function [w1, w2, varargout] = split(v)
            %
            % takes as input a single PDE_vector, returns n PDE_vectors of
            % vector dimension 1, with n = v.size_vector
            if nargout>v.size_vector
                error('Too many outputs requested')
            end
            length_index = v.size_real/nargout;
            if mod(length_index,1)~=0
                error('Numebr of parameters must be compatible with requested outputs')
            end
            for i = 1:nargout
                index_i = (i-1)*length_index+(1:length_index);
                varargout{i} =PDE_vector(v.parameters(index_i),v.vector1D(i,:),v.vector2D(i,:,:));
            end
        end
        % SPLIT
        
        
        
        % PLOT
        function plot(X, varargin)
            y = X.vector1D;
            z = X.vector2D;
            npoints_y=200*ceil(2*pi);
            xt1=zeros(npoints_y,1);
            t1=linspace(0,2*pi,npoints_y);
            
            m=X.node_space;
            for i=1:npoints_y
                t=t1(i);
                xt1(i)=sum(y.*exp(1i*[-m:m]*t));
            end
            
            %plot_y = ifft(y);
            %f1 = figure;
            figure('DefaultAxesFontSize',14)
            plot(t1,real(xt1),varargin{:});
            
            npoints_y = 200;
            t1 = linspace(0,2*pi,npoints_y);
            npoints_z = 20;
            t2 = linspace(0,1/real(X.parameters(1)),npoints_z);
            xt2 = zeros(length(t2),length(t1));
            
            x_vec = -m:m;
            t_vec = -X.node_time:X.node_time;
            index_t = repmat(t_vec.',1,length(x_vec));
            index_x = repmat(x_vec,length(t_vec),1);
            
            for i=1:npoints_z
                t=t2(i);
                for j = 1:npoints_y
                    x = t1(j);
                    xt2(i,j)=sum(sum(squeeze(z(1,:,:)).*exp(1i*index_t*t+1i*index_x*x)));
                end
            end
            
            %npoints=200*ceil(1/real(X.parameters(1)));
            %xt2=zeros(npoints,1);
            f2=figure;
            surf(t2,t1,real(xt2).',varargin{:})
            shading interp
            xlabel('time','Interpreter','Latex','FontSize',20)
            ylabel('space','Interpreter','Latex','FontSize',20)
        end
        % end PLOT
        
        
        
        
        %
        %         %PLOT
        %         function plot(X,varargin)
        %             % default plotting option for Xi_vector
        %             if abs(X.parameters)>0.1
        %             npoints=200*ceil(1/real(X.parameters(1)));
        %             t1=linspace(0,1/X.parameters(1),npoints);
        %             else
        %             npoints=200*ceil(1/0.1);
        %             t1=linspace(0,1/0.1,npoints);
        %             end
        %             xt1=zeros(npoints,X.size_vector);
        %             %xt2=zeros(200,1);
        %             m=X.nodes;
        %             for i=1:npoints
        %                 t=t1(i);
        %                 for j=1:X.size_vector
        %                     xt1(i,j)=sum(X.vector(j,:).*exp(1i*[-m:m]*2*pi*t*X.parameters(1)));
        %                 end
        %                 %xt2(i)=sum(xBar{3}.*exp(1i*[-m:m]'*2*pi*t*xBar{1}));
        %             end
        %
        %             %scrsz = get(0,'ScreenSize');
        %             %figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
        %             for j=1:X.size_vector
        %                 subplot(1,X.size_vector,j);
        %                 plot(t1,real(xt1(:,j)),varargin{:})
        %                 xlabel('t')
        %                 ylabel(sprintf('component %d',j));
        %                 xlim([0,1/X.parameters(1)])
        %             end
        %         end
        %         %PLOT
        %
        %         %PLOT2
        %         function plot2(X,varargin)
        %             % 2D plotting option for Xi_vector
        %             if X.size_vector ~= 2
        %                 error('This function is suppose to bejust for 2D vectors');
        %             end
        %             if abs(X.parameters(1))>0.1
        %             npoints=200*ceil(1/real(X.parameters(1)));
        %             t1=linspace(0,1/X.parameters(1),npoints);
        %             else
        %             npoints=200*ceil(1/0.1);
        %             t1=linspace(0,1/0.1,npoints);
        %             end
        %             xt1=zeros(npoints,X.size_vector);
        %             %xt2=zeros(200,1);
        %             m=X.nodes;
        %             for i=1:npoints
        %                 t=t1(i);
        %                 for j=1:X.size_vector
        %                     xt1(i,j)=sum(X.vector(j,:).*exp(1i*[-m:m]*2*pi*t*X.parameters(1)));
        %                 end
        %                 %xt2(i)=sum(xBar{3}.*exp(1i*[-m:m]'*2*pi*t*xBar{1}));
        %             end
        %
        %             %scrsz = get(0,'ScreenSize');
        %             %figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
        %
        %             plot(real(xt1(:,1)),real(xt1(:,2)),varargin{:})
        %             xlabel('x')
        %             ylabel('y');
        %             xlim([-0.6/X.parameters(1),0.6/X.parameters(1)])
        %             axis equal
        %             if max(abs(imag(xt1)))>max(abs(real(xt1)))
        %                 warning('The solution plotted has an iportant imaginary part')
        %             end
        %         end
        %         %PLOT2
       
        
    end
end