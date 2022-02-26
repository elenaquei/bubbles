function alpha = interpolation(alpha0, varargin)
    function varargout = comma_separated_ij(i,j,varargin)
        varargout = cell(length(varargin),1);
        for k = 1:length(varargin)
            varargout{k} = varargin{k}{i,j};
        end
    end

if isa(alpha0, 'float') || isa(alpha0, 'intval')
    alpha = interpolation_float(alpha0, varargin{:});
elseif isa(alpha0, 'Xi_vector')
    alpha = interpolation_Xivec(alpha0, varargin{:});
elseif isa(alpha0, 'scalar_eqs')
    alpha = interpolation_scal(alpha0, varargin{:});
elseif isa(alpha0, 'full_problem')
    alpha = interpolation_full_problem(alpha0, varargin{:});
elseif isa(alpha0, 'cell')
    alpha = cell(size(alpha0));
    for i = 1:size(alpha0,1)
        for j = 1: size(alpha0,2)
            alpha{i,j} = interpolation(alpha0{i,j}, comma_separated_ij(i,j,varargin{:}));
        end
    end
else
    error('Not valid input')
end

end

function alpha = interpolation_full_problem(alpha0, varargin)
    function varargout = scalar_eqs_var(varargin)
        varargout = cell(nargin, 1);
        for k = 1:nargin
            varargout{k} = varargin{k}.scalar_equations;
        end
    end
alpha = alpha0;
alpha.scalar_equations = interpolation_scal(alpha0.scalar_equations, scalar_eqs_var(varargin{:}));
end

function alpha = interpolation_scal(alpha0, varargin)
    function varargout = lin_coef_var(i,varargin)
        varargout = cell(nargin, 1);
        for k = 1:nargin-1
            varargout{k} = varargin{k}.linear_coef{i};
        end
    end
alpha = alpha0;
for i =1:3
    alpha.linear_coef{i} = interpolation_float(alpha0.linear_coef{i},...
        lin_coef_var(i,varargin{:}));
end
end

function alpha = interpolation_Xivec(alpha0, varargin)
    function varargout = Xi_vec2vec_var(varargin)
        varargout = cell(nargin, 1);
        for k = 1:nargin
            varargout{k} = Xi_vec2vec(varargin{k});
        end
    end
vec = interpolation_float( Xi_vec2vec(alpha0),Xi_vec2vec_var(varargin{:}));
alpha= vec2Xi_vec(vec, alpha0);
end

function alpha = interpolation_float(varargin)
    function sum_r = sum_real_intval(varargin)
        sum_r = 0;
        for k = 1:nargin
            sum_r = sum_r + real(intval(varargin{k}));
        end
    end
    function sum_r = sum_imag_intval(varargin)
        sum_r = 0;
        for k = 1:nargin
            sum_r = sum_r + imag(intval(varargin{k}));
        end
    end
    function bool = all_imag(varargin)
        bool = 1;
        for k = 1:nargin
            if any(any(imag(varargin{k})))
                bool = 0;
                return
            end
        end
    end
    function distance = max_distance(point, varargin)
        distance = 0*point;
        for k = 1: nargin - 1
            distance = max(distance, sup(abs(intval(point) - intval(varargin{k}))));
        end
    end

if nargin == 1
    alpha = intval(varargin{1});
    return
end
if  all_imag(varargin{:})
    mid_point = ( sum_real_intval(varargin{:}))/nargin;
else
    mid_point = (sum_real_intval(varargin{:}))/nargin+ ...
        1i*(sum_imag_intval(varargin{:}))/nargin;
end
alpha = midrad(mid_point.mid,max_distance(mid_point.mid, varargin{:}));
end

function val_intval = three_intval_real(a,b)
val_min = min(a, b);
val_max = max(a, b);
if isintval(val_min)
    val_min = inf(val_min);
end
if isintval(val_max)
    val_max = sup(val_max);
end
val_intval = infsup(val_min, val_max);
end