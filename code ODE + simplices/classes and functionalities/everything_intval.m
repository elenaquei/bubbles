function [varargout] = everything_intval(varargin)

nOutputs = nargout;
varargout = cell(1,nOutputs);

for k = 1:nOutputs
    varargout{k} = intval(varargin{k});
end