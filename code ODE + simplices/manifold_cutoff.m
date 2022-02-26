function bool_growth = manifold_cutoff(node)
% function bool_growth = manifold_cutoff(node)
%
% this function SHOULD BE CHANGED BY THE USER for the specific case at hand
%
% this function returns a boolean indicating if the 2d-manifold should grow
% from the given node
% this function is not called for the starting node.
%
% INPUT
% node
% OUTPUT 
% bool:    should the manifold grow from this node
%
% notice that an additional prioritisation of nodes is also taking place,
% so it might be possible that some nodes are never grown even if they
% return True

%% amplitude cutoff
% in cases of bubbles, this allows to only grow the manifold for positive
% amplitudes. A small negative amplitude is accepted for completeness

amplitude = node.solution.scalar(4);
if amplitude>=-2E-2
    bool_growth = 1;
else
    bool_growth = 0;
end

%% additional cutoff for lorenz 84 - comment or decomment
%
% if amplitude>=-0.1 && amplitude<=0.3
%     bool_growth = 1;
% else
%     bool_growth = 0;
% end
