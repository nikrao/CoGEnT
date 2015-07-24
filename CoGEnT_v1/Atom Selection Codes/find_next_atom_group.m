function [next_atom] = find_next_atom_group(gradf,groups)
% INPUTS:
% gradf = gradient of f at current iterate
% groups= sparse matrix . rows = # elements, cols = # groups
% OUTPUT:
% next_atom = next atom

gradf = -gradf;
% the atomic dual norm is the max over the norms of the groups

val = repmat(gradf,1,size(groups,2)).*groups;
val = sum(val.^2);

[~, bestind] = max(val);

% we need to maximally correlate the negative gradient
next_atom = zeros(length(gradf),1);
active_idx = find(groups(:,bestind)==1);
next_atom(active_idx) = gradf(active_idx);
next_atom = next_atom/norm(next_atom);

end