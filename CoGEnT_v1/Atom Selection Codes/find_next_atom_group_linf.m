function [next_atom] = find_next_atom_group_linf(gradf,groups)
% INPUTS:
% gradf = gradient of f at current iterate
% groups= groupings over the variables 
% OUTPUT:
% next_atom = next atom

gradf = -gradf;
% the atomic dual norm is the max over the L1 norms of the groups
mult = sqrt(sum(groups));
val = repmat(gradf,1,size(groups,2)).*groups;
val = sum(abs(val))./mult;

[~, bestind] = max(val); %this is the group that corresponds to the max norm

% we need to maximally correlate the negative gradient
next_atom = zeros(length(gradf),1);
active_idx = find(groups(:,bestind)==1);
next_atom(active_idx) = sign(gradf(active_idx));

% remove zeros from the 'active_idx' indices
idx = find(next_atom(active_idx) == 0);
vec = 2*binornd(1,0.5,length(idx),1)-1;
temp = next_atom(active_idx);
temp(idx) = vec;
next_atom(active_idx) = temp/norm(temp);

end