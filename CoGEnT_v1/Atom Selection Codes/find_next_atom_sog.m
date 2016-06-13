function [next_atom] = find_next_atom_sog(gradf,groups,k)
% INPUTS:
% gradf = gradient of f at current iterate
% groups= sparse matrix . rows = # elements, cols = # groups
% OUTPUT:
% next_atom = next atom

gradf = -gradf;
% the atomic dual norm is the max over the norms of the groups

val = repmat(gradf,1,size(groups,2)).*groups;

% retain only the top k of each column
val = abs(val);
[val,indx] = sort(val,1,'descend');
val = val(1:k,:);
val = sum(val.^2);

[~, bestind] = max(val);
indx = indx(1:k,bestind);

% we need to maximally correlate the negative gradient
next_atom = zeros(length(gradf),1);
%active_idx = find(groups(:,bestind)==1);
%active_idx = groups(indx,bestind);
next_atom(indx) = gradf(indx);


next_atom = next_atom/norm(next_atom);

end