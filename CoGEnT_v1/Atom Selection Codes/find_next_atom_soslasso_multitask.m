function [next_atom] = find_next_atom_soslasso_multitask(gradf,groups)
% INPUTS:
% gradf  = gradient of f at current iterate
% groups = groupings over the variables
% NORM : \sum ||xg||_2 * sqrt(|G|) + ||xg||_1
% OUTPUT:
% next_atom = next atom

[nr nc] = size(gradf);
gradf = -gradf(:);


temp  = (sqrt(sum((repmat(gradf,1,size(groups,2)).*groups).^2)))./(sqrt(sum(groups)));
[~,bestgrp] = max(temp);

inds = find(groups(:,bestgrp)==1); % this is the group that maximizes norm
% now pick an element in this group that is the "best"
temp = zeros(length(gradf),1);
temp(inds) = gradf(inds)/norm(gradf(inds));
temp = abs(temp);
[~,inds] = max(temp);
next_atom = zeros(length(gradf),1);
next_atom(inds) = sign(gradf(inds));

next_atom = reshape(next_atom,nr,nc);
end