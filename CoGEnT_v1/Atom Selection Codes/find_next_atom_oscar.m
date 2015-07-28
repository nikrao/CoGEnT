function [next_atom] = find_next_atom_oscar(gradf,Li)
% INPUTS:
% gradf = gradient of f at current iterate
% Li = multiplier for the Linf norm
% OUTPUT:
% next_atom = next atom

gradf = -gradf;

n = length(gradf);
w = ones(n,1) + Li*(n*ones(n,1) - (1:n)');
t = cumsum(w);
t = 1./t;

[g,inds] = sort(abs(gradf),'descend');
g = cumsum(g);
tg = t.*g;
[~,kstar] = max(tg);

next_atom = zeros(n,1);
next_atom(inds(1:kstar)) = gradf(inds(1:kstar));
next_atom = next_atom/norm(next_atom);

% 
% p = length(gradf);
% % inds = nchoosek(1:p,2)';
% val = sum(abs(gradf(pairs)));
% val = [val abs(gradf)'];
% 
% 
% [~, bestind] = max(val); %this is the group that corresponds to the max norm
% 
% % we need to maximally correlate the negative gradient
% next_atom = zeros(length(gradf),1);
% if bestind <= cols(pairs)
%     active_idx = pairs(:,bestind);
% else
%     active_idx = bestind-cols(pairs);
% end
% next_atom(active_idx) = sign(gradf(active_idx));
% if (length(active_idx) == 2)
%     next_atom = next_atom*c;
% end
% 
% % remove zeros from the 'active_idx' indices
% idx = find(next_atom(active_idx) == 0);
% vec = 2*binornd(1,0.5,length(idx),1)-1;
% temp = next_atom(active_idx);
% temp(idx) = vec;
% next_atom(active_idx) = temp;

end