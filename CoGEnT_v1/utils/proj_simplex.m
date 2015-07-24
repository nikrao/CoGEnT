function t = proj_simplex(v,tau)

% function to perform projection on the simplex (tau) on O(n)

n = length(v);
% intializations
U = 1:n;
s = 0; rho = 0;

while ~isempty(U)
    k = randsample(U,1);
    % pick an entry in U at random
    
    G = find(v(U)>=v(k));
    L = find(v(U)<v(k));
    
    delrho = numel(G);
    dels   = sum(v(G));
    
    if (s+dels-(rho+delrho)*v(k)<tau)
        s = s+dels;
        rho = rho+delrho;
        U = L;
    else
        U = setdiff(G,k);
    end
    
    theta = (s - tau)/rho;
end

t = max(0,v-theta);