function next_atom = find_next_atom_DCT(gradf)

% function to find most correlated DCT atom
gradf = -gradf;
n = length(gradf);
vals = zeros(n-1,1);
for k = 2:n
    v = cos((0:n-1)*pi*(k-1)/n)';
    v = v-mean(v);
    v = v/norm(v);
    vals(k-1) = gradf'*v;
end
[~,k] = max(abs(vals));
k = k+1;
v = cos((0:n-1)*pi*(k-1)/n)';
v = v-mean(v);
next_atom = (v/norm(v))*sign(vals(k-1));
end
