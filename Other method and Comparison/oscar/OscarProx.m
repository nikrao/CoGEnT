function y = OscarProx(W,lambda,l2)


[n t] = size(W);
w = lambda + lambda*l2*(n*ones(n,1) - (1:n)');
y = zeros(n,t);


[tv, Pv] = sort(abs(W), 'descend');

u = sign(tv).*max(0,tv - w);

y = u(Pv);

end
