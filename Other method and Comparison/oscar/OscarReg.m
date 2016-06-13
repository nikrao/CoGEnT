function y  = OscarReg(W,l2)


[n t] = size(W);
w = ones(n,1) + l2*(n*ones(n,1) - (1:n)');
   
y = w.*sort(abs(W),'descend');
y = sum(y);


end
