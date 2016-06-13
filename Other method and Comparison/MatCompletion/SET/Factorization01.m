function VS = Factorization01(X0,IndexM,U,VS0)
% fix U, calculate the best match V according to VS0

[m,n] = size(X0);
r = size(U,2);

if isempty(VS0) == 1
    VS0 = zeros(n,r);
end

Vt = zeros(r,n);
for cn=1:n
    cnz_p = IndexM(:,cn)~=0;
    
    [Vt(:,cn),flag] = ...
        qmr(U(cnz_p,:)'*U(cnz_p,:),U(cnz_p,:)'*X0(cnz_p,cn),...
        1e-10,max(m,n)^2,[],[],VS0(cn,:)');
end
VS = Vt'; 