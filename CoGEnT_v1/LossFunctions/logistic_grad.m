function f = logistic_grad(y,Phi,wt,x)

if wt == 0

    tmp = 1 + exp(-y.*(Phi*x));
    f = Phi'*(y./tmp);

else
    
    alpha = sum(y==1)/numel(y);
    tmp = 1 + exp(-y.*(Phi*x));
    avec = (1-alpha)*ones(length(y),1);
    avec(y==1) = alpha;
    
    tmp = avec.*y./tmp;
    f = Phi'*tmp;

end