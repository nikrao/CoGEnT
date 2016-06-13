function f = sqhinge_grad(y,Phi,wt,x)

if wt == 0
    
    f = -2*y.*max(0,1 - y.*(Phi*x));
    f = Phi'*f;

else
    
    alpha = sum(y==1)/numel(y);
    
    mul =max(0,1 - y.*(Phi*x));
    avec = (1-alpha)*ones(length(y),1);
    avec(y==1) = alpha;
    
    f = -2*avec.*y.*mul;
    f = Phi'*f;
    
end