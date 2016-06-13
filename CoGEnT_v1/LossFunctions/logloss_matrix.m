function f = logloss_matrix(y,Phi,wt,x)

if wt == 0
   
    x = x(Phi);
    f = sum(log(1 + exp(-y.*x)));
    
else
    
    a = sum(y==1)/numel(y);
    avec = (1-a)*ones(length(y),1);
    avec(y==1) = a;
    f = sum(avec.*log(1 + exp(-y.*x)));

end

end