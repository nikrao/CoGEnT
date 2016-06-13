function f = sqhinge_loss(y,Phi,wt,x)

if wt == 0
    
    f = sum(max(0,1 - y.*Phi*x).^2);
    
   
else
    
    alpha = sum(y==1)/numel(y);
    
    mul =max(0,1 - y.*Phi*x).^2;
    avec = (1-alpha)*ones(length(y),1);
    avec(y==1) = alpha;
    
    f = sum(avec.*mul);
    

end