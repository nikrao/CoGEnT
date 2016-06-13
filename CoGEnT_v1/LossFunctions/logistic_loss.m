function f = logistic_loss(y,Phi,wt,x)

if wt == 0

    f = sum(log(1 + exp(-y.*(Phi*x))));

else
    
    alpha = sum(y==1)/numel(y);
    tmp = log(1 + exp(-y.*(Phi*x)));
    avec = (1-alpha)*ones(length(y),1);
    avec(y==1) = alpha;
    f = sum(tmp.*avec);

end