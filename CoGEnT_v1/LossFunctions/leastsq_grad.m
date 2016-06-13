function f = leastsq_grad(y,Phi,x)
    f = -Phi'*(y - Phi*x);

end