function [ alpha1 ] = linesearch( y,x )
cvx_solver sdpt3
cvx_begin quiet
variable alpha1(1,1)
minimize norm(y-x*alpha1)
subject to 
alpha1>=0
alpha1<=1
cvx_end

end

