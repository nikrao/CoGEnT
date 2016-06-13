function  [coeft_n]=grad_proj_pari(At1,y,tau);

cvx_solver sdpt3
cvx_begin quiet
variable coeft_n(size(At1,2))
minimize sum_square(y-At1*coeft_n)
subject to 
norm(coeft_n,1)<=tau
cvx_end
end

