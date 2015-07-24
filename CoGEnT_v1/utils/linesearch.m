function a = linesearch(lossfun,x,p)
tau = 0.5;
c = 0.5;

m = -norm(p)^2;
t = -c*m;j = 0;
j = 0;
keep_going = true;



while keep_going
	LHS = lossfun(x)-lossfun(x + a*p); RHS = a*t;
	if LHS<RHS
	j = j+1;
	a = tau*a;
	else
	keep_going = false;
	end
	if j >= maxit
	keep_going = false;
	end
	
end


