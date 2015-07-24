function [coeft] = grad_proj_gen(A,y,coeft,tau,gpiter,gptol,lossfun,gradfun)
	
xz= coeft;
t = 1;
t_old = 0;
iter = 0;
gamma = 1;
gamma_inc = 2;
xold = 0;

while iter < gpiter
    alpha = (t_old - 1) /t;
    
	
    xs = (1 + alpha) * xz - alpha * xold;
    
    % compute function value and gradients of the search point
	gxs = gradfun(xs,A,y)
	
    
    % the Armijo Goldstein line search scheme
    while true
		
		fs = lossfun(xs,A,y);
		% project temp onto L1 ball of tau
		xzp = x = projsplx(xs - gxs/gamma,tau);		
		fzp = lossfun(xzp,A,y);
        
        delta_Wzp = xzp - xs;
        nrm_delta_Wzp = norm(delta_Wzp, 'fro')^2;
        r_sum = (nrm_delta_Wzp)/2;
        
        fzp_gamma = fs + sum(sum(delta_Wzp .* gWs))...
            + gamma/2 * nrm_delta_Wzp;
        
        if (r_sum <=1e-20)
            break;
        end
        
        if (fzp <= fzp_gamma)
            break;
        else
            gamma = gamma * gamma_inc;
        end
    end
    
    xold = xz;
    xz = xzp;
        
    if (bFlag)
        break;
    end
    
	if (lossfun(xz,A,y) - fold <= gptol)
		break;
	end
	
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
end

coeft = xz;	

end		