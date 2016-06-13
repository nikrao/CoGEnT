function [f_final, coeft_final] = CoGEnT_offgrid(x,S,maxiter,epsilon, tau)


% x       : observations
% maxiter : max. iterations
% epsilon : merging window
% tau     : atomic norm parameter

%% Initialization for the Algorithm
n = length(x);
init_grid = linspace(0,1,10000)';  % initial grid, uniformly spaced
f0=rand;  % random initial starting point

At = 1/sqrt(n)*exp(1i*2*pi*f0*S);
% ^ initial atom
coeft=At\x;
coeft = real(coeft);
xp=At*coeft;

fbins= f0;


%% ALGORITHM PARAMETERS
maxiter_gp = 20;  
tol_gp = 1e-3;
eta = 0.8;
backcount = 20;  % to check convergence
tol = 1e-10 ; % to check convergence
grid=init_grid;
drops = 0;
%% begin iterations

for ii = 2:maxiter
    
    oldf = 0.5*norm(xp - x)^2;
    gradf=(x-xp);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Adaptive Grid Refinement %%%%%%%%
    for jj = 1:length(fbins)-1
        grid = [grid; 0.5*(fbins(jj) + fbins(jj+1))]; 
    end
    
    grid=grid.*(grid>=0).*(grid<=1);
    grid = sort(grid);
    zergrid = find(grid==0);
    grid(zergrid) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fval = evaluate_trig(gradf,grid,S, 0);
    [~,f] = max(fval);

    
    fbins = [fbins; grid(f)]; %append the new selected frequency
    
    
    a = 1/sqrt(n)*exp(1i*2*pi*grid(f)*S);  % new atom
    
    At = [At a];   % the current atomic set
    coeft_old = coeft;
    
    phinext=(tau*a-xp);
    gamma=-((x-xp)'*phinext)/(phinext'*phinext);
    coeft = [(1-gamma)*coeft;gamma*tau];

    [coeft,~]=gp_reduced_complex(At,x,coeft,tau,maxiter_gp,tol_gp);
    
    xp=At*coeft;
    newf = 0.5*norm(xp - x)^2;
    err(ii) = newf;
    
    
    %%%%%%%% BACKWARD STEP %%%%%%%%%%%%%
        
        gradfb = (xp-x);
        [atom_rem, idx_rem] = min(-coeft.*(At'*gradfb)+0.5*coeft.*coeft.* diag(At'*At));
        A_new = At;
        A_new(:,idx_rem) = [];
        fbin_new = fbins;
        fbin_new(idx_rem) = [];
        coeft_new = coeft;coeft_new(idx_rem) = [];
        
        [fnew,coeft_new] = refine_frequencies(fbin_new,coeft_new,epsilon);
        
        % add new frequencies to grid
        idx_rem = find(grid == fbins(idx_rem));
        grid(idx_rem) = [];
        grid = [grid;fnew];
        
        % form corresponding atoms
        Anew = [];
        for ff = 1:length(fnew)
            Anew = [Anew 1/sqrt(n)*exp(1i*2*pi*fnew(ff)*S)];
        end
        [coeft_new,~]=gp_reduced_complex(Anew,x,zeros(length(coeft_new),1),tau,maxiter_gp,tol_gp);
        xnew = Anew*coeft_new;
        
        newf_rem = 0.5*norm(x - xnew)^2;
        
         if ((newf_rem - newf) <= ((oldf-newf)*eta))
            At = Anew;
            xp = xnew;
            coeft = coeft_new;
            fbins = fnew;
            err(ii) = newf_rem;
            drops = 1+drops;

         end
    
    

    
    % convergence check
    if ((err(ii-1)-err(ii) < tol) && (ii > 2) && (err(ii-1)>= err(ii)))
        break;
    end
    if drops >= backcount
        break;
    end
        
end

%% merge the coefficients corresponding to a single frequency

f_final = fbins;
coeft_final = coeft;
% if both 1 and 0 are present, merge them into 1
is1 = find(f_final == 1);
is0 = find(f_final == 0);

if ~isempty(is1) && ~isempty(is0)
    f_final(is0) = [];
    coeft_final(is1) = coeft_final(is1) + coeft_final(is0);
    coeft_final(is0) = [];
end


% merge coefficients epsilon close to each other
[f_final,coeft_final] = refine_frequencies(f_final,coeft_final,epsilon);
A_final=[];
for f = 1:length(f_final)
   A_final = [A_final 1/sqrt(n)*exp(1i*2*pi*f_final(f)*S)]; 
end
[coeft_final,~]=gp_reduced_complex(A_final,x,zeros(length(coeft_final),1),tau,maxiter_gp,tol_gp);


% remove zeros
inds = find(coeft_final~=0);
f_final = f_final(inds);
coeft_final = coeft_final(inds);

end

