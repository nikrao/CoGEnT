% compressed sensing off the grid using atomic norms

clear;
clc;
close all;

N =100;  % signal length
n =40;   % number of random measurements obtained
S=randsample(N,n);
S=S-1;    % because the indices go from 0,1,...,N-1
noiserange = linspace(0,1,100);
ntests = 100;
nf = 0.16/n;

nind = 0;
for sig = noiserange
    nind = 1:nind;
    
    for test = 1:tests
        
        % act_freq = linspace(0.05,1,10);  % actual frequencies on the grid
        act_freq = rand(5,1);            % arbitrary random frequencies
        a0=rand(length(act_freq),1);     % random weights for the actual frequencies
        
        keepitreal = 0;                  % if ON, only deals with DCT instead of DFT
        
        [act_freq,afind] = sort(act_freq);
        a0 = a0(afind);
        NEAR_FIELDl = act_freq-nf;NEAR_FIELDr = act_freq+nf;
        NEAR_FIELD = [NEAR_FIELDl NEAR_FIELDr];
        % this is the near field region
        ffl = 0;
        FAR_FIELD = [];
        for ii = 1:size(NEAR_FIELD,1)
            vec = [ffl NEAR_FIELD(ii,1)];
            ffl = NEAR_FIELD(ii,2);
            FAR_FIELD = [FAR_FIELD;vec];
        end
        FAR_FIELD = [FAR_FIELD;[FAR_FIELD(end,2) 1]];
        
        Atrue = [];
        % true  signal
        x=zeros(n,1);
        for jj = 1:length(act_freq)
            
            if ~keepitreal
                x = x + 1/sqrt(n)*a0(jj)*exp(1i*2*pi*act_freq(jj)*S);
                Atrue = [Atrue exp(1i*2*pi*act_freq(jj)*S)];
            else
                x = x + 1/sqrt(n)*a0(jj)*cos(2*pi*act_freq(jj)*S);
                Atrue = [Atrue cos(2*pi*act_freq(jj)*S)];
            end
        end
        xtrue = x;
        
        x = x+sig*randn(n,1);  % observations
        xp=zeros(n,1);
        
        f0=rand;  % random initial starting point
        
        init_grid = linspace(0,1,10000)';  % initial grid, uniformly spaced
        tau = norm(a0,1) + 0.05;   % regularization parameter
        
        
        if keepitreal
            At = 1/sqrt(n)*cos(2*pi*f0*S);
        else
            At = 1/sqrt(n)*exp(1i*2*pi*f0*S);
        end
        % ^ initial atom
        coeft=At\x;
        xp=At*coeft;
        
        maxiter = 500;
        bcount = 0;
        fbins=[f0];
        
        
        % ALGORITHM PARAMETERS
        gridding=0;
        maxiter_gp = 200;  % as of now we are not using gradient projections
        tol_gp = 1e-4;
        eta = 0.5;
        usecvx = 1;   % if ON, use CVX for the optimization step
        verbose = 0;  % print messages and diagnostics
        backcount = 20;  % to check convergence
        tol = 1e-10 ; % to check convergence
        
        
        grid=init_grid;
        drops = 0;
        for ii = 2:maxiter
            
            if ii>20
                gridding=1;
            end
            
            if verbose
                disp(ii)
            end
            
            oldf = 0.5*norm(At*coeft - x)^2;
            if ii==1
                foba_err(ii) = oldf;
            end
            gradf=(x-At*coeft);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Adaptive Grid Refinement %%%%%%%%
            if gridding==1;
                Delta_grid1=sort(fbins);
                Delta_grid1=circshift(Delta_grid1',1)';
                Delta_grid=min(abs(Delta_grid1-sort(fbins)));
                Delta_grid=min(Delta_grid,1-max(fbins));
                
                R=Delta_grid/10; %Set Refinement to be a tenth of the signal resolution  Delta
                N_R=20; %Discretize around each point in the support +-R at a discretization level of N_R
                
                for i=1:length(fbins)
                    grid_temp=linspace(fbins(i)-R,fbins(i)+R,N_R);
                    grid=[grid; grid_temp'];
                end
                
            end
            
            grid=grid.*(grid>=0).*(grid<=1);
            grid = sort(grid);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fval = evaluate_trig_pari(gradf,grid,S, keepitreal);
            [~,f] = max(fval);
            
            fbins_old = fbins;
            
            fbins = [fbins; grid(f)]; %append the new selected frequency
            
            
            a = 1/sqrt(n)*exp(1i*2*pi*grid(f)*S);  % new atom
            if keepitreal
                a =  1/sqrt(n)*cos(2*pi*grid(f)*S);
            end
            
            At_old = At;
            At = [At a];   % the current atomic set
            coeft_old = coeft;
            coeft = [coeft;0];
            
            if ~usecvx
                [coeft1,err]=gp_reduced(At,x,coeft,tau,maxiter_gp,tol_gp);
            else
                cvx_begin quiet
                variable coeft(size(At,2))
                minimize sum_square(real(x)-real(At)*coeft)+sum_square(imag(x)-imag(At)*coeft)
                %minimize sum_square(real(x)-real(At)*coeft + imag(x)-imag(At)*coeft)
                subject to
                coeft >= 0;
                sum(coeft) <= tau;
                cvx_end
                
                err = 1;
            end
            if err==0
                coeft = coeft1;
            end
            
            xp=At*coeft;
            newf = 0.5*norm(xp - x)^2;
            foba_err(ii) = newf;
            
            do_multiple = 1;
            while do_multiple
                gradfb = (xp-x);
                [atom_rem, idx_rem] = min(-coeft.*(At'*gradfb)+0.5*coeft.*coeft.* diag(At'*At));
                
                if idx_rem == size(At,2)
                    drops = drops + 1;
                end
                
                % no backward steps initially
                if ii<=10
                    break
                end
                
                A_new = At;
                A_new(:,idx_rem) = [];
                fbin_new = fbins;
                if verbose
                    fprintf(1,'\n deleted coefficient value %10.5e \n', fbin_new(idx_rem))
                end
                
                fbin_new(idx_rem) = [];
                
                if ~usecvx
                    [beta_new,err]=gp_reduced(A_new,x,fbin_new,tau,maxiter_gp,tol_gp);
                else
                    cvx_begin quiet % 1 D line search
                    variable beta_new(size(A_new,2))
                    minimize sum_square(real(x)-real(A_new)*beta_new)+sum_square(imag(x)-imag(A_new)*beta_new)
                    subject to
                    beta_new >= 0;
                    sum(beta_new) <= tau;
                    cvx_end
                end
                
                if err==1
                    beta_new = coeft;
                    beta_new(idx_rem) = [];
                end
                
                % reoptimize xp based on the "new" atomic set.
                x_rem = A_new*beta_new;
                
                newf_rem = 0.5* norm(x - x_rem)^2;
                
                if ((newf_rem - newf < (oldf-newf)*eta))
                    bcount = bcount+1;
                    At = A_new;
                    xp = x_rem;
                    coeft = beta_new;
                    fbins = fbin_new;
                    foba_err(ii) = newf_rem;
                    if size(At,2)<=1
                        do_multiple = 0;
                    end
                    
                    if verbose
                        fprintf(1,'\n iter after backward step: %4d, f: %10.5e \n', ii,newf_rem);
                    end
                else
                    do_multiple = 0;
                end
                
            end
            
            if verbose
                norm(x - At*coeft)^2/N
                stem(act_freq,a0); hold on; stem(fbins',abs(coeft),'r .');
                hold off
                pause(.1);
            end
            
            % convergence check
            if ((foba_err(ii-1) - foba_err(ii) < tol) && (ii > 2) && (foba_err(ii-1)>foba_err(ii)))
                if verbose
                    fprintf(2,'\n Convergence achieved \n');
                end
                break;
            end
            if drops >= backcount
                if verbose
                    fprintf(2,'\n Dropcount achieved \n');
                end
                break;
            end
            
            
        end
        
        
        % merge the coefficients corresponding to a single frequency
        U = unique(fbins);
        for jj = 1:length(U);
            
            fi = U(jj);
            idx = find(fbins == fi);
            f_final(jj) = fi;
            coeft_final(jj) = sum(coeft(idx));
            
        end
        
        % if both 1 and 0 are present, merge them into 1
        is1 = find(f_final == 1);
        is0 = find(f_final == 0);
        if ~isempty(is1) && ~isempty(is0)
            f_final(is0) = [];
            coeft_final(is1) = coeft_final(is1) + coeft_final(is0);
            coeft_final(is0) = [];
        end
        is1 = find(act_freq == 1);
        is0 = find(act_freq == 0);
        if ~isempty(is1) && ~isempty(is0)
            act_freq(is0) = [];
            a0(is1) = a0(is1) + a0(is0);
            a0(is0) = [];
        end
        
        % CALCULATE ERRORS
        % spurious amplitudes 
        Flength = length(f_final);
        farlength  =size(FAR_FIELD,1);
        farfieldsum = 0;
        isfarfield = 0;
        for ii = 1:Flength
            thisfreq = f_final(ii);
            isfarfield = 0;
            for jj = 1:farlength
                lb = FAR_FIELD(jj,1);ub = FAR_FIELD(jj,2);
                if (thisfreq>lb)&&(thisfreq<ub)
                    isfarfield = 1;
                    break
                end
            end
            
            if isfarfield==1
                farfieldsum = farfieldsum+abs(coeft_final(ii));
            end
        end
        FF_SUM(test) = farfieldsum;
        
        % near region approximation
        nearfieldsum = 0;
        nearlength = size(NEAR_FIELD,1);
        
        
        % weighted frequency deviation
            
                
            
        
        
                
                
                
                
                
                
    end
    
    stem(act_freq,a0); hold on; stem(f_final,coeft_final,'r .');
    % find all freqs in the far region

        
    
    
