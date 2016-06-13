% compressed sensing off the grid using atomic norms 

clear;
clc;
close all;

N =1000;  % signal length
n =250;   % number of random measurements obtained
S=randsample(N,n);
S=S-1;    % because the indices go from 0,1,...,N-1
minfreq = 0.05;
maxfreq = 0.95;
act_freq = rand(10,1)+minfreq;            % arbitrary random frequencies
goodinds = find(act_freq<=maxfreq);
act_freq = act_freq(goodinds);

a0=rand(length(act_freq),1);     % random weights for the actual frequencies
Atrue = [];
% true  signal
x=zeros(n,1);
for jj = 1:length(act_freq)
    
    x = x + 1/sqrt(n)*a0(jj)*exp(1i*2*pi*act_freq(jj)*S);
    Atrue = [Atrue exp(1i*2*pi*act_freq(jj)*S)];
    
end
sig = .01;  % noise std deviation
xtrue = x;

x = x+sig*randn(n,1);  % observations

epsilon = 1e-3; % "merge frequencies in this window"

%% Initialization for the Algorithm

init_grid = linspace(minfreq,maxfreq,10000)';  % initial grid, uniformly spaced
f0=init_grid(1);  % random initial starting point
tau = norm(a0,1) + sig;   % regularization parameter (clairvoyant)

At = 1/sqrt(n)*exp(1i*2*pi*f0*S);
% ^ initial atom
coeft=At\x;
coeft = real(coeft);
xp=At*coeft;

maxiter = 100;
bcount = 0;
fbins= f0;


%% ALGORITHM PARAMETERS
maxiter_gp = 10;  
tol_gp = 1e-2;
eta = 0.8;
verbose = 1;  % print messages and diagnostics
backcount = 20;  % to check convergence
tol = 1e-10 ; % to check convergence
do_back = 1;
cg = false;       % if on, run only conditional gradient
if cg==true
    do_back = 0;
end
grid=init_grid;
drops = 0;
%% begin iterations

for ii = 2:maxiter
    
    if verbose==true
        disp(ii)
    end
    
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

    fbins_old = fbins;
    
    fbins = [fbins; grid(f)]; %append the new selected frequency
    
    
    a = 1/sqrt(n)*exp(1i*2*pi*grid(f)*S);  % new atom
    
    At = [At a];   % the current atomic set
    coeft_old = coeft;
    
    if cg==true
        phinext=(tau*a-xp);
        gamma=-((x-xp)'*phinext)/(phinext'*phinext);
        coeft = [(1-gamma)*coeft;gamma*tau];
    else
        coeft = [coeft; 0];
    end

    [coeft,~]=gp_reduced_complex(At,x,coeft,tau,maxiter_gp,tol_gp);
    
    xp=At*coeft;
    newf = 0.5*norm(xp - x)^2;
    err(ii) = newf;
    
    
    %%%%%%%% BACKWARD STEP %%%%%%%%%%%%%
    if do_back == 1
        
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
            bcount = bcount+1;
            At = Anew;
            xp = xnew;
            coeft = coeft_new;
            fbins = fnew;
            err(ii) = newf_rem;
            drops = 1+drops;
            if verbose==true
                fprintf(1,'\n iter after backward step: %4d, f: %10.5e \n', ii,newf_rem);
            end
         end
    end
    
    if verbose==true
        stem(act_freq,a0); hold on; stem(fbins',(coeft),'r .');
        hold off
        pause(.05);
    end
    
    % convergence check
    if ((err(ii-1)-err(ii) < tol) && (ii > 2) && (err(ii-1)>= err(ii)))
        if verbose==true
            fprintf(2,'\n Convergence achieved \n');
        end
        break;
    end
    if drops >= backcount
        if verbose==true
            fprintf(2,'\n Dropcount achieved \n');
        end 
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

is1 = find(act_freq == 1);
is0 = find(act_freq == 0);

if ~isempty(is1) && ~isempty(is0)
    act_freq(is0) = [];
    a0(is1) = a0(is1) + a0(is0);
    a0(is0) = [];
end

% merge coefficients epsilon close to each other
[f_final,coeft_final] = refine_frequencies(f_final,coeft_final,epsilon);
A_final=[];
for f = 1:length(f_final)
   A_final = [A_final 1/sqrt(n)*exp(1i*2*pi*f_final(f)*S)]; 
end
[coeft_final,~]=gp_reduced_complex(A_final,x,zeros(length(coeft_final),1),tau,maxiter_gp,tol_gp);

% display result
close all
stem(act_freq,a0); hold on; plot(f_final,coeft_final,'r .');
xlabel('Frequency','FontSize',22);
ylabel('a0','FontSize',22)
legend('True','Recovered')

