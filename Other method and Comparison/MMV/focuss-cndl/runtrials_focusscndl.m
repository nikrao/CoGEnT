function [paramarray, perfarray] = runtrials_focusscndl
% runtrials_focusscndl -- Runs trials of FOCUSS-CNDL learning algorithm (Kreutz:2003)
%
% [paramarray, perfarray] = runtrials_focusscndl
%
% JFM   8/28/2003
% Rev:  11/29/2004


% Experimental setup
numtrials = 1;
writeoutput = 1;
filenamebase = 'trial';
params.trialnum = 1;        % Trial number
params.outputfilename = sprintf('%s%d', filenamebase, params.trialnum);


% Data generation paramters
params.generatedata = 1;    % Generate random data of size given below
%params.workspace = 'icadata20x20div8.mat';  % If generatedata = 0, load from from .mat workspace : A Aorig xorig btrain x
%params.workspace = 'learned67.mat';
params.rows = 20;           % Size of A
params.cols = 30;
params.numtrain = 1000;     % Number of training vectors (x_k)
params.div = 7;             % Number of active elements in each x_k
params.divrand = 0;         % Randomness in div (can choose fewer than div)
params.SNRdB = 0;           % SNR of generated data (Btrain), 0 for no noise
params.whiten = 0;          % Whiten/sphere data 


% ------ Algorithm parameters ------
params.algorithm = 'focuss-cndl';  % Algorithm
params.maxiter = 500;       % Not pased to ext_ica.m for now
params.itersize = 1000;      % Number of vectors to update at each iteration
params.saveatiter = 10;
params.p = 1.0;             % Choice of diversity measure
params.gammaA = 1.0;        % Learning rate for A matrix (see eq. 43) (was 1.0 in Kreutz:2003) (in Palmer:2003, this is alpha)
params.alpha = 1.0;         % Regularization for x's (see eq. 28)
params.mu = 1.0;            % Balance between x_k and x_{k+1} updates
params.lambdamax = 2.0e-3;  % Regularization parameter limit

params.Ainit = 'random';    % Initial condition for A, 'Aorig', 'random', 'btrain', 'randblob'
params.Ainitnorm = 'col_Frob'; % Normalization for A, 'unit_col', 'col_Frob'
params.Asavehistory = 0;    % Save A at iterations (0 not to save)
params.Aupdateiter = 200;   % Number of x_k's to FOCUSS before updating A (if Aintupdate=1)

params.xinit = 'pinv';      % Pseudoinverse solution for first iteration, 'pinv', 'random'
params.xiter = 1;           % Number of times to update x vectors after each update of A matrix, for algorithm in Kreutz 2000, this should be 1.
params.focussconverge = 1e-6;  % Focuss convergence limit for update of x_k

params.reinit = 1;          % Reinitialize x_k's if the diversity is too high after many iterations
params.reinititer = 175;
params.reinit_type = 'focuss';  % Reinit with x = 'pinv',  'rand', 'focuss'

params.sparsifyiter = 1;    % Select the highest valued elements every # iterations
params.positive_code = 0; 
params.limit_x = 0;          % 1 (true) to restrict x_k to [limit_x_min, limit_x_max]
params.limit_x_min = -10.0;
params.limit_x_max = 10.0;


% Internal algorithm performance display preferences
params.compareA = 1;        % Compare the columns of the learned A matrix to original
params.timeticks = 0;       % 1 to display '.' after a certain amount of time
params.display = 0;         % Display plots during training (works best on UNIX)
% ------ End algorithm parameters ------


% Performance measurement parameters
params.Atolerance = 0.01;   % Tolerance of matched columns of A 
params.Xtolerance = 0.05;   % Tolerance of learned solutions X


% Setup repeated trials
 
paramarray(1) = params;

for i = 2:20
    params.trialnum = params.trialnum + 1;
    params.outputfilename = sprintf('%s%d', filenamebase, params.trialnum); 
    paramarray(i) = params;
end

%params.p = 0.8;
%params.reinit = 1;
for i = 21:40
    params.trialnum = params.trialnum + 1;
    params.outputfilename = sprintf('%s%d', filenamebase, params.trialnum); 
    paramarray(i) = params;
end



% Initialize performance struct
perf.nummatch = 0;
perf.rmse = 0;
perf.numxmatch = 0;
perf.sparsity = 0;
perf.time = 0;
perf.cputime = 0;
perf.x100norm = 0;
perf.x100dist = 0;
perf.xsolved = 0;
perf.xmaxdist = 0;
perf.snr = 0;
perf.snrstdev = 0;
perfarray(1:numtrials) = perf;


% Run the trials
for i = 1:numtrials
    
    c1 = fix(clock);    cpu1 = cputime;
    % test_generate ?? maybe on cairo?
    %[rmse, A, Aorig, x, xorig, diversity, btrain] = test_generate(paramarray(i));
    [rmse, A, Aorig, x, xorig, diversity, btrain, x100norm] = learnrd(paramarray(i));
    
    c2 = fix(clock);     cpu2 = cputime;
    
    perfarray(i).cputime = cpu2 - cpu1;
    perfarray(i).time = etime(c2, c1);
    perfarray(i).rmse = rmse(paramarray(i).maxiter);
    perfarray(i).sparsity = mean( diversity(paramarray(i).maxiter, :) );
   % perfarray(i).x100norm = x100norm(paramarray(i).maxiter);
    
    [num, mindist, normf1, matchedframe, index] = compareframe(Aorig, A, paramarray(i).Atolerance);
    perfarray(i).nummatch = num;
    [dist, nummatch] = comparesolution(xorig, x, index, paramarray(i).Xtolerance);
    perfarray(i).xsolved = nummatch;    
    perfarray(i).x100dist = dist(100);
    perfarray(i).xmaxdist = max(dist);
    
    % Calc SNR on square matricies
    if(size(A,1) == size(A,2))
        [xreorder] = reorderica(x, A*Aorig);
        [perfarray(i).snr, perfarray(i).snrstdev] = calcsnr(xreorder', xorig');
    else
        perfarray(i).snr = -1;
        perfarray(i).snrstdev = -1;
    end
        
    disp(sprintf('CPU time:   %d', perfarray(i).cputime) );
    
    %  plotresults(paramarray(i).outputfilename, paramarray(i).outputfilename);
    
    % Write results to file
    if(writeoutput == 1)
        [fid, message] = fopen(sprintf('%sresults.out', filenamebase), 'a');
        if(fid == -1)
            disp(message);
            return;
        end
    
        fprintf(fid, '%d\t%d\t%f\t%f\t%d\t%f\t%f\t%f\t', ...
            paramarray(i).trialnum, perfarray(i).nummatch, perfarray(i).sparsity, ...
            perfarray(i).rmse, perfarray(i).xsolved, perfarray(i).xmaxdist, ...
            perfarray(i).x100dist, perfarray(i).x100norm );
        
        fprintf(fid, '%8.2f\t%d\t%d\t%d\t%d\t%d\t', ...
            perfarray(i).cputime / 60,     ...
            paramarray(i).rows, paramarray(i).cols, paramarray(i).numtrain, ...
            paramarray(i).div, paramarray(i).maxiter);
    
        fprintf(fid, '%8.2f\t%8.2f\t\n', ...
            perfarray(i).snr, perfarray(i).snrstdev);
        
        fclose(fid);
    end
end


% Print results to screen 
for i = 1:numtrials
    disp(sprintf('----- Output file:  %s -----', paramarray(i).outputfilename));
    disp(sprintf('Time used:  %8.2f minutes', perfarray(i).time / 60) );
    disp(sprintf('CPU time:   %8.2f minutes', perfarray(i).cputime / 60) );
    disp(sprintf('Num match:  %d', perfarray(i).nummatch) );
    disp(sprintf('Sparsity:   %f', perfarray(i).sparsity) );
    disp(sprintf('Rmse:       %f', perfarray(i).rmse) );
    
end

            
            
        