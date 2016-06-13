function [paramarray, perfarray] = runtrials_image
% runtrails_image -- Runs trials of dictionary learning algorithm
%
% [paramarray, perfarray] = runtrials_image
%
% JFM   1/12/2001
% Rev:  2/19/2004

% ------ Experimental setup ------
numtrials = 1;
writeoutput = 1;
filenamebase = 'd:\project\temp5\image';
params.trialnum = 40;       % Trial number (of first trial)
params.outputfilename = sprintf('%s%d', filenamebase, params.trialnum); 

% ------ Data generation parameters ------
params.generatedata = 0;    % Generate random data of size given below
params.workspace = 'imagedata_lw-exp1_1'; %'imagedata_coil20small_3.mat'; % % If generatedata = 0, load from
                                % from .mat workspace : A Aorig xorig btrain x 
workspace = load(params.workspace, 'patch_width', 'patch_height', 'numsamples');

params.patchwidth = workspace.patch_width;      % Width of image patch
params.patchheight = workspace.patch_height;     % Height of image patch
params.rows = params.patchwidth * params.patchheight;           % Size of A
params.cols = 128;
params.outwidth = 8;        % Number of patches wide the output images are
params.outheight = 8;       % Number of patches height the output images are
params.numtrain = workspace.numsamples;    % Number of training vectors (x_k)
params.div = 20;            % Number of active elements in each x_k
params.whiten = 0;          % Whiten/sphere data (may not be implemented anymore)


% ------ Algorithm parameters ------
params.algorithm = 'focuss-cndl';  % Algorithm (see learnrd.m)
params.maxiter = 499; %150; %500;
params.itersize = params.numtrain; %100;       % Number of vectors to update at each iteration
params.xiter = 1; %50;           % Number of times to update x vectors after each update of A matrix, for algorithm in Kreutz 2000, this should be 1.
params.saveatiter = 5;
params.p = 1.0;             % Choice of diversity measure
params.gammaA = 0.01;        % Learning rate for A matrix (see eq. 35) (in Palmer:2003, this is alpha)
params.alpha = 1.0;         % Regularization for x's (see eq. 28)
params.mu = 1.0;            % Balance between x_k and x_{k+1} updates
params.lambdamax = 2.0e-4;  % Regularization parameter limit
params.Ainit = 'random';    % Initial condition for A, 'Aorig', 'random', 'btrain', 'randblob'
params.Ainitnorm = 'col_Frob'; % Normalization for A, 'unit_col', 'col_Frob'
params.Aupdateiter = 200;   % Number of x_k's to FOCUSS before updating A (if Aintupdate=1)
params.Asavehistory = 1;    % Save A at iterations (0 not to save)
params.xinit = 'pinv';      % Pseudoinverse solution for first FOCUSS iteration, 'pinv', 'random', 'rand01'
params.focussconverge = 1e-6;  % Focuss convergence limit for update of x_k
%params.mmpconverge = 1e-2;  % MMP convergence limit for selection of x_k
params.reinit = 1;          % Reinitialize x_k's if the diversity is too high after many iterations
params.reinititer = 50;
params.reinit_type = 'focuss';  % Reinit with x = 'pinv',  'rand', 'focuss'
params.sparsifyiter = 1;    % Select the highest valued elements every # iterations
params.positive_code = 0;    % Only allow x_k >= 0
params.limit_x = 0;          % Restrict x_k to [limit_x_min, limit_x_max]
params.limit_x_min = -10.0;
params.limit_x_max = 10.0;

% Internal algorithm performance display preferences
params.compareA = 0;        % Compare the columns of the learned A matrix to original
params.display = 0;         % Display plots during training (works best on UNIX)
params.timeticks = 0;       % 1 to display '.' after a certain amount of time
% ------ End algorithm parameters ------


% Performance measurement parameters
params.Atolerance = 0.01;    % Tolerance of matching columns in A
params.Xtolerance = 0.05;

% Setup repeated trials
paramarray(1) = params;

for i = 2:20
    params.p = params.p - 0.1;
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

perfarray(1:numtrials) = perf;


% Run the trials
for i = 1:numtrials
    
    c1 = fix(clock);
    cpu1 = cputime;
    [rmse, A, Aorig, x, xorig, diversity, btrain, x100norm] = learnrd(paramarray(i));
    c2 = fix(clock);
    cpu2 = cputime;
    
    perfarray(i).cputime = cpu2 - cpu1;
    perfarray(i).time = etime(c2, c1);
    perfarray(i).rmse = rmse(paramarray(i).maxiter);
    perfarray(i).sparsity = mean( diversity(paramarray(i).maxiter, :) );
    perfarray(i).x100norm = x100norm(paramarray(i).maxiter);
    
    [num, mindist, normf1, matchedframe, index] = compareframe(Aorig, A, paramarray(i).Atolerance);
    perfarray(i).nummatch = num;
    [dist, nummatch] = comparesolution(xorig, x, index, paramarray(i).Xtolerance);
    perfarray(i).xsolved = nummatch;    
    perfarray(i).x100dist = dist(100);
    perfarray(i).xmaxdist = max(dist);
    
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
        fprintf(fid, '%8.2f\t%d\t%d\t%d\t%d\t%d\t%6.2f\n', ...
        perfarray(i).cputime / 60,     ...
        paramarray(i).rows, paramarray(i).cols, paramarray(i).numtrain, ...
        paramarray(i).div, paramarray(i).maxiter, paramarray(i).p  );
    
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

            
            
        