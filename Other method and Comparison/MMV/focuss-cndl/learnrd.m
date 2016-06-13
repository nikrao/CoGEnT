function [rmse, A, Aorig, x, xorig, diversity, btrain, x100norm] = learnrd(params)
% learnrd -- Learn a random dictionary, see Section 3.2 Kreutz 2000 July
%
%
% [rmse, A, Aorig, x, xorig, diversity, btrain, x100norm] = learnrd(params)
%
% JFM   7/20/2000
% Rev:  1/7/2004

%!rm block

rows = params.rows;      % Size of A
cols = params.cols;
numtrain = params.numtrain;
div = params.div;        % Number of active elements in each x

% Seed random number generators
randn('state',sum(100*clock));
rand('state',sum(100*clock));

A0 = randn(rows, cols);

if(params.generatedata == 0)
    % Load the inital conditions and data from a matlab workspace
    load(params.workspace);
    
    % Initial conditions for A
    if(strcmp(params.Ainit,'random') == 1)
        Ain = A0;                 % Random
    elseif(strcmp(params.Ainit,'Aorig') == 1)        
        Ain = Aorig;              % Original A, no dictionary learning
    elseif(strcmp(params.Ainit,'btrain') == 1)
        Ain = btrain(:, 1:cols);  % First columns of data vector    
    elseif(strcmp(params.Ainit,'randblob') == 1)
        % Use Lewicki's blob initialization for A
        Ain = initApatch('randblob',rows, cols); 
    else
        disp('Invalid initial condition for A matrix (params.Ainit)');
        return;
    end
    
    % Initial condition for x vectors
    % !! Some of these won't work with overcomplete dictionaries, try pinv
    % for overcomplete
    [xrows, xcols] = size(xorig);
    if(strcmp(params.xinit,'random') == 1)
        x = randn(xrows, xcols);
    elseif(strcmp(params.xinit,'rand01') == 1)
        x = rand(xrows, xcols);        
    elseif(strcmp(params.xinit,'pinv') == 1)
        x = [];
    elseif(strcmp(params.xinit,'zeros') == 1)
        x = zeros(xrows, xcols);
    else
        disp('Invalid initial condition for x matrix (params.xinit)');
        return;       
    end
    
else        
    % Generate data
    [Aorig btrain xorig] = createdata(numtrain, div, rows, cols, params.SNRdB, params.divrand);

    % Initial conditions for A
    if(strcmp(params.Ainit,'random') == 1)
        Ain = A0;                   % Random
    elseif(strcmp(params.Ainit,'Aorig') == 1)        
        Ain = Aorig;                % Original A, no dictionary learning
    elseif(strcmp(params.Ainit,'btrain') == 1)
       % Ain = btrain(:, 1:cols);    % First columns of data vector    
       Ain = btrain(:, 1:cols) + 0.1 * A0;
    elseif(strcmp(params.Ainit,'blob') == 1)
        % Use Lewicki's blob initialization for A
        Ain = initApatch('blob',rows, cols);     
    else
        disp('Invalid initial condition for A matrix');
        return;
    end
    
    % Initial condition for x vectors
    [xrows, xcols] = size(xorig);
    if(strcmp(params.xinit,'random') == 1)
        x = randn(xrows, xcols);
    elseif(strcmp(params.xinit,'rand01') == 1)
        x = rand(xrows, xcols);
    elseif(strcmp(params.xinit,'pinv') == 1)
        x = [];    
    elseif(strcmp(params.xinit,'zeros') == 1)
        x = zeros(xrows, xcols);
    else
        disp('Invalid initial condition for x matrix (params.xinit)');
        return;       
    end
    
    % Linear combination of training vectors
    %Ain = zeros(rows, cols);
    %for i = 1:cols
    %    for j = 1:div
    %        index = fix(rand * numtrain) + 1;
    %        Ain(:, i) = Ain(:, i) + btrain(:, index);
    %    end
    %end
    
    %Ain = (Aorig + .2*A0)/(norm(Aorig + .2*A0, 'fro'));    
    
    %load learned113; %Ain = Aorig;
    %save learned113 Ain Aorig xorig btrain x   
    %save focuss5 Aorig btrain xorig    
    %save icadata20x20div8 Aorig btrain xorig
end

if(strcmp(params.Ainitnorm,'Frob'))
    % Normalize by Frobenius norm
    Ain = Ain / norm(Ain, 'fro');
    
    disp('Ainit norm: Frob');
elseif(strcmp(params.Ainitnorm,'col_Frob'))
    % Normalize column-wise to unit Frob norm(like Engan)
    for i = 1:cols
        nrm = norm(Ain(:, i) );
     %   disp(sprintf('Col %d, norm = %f', i, nrm));
        Ain(:, i) = Ain(:, i) / (sqrt(cols) * nrm);
    end
    
    disp('Ainit norm: Column norm to ||A||_Frob = 1');
elseif(strcmp(params.Ainitnorm,'unit_col'))
    % Normalize to unit column-norm 
    for i = 1:cols
        nrm = norm(Ain(:, i) );
        Ain(:, i) = Ain(:, i) / nrm;
    end
    
    disp('Ainit norm: Unit column norm');
else
    disp('Ainit norm: No normalization');
end


if(strcmp(params.algorithm,'focuss-dl') == 1)
    [rmse, A, x, diversity, x100norm] = trainrd_focussdl(btrain, x, params, Ain, div, Aorig, xorig); 
elseif(strcmp(params.algorithm,'focuss-cndl') == 1)
    [rmse, A, x, diversity, x100norm] = trainrd_focusscndl(btrain, x, params, Ain, div, Aorig, xorig); 
elseif(strcmp(params.algorithm,'focuss-cndl2') == 1)
    [rmse, A, x, diversity, x100norm] = trainrd_focusscndl2(btrain, x, params, Ain, div, Aorig, xorig); 
elseif(strcmp(params.algorithm,'palmer') == 1)
    [rmse, A, x, diversity, x100norm] = trainrd_palmer(btrain, x, params, Ain, div, Aorig, xorig); 
elseif(strcmp(params.algorithm,'mmp') == 1)
    [rmse, A, x, diversity, x100norm] = trainrd_mmp(btrain, x, params, Ain, div, Aorig, xorig);
elseif(strcmp(params.algorithm,'extica') == 1)
    [rmse, A, x, diversity, x100norm] = trainrd_extica(btrain, x, params, Ain, div, Aorig, xorig);
elseif(strcmp(params.algorithm,'fastica') == 1)
    [rmse, A, x, diversity, x100norm] = trainrd_fastica(btrain, x, params, Ain, div, Aorig, xorig);    
else
    disp('Invalid vector selection method (params.vectorselection)');
end


return;
