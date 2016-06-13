% test CoGEnT for demixing astronomical images
clear;
clc;
close all;
warning off
dct = false;
slr = true;
%% DCT + SPARSE
if dct==true
    I = imread('galstar.jpg');
    I = rgb2gray(I);
    I = im2double(I);
    I = imresize(I,[128,128]);
    I = I-mean(mean(I));

    y = I(:);
    p = length(y);
    fprintf('... image loaded and resized ... \n');
    
    fprintf(' making DCT dictionary...')
    D = overDCTdict(p,p);fprintf(' done \n')
    
    
    selfun1 = @(gradf) find_next_atom_DCT_D(gradf,D);
    selfun2 = @(gradf) find_next_atom_l1(gradf);
    Ainit2 = [1; zeros(p-1,1)];
    Ainit1 = D(:,1);
    tauset1 = linspace(1,5000,200); % DCT
    tauset2 = 1:50; %linspace(500,5000,100);
    
    Phi = sparse(p,p);
    for jj = 1:p
        Phi(jj,jj) = 1;
    end
    

    fprintf('... beginning iterations ...  \n')
    errtau = cell(length(tauset1),1);
    if matlabpool('size')==0
        matlabpool open local 10
    end
    parfor t1 = 1:length(tauset1)
        tau1 = tauset1(t1);
        temperr = zeros(1,length(tauset2));
        for t2 = 1:length(tauset2)
            tau2 = tauset2(t2);
            [x1,x2,At1, At2 iter, obj, time, back_count] = CoGEnT_Demix(y, Phi, tau1, tau2,...
                Ainit1, Ainit2, selfun1, selfun2,...
                'maxiter',500,...
                'debias',1,...
                'verbose',0,...
                'tol',1e-4);
           
             temperr(t2) = norm(x1+x2 - y)^2/numel(y);
        end
        fprintf('tau1 = %f \n', tau1);
        errtau{t1} = temperr;
    end
    matlabpool close
    fprintf('done \n');
    
    save demix_DCT errtau tauset1 tauset2 y I

end
%% SPARSE + LOW RANK

if slr == true
    matsize = [64,64];
    img = imread('chess2.jpg');
    img = img(12:end-10,14:end-12,:);
    img = rgb2gray(img);
    img = im2double(img);
    img = imresize(img,[matsize(1),matsize(2)]);
    y = img(:);
    fprintf('... image loaded ...');
    
    % y has to be decomposed as sparse + low rank
    p = length(y);
    selfun1 = @(gradf) find_next_atom_l1(gradf);
    selfun2 = @(gradf) find_next_atom_nucnorm(gradf,matsize(1),matsize(2));
    Ainit1 = [1; zeros(p-1,1)];
    Ainit2 = randn(matsize(1),1)*randn(1,matsize(2));
    Ainit2 = Ainit2(:)/norm(Ainit2(:));
    tauset1 = linspace(0.01,10,50);
    tauset2 = linspace(0.01,10,50);
    Phi = sparse(p,p);
    for jj = 1:p
        Phi(jj,jj) = 1;
    end
    
    fprintf('... beginning iterations ...  \n')
    i1 = 0;
    for tau1 = tauset1
        i1 = 1+i1;
        i2 = 0;
        for tau2 = tauset2
            i2 = 1+i2;
            [x1,x2,At1, At2 iter, obj, time, back_count] = CoGEnT_Demix_lowrank(y, Phi, tau1, tau2,...
                Ainit1, Ainit2, selfun1, [matsize],...
                'tol',1e-4,...
                'verbose',0,...
                'maxiter',1000,...
                'dropcount',100);
            
            fprintf('.');
            
            errtau(i1,i2) = norm(x1+x2 - y)^2/numel(y);
        end
        fprintf('tau1 = %f \n', tau1);
    end
    fprintf('done \n');
end
