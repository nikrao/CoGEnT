%%% TEST CoGEnT for sparse + low rank demixing 
clear;
clc;

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

Ainit1 = [1; zeros(p-1,1)];
selfun1 = @(gradf) find_next_atom_l1(gradf);
tauset1 = linspace(10,100,10);


Ainit2 = randn(matsize(1),1)*randn(1,matsize(2));
Ainit2 = Ainit2(:)/norm(Ainit2(:));
selfun2 = @(gradf) find_next_atom_nucnorm(gradf,matsize(1),matsize(2));
tauset2 = linspace(1,100,10);


Phi = speye(p);



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = 0;
for tau1 = tauset1
    t1 = 1+t1;
    t2 = 0;
    for tau2 = tauset2
        t2 = 1+t2;
        
        [x1,x2,At1, At2 iter, obj, time, back_count] = CoGEnT_Demix_lowrank4(y, Phi, tau1, tau2,...
            Ainit1, Ainit2, selfun1, [matsize],...
            'tol',1e-4,...
            'verbose',0,...
            'maxiter',200,...
            'dropcount',100,...
            'debias',1);

        
        errs(t1,t2) = norm(y - x1 - x2)^2;
        fprintf('.');
    end
    
    fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
          