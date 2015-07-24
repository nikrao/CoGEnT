% OSCAR using CoGEnT

clear;
clc;
close all;

% Measurements and true signal
n = 300;
p = 1000;
Phi = randn(n,p);
beta = zeros(p,1);
a = 7*ones(20,1)+ randn(20,1);
b = 9*ones(20,1)+ randn(20,1);
c = -8*ones(20,1)+ randn(20,1);
beta(1:20) = a; beta(301:320)=b; beta(601:620) = c;
noise = 0.05;
y = Phi*beta + noise*randn(n,1);

% cogent parameters
maxiter = 500;
eta = 0.5;
gptol = 1e-3;
gpiter = 50;
gp_forward = 1;


fprintf('... Measurements Obtained. Begin CoGEnT ... \n')

Ainit = [1; zeros(p-1,1)];
tauset = 250:25:450;
cset = 2.^(-10:2:10);

for tind = 1:length(tauset)
    tau = tauset(tind);
    for cind = 1:length(cset);
        c = cset(cind);
        
        
        selfun =@(gradf) find_next_atom_oscar(gradf,c);
        [x, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun, ...
            'maxiter',maxiter,...
            'eta', eta,...
            'gptol', gptol,...
            'gpiter', gpiter,...
            'gp_forward',gp_forward,...
            'stopcrit',1,...
            'debias',1);
        
        fprintf('.');
        
        errt(tind,cind) = norm(x-beta)^2/p;
        
    end
    fprintf('\n');
end

fprintf('done \n')


[t c]  = find(errt == min(min(errt)));

tau = tauset(t(1));
c = cset(c(1));
selfun =@(gradf) find_next_atom_oscar(gradf,c);
        [x, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun, ...
            'maxiter',maxiter,...
            'eta', eta,...
            'gptol', gptol,...
            'gpiter', gpiter,...
            'gp_forward',gp_forward,...
            'stopcrit',1,...
            'debias',1);
 
 stem(beta);
 hold on
 plot(x,'r .');
 legend('Original','Recovered')
 fprintf('\n MSE = %f \n',norm(x-beta)^2/p);
 fprintf('\n Mean Hamming Error = %f \n',sum((x~=0).*(beta~=0))/p);