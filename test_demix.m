
%%
clear;
clc;
p = 1000;
D = overDCTdict(p,p);
fprintf('\n')
k = 10;
actdct = randsample(p,k);
sdct = zeros(p,1);
sdct(actdct) = 2*rand(k,1)-1;
xdct = inv(D)*sdct;
actl1 = randsample(p,k);
xl1 = zeros(p,1);
xl1(actl1) = 2*rand(k,1)-1;


%  observations
Phi = sparse(eye(p));
% m = 500;Phi = randn(m,p)/sqrt(m);
y = Phi*(xl1 + xdct);

selfun2 = @(gradf) find_next_atom_DCT_D(gradf,D);
selfun1 = @(gradf) find_next_atom_l1(gradf);
Ainit2 = D(:,1);
% Ainit2 = [1; zeros(p-1,1)];
Ainit1 = [1; zeros(p-1,1)];
tauset2 = norm(sdct,1)*linspace(0.8,1.2,10);
tauset1 = norm(xl1,1)*linspace(0.8,1.2,10);
maxiter = 100;

for t1 = 1:length(tauset1);
    tau1 = tauset1(t1);
    
    for t2 = 1:length(tauset2);
        tau2 = tauset2(t2);
        
        [x1,x2,~, ~] = CoGEnT_Demix(y, Phi, tau1, tau2,...
    Ainit1, Ainit2, selfun1, selfun2,...
            'maxiter',maxiter,...
            'debias',1,...
            'verbose',0,...
            'backward',1);
        
        errtau(t1,t2) = norm(x1-xl1) + norm(D*x2-sdct);
        fprintf('.')
    end
    fprintf('\n')
end

[t1 t2] = find(errtau == min(min(errtau)));
tau1 = tauset1(t1(1));tau2 = tauset2(t1(1));
[x1,x2,At1, At2] = CoGEnT_Demix(y, Phi, tau1, tau2,...
            Ainit1, Ainit2, selfun1, selfun2,...
            'maxiter',maxiter,...
            'debias',1,...
            'verbose',0,...
            'backward',1);

        
        %%
close all
subplot(211)
inds = find(sdct~=0);
stem(inds,sdct(inds))
hold on
dd = D*x2;
inds = find(abs(dd)>=0.002);
plot(inds,dd(inds),'r .')
title('DCT recovery')
% legend('True','Recovered')
xlim([1,1000]);
subplot(212)
inds = find(xl1~=0);
stem(inds,xl1(inds))
hold on
inds = find(x1~=0);
plot(inds,x1(inds),'r .')
title('Sparse recovery')
% legend('True','Recovered')
xlim([1,1000]);


