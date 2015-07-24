function [next_atom] = find_next_atom_soslasso(gradf,groups,L1)
% INPUTS:
% gradf  = gradient of f at current iterate
% groups = groupings over the variables
% NORM : \sum ||xg||_2 * sqrt(|G|) + L1*||xg||_1
% OUTPUT:
% next_atom = next atom

gradf = -gradf;
p = length(gradf);


temp  = (sqrt(sum((repmat(gradf,1,size(groups,2)).*groups).^2)))./(sqrt(sum(groups)));
[~,bestgrp] = max(temp);

inds = find(groups(:,bestgrp)==1); % this is the group that maximizes norm

next_atom = zeros(p,1);
g = gradf(inds);
sg = length(g);
cvx_begin quiet
variable a(sg,1)
maximize sum(a.*g)
subject to
L1*sum(abs(a)) + sqrt(sg)*norm(a) <= 1;
cvx_end

% a = sos_atompick(-g,L1);


next_atom(inds) = a;
next_atom = next_atom/norm(next_atom);

end

% function W = sos_atompick(x,T)
% 
% sg = length(x);
% Wz = zeros(sg,1);
% t = 1;
% t_old = 0;
% iter = 0;
% gamma = 1;
% gamma_inc = 2;
% Wz_old = Wz;
% bFlag = 0;
% 
% L = Wz'*x;
% R = norm(Wz)*sqrt(sg) + T*norm(Wz,1);
% Fnew = L+R;
% 
% maxiter =100;
% tol = 1e-4;
% while iter<maxiter
%     %tic;
%     Fold = Fnew;
%     alpha = (t_old - 1) /t;
%     
%     Ws = (1 + alpha) * Wz - alpha * Wz_old;
%     
%     % compute current gradient and function value
%     Fs = Ws'*x;
%     gWs = Ws;
%     
%     % the Armijo Goldstein line search scheme
%     numtries = 1;
%     while true
%         
%         %%%%%%%%% PROX PART
%         temp = Ws - gWs/gamma;
%         temp = sign(temp).*max(abs(temp)-T/gamma,0);
%         if norm(temp)>0
%            temp = (temp/norm(temp))*max(0,norm(temp)-sqrt(sg)/gamma);
%         end
%         Wzp = temp;
%         %%%%%%%%%%%
%         
%         
%         Fzp = Wzp'*x;
%         
%         delta_Wzp = Wzp - Ws;
%         nrm_delta_Wzp = norm(delta_Wzp)^2;
%         r_sum = (nrm_delta_Wzp)/2;
%         
%         Fzp_gamma = Fs + sum(sum(delta_Wzp .* gWs))...
%             + gamma/2 * nrm_delta_Wzp;
%         
%         if (r_sum <=1e-20)
%             bFlag=1; % this shows that, the gradient step makes little improvement
%             break;
%         end
%         
%         if (Fzp <= Fzp_gamma)
%             if numtries > 1
%                 break;
%             else
%                 numtries = numtries + 1;
%                 gamma = gamma/gamma_inc;
%             end
%         else
%             gamma = gamma * gamma_inc;
%             numtries = numtries + 1;
%         end
%     end
%     
%     % update the variables
%     Wz_old = Wz;
%     Wz = Wzp;
%     
%     Fnew = Fzp + norm(Wzp)*sqrt(sg) + T*norm(Wzp,1);
%     
%     if (bFlag==true)
%         break;
%         % terminate if gradient change is small
%     end
%     
%     if iter>=2
%         if ( abs(Fold - Fnew)  <=...
%                 tol* Fold)
%             break;
%         end
%     end
%     
%     % update parameters
%     iter = iter + 1;
%     t_old = t;
%     t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
%     
%     
% end
% W = Wzp;
% 
% end
