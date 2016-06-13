function [U,V,t,e,y,stats, err2] = sdplrmc(m, n, I, J, ent, r, translate, noise, sfm, ...
pars_in, U0, V0, t0, init_penalty, ent_gt)

% Last Change: Sun Mar 07 12:00 PM 2010 C
% Author: Sameer Sheorey, TTI-C
%
% [U,V,t,e] = SDPLRMC(m, n, I, J, ent, r, translate, noise)
%
% Low rank matrix completion with semidefinite programming. 
% Uses sdplr (Burer et al) to solve SDP with low rank constraints.
% http://dollar.biz.uiowa.edu/~sburer/www/doku.php?id=software#sdplr
%
% Given: Entries (I,J) \in \Omega of the matrix M specified in V. (This is the
% output of find if unknown entries are zero). M is of size m rows and n
% columns and is known to be of rank r or less. If noise is expected in the
% entries, the parameter noise is nonzero and specifies a weight for errors.
% If M = M_ + t*one', where t is a translation vector, then the factors of M_
% and t can be recovered by setting translate=true. r is now the rank of M_.
%
% The missing entries are found by minimizing ||M_||_*, the nuclear norm of M_.
% Suppose M_ = U*V' be the SVD of M_, with the singular values absorbed
% equally into the two factors. Let R = [U; V]. U, V and R all have r columns. 
% Then X = R*R' is symmetric positive semidefinite (PSD) and \tr X = ||M_||_*
% Further, X can be written as a 2x2 block matrix [U*U' M_; M_' V*V'].
%
% No translation, no noise :
%
% min \tr X  s.t. X is symmetric PSD,
%                 X_{i,j+n} = M_{i,j} for (i,j) \in \Omega
%
% Translation, but no noise :
%
% min \tr X  s.t. X is symmetric PSD,
%                 X_{i,j+n} = M_{i,j} - t1_i + t2_i  for (i,j) \in \Omega
%                 t1, t2 >= 0     (Limitation of sdplr)
% Note that there is no unique way of expressing a rank r+1 matrix as a sum of a
% rank r matrix and a translation. This formulation yields the min nuclear
% norm solution. Other solutions are possible by adding a penalty on the L_1
% norm of the translation: tweight * \sum_i t1_i + t2_i
%                 
% Both translation and noise :
%
% L_1 error:
% min \tr X + \lambda \sum_\Omega e1_{i,j} + e2_{i,j}                
%            s.t. X is symmetric PSD,
%                 X_{i,j+n} - e1_{i,j} + e2_{i,j} = M_{i,j} - t1_i + t2_i  
%                 for (i,j) \in \Omega
%                 t1, t2, e1, e2 >= 0
%
% In this case, e_{i,j} = e1_{i,j} - e2_{i,j} is the error in the low rank +
% translation model. min e1+e2 for a fixed e = e1-e2 with e1,e2>=0 is attained
% when e1+e2 = |e|. So the second term in the objective is the L_1 error. This
% should be more robust than an L_2 error.
%
% L_2 error:
% min \tr X + \lambda (\tr E - 1)
%           s.t. X is symmetric PSD,
%                X_{i,j+n} + e(l) = M_{i,j} - t_i
%                for (i,j) \in \Omega
%                E is rank 1 symmetric PSD,
%                E = [1; e] * [1; e]'
%
tweight = 0; % Weight for t1+t2 in objective when translates are 
        % given by t = t1-t2; t1, t2>=0
        % tweight = 0 gives min nuclear norm untranslated matrix. Other values of tweight
        % will give different solutions.

if ~exist('translate', 'var') || isempty(translate)
  translate = false;
end 
if ~exist('noise', 'var') || isempty(noise)
  noise = {0,0};
elseif ~iscell(noise),  % Default L_2 error
  noise = {2, noise};   % / sqrt(m+n);    % Normalization
end 

assert(any(noise{1} == [0 1 2]), ...
'Error type should be 0 (none), 1 (L1) or 2 (L2).\n');
assert(~islogical(noise) || noise<0, ['Noise should be a positive number, '...
'not logical\n']);

k = length(ent);
if(exist('r', 'var') && ~isempty(r))
  dof = r*(m+n-r);
  if(translate),      % This is probably an overestimate
    dof = dof + m;
  end
  if(k<dof)
    warning('Too few (%d) known entries. Atleast %d needed.', k, dof);
  end
  %assert(k>=dof, 'Too few (%d) known entries. Atleast %d needed.\n', k, dof);
else
  r = [];
  warning('No rank specified.');
end

if ~exist('sfm', 'var')
  sfm = [];
else
  assert(mod(m,2)==0, ['Joint measurement matrix must have an even number'...
 ' of rows (2D point coordinates)']);  
 assert(m>=4, 'Atleast 2 images needed');
end

I = I(:);
J = J(:);
ent = ent(:);

fprintf('%d x %d matrix of rank %d with %d known entries.', m, n, r, k);

%global Atr c K pars

K.s = m+n;
K.l = 0;
if isempty(sfm)
  c = reshape(speye(m+n), [], 1);
else
  c = sparse((m+n)^2, 1);       % All zeros
end
cols = sub2ind([m+n m+n], J+m, I); % sdplr expects lower triangular entires only
rows = (1:k)';
vals = 0.5*ones(k,1);
ncons = k;
% [X e1 e2 t1 t2]


% if(translate)
%     rows = [rows' 1:k ncons+1]';
%     ncons = ncons + 1;
%     cols = [cols; sum(K.s.^2) + K.l + 1 + I; sum(K.s.^2) + K.l + 1];
%     vals = [vals; ones(k+1,1)];
%     ent = [ent; 1];
%     K.s = [K.s m+1];
%     if (tweight==0)
%       c = [c; sparse((m+1)^2,1)];
%     else
%      c = [c; tweight *reshape(speye(m+1)-sparse(1,1,1,m+1,m+1), [], 1)];
%    end
% end

if(noise{1}==1)         % L_1 error
  c = [c; noise{2}*ones(2*k,1)];
  rows = [rows' 1:k 1:k]';
  cols = [cols; sum(K.s.^2) + K.l + (1:2*k)'];
  vals = [vals; -ones(k,1); ones(k,1)];
  K.l =  K.l + 2*k;
end

if(noise{1}==2)         % L_2 error
  c = [c; noise{2}*(reshape(speye(k+1), [], 1))];
  rows = [rows' 1:k ncons+1]';
  ncons = ncons + 1;
  cols = [cols; sum(K.s.^2) + K.l + (1:k)' + 1; sum(K.s.^2) + K.l + 1];
  vals = [vals; ones(k+1,1)];
  ent = [ent; 1];
  K.s = [K.s k+1];
end

if(translate)
  rows = [rows' 1:k 1:k]';
  vals = [vals; -ones(k,1); ones(k,1)];
  if (tweight==0)
    c = [c; sparse(2*m,1)];
  else
    c = [c; tweight * ones(2*m,1)];
  end
  cols = [cols; sum(K.s.^2) + K.l + [2*I-1; 2*I]];
  K.l =  K.l + 2*m;
end


% SfM (Orthographic projection) (rank 3):
% No longer find SVD by minimizing \tr X.
% Extra constraints: U is camera matrix with orthogonal rows, so upper left
% block of X looks like this :
% [1      ...
% [0 1    ...
% [? ? 1  ...
% [? ? 0 1...]
%

sfm_pen_fac = 1;
if(strcmpi(sfm,'ortho'))
  fprintf('Orthographic projection\n');
  rows = [rows' ncons+(1:3*m/2)]';
  ncons = ncons + 3*m/2;
  cols = [cols; sub2ind([m+n m+n], [1:m 2:2:m], [1:m 1:2:m])'];
  % Diagonals = 1, off diagonal = 0
  vals = [vals; sfm_pen_fac*ones(3*m/2,1)];
  ent = [ent; sfm_pen_fac*ones(m, 1); zeros(m/2, 1)];   % More matrix entries known
end

% SfM (Scaled Orthographic, or weak perspective) (rank 3)
% Upper left Block of X looks like this :
% [1        ...
% [0 1      ...
% [? ? s2   ...
% [? ? 0 s2 ...
% This means X_(2i-1,2i-1) = X(2i,2i)
% Note that the scale of the first camera is fixed to 1.

if(strcmpi(sfm,'sortho'))   % Diagonals (equal) then off-diags (zero)
  fprintf('Scaled Orthographic projection\n');
  nextrows = [1 2 2+(1:m/2-1) 2+(1:m/2-1) m/2+1+(1:m/2)]';
  rows = [rows; nextrows + ncons];
  ncons = ncons + m+1;
  % diags (all 1s, then all -1s) then off diags (0)
  cols = [cols; sub2ind([m+n m+n], [1 2 3:2:m 4:2:m 2:2:m], ...
  [1 2 3:2:m 4:2:m 1:2:m])'];
  % Diagonals = 1, off diagonal = 0
  vals = [vals; 1; 1; ones(m/2-1,1); -ones(m/2-1,1); ones(m/2,1)];
  ent = [ent; 1; 1; zeros(m, 1)]; 
end

if(strcmpi(sfm, 'affine'))
  % No other constraints
end

% A is transposed to save space

Atr = sparse(cols, rows, vals, sum(K.s.^2) + K.l, ncons, length(vals));
whos c
whos Atr
clear cols rows vals
%figure;imshow(full(Atr)', [-1 1], 'InitialMagnification', 'fit');
%keyboard

pars.feastol = 1.0e-5; % default = 1.0e-5
pars.centol = 0.1; %Desired level of accuracy in intermediate calculations (default = 1.0e-1)
pars.dir = 1;               % Truncated Newton
pars.penfac = 2;        % Default = 2
pars.numlbfgs = 4;          % Default = 4
if isempty(r)
  pars.ranktol = 1e-3;        % Drop SVs with lower values
else
  pars.forcerank = r;
%   if (translate)
%     pars.forcerank = [pars.forcerank 1];
%   end
  if noise{1}==2
    pars.forcerank = [pars.forcerank 1];
  end
  if(K.l>0)
    pars.forcerank = [pars.forcerank 1];
  else
    K = rmfield(K, 'l');
  end
end
pars.soln_factored = 1;     % Only want factors, not whole matrix
pars.limit = 120;          % Time limit (seconds)
pars.reduce = 0;            % Reduce problem dim with iterations. Default = 0
pars.writeraw = 0;

% Override with supplied parameters, if any
if exist('pars_in', 'var') && isstruct(pars_in)
  pt = fieldnames(pars_in);
  for a=1:length(pt)
    pars.(pt{a}) = pars_in.(pt{a});
  end
end


% Initial solutions
if exist('init_penalty', 'var') && ~isempty(init_penalty)
  info0.penalty = init_penalty;
end
r0 = {};
if exist('U0', 'var') && ~isempty(U0)
  r0{1} = [U0; V0];
  next = 1; nextcell = 2;
  if(noise{1}~=0)
      if(noise{1}==1)
          r0{nextcell}(1:2*k,1) = randn(2*k,1);   % Small errors =? initialize with zeros
          next = 2*k+1;
      else
          r0{nextcell}(1,1) = 1;
          r0{nextcell}(2:k+1,1) = randn(k,1);
          nextcell = nextcell + 1;
      end
  end
  if exist('t0', 'var') && ~isempty(t0)
    r0{nextcell}(next:next+2*m-1,1) = t0;
  end
  if isempty(init_penalty)
    [x,y,stats] = sdplr(Atr, ent, c, K, pars, [], [], [], r0);
  else
    [x,y,stats] = sdplr(Atr, ent, c, K, pars, [], [], [], info0, r0);
  end
else
   %keyboard();
%  x = sdplrlbls(Atr, ent, c, K, [], pars); x = {x};

  if isempty(init_penalty)
    [x,y,stats] = sdplr(Atr, ent, c, K, pars, [], [], []);
  else
    [x,y,stats] = sdplr(Atr, ent, c, K, pars, [], [], [], info0);
  end
end

%keyboard;
U = x{1}(1:m,:);
V = x{1}(m+1:m+n,:);
next = 1; nextcell=2;
% if(translate)
%   t = x{nextcell}(2:m+1);
%   nextcell = nextcell + 1;
% end
if(noise{1}~=0)
  if(noise{1}==1)
    e = x{nextcell}(next:next+2*k-1).^2;
    e = e(k+1:2*k) - e(1:k);
    next = 2*k+1;
  else
    e = x{nextcell}(2:k+1);
    nextcell = nextcell + 1;
  end
%  e = sparse(I, J, e, m, n);
else
    e=0;
%  e = sparse(m,n);  % All zero sparse
end
if(translate)
  t = x{nextcell}(next:next+2*m-1).^2;
  t = t(2:2:end)-t(1:2:end);
else
  t = sparse(m,1);  % All zero sparse
end

% Compute error measures (rms, mean, max)
err2=0; err1=0; errm=0;
for a=1:k
  d = abs( U(I(a),:) * V(J(a), :)' + t(I(a)) - ent_gt(a) );
  err1 = err1 + d;
  if(errm < d)
    errm = d;
  end
  err2 = err2 + d^2;
end
err2 = sqrt(err2/k);
err1 = err1/k;
fprintf('Reprojection error on Omega (rms): mean = %g, max = %g (pixels)\n', err2, errm);
fprintf('Reprojection error on Omega (rms) using error variable e: mean = %g, max = %g (pixels)\n', norm(e)/sqrt(k), max(abs(e)));
fprintf('Relative error on Omega is %g\n', sqrt(k)*err2/norm(ent,2));
fprintf('Relative error on Omega, based on e, is %g\n', norm(e,2)/norm(ent,2));

% DEBUG
%X = U*V' + e + repmat(t, 1, n);
%rms = norm(X(sub2ind([m, n], I, J))-ent)/sqrt(k);
%x{1}
%full(e)
%if(rms > 1e-3)
%  keyboard;
%end
%keyboard