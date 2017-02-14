function [x,model_err] = minfunc_frugal(fh,x,LB,UB,options)
% minfunc_frugal - Stochastic spectral gradient with coarse
% wavefield solves
%
% use:
%   [xn,info] = minfunc_frugal(fh,x0,options)
%
% input:
%   fh - function handle to misfit of the form [f,g] = fh(x,I,eps)
%        where I is an indexset of size min(b0 + k*beta,bmax) at iteration
%        k and eps is a tolerance alpha^k*eps0.
%        f is the function value, g is the gradient of the same size
%        as the input vector x. 
%   x0 - initial guess
%
%   options.itermax - max iterations [default 10]
%   options.tol     - tolerance on 2-norm of gradient [1e-6]
%   options.M       - history size [5]
%   options.fid     - file id for output [1]
%   options.write   - save iterates to disk [0]
% 
%
%   options.eps0    - initial tolerance for misfit [1e-4]
%   options.epsmin  - smallest tolerance [1e-4]
%   options.alpha   - controls decrease of tolerance max(epsk = alpha^k*eps0,epsmin) [1]
%   options.b0      - initial batchsize [bmax]
%   options.bmax    - max. batch size   
%   options.beta    - controls growth of batchsize bk = min(b0 + k*beta,bmax) [0]
%   options.seed    - random seed for picking indexset I [1]
%   options.redraw  - redraw indexset at each iteration. [0]
%
% output:
%   xn - final estimate
%
%
% Author: Curt Da Silva, 2016
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if nargin<3
    options = [];
end

% various parameters
if isfield(options,'bmax')
   bmax = options.bmax; 
   
else
    error('options.bmax must be provided.');
end
nlabs = parpool_size();
assert(length(bmax)==nlabs,'bmax should be per-worker max');
maxIter    = getoption(options,'maxIter',5);
innerit    = getoption(options,'innerit',1);
b0         = getoption(options,'b0',ones(size(bmax)));
binc       = getoption(options,'binc',zeros(size(bmax)));
seed       = getoption(options,'seed',1);
rand_draw  = getoption(options,'rand_draw',@(x) rand_set(x,bmax));
solver     = getoption(options,'solver','lbfgs');
itr_save   = getoption(options,'itr_save',[]);
true_sol   = getoption(options,'true_sol',[]);
resume_previous = getoption(options,'resume_previous',false);



% initialization
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);
          
converged = 0;
iter = 0;
bk   = b0;

% initial evaluation
x = vec(x);

proj = @(x) projectBounds(x,LB,UB);

opts = struct;
if strcmp(solver,'gn')
    opts.maxIter = innerit;
else
    opts.maxIter = ceil(innerit*max(bmax./bk)-1);
end
opts.verbose = 0;
opts.true_sol = true_sol;

maxevals = innerit*maxIter;

xhist = cell(maxIter,1);

if resume_previous && ~isempty(itr_save) && exist(itr_save,'file')
    load(itr_save);
    iter = find(~cellfun('isempty',xhist),1,'last');
    disp('Resuming previous progress');
    x = xhist{iter};
end

fevals = 0;
num_failures = 0;
model_err = [];
% main loop
while ~converged
    Ik = rand_draw(bk);
    func = @(x) fh(x,Ik);
    switch solver
      case 'spg'
        [x,fsave,funEvals,projects,iter_save] = minConf_SPG(func,x,proj,opts);
      case 'lbfgs'
        [x,fsave,funEvals,err,iter_save,misfit_aux] = minConf_TMP(func,x,LB,UB,opts);
      case 'gn'
        x = gauss_newton(func,x,LB,UB,opts);
    end     
    
    fevals = fevals + funEvals*(max(bk./bmax));
    disp(['iter : ' num2str(iter) ' fstart  ' num2str(fsave(1),'%3.3e') ' fend ' num2str(fsave(end),'%3.3e') ' feval budget ' num2str(fevals/maxevals*100) '%']);
    
    % If the function hasn't decreased enough, increase the points used
    if sum(binc) > 0 && (length(fsave)<=2 || (max(fsave)-min(fsave)) < 0.05*min(fsave)   )
        bk = min(bk + binc,bmax);     
        disp('Stagnation, increasing number of summands');
        num_failures = num_failures + 1;
    else
        num_failures = 0;
    end
    
    iter = iter+1;
    if ~isempty(itr_save)
        xhist{iter} = x;
        if exist(itr_save,'file')
            save(itr_save,'-append','xhist');    
        else
            save(itr_save,'-v7.3','xhist');
        end        
    end
    converged = iter >= maxIter || (fevals >= maxevals) || (num_failures >= 2);
    if not(isempty(true_sol))
        model_err = [model_err err];
    end
end

end



