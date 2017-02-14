classdef LinSolveOpts < dynamicprops & matlab.mixin.Copyable
%LINSOLVEOPTS Parameter object for solving linear systems
%
%  
% Options
%   .tol          - relative residual tolerance (default: 1e-6)
%
%   .maxit        - maximum (outer) iterations (default: 10000)
%
%   .maxinnerit   - maximum (inner) iterations, for certain solvers (default: 10)
%
%   .solver       - linear solver to use, one of
%      LinSolveOpts.SOLVE_GMRES     - GMRES
%      LinSolveOpts.SOLVE_FGMRES    - GMRES w/ a flexible preconditioning option
%      LinSolveOpts.SOLVE_LU        - LU decomposition (explicit matrices, 2D problems)
%      LinSolveOpts.SOLVE_BACKSLASH - Matlab's backslash
% 
%   .precond      - preconditioner to use, one of
%      LinSolveOpts.PREC_IDENTITY - identity preconditioner (default)
%      LinSolveOpts.PREC_MLGMRES  - multi-level GMRES preconditioner
%  
%      OR 
%      a LinSolveOpts object, which specifies an iterative solver to use 
%     
    
    properties
        tol, maxit, maxinnerit, solver, precond,params;
    end
    
    properties (Constant)
        SOLVE_GMRES    = 'solve_gmres';        
        SOLVE_FGMRES   = 'solve_fgmres';
        SOLVE_LU       = 'solve_lu';
        SOLVE_BACKSLASH= 'solve_backslash';
        PREC_IDENTITY  = 'prec_identity';        
        PREC_MLGMRES   = 'prec_mlgmres';
        PREC_SHIFTLAP  = 'prec_shiftlap';
        PREC_CALANDRA12 = 'prec_calandra12';
    end
    
    methods
        function opts = LinSolveOpts()
            opts.tol        = 1e-6;
            opts.maxit      = 10000;
            opts.maxinnerit = 5;
            opts.solver     = LinSolveOpts.SOLVE_FGMRES;
            opts.precond    = LinSolveOpts.PREC_IDENTITY;
            opts.params     = struct;
        end
    end
    
end

