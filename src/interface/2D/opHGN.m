classdef opHGN < opSpot
% SPOT wrapper for the Gauss-Newton Hessian of F.m
%
% use:
%   H = opHGN(m,Q,model)
%
% see PDEfunc.m for further documentation
%
    properties
        mt,Q,model,params;
    end
        
    methods

       function op = opHGN(mt,Q,model,params)            
           n = length(mt);                                       
           op = op@opSpot('oppHGN', n,n);
           op.cflag     = 0;  
           op.linear    = 1;
           op.children  = []; 
           op.sweepflag = 0;
           op.mt        = mt;
           op.Q         = Q;
           op.model     = model;           
           if exist('params','var')==0||isempty(params)
               params = struct; 
           end
           params = default_fwi_params2d(params);
                      
           params.wri = false;
           op.params = params;
       end                             
    end
        
    methods ( Access = protected )
       function y = multiply(op,x,~)           
           y = PDEfunc(PDEopts.HESS_GN,op.mt,op.Q,x,[],op.model,op.params);                  
       end        
    end     
end 

    
    
    
