classdef opDF < opSpot
% SPOT wrapper for the Jacobian of F.m
%
% use:
%   J = opDF(m,Q,model)
%
% see PDEfunc.m for further documentation
%
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

    properties
        mt,Q,model,nfreq,nt,params;
    end    
    
    methods
       function op = opDF(mt,Q,model,params)
           nsrc  = size(Q,2);
           nrec  = length(model.xrec)*length(model.zrec);
           nfreq = length(model.freq);
           m = nsrc*nrec*nfreq;
           n = length(mt);           
           
           op = op@opSpot('opDF', m, n);
           op.cflag     = 1;
           op.linear    = 1;
           op.children  = [];
           op.sweepflag = 0;
           op.mt        = mt;
           op.Q         = Q;
           op.model     = model;
           op.nfreq     = nfreq;
           op.nt        = nsrc*nrec;
           if exist('params','var')==0 || isempty(params)
               params = struct;
           end
           params.wri = false;
           params = default_fwi_params2d(params);
           op.params = params;
       end        
    end
        
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
           if mode == 1
               y = PDEfunc(PDEopts.JACOB_FORW,op.mt,op.Q,x,[],op.model,op.params);
           else
               y = PDEfunc(PDEopts.JACOB_ADJ,op.mt,op.Q,x,[],op.model,op.params);
           end   
           y = vec(y);
       end        
    end     
end 

    
    
    
