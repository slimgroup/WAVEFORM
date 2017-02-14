function [H,dH_forw,dH_adj,DdH_adj] = Poiss2D(epsilon,n)
% Finite volume discretization of the Poisson equation
%
%   div(epsilon grad(u)) = q
% 
%   
% From http://www.ece.utah.edu/~ece6340/LECTURES/Feb1/Nagel%202012%20-%20Solving%20the%20Generalized%20Poisson%20Equation%20using%20FDM.pdf
%
% Usage:
%   [H,dH_forw,dH_adj] = Poiss2D(epsilon,n);
%
% Input:
%   epsilon - nz x nx
% 
% Output:
%   H       - nz*nx x nz*nx system matrix
%   dH_forw - function of (deps,u), derivative mapping deps -> DH(eps)[deps]u
%   dH_adj  - function of (z,u), adjoint of dH_forw
    

    assert(numel(epsilon)==prod(n(1:2)),'epsilon must be nz x nx');
    nz = n(1); nx = n(2); 
    
    int_z = ones(nz,1); int_x = ones(nx,1);
    [IZ,IX] = ndgrid(int_z,int_x);
    
    Rext = opKron(opExtension(nx,[1 1],1),opExtension(nz,[1 1],1));
    Rext_adj = opKron(opRestriction(nx+2,2:(nx+1)),opRestriction(nz+2,2:(nz+1)));
    e = Rext*vec(epsilon);
    Rz = @(offset) opRestriction(nz+2,(2+offset):(nz+1+offset));
    Rx = @(offset) opRestriction(nx+2,(2+offset):(nx+1+offset));
    R = @(oz,ox) opKron(Rx(ox),Rz(oz));        
    
    b0 = -(R(0,0)+R(-1,0)+R(0,-1)+R(-1,-1));
    b1 = opDiag_swp(vec(IZ))*(1/2*(R(0,0)+R(0,-1)));
    b2 = opDiag_swp(vec(IX))*(1/2*(R(0,0)+R(-1,0)));
    b3 = opDiag_swp(vec(IZ))*(1/2*(R(-1,-1)+R(-1,0)));
    b4 = opDiag_swp(vec(IX))*(1/2*(R(0,-1)+R(-1,-1)));
    e1 = ones(prod(n),1);    
    
    A = @(oz) spdiags(e1,oz,prod(n),prod(n));
    sdiag = @(x) spdiags(x,0,length(x),length(x));
    H = A(0)*sdiag(b0*e) + A(1)*sdiag(b1*e) + A(nz)*sdiag(b2*e) + A(-1)*sdiag(b3*e) + A(-nz)*sdiag(b4*e);
    dH_forw = @(de,u) (A(0)*sdiag(b0*(Rext*de)) + A(1)*sdiag(b1*(Rext*de)) + A(nz)*sdiag(b2*(Rext*de)) + A(-1)*sdiag(b3*(Rext*de)) + A(-nz)*sdiag(b4*(Rext*de)))*u;
    dH_adj = @(z,u) Rext_adj*(b0'*(conj(u) .* (A(0)'*z)) + b1'*(conj(u) .* (A(1)'*z)) + b2'*(conj(u) .* (A(nz)'*z)) + b3'*(conj(u) .* (A(-1)'*z)) + b4'*(conj(u) .* (A(-nz)'*z)));
    DdH_adj = @(u,du,dm,z) Rext_adj*(b0'*(conj(du) .* (A(0)'*z)) + b1'*(conj(du) .* (A(1)'*z)) + b2'*(conj(du) .* (A(nz)'*z)) + b3'*(conj(du) .* (A(-1)'*z)) + b4'*(conj(du) .* (A(-nz)'*z))); 
    
end
