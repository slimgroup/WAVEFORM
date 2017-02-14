function [f,g] = quad_model(gk,hk,x,xk)
    d = x-xk;
    hd = hk*d;
    f = gk'*d + 0.5*d'*hd;
    g = gk + hd;
end