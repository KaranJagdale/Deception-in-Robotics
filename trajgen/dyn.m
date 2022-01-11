function xdot = dyn(t,x, K, l)
    
    xdot = K*x + l;
end