function x = state(i, K, phi, Gamma, L, x0)
    f_term = (phi + Gamma*K)^(i+1)*x0;
    s_term = 0; 
    for j = 0:i
        s_term = s_term + (phi + Gamma*K)^i;
    end
    s_term = s_term*Gamma*L;
    x = f_term + s_term;
end