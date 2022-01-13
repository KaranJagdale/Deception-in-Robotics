function Kmat = vectom(x,n)
    
    Kmat = zeros(n,n);
    for i = 1:n
       Kmat(i,:) = x(1+n*(i-1):n*(i-1)+n); 
    end    
end