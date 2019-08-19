function [ rho ] = comp_rho( n, X )
%COMP_RHO Computes the bound rho from Stewart '89 and O'Leary '90 for
% Projections: W  = X(X'DX)^{-1} X'D, where D is p.d diagonal.

%{

    INPUTS:
    n : Problem dimension
    X : Range of W.

    OUTPUT:
    rho: Bound

%}

Q       = orth(X);

comb    = 2^n;

col     = zeros(n,1);
SM      = zeros(n,comb);
SM      = select_oblique(n,0,0,col,SM);

rho     = 1;

for i = 1:comb
    
    Sel     = diag(SM(:,i));
    tr      = trace(Sel);     
    
    if tr > 0
        sig     = svd((Sel*Q));
    
        sig_p   = sig(abs(sig)>1e-10);
    
        sig_c = min(sig_p);
        
        if sig_c <= rho
            rho = sig_c;
        end
        
    end
    
end

end

