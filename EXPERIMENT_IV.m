%-------------------------- AOP EXPERIMENT IV ----------------------------%
% 08/16/18, J.B.,
% Code accompanying the manuscript "On the Eigendecomposition and Singular 
% Value Decomposition of Oblique Projection Matrices", J.J.Brust,
% R.F.Marcia, C.G.Petra, 2018.

% This experiment uses Algorithm 1 to compute the spectral norms of the
% oblique projections 
%
% W = X (X' D X)^{-1} X' D,
%
% as in as in (Stewart,1989) or (O'Leary 1990). Here 
% X (n x m) and D (m x m) p.d. A comparison to a bound proposed in the
% two articles is made.
%
%  || W ||_2 <= 1/rho = 1/ min_I inf_+ (U_I),
%
% where U is an orthonormal basis of X.
%
% Versions:
% 08/14/18, Tests for manuscript. Run over specified values of n with
% fixed ratio to m, in particular m = floor(n/scale)
% 08/29/18, J.B. Reruns for multiple repetitions, and storing of data.
% 09/09/18, J.B., Preparations for release.
% 10/28/18, J.B., Renaming: EXPERIMENT III <- EXPERIMENT IV

clc;
clear;

fprintf('------- AOP: Algorithms For Oblique Projection Matrices ------- \n');
fprintf('-------      J.J.Brust, R.F.Marica, C.G.Petra, 2018     ------- \n \n');

fprintf('EXPERIMENT IV: Spectral Norm Computations \n');
fprintf('ALGORITHMS: ALG. 1 and Bound (Stewart,O''Leary) \n');

fprintf(' n \t m \t ALG. 1 \t MTLB. \t \t BND. \n');

ns      = [6,7,8,9,10,11,12,20];
numns   = length(ns);

sigs1    = zeros(numns,1);
sigsm    = zeros(numns,1);
bnds     = zeros(numns,1);
ms       = zeros(numns,1);

nreps    = 10; % 1
data     = zeros(numns,3,nreps);

idxrep  = 5; % 1
scale   = 3;


tic;
for i_rep = 1:nreps
    
    for i = 1:numns
        
        n = ns(i);
        m = floor(n/scale);
        
        Im  = eye(m);
        In  = eye(n);
        X   = randn(n,m);
        D   = abs(randn(n,1));
        Y   = diag(D)*X;
        
        Q   = orth(X);
        
        W   = X*((X'*Y)\Y');
        QW  = Q'*(In-W);
        
        gamm2   = eigs((QW*QW'),1);
        sig1    = sqrt(1+gamm2);
        
        sigm    = max(svd(W));
        
        rho     = comp_rho(n,X);
        rho_i   = 1/rho;
        
        sigs1(i)    = sig1;
        sigsm(i)    = sigm;
        bnds(i)     = rho_i;
        ms(i)       = m;
        
        if i_rep == idxrep
       
            fprintf(' %i \t %i \t %1.4f \t %1.4f \t %4.4f \n',n, m, sig1, sigm, rho_i);
 
        end
        
    end
    
    data(:,1,i_rep) = sigs1(:);
    data(:,2,i_rep) = sigsm(:);
    data(:,3,i_rep) = bnds(:);
    
end

t_ex3 = toc;
fprintf('\n Time to run EXPERIMENT IV: %1.4f s. \n', t_ex3);

sigs1_idx   = data(:,1,idxrep);
sigsm_idx   = data(:,2,idxrep);
bnds_idx    = data(:,3,idxrep);

save('EXPERIMENT_IV','sigs1_idx','sigsm_idx','bnds_idx','data','ms');